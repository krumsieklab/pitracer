library(tidyverse)
library(visNetwork)
library(readxl)
library(data.table)

work_dir <- getwd()

source(glue::glue(work_dir, "/src/tracing/tr_load_vis_helpers.R"))
source(glue::glue(work_dir, "/src/tracing/tr_load_vis_objects.R"))

selectize_style <- "z-index:1;
                    min-width:160px;
                    width:200px;
                    height:26px;
                    margin:2px 0 0;
                    list-style:none;
                    font-size:13px;
                    text-align:left;
                    background-color:#ffffff;
                    border:1px solid #cccccc;
                    border:1px solid rgba(0,0,0,0.d5);
                    border-radius:3px;
                    -webkit-background-clip:padding-box;
                    background-clip:padding-box;
                    display:block;"


# Visualization functions --------------------------------------------------


visualize_trace_gene_to_gene <- function(path_row, show_clusters = FALSE) {


  # Check if cluster_num column exists, if not add it ---------------------------

  if (!"cluster_num" %in% colnames(path_row)) {
    stop("cluster_num column missing. Consider adding it, even as a dummy column")
  }


  # Proceed with creating the visualization ---------------------------------

  sprintf("visualizing trace: %s",
          tr_gen_get_endpoints(path_row)) %>%
    print()

  source_nodes <- path_row$gene %>% unique()
  target_nodes <- path_row$met_ids %>% unique()
  shortest_path <- path_row$association_path %>%
    sapply(length) %>%
    min()
  main_nodes <<-
    path_row$association_path %>%
    unlist() %>%
    unique()

  # Do this for interactive blacklisting
  main_nodes_dict <-
    tibble(original_name = main_nodes,
           translated_name = main_nodes,
           endpoint_node = if_else(original_name %in% c(source_nodes, target_nodes), TRUE, FALSE))

  edges <-
    path_row %>%
    unnest_path() %>%
    add_moa() %>%
    translate_generic_moa() %>%
    dplyr::mutate(edge_type = "main") %>%
    ungroup() %>%
    group_by(rxn, from, to, edge_type) %>%
    summarise(cluster_nums = str_c(cluster_num, collapse = ", ")) %>%
    count_rxn_occ()

  omnipath_refs <- edges %>%
    filter(grepl("omnipathdb", edges$rxn)) %>%
    select(rxn, from, to) %>%
    left_join(omnipathdb_edges, by = c("from" = "source_genesymbol", "to" = "target_genesymbol")) %>%
    rowwise() %>%
    mutate(sources = str_split(sources, ";") %>% unlist() %>% list()) %>%
    unnest(c(sources)) %>%
    rowwise() %>%
    mutate(references = str_split(references, ";") %>% unlist() %>% list()) %>%
    unnest(c(references)) %>%
    separate(references, c("reference", "reference_id"), ":")

  omnipath_refs$reference <- replace(omnipath_refs$reference,
                                     omnipath_refs$reference != omnipath_refs$sources,
                                     NA)

  omnipath_refs$reference_id <- replace(omnipath_refs$reference_id,
                                        is.na(omnipath_refs$reference),
                                        NA)

  omnipath_refs <- omnipath_refs %>%
    distinct(rxn, from, to, sources, reference, .keep_all = TRUE)

  omnipath_refs$url <- apply(omnipath_refs, 1, function(row) {
    if(is.na(row["reference_id"])) {
      omnipathdb_source_list %>%
        .[row["sources"]] %>%
        unlist() %>%
        paste0(collapse = ";")
    } else {
      sprintf("https://pubmed.ncbi.nlm.nih.gov/%s", row["reference_id"])
    }
  })

  main_reactions <- edges %>% distinct(rxn)


  # Convert to final edge dataframe ---------------------------------

  vis_edges <-
    bind_rows(edges %>% transmute(from, to = rxn, edge_type, cluster_nums),
              edges %>% transmute(from = rxn, to, edge_type, cluster_nums)) %>%
    dplyr::mutate(ggm_edge = FALSE)


  # Convert to final nodes dataframe ----------------------------------------

  vis_nodes <-
    edges %>%
    gather(key = "tmp", value = "id",
           -c(edge_type, cluster_nums, rxn_rep)) %>%
    transmute(id,
              reversible = FALSE,
              node_type = if_else(tmp == "rxn", "rxn", "gene"),
              node_priority = edge_type,
              crowded = FALSE,
              cluster_nums) %>%
    unnest(cluster_nums = str_split(cluster_nums, ", ")) %>%
    distinct() %>%
    group_by(id) %>%
    mutate(cluster_nums = str_c(cluster_nums, collapse = ", ")) %>%
    distinct()

  vis_nodes$title <- vis_nodes %>%
    apply(1, function(row) {
      if (row["node_type"] == "gene") {
        sprintf("<h4>%s</h4><br><b>References:</b><br>GeneCards: <a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s' target='_blank'>[1]</a>",
                row["id"],
                row["id"])
      } else if (row["node_type"] == "rxn" && !(row["id"] %in% omnipath_refs$rxn)) {
        sprintf("<b>%s</b><br>No references", row["id"])
      } else {
        omnipath_refs %>%
          filter(rxn == row[["id"]]) %>%
          rowwise() %>%
          mutate(url = str_split(url, ";") %>% unlist() %>% list()) %>%
          unnest(c(url)) %>%
          split(.$sources) %>%
          sapply(function(df) {
            if (!all(is.na(df$reference))) {
              df <- df[which(!is.na(df$reference)), ]
            }
            df$url %>%
              sapply(function(url) {
                sprintf("<a href='%s' target='_blank'>%s</a>",
                        url,
                        paste0("[", match(url, df$url), "]")
                )
              }) %>% unname()
          }, simplify = FALSE) %>%
          lapply(function(source_urls) {
            paste(source_urls, collapse = ", ")
          }) %>%
          paste(names(.), ., sep = ": ", collapse = "<br>") %>%
          paste(paste0("<b>", row[["id"]], "</b>"), "References:", ., sep = "<br>")
      }
    })



  # Format edges -----------------------------------------------------------


  vis_edges <-
    vis_edges %>%
    dplyr::mutate(
      width = case_when(
        edge_type == "main" ~ edge_main_width,
        edge_type == "main_secondary" ~ edge_main_secodary_width,
        TRUE ~ edge_minor_width
      ),

      color = case_when(
        edge_type == "main" ~ edge_main_color,
        edge_type == "main_secondary" ~ edge_main_secondary_color,
        TRUE ~ edge_minor_color
      ),

      dashes = if_else(str_detect(from, "_ggm") | str_detect(to, "_ggm"),
                       TRUE,
                       FALSE)
    )



  # Format nodes ------------------------------------------------------------

  vis_nodes <-
    vis_nodes %>%
    dplyr::mutate(
      color.background = case_when(
        id %in% source_nodes ~ node_source_color,
        id %in% target_nodes ~ node_target_color,

        node_type == "recon" & node_priority == "main" ~ node_metabolite_recon_color,
        node_type == "ggm"   & node_priority == "main" ~ node_metabolite_ggm_color,
        node_type == "rxn"   & node_priority == "main" ~ node_action_color,
        node_type == "gene"  & node_priority == "main" ~ node_gene_color,

        node_type == "recon" & node_priority == "main_secondary" ~ node_metabolite_recon_color,
        node_type == "gene"  & node_priority == "main_secondary" ~ node_gene_color,

        node_type == "recon" & node_priority == "minor" ~ node_minor_metabolite_color,
        node_type == "ggm"   & node_priority == "minor" ~ node_minor_metabolite_color,
        node_type == "rxn"   & node_priority == "minor" ~ node_minor_action_color,
        node_type == "gene"  & node_priority == "minor" ~ node_minor_gene_color,


        TRUE ~ "black"
      ),
      borderWidth = if_else(crowded | reversible, 2.5, 0),
      color.border = case_when(
        crowded ~ node_crowded_border_color,
        node_priority == "main"  & reversible ~ node_reversible_border_color,
        node_priority == "minor" & reversible ~ node_minor_reversible_border_color,
        TRUE ~ color.background
      ),
      color.highlight.border  = color.border,
      color.highlight.background  = color.background,
      size = case_when(
        id %in% c(source_nodes, target_nodes) ~ node_endpoint_size,
        node_priority == "main" ~ node_main_size,

        node_type == "recon" & node_priority == "main_secondary" ~ node_main_secondary_size,
        node_type == "gene"  & node_priority == "main_secondary" ~ node_main_secondary_size,

        node_type == "rxn" ~ node_rxn_size,
        TRUE ~ node_minor_size
      ),
      font.size = case_when(
        id %in% c(source_nodes, target_nodes) ~ node_endpoint_font_size,
        node_priority == "main" ~ node_main_font_size,
        node_priority == "main_secondary" ~ node_main_secondary_font_size,
        TRUE ~ node_minor_font_size),
      font.color = case_when(
        node_priority == "main" ~ node_main_font_color,
        node_priority == "main_secondary" ~ node_main_secondary_font_color,
        TRUE ~ node_minor_font_color),
      shape = case_when(
        node_type == "rxn" ~ node_rxn_shape,
        id %in% c(source_nodes, target_nodes) ~ node_endpoint_shape,
        TRUE ~ node_all_others_shape
      )
    )




  # Display clusters setting -------------------------------------------------

  if (show_clusters) {

    selected_by <- list(variable = "cluster_nums",
                        main = "Select path by cluster",
                        multiple = TRUE,
                        style = selectize_style,
                        highlight = TRUE)

  } else {
    selected_by <- NULL
  }



  # Visualize network -------------------------------------------------------

  recon_vis <-
    visNetwork(nodes = vis_nodes,
               edges = vis_edges,
               height = "800px",
               width = "1200px") %>%
    visIgraphLayout(layout = "layout_with_fr", randomSeed = 123) %>%
    visEdges(arrows = "to",
             color = "gray") %>%
    visOptions(
      nodesIdSelection = list(enabled = TRUE,
                              main = "Select node by name",
                              style = selectize_style
      ),

      selectedBy = selected_by,

      highlightNearest = list(enabled = TRUE,
                              labelOnly = FALSE,
                              algorithm = "hierarchical",
                              degree = 0)) %>%
    visLayout(randomSeed = preset_seed)

  list(visualization = recon_vis,
       trace_name = sprintf("%s ---> %s (%s)", source_nodes, target_nodes, shortest_path),
       main_nodes_dict = main_nodes_dict)

}


visualize_trace_non_gene_to_gene <- function(path_row, show_clusters = FALSE) {

  # Check if cluster_num column exists, if not add it ---------------------------


  if (!"cluster_num" %in% colnames(path_row)) {
    stop("cluster_num column missing. Consider adding it, even as a dummy column")
  }



  # Proceed with creating the visualization ---------------------------------

  sprintf("visualizing trace: %s",
          tr_gen_get_endpoints(path_row)) %>%
    print()
  source_nodes <- path_row$gene %>% unique()
  target_nodes <- path_row$met_ids %>% unique()
  shortest_path <- path_row$association_path %>%
    sapply(length) %>%
    min()
  main_nodes <<-
    path_row$association_path %>%
    unlist() %>%
    unique()

  # missing reversibility information
  edges <-
    path_row %>%
    unnest_path() %>%
    add_moa() %>%
    translate_generic_moa() %>%
    dplyr::mutate(edge_type = "main") %>%
    ungroup() %>%
    group_by(rxn, from, to, edge_type) %>%
    summarise(cluster_nums = str_c(cluster_num, collapse = ", ")) %>%
    count_rxn_occ()


  main_reactions <- edges %>% distinct(rxn)


  edges_expanded <-
    edges %>%
    fill_recon() %>%
    # some reversible reactions might have been traversed in
    # reverse, so must readjust the recon reactions that are
    # added by fill_recon()
    orient_rxns() %>%
    expand_reactions()

  # secondary here means all metabolites/genes participating
  # in the main reactions
  edges_2ndary_met <- get_secondary_mets(edges_expanded)


  edges_2ndary_genes <-
    edges_expanded %>%
    dplyr::transmute(rxn, rxn_rep, cluster_nums, genes, edge_type = "main_secondary") %>%
    unnest(genes = str_split(genes, ", ")) %>%
    distinct() %>%
    filter(genes != "",
           !genes %in% source_nodes) %>%
    dplyr::transmute(rxn, rxn_rep, cluster_nums, from = genes, to = rxn,
                     edge_type, ggm_edge = FALSE)



  edges_neighbor_genes <-
    recon_rxn_genes %>%
    unnest(genes) %>%
    # unnest(c(genes)) %>%
    filter(!genes %in% main_nodes,
           !rxn %in% main_reactions) %>%
    dplyr::transmute(rxn, rxn_rep = 1,
                     from = genes, to = rxn,
                     edge_type = "minor", ggm_edge = FALSE) %>%
    distinct()

  edges_neighbor_mets <-
    rxn_steps %>%
    dplyr::mutate(from = reactant,
                  to = product,
                  from_has_main = from %in% main_nodes,
                  to_has_main = to %in% main_nodes) %>%
    # this step already ensures that main edges are excluded
    # since those are the only reactions that have the main
    # reaction steps. Double checked with anti_join below
    anti_join(main_reactions, by = "rxn") %>%
    dplyr::group_by(rxn) %>%
    filter(any(from_has_main) | any(to_has_main)) %>%
    ungroup() %>%
    dplyr::transmute(rxn, rxn_rep = 1, from, to, edge_type = "minor") %>%
    distinct()


  # Combine metabolite edges ---------------------------------------------------

  edges_all <- bind_rows(edges, edges_2ndary_met)

  if (show_neighboring_rxns) {

    edges_all <- bind_rows(edges_all, edges_neighbor_mets)

  }


  # Combine gene edges ------------------------------------------------------

  edges_genes <- bind_rows(edges_2ndary_genes, edges_neighbor_genes)

  # Get maximum 5 reactions per main node -----------------------------------

  n_show <- 5

  rxn_nodes_main <-
    edges_all %>%
    gather(directionality, molname, -c(rxn, rxn_rep, edge_type)) %>%
    filter(molname %in% main_nodes) %>%

    # select only metabolites
    dplyr::arrange(molname, edge_type) %>%

    # make sure that only a single reaction is associated to
    # a metabolite. Can do debugging using previous lines, if
    # you want to see whether the nodes were sorted properly
    distinct(rxn, molname)

  rxn_node_count <- dplyr::count(rxn_nodes_main, molname)

  rxn_nodes_selected <-
    rxn_nodes_main %>%
    distinct(molname, rxn) %>%
    group_by(molname) %>%
    dplyr::slice(1:n_show) %>%
    # NOTE that some metabolites might share a rxn,
    # come back here if there is an error downstream in
    # the visualization
    .$rxn %>%
    unique()

  ### Do the shuffling
  ### shuffle first and then sort by edge_type (ignoring from and to)


  # Narrow down number of edges for visualization ---------------------------

  edges_all <-
    edges_all %>%
    filter(rxn %in% rxn_nodes_selected) %>%
    add_reversibility()

  ## select gene edges based on rxn_nodes_selected here

  edges_genes <-
    edges_genes %>%
    filter(rxn %in% rxn_nodes_selected) %>%
    distinct() %>%
    # transform into gathered_edge format
    dplyr::transmute(rxn, rxn_rep, cluster_nums,
                     edge_type, directionality = "from",
                     molname = from, ggm_edge) %>%
    add_reversibility()



  # Get interface metabolites so can orient them properly -------------------

  interface_mets <-
    edges_all %>%
    filter(!str_detect(from, "\\[|__ggm") & str_detect(to, "\\[|__ggm") ) %>%
    left_join(edges_expanded %>% dplyr::distinct(rxn, reactant), by = "rxn") %>%
    group_by(rxn) %>%
    dplyr::mutate(direction_sub = if_else(to %in% reactant, "from", "to")) %>%
    ungroup() %>%
    dplyr::distinct(rxn, molname = to, direction_sub)


  # get readable names for mets and rxns ------------------------------------


  gathered_edges <-
    gather_edges(edges_all) %>%
    fix_interface_direction(interface_mets) %>%

    # NOTE: THE FOLLOWING IS PROBLEMATIC FOR INTERFACE METABOLITES
    # have to slice, since gathered_edges will have redundant
    # rows that only differ in edge_type
    arrange(molname, rxn, rxn_rep, edge_type) %>%
    group_by(molname, rxn, rxn_rep) %>%
    dplyr::slice(1) %>%
    ungroup()

  # Do this for interactive blacklisting
  main_nodes_dict <-
    tibble(original_name = main_nodes,
           translated_name = main_nodes,
           endpoint_node = if_else(original_name %in% c(source_nodes, target_nodes), TRUE, FALSE))

  # completely translate gathered_edges
  # can have a flag for this later

  if (translate_names) {

    gathered_edges <-
      gathered_edges %>%
      translate_metabolite_name(recon_lookup, "molname") %>%
      translate_rxn_name()

    edges_genes <-
      edges_genes %>%
      translate_rxn_name()

    main_nodes_dict <-
      tibble(original_name = main_nodes,
             translated_name =
               main_nodes %>%
               translate_metabolite_name_list(recon_lookup),
             endpoint_node = if_else(original_name %in% c(source_nodes, target_nodes), TRUE, FALSE))

    # translate main node name if translating edge names
    source_nodes <- translate_metabolite_name_list(source_nodes, recon_lookup)
    target_nodes <- translate_metabolite_name_list(target_nodes, recon_lookup)
    main_nodes <-  translate_metabolite_name_list(main_nodes, recon_lookup)

    main_reactions <- translate_rxn_name(main_reactions)
    rxn_node_count <-
      rxn_node_count %>%
      translate_metabolite_name(recon_lookup, "molname")

  }


  crowded_mets <-
    rxn_node_count %>%
    filter(n > 5) %>%
    .$molname

  # annotate edges ----------------------------------------------------------

  gathered_edges <-
    gathered_edges %>%
    add_metabolite_numbering(exclude_nodes = main_nodes) %>%
    add_ggm_edge_anno()

  if (remove_minor_met_nodes) {

    gathered_edges <- filter(gathered_edges, edge_type != "minor")

  }

  if (remove_main_secondary_met_nodes) {

    gathered_edges <- filter(gathered_edges, edge_type != "main_secondary")

  }

  if (remove_minor_gene_nodes) {

    edges_genes <- filter(edges_genes, edge_type != "minor")

  }

  if (remove_main_secondary_gene_nodes) {

    edges_genes <- filter(edges_genes, edge_type != "main_secondary")

  }

  gathered_edges <- bind_rows(gathered_edges, edges_genes)


  # Convert to final edge dataframe ---------------------------------


  vis_edges <-
    gathered_edges %>%
    dplyr::mutate(rxn = if_else(rxn_rep > 1,
                                str_c(rxn, " (#", rxn_rep, ")"),
                                rxn)) %>%
    dplyr::transmute(from = if_else(directionality == "from",
                                    molname,
                                    rxn),
                     to = if_else(directionality == "to",
                                  molname,
                                  rxn),
                     edge_type,
                     ggm_edge,
                     cluster_nums) %>%
    # remove double edges with differet edge type
    arrange(from, to, ggm_edge, edge_type) %>%
    group_by(from, to, ggm_edge) %>%
    dplyr::slice(1)

  # Convert to final nodes dataframe ----------------------------------------


  nodes_rxn <-
    gathered_edges %>%
    dplyr::transmute(id = rxn,
                     rxn_rep,
                     reversible,
                     node_type = "rxn",
                     node_priority = if_else(id %in% main_reactions$rxn,
                                             "main",
                                             "minor"),
                     crowded = FALSE,
                     cluster_nums) %>%
  dplyr::mutate(title = edges_expanded[match(id, edges_expanded$rxn_name_short), c("rxn", "rxn_name")] %>%
  {
    ifelse(is.na(.$rxn),
           "No references available",
           sprintf("<h5>%s</h5><br><b>References:</b><br>Virtual Metabolic Human: <a href='https://www.vmh.life/#reaction/%s' target='_blank'>[1]</a>",
                   .$rxn_name,
                   .$rxn)
    )
  },
  id = if_else(rxn_rep > 1,
               str_c(id, " (#", rxn_rep, ")"),
               id)
  ) %>%
    distinct() %>%
    dplyr::select(-rxn_rep)

  nodes_mol <-
    gathered_edges %>%
    dplyr::transmute(id = molname,
                     node_priority = edge_type,
                     node_type = case_when(
                       str_detect(molname, "\\[[a-z]\\]") ~ "recon",
                       str_detect(molname, "__ggm") ~ "ggm",
                       TRUE ~ "gene"
                     ),
                     reversible = FALSE,
                     crowded = if_else(id %in% crowded_mets,
                                       TRUE,
                                       FALSE),
                     cluster_nums
    ) %>%
    dplyr::mutate(title = id %>%
                    sapply(function(lst) {
                      lst %>%
                        str_split(" \\(#") %>%
                        lapply("[[", 1) %>%
                        unlist() %>%
                        str_split(" \\[[a-z]\\]+") %>%
                        lapply("[[", 1) %>%
                        sapply(
                          function(node_name) {
                            ifelse(node_name %in% genes,
                                   node_name,
                                   full_metabolite_database$recon_id[match(node_name, full_metabolite_database$name)])
                          }
                        )
                    }) %>%
                    data.frame(name = names(.),
                               id = .,
                               stringsAsFactors = FALSE) %>%
                    {
                      ifelse(is.na(.$id),
                             "No references available",
                             sprintf("<h5>%s</h5><br><b>References:</b><br>%s: <a href='%s%s' target='_blank'>[1]</a>",
                                     .$name,
                                     ifelse(node_type == "recon",
                                            "Virtual Metabolic Human",
                                            "GeneCards"),
                                     ifelse(node_type == "recon",
                                            "https://www.vmh.life/#metabolite/",
                                            "https://www.genecards.org/cgi-bin/carddisp.pl?gene="),
                                     .$id)
                      )
                    }
    ) %>%
    distinct() %>%
    unnest(cluster_nums = str_split(cluster_nums, ", ")) %>%
    group_by_at(vars(-c(cluster_nums, node_priority))) %>%
    # Doing the following will ensure that cluster_nums are selected properly,
    # even if they assign main_secondary to certain nodes
    summarise(node_priority = unique(node_priority) %>% str_c(collapse = ", "),
              cluster_nums = unique(cluster_nums) %>% str_c(collapse = ", ")) %>%
    unnest(node_priority = str_split(node_priority, ", ")) %>%
    ungroup() %>%
    arrange(id, node_priority) %>%
    dplyr::group_by(id) %>%
    dplyr::slice(1)


  vis_nodes <- bind_rows(nodes_rxn, nodes_mol)

  # Format edges -----------------------------------------------------------


  vis_edges <-
    vis_edges %>%
    dplyr::mutate(
      width = case_when(
        edge_type == "main" ~ edge_main_width,
        edge_type == "main_secondary" ~ edge_main_secodary_width,
        TRUE ~ edge_minor_width
      ),

      color = case_when(
        edge_type == "main" ~ edge_main_color,
        edge_type == "main_secondary" ~ edge_main_secondary_color,
        TRUE ~ edge_minor_color
      ),

      dashes = if_else(str_detect(from, "_ggm") | str_detect(to, "_ggm"),
                       TRUE,
                       FALSE)
    )



  # Format nodes ------------------------------------------------------------

  vis_nodes <-
    vis_nodes %>%
    dplyr::mutate(
      color.background = case_when(
        id %in% source_nodes ~ node_source_color,
        id %in% target_nodes ~ node_target_color,

        node_type == "recon" & node_priority == "main" ~ node_metabolite_recon_color,
        node_type == "ggm"   & node_priority == "main" ~ node_metabolite_ggm_color,
        node_type == "rxn"   & node_priority == "main" ~ node_action_color,
        node_type == "gene"  & node_priority == "main" ~ node_gene_color,

        node_type == "recon" & node_priority == "main_secondary" ~ node_metabolite_recon_color,
        node_type == "gene"  & node_priority == "main_secondary" ~ node_gene_color,

        node_type == "recon" & node_priority == "minor" ~ node_minor_metabolite_color,
        node_type == "ggm"   & node_priority == "minor" ~ node_minor_metabolite_color,
        node_type == "rxn"   & node_priority == "minor" ~ node_minor_action_color,
        node_type == "gene"  & node_priority == "minor" ~ node_minor_gene_color,


        TRUE ~ "black"
      ),
      borderWidth = if_else(crowded | reversible, 2.5, 0),
      color.border = case_when(
        crowded ~ node_crowded_border_color,
        node_priority == "main"  & reversible ~ node_reversible_border_color,
        node_priority == "minor" & reversible ~ node_minor_reversible_border_color,
        TRUE ~ color.background
      ),
      color.highlight.border  = color.border,
      color.highlight.background  = color.background,
      size = case_when(
        id %in% c(source_nodes, target_nodes) ~ node_endpoint_size,
        node_priority == "main" ~ node_main_size,

        node_type == "recon" & node_priority == "main_secondary" ~ node_main_secondary_size,
        node_type == "gene"  & node_priority == "main_secondary" ~ node_main_secondary_size,

        node_type == "rxn" ~ node_rxn_size,
        TRUE ~ node_minor_size
      ),
      font.size = case_when(
        id %in% c(source_nodes, target_nodes) ~ node_endpoint_font_size,
        node_priority == "main" ~ node_main_font_size,
        node_priority == "main_secondary" ~ node_main_secondary_font_size,
        TRUE ~ node_minor_font_size),
      font.color = case_when(
        node_priority == "main" ~ node_main_font_color,
        node_priority == "main_secondary" ~ node_main_secondary_font_color,
        TRUE ~ node_minor_font_color),
      shape = case_when(
        node_type == "rxn" ~ node_rxn_shape,
        id %in% c(source_nodes, target_nodes) ~ node_endpoint_shape,
        TRUE ~ node_all_others_shape
      )
    )




  # Display clusters setting -------------------------------------------------

  if (show_clusters) {

    selected_by <- list(variable = "cluster_nums",
                        main = "Select path by cluster",
                        multiple = TRUE,
                        style = selectize_style,
                        highlight = TRUE)

  } else {

    selected_by <- NULL


  }



  # Visualize network -------------------------------------------------------

  recon_vis <-
    visNetwork(nodes = vis_nodes,
               edges = vis_edges,
               height = "800px",
               width = "1200px") %>%
    visIgraphLayout(layout = "layout_with_fr", randomSeed = 123) %>%
    visEdges(arrows = "to",
             color = "gray") %>%
    visOptions(
      nodesIdSelection = list(enabled = TRUE,
                              main = "Select node by name",
                              style = selectize_style
      ),

      selectedBy = selected_by,

      highlightNearest = list(enabled = TRUE,
                              labelOnly = FALSE,
                              algorithm = "hierarchical",
                              degree = 0)) %>%
    visLayout(randomSeed = preset_seed)

  list(visualization = recon_vis,
       trace_name = sprintf("%s ---> %s (%s)", source_nodes, target_nodes, shortest_path),
       main_nodes_dict = main_nodes_dict)

}


tr_vis_is_gene_to_gene <- function(path_row) {

  path_row %>%
    mutate(gene_to_gene = if_else(!str_detect(gene, "\\[|__ggm") &
                                    !str_detect(metabolite, "\\[|__ggm"),
                                  TRUE, FALSE)) %>%
    .$gene_to_gene %>%
    unique()

}




tr_visualize_trace <- function(path_row, show_clusters = FALSE) {

  if (tr_vis_is_gene_to_gene(path_row)) {
    res <- visualize_trace_gene_to_gene(path_row, show_clusters)

  } else {
    res <- visualize_trace_non_gene_to_gene(path_row, show_clusters)

  }

  res

}