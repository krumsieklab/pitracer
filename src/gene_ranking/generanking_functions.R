# Contains all gene ranking functionality

library(tidyverse)
library(data.table)
library(igraph)
library(visNetwork)
library(svMisc)

# add_upstream should be set to TRUE, if you want to add upstream genes
# remove_transport should be set to TRUE, if you want to remove transport reactions
# Warning, if you remove transport reactions, you might not be able to compare
# it to the complete distance matrix, since the compartment information is removed
get_distances <- function(mets,
                          trace_objects,
                          k,
                          add_upstream,
                          remove_transport = FALSE,
                          use_dist = in_distances,
                          directed = TRUE, # directed, meaning if distance to be calculated from trace is directed or not
                          shinyProgress = NULL # for Shiny interaction only
) {
  
  if (missing(add_upstream)) stop("Please set add_upstream as either TRUE or FALSE")
  
  readable_mets <-
    recon_lookup %>%
    filter(id %in% mets) %>%
    .$name %>%
    unique()
  
  
  # Get trace edges
  trace_edges <-
    trace_objects[mets] %>%
    bind_rows() %>%
    filter(has_trace) %>%
    unnest(association_path) %>%
    dplyr::rename(start = from, end = to) %>%
    group_by(start, end) %>%
    
    
    # adding the ff so I can get the right edges per path,
    # will be irrelevant later, since there is one gene
    # going to all those metabolites
    mutate(path_number = row_number()) %>%
    unnest(from = association_path) %>%
    mutate(to = lead(from)) %>%
    
    # remove last step, which basically circles on in the path
    ungroup() %>%
    group_by(start, end, path_number) %>%
    dplyr::mutate(step_ = row_number()) %>%
    filter(step_ != max(step_)) %>%
    ungroup() %>%
    
    filter(path_number <= k) %>%
    
    dplyr::select(-c(has_trace, path_number, step_)) %>%
    
    # distinct since it doesn't matter how many times it
    # occurs within for a gene
    distinct() %>%
    ungroup() %>%
    filter(start %in% readable_mets & end %in% readable_mets)
  
  # Get a df of all start and end nodes
  terminal_nodes <-
    trace_edges %>%
    distinct(start, end)
  
  
  # Score edges based on both start and end ---------------------------------
  
  
  edge_scores <-
    trace_edges %>%
    distinct(from, to)
  
  
  
  # Add auxiliary genes ----------------------------------------------
  
  # genes connected to interaction nodes
  addon_edges <-
    left_join(edge_scores, met_network, by = c("from", "to")) %>%
    distinct() %>%
    left_join(generic_lookup, by = c("moa" = "gen_rxn")) %>%
    left_join(recon_genes, by = "rxn") %>%
    dplyr::select(-c(moa, weight, rxn)) %>%
    filter(!is.na(genes)) %>%
    ungroup() %>%
    
    # do this step, since generic_lookup was made in a way to have
    # interaction genes. I created recon_genes here to make it more
    # intuitive if I read this in the future
    filter(from != genes,
           to != genes,
           !genes %in% terminal_nodes$start,
           !genes %in% terminal_nodes$start) %>%
    distinct()
  
  addon_edges_from <-
    addon_edges %>%
    dplyr::transmute(to = from,
                     from = genes) %>%
    dplyr::select(from, to)
  
  addon_edges_to <-
    addon_edges %>%
    dplyr::transmute(from = genes,
                     to)
  
  addon_edges_all <-
    bind_rows(addon_edges_from, addon_edges_to) %>%
    distinct()
  
  
  addon_nodes <-
    addon_edges %>%
    dplyr::select(gene = genes) %>%
    distinct() %>%
    group_by(gene) %>%
    ungroup() %>%
    dplyr::mutate(auxiliary = TRUE)
  
  
  # Get distance matrix without upstream genes
  
  if(directed) {
    
    distance_raw <-
      bind_rows(addon_edges_all, edge_scores) %>%
      distinct() %>%
      graph_from_data_frame() %>%
      distances(mode = "out")
    
  } else {
    
    distance_raw <-
      bind_rows(addon_edges_all, edge_scores) %>%
      distinct() %>%
      graph_from_data_frame() %>%
      distances(mode = "all")
    
  }
  
  
  
  # Keep only genes-metabolite distances
  dist_genes <-
    rownames(distance_raw) %>%
    .[!str_detect(., "\\[")] %>%
    unique()
  
  dist_mets <-
    rownames(distance_raw) %>%
    .[str_detect(., "\\[")] %>%
    unique()
  
  distance_raw <- distance_raw[dist_genes, dist_mets]
  
  if (add_upstream) {
    # get edge list with distances
    trace_distances <-
      distance_raw %>%
      as.data.frame() %>%
      rownames_to_column("from") %>%
      gather(to, dist, -from) %>%
      filter(dist != Inf) %>%
      distinct() %>%
      dplyr::select(from, to, weight = dist) %>%
      # only keep gene-metabolite distances
      dplyr::filter(!str_detect(from, "\\["))
    
    # construct network of traces with genes
    # take those genes, and see what other genes are reachable
    # Check reachability of upstream genes to trace genes --------------------
    
    reachable_genes <-
      addon_nodes$gene %>%
      unique() %>%
      use_dist[., ] %>%
      as.data.frame() %>%
      rownames_to_column("to") %>%
      gather(from, dist, -to) %>%
      filter(dist != Inf) %>%
      filter(!from %in% trace_distances$from) %>%
      distinct() %>%
      dplyr::select(from, to, weight = dist) %>%
      # ensure that only genes are here
      dplyr::filter(!str_detect(from, "\\[") & !str_detect(to, "\\["))
    
    counter <<- 0
    total_count <<- reachable_genes$from %>% unique() %>% length()
    print ("Screening for upstream genes...")
    
    upstream_distances <-
      reachable_genes$from %>% unique() %>%
      lapply(function(gene) {
        counter <<- counter + 1
        progress(counter, max.value = total_count)
        
        if (!is.null(shinyProgress)) {
          shinyProgress$set(
            value = 15 + floor(counter/total_count*70),
            message = "Screening for upstream genes...",
            detail = sprintf("%s/%s", counter, total_count)
          )
        }
        
        reachable_genes %>%
          filter(from == gene) %>%
          inner_join(trace_distances, by = c("to" = "from")) %>%
          dplyr::transmute(from,
                           to = to.y,
                           dist = weight.x + weight.y) %>%
          group_by(to) %>%
          filter(dist == min(dist))
        
      }) %>%
      bind_rows()
    
    distance_raw <-
      upstream_distances %>%
      bind_rows(trace_distances %>% dplyr::rename(dist = weight)) %>%
      ungroup() %>%
      distinct()
    
    if (remove_transport) {
      
      distance_raw <-
        distance_raw %>%
        mutate(to = str_remove(to, "\\[.\\]")) %>%
        group_by(from, to) %>%
        arrange(dist) %>%
        
        # select the shortest distance for each metabolite
        slice(1) %>%
        ungroup() %>%
        
        spread(to, dist, fill = Inf) %>%
        column_to_rownames("from") %>%
        as.matrix()
      
    } else {
      
      distance_raw <-
        distance_raw %>%
        spread(to, dist, fill = Inf) %>%
        column_to_rownames("from") %>%
        as.matrix()
    }
    
  }
  
  return(distance_raw)
  
}


get_gene_weights <- function(input_mat,
                             input_trace_genes,
                             input_indistance,
                             sort_order) {
  
  gene_weights <-
    input_mat %>%
    rownames() %>%
    unique() %>%
    input_indistance[input_trace_genes, .] %>%
    as.data.frame() %>%
    rownames_to_column("to") %>%
    gather(from, distance, -to) %>%
    filter(!distance %in% c(0, Inf)) %>%
    group_by(from) %>%
    dplyr::summarise(gene_weight = 1 - (n()/length(input_trace_genes)))
  
  
  gene_weights_diag <-
    sort_order %>%
    data.frame(from = .) %>%
    left_join(gene_weights, by = "from") %>%
    dplyr::transmute(from,
                     gene_weight = if_else(is.na(gene_weight), 1, gene_weight)) %>%
    column_to_rownames("from") %>%
    as.matrix() %>%
    .[sort_order, ] %>%
    diag()
  
  rownames(gene_weights_diag) <- sort_order
  colnames(gene_weights_diag) <- sort_order
  
  gene_weights_diag
  
}


get_met_pw <- function(result_obj, experiment_metabolites) {
  
  met_scores <-
    result_obj$S %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    gather(met, score, -gene) %>%
    left_join(experiment_metabolites, by = c("met" = "recon_id")) %>%
    dplyr::transmute(gene,
                     met,
                     met_name = name,
                     score,
                     weight = if_else(as.logical(is_significant), 1, 0),
                     weighted_score = score * weight) %>%
    inner_join(pathways, by = "met") %>%
    dplyr::select(-weight) %>%
    arrange(desc(weighted_score))
  
  
  # just adding the weighted metabolite scores of all metabolites in a pathway
  pw_sum <-
    met_scores %>%
    dplyr::group_by(gene, subsystem) %>%
    dplyr::summarise(pw_score = sum(weighted_score)) %>%
    ungroup() %>%
    group_by(gene) %>%
    # filter(pw_score > 0) %>%
    arrange(pw_score) %>%
    left_join(result_obj$scores, by = "gene") %>%
    ungroup() %>%
    dplyr::transmute(gene,
                     gene_rank,
                     pathway = subsystem,
                     pathway_score = pw_score) %>%
    arrange(desc(pathway_score))
  
  
  list(metabolite_scores = met_scores,
       pathway_scores    = pw_sum)
  
}


get_gene_scores <- function(D, # this is the pre-calculated distance matrix from the traces, to be distinguished from the global piTracer distance matrix
                            metabolites,
                            w, # target metabolite set, 1 = target, 0 = no target
                            tracer_dist, # this is the pre-calculated global piTracer distance matrix, not used when no upstream
                            with_upstream = FALSE) {
  
  
  #### scoring ----
  
  # transform weights?
  wtrans <- function(x)x # leave them as 1/0
  
  # set distance matrix
  D <- D %>% as.matrix()
  
  # # transform distances to scores
  S <- 1/D
  
  # subset everything to the common metabolites
  metcommon <- intersect(metabolites, colnames(S))
  w <- w[match(metcommon, metabolites)]
  S <- S[, metcommon]
  
  enzyme_names <- rownames(S)
  
  if (with_upstream) {
    
    # get all trace enzymes
    trace_enzymes <-
      dir_without_upstream %>%
      rownames()
    
    diag_g <- get_gene_weights(D, trace_enzymes, tracer_dist, enzyme_names)
    
  } else {
    
    diag_g <- rep(1, nrow(D)) %>% diag()
    rownames(diag_g) <- enzyme_names
    colnames(diag_g) <- enzyme_names
    
  }
  
  # score all genes
  G <- diag_g %*% S %*% wtrans(w)
  
  # make dplyr friendly
  scores <- G %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    dplyr::rename(score = V1, gene = rowname) %>%
    dplyr::mutate(gene_rank = dense_rank(dplyr::desc(score))) %>%
    dplyr::arrange(gene_rank)
  
  # return dplyr-friendly version as well as individual components
  list(
    scores = scores,
    G = G,
    diag_g = diag_g,
    S = S,
    w = w,
    wtrans = wtrans(w)
  )
  
}

# main ranking function that calculates the gene ranking
calculate_generanks <- function(user_input, comp_priority, top_n_results = 30, shinyProgress = NULL) {
  
  # Process input into a format the algorithm understands ---------------------
  if (!is.null(shinyProgress)) {
    shinyProgress$set(value = 0, message = "Processing input...")
  }
  
  names(user_input) <- c("metabolite_name", "is_significant")
  
  # assemble data frame of metabolites, fix HMDB and map to recon
  experiment_metabolites <- data.frame(
    name = user_input$metabolite_name,
    recon_id = user_input$metabolite_name %>% hmdb5_fixer() %>% to_recon_mapper(lst = ., comp_priority = comp_priority),
    is_significant = user_input$is_significant %>%
      toupper() %>%
      as.logical()
  ) %>% tibble::remove_rownames()
  
  
  # Create comparison table so the user understands what became of the input ----
  user_input_mapping <- data.table(
    input_name = user_input$metabolite_name,
    mapped_metabolites = experiment_metabolites$recon_id,
    is_significant = experiment_metabolites$is_significant
  )

  # throw error if there are unmapped metabolites
  if (any(is.na(experiment_metabolites$recon_id))) {
    stop(sprintf("Could not find the following input identifiers in the database: %s",
                 paste0(experiment_metabolites$name[is.na(experiment_metabolites$recon_id)], collapse = ", ")))
  }
  
  # this should not occur anymore, but leaving it in as a safety measure
  if (nrow(experiment_metabolites) == 0) {
    stop("No valid metabolites entered")
  }
  
  # Get traces ---------------------------------------------------------------
  if (!is.null(shinyProgress)) {
    shinyProgress$set(value = 5, message = "Fetching traces...")
  }
  
  # calculate traces
  trace_objects <- fast_trace_loop(
    start_nodes = experiment_metabolites$recon_id,
    end_nodes = experiment_metabolites$recon_id,
    network = trace_net,
    k = 10
  )
  
  # no distances found?
  if (trace_objects %>% sapply(is.null) %>% all()) {
    stop("No distances found")
  }
  
  # Remove metabolites that have no traces
  experiment_metabolites <- experiment_metabolites %>%
    filter(recon_id %in% names(trace_objects))
  
  # Calculate outward distances
  if (!is.null(shinyProgress)) {
    shinyProgress$set(value = 10, message = "Calculating outward distances...")
  }
  
  # directed distances without upstream effects
  dir_without_upstream <<- get_distances(
    experiment_metabolites$recon_id,
    trace_objects,
    k = 10,
    add_upstream = FALSE,
    remove_transport = FALSE
  )
  
  # Calculate upstream distances
  if (!is.null(shinyProgress)) {
    shinyProgress$set(value = 15, message = "Calculating upstream distances...")
  }
  
  # undirected distances, with upstream effects
  met_undir_dir_noppi_upstream <- get_distances(
    experiment_metabolites$recon_id,
    trace_objects,
    k = 10,
    add_upstream = TRUE,
    remove_transport = FALSE,
    directed = FALSE,
    use_dist = no_ppi_met_undir_ins,
    shinyProgress = shinyProgress
  )
  
  # Calculate gene scores
  if (!is.null(shinyProgress)) {
    shinyProgress$set(value = 85, message = "Calculating gene scores...", detail = NULL)
  }
  
  # generate main result table, with one single score per gene
  Gres_met_undir_noppi_upstream <- get_gene_scores(
    met_undir_dir_noppi_upstream,
    experiment_metabolites$recon_id,
    experiment_metabolites$is_significant %>% as.logical() %>% as.numeric(),
    in_distances_noppi_met_undir,
    with_upstream = TRUE)
  
  # filter out top n results and sort
  main_results <<-
    Gres_met_undir_noppi_upstream$scores %>%
    filter(gene_rank <= top_n_results) %>%
    arrange(gene_rank)
  
  # gather pathway-to-recon mapping
  pathways <<-
    data.reactions %>%
    dplyr::select(reactants, products, subsystem) %>%
    gather(type, met, -subsystem) %>%
    unnest(met) %>%
    dplyr::mutate(met = str_remove(met, "M_"),
                  met = str_replace(met, "__91__", "\\["),
                  met = str_replace(met, "__93__", "\\]")) %>%
    filter(met %in% colnames(dir_without_upstream)) %>%
    filter(!str_detect(subsystem, "Transport")) %>%
    filter(subsystem != "Miscellaneous") %>%
    distinct(subsystem, met)
  
  

  # Create metabolite contributions
  if (!is.null(shinyProgress)) {
    shinyProgress$set(value = 95, message = "Calculating metabolite contributions...")
  }
  
  met_contributions <<- get_met_pw(Gres_met_undir_noppi_upstream, experiment_metabolites) %>%
    .$metabolite_scores %>%
    filter(gene %in% main_results$gene) %>%
    group_by(gene, met, met_name, score, weighted_score) %>%
    summarise(pathways = str_c(subsystem, collapse = ", ")) %>%
    arrange(match(gene, main_results$gene), desc(weighted_score)) %>%
    as.data.frame()
  
  if (!is.null(shinyProgress)) {
    shinyProgress$set(value = 100, message = "Done")
  }
  
  # generate pathway scores, first ranked by genes as in the main results, then by pathway score
  pathway_scores <- 
    met_contributions %>% 
    # group by gene and pathway pair
    group_by(gene, pathways) %>% 
    # sum up all scores
    summarize(pathway_score = sum(weighted_score)) %>% 
    # sort by the original gene order, the by descending pathway score
    arrange(factor(gene, levels=main_results$gene), desc(pathway_score)) 
  
  
  # assemble final results
  list(
    # list with genes and scores
    main_results = main_results,
    # individual contributions of metabolites
    met_contributions = met_contributions,
    # pathway scores
    pathway_scores = pathway_scores,
    # user input to recon matching
    user_input_mapping = user_input_mapping
  )
  
}


