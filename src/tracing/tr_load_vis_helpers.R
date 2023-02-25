library(tidyverse)
library(readxl)
library(data.table)

# Helper functions ---------------------------------------------------------------

get_readable_names <- function(edge_df, dict, col_name) {

  left_index1 <- enquo(col_name)
  right_index1 <- "id"
  by1 <-  set_names(quo_name(right_index1), quo_name(left_index1))

  trans_df <-
    edge_df %>%
    left_join(dict, by = by1)

  edge_df[[col_name]] = if_else(!is.na(trans_df[["name"]]),
                                trans_df[["name"]],
                                edge_df[[col_name]] )

  edge_df
}

unnest_path <- function(path_row) {

  path_row %>%
    ungroup() %>%

    # number each k-path
    group_by(gene, metabolite, met_ids) %>%
    dplyr::mutate(path_number = row_number()) %>%

    # unroll paths to make it into a from-to format
    unnest(association_path) %>%
    group_by(gene, metabolite, met_ids, path_number) %>%
    mutate(to = lead(association_path, 1)) %>%
    dplyr::slice(1:(n()-1)) %>%
    ungroup() %>%

    transmute(from = association_path,
              to,
              source_node = gene,
              target = met_ids,
              cluster_num)

}

add_moa <- function(unnested_paths) {

  # do this for speedup
  ppitf_network_mini <-
    ppitf_network %>%
    filter(from %in% main_nodes, to %in% main_nodes)

  # separate this function to met vs ppitf networks
  unnested_paths %>%
    left_join(met_network, by = c("from", "to")) %>%
    left_join(ppitf_network_mini, by = c("from", "to")) %>%
    mutate(moa = if_else(!is.na(moa.x), moa.x, moa.y)) %>%
    dplyr::select(-c(weight.x, weight.y, moa.x, moa.y))

}

translate_generic_moa <- function(paths_with_moa) {

  paths_with_moa %>%
    left_join(generic_step_lookup, by = c("moa" = "gen_rxn")) %>%
    dplyr::mutate(rxn = if_else(is.na(rxn), moa, rxn)) %>%
    dplyr::select(-moa) %>%
    distinct()

}

fill_recon <- function(translated_paths) {

  # This adds recon reactants, products, and reversibility info
  left_join(translated_paths, recon_rxns, by = "rxn")  %>%
    rowwise() %>%
    dplyr::mutate(reactants = paste(unlist(reactants), collapse=", "),
                  products = paste(unlist(products), collapse=", "),
                  genes = paste(unlist(genes), collapse=", "))

}


count_rxn_occ <- function(df) {
  # take custon column
  # number repeated reactions
  df %>%
    dplyr::group_by(rxn) %>%
    dplyr::mutate(rxn_rep = 1:n()) %>%
    ungroup()

}

expand_reactions <- function(df_with_rxns) {

  df_with_rxns %>%
    unnest(reactant = str_split(reactants, ", ")) %>%
    unnest(product = str_split(products, ", "))

}

regexify <- function(recon_string) {
  str_replace(recon_string, "\\[", "\\\\[") %>%
    str_replace("\\]", "\\\\]")
}

orient_rxns <- function(df_filled) {

  df_filled %>%
    rowwise() %>%
    dplyr::mutate(to_in_products = if_else(str_detect(products, regexify(to)),
                                           TRUE, FALSE),
                  reactants_tmp = reactants,
                  reactants = if_else(reversible & !to_in_products & !str_detect(rxn, "_ggm_"),
                                      products, reactants),
                  products = if_else(reversible & !to_in_products & !str_detect(rxn, "_ggm_"),
                                     reactants_tmp, products))

}


get_secondary_mets <- function(df) {

  df %>%
    # remove ggm edges
    filter(!str_detect(rxn, "_ggm_")) %>%

    # remove the main edges
    filter((from != reactant | to != product)) %>%
    filter((from != product | to != reactant)) %>%

    # annotate which metabolites are connected
    left_join(rxn_steps,
              by = c("rxn", "reactant", "product")) %>%
    distinct() %>%

    # label the edges
    dplyr::mutate(edge_type = case_when(
      (from == reactant | to == product |
         from == product | to == reactant) & connected == TRUE ~ "main_secondary",
      TRUE ~ "minor"
    )) %>%

    dplyr::select(rxn, rxn_rep, cluster_nums,
                  from = reactant, to = product,
                  edge_type)

}

remove_main_edges <- function(df, main_from, main_to) {

  df %>%
    filter((!reactant %in% main_from) | (!product %in% main_to)) %>%
    filter((!product %in% main_from) | (!reactant %in% main_to))

}


translate_rxn_name <- function(df) {

  df %>%
    left_join(rxn_readable_names, by = "rxn") %>%
    dplyr::mutate(rxn = if_else(!is.na(rxn_name), rxn_name, rxn)) %>%
    dplyr::select(-rxn_name)


}


add_reversibility <- function(df) {

  df %>%
    left_join(recon_rxns %>% dplyr::select(rxn, reversible),
              by = "rxn") %>%
    dplyr::mutate(reversible = case_when(
      str_detect(rxn, "_ggm_") ~ TRUE,
      TRUE ~ reversible))


}

gather_edges <- function(df) {

  df %>%
    gather(directionality,
           molname, -c(rxn, rxn_rep, cluster_nums, edge_type, reversible)) %>%
    distinct()

}

translate_metabolite_name <- function(gathered_df, dict, col_name) {

  left_index1 <- enquo(col_name)
  right_index1 <- "id"
  # have to do this to be able to pass custom colnames in
  # function
  by1 <-  set_names(quo_name(right_index1), quo_name(left_index1))

  translated_df <- left_join(gathered_df, dict, by = by1)


  gathered_df[[col_name]] = if_else(!is.na(translated_df[["name"]]),
                                    translated_df[["name"]],
                                    gathered_df[[col_name]] )

  gathered_df
}


translate_metabolite_name_list <- function(met_list, dict) {

  met_df <- tibble(molname = met_list)
  translate_metabolite_name(met_df, dict, "molname") %>%
    .$molname

}


add_metabolite_numbering <- function(gathered_df, exclude_nodes) {

  gathered_df %>%
    group_by(molname) %>%
    dplyr::mutate(molname_rep = 1:n()) %>%
    ungroup() %>%
    dplyr::mutate(molname = ifelse(molname_rep!=1 & !molname %in% exclude_nodes,
                                   str_c(molname, " (#", molname_rep, ")"),
                                   molname)) %>%
    dplyr::select(-molname_rep)

}

add_ggm_edge_anno <- function(gathered_df) {

  gathered_df %>%
    dplyr::mutate(ggm_edge = if_else(str_detect(rxn, "_ggm_"),
                                     TRUE, FALSE))

}

fix_interface_direction <- function(gathered_df, interface_df) {

  gathered_df %>%
    left_join(interface_df, by = c("rxn", "molname")) %>%
    dplyr::mutate(directionality = if_else(!is.na(direction_sub),
                                           direction_sub,
                                           directionality)) %>%
    dplyr::select(-direction_sub)


}