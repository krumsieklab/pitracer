library(tidyverse)
library(readxl)
library(data.table)

# Load objects --------------------------------------------------------------

work_dir <- getwd()

load(glue::glue(work_dir, "/precalculated_data/parsed_recon_objects.RData"))
#rxn_scores <- fread(glue::glue(work_dir, "/precalculated_data/scores.csv")) # debug: commented out, not needed
connected_mets <- fread(glue::glue(work_dir, "/precalculated_data/recon_reactions_rescored.csv"))

recon_lookup <- read_excel(glue::glue(work_dir, "/precalculated_data/recon_id_name_lookup.xlsx"))
generic_step_lookup <- readRDS(glue::glue(work_dir, "/precalculated_data/recon_to_generic_reactions_lookup.rds"))

recon_rxns <- read_rds(glue::glue(work_dir, "/precalculated_data/recon_reactions.rds"))


# Assign recon helper dataframes ----------------------------------------

rxn_readable_names <-
  recon_rxns %>%
  rowwise() %>%
  dplyr::transmute(rxn,
                   rxn_name = if_else(use_short_rxn_names,
                                      rxn_name_short,
                                      rxn_name)) %>%
  ungroup()

recon_rxn_genes <- filter(recon_rxns, length(genes) > 0)


# Create edges dataframe with oriented recon reactions --------------------

tmp_recon_forward <-
  recon_rxns %>%
  rowwise() %>%
  dplyr::mutate(reactants = paste(unlist(reactants), collapse=", "),
                products = paste(unlist(products), collapse=", ")) %>%
  expand_reactions() %>%
  left_join(connected_mets %>% dplyr::mutate(weight = 1),
            by = c("rxn", "reactant", "product")) %>%
  dplyr::transmute(rxn, reactant, product,
                   connected = if_else(weight != 1 | is.na(weight), FALSE, TRUE),
                   reversible)

tmp_recon_backward <-
  tmp_recon_forward %>%
  filter(reversible) %>%
  dplyr::mutate(reactant_tmp = reactant,
                reactant = product,
                product = reactant_tmp) %>%
  dplyr::select(-reactant_tmp)

tmp_ggm <-
  translate_generic_moa(met_network) %>%
  filter(str_detect(rxn, "_ggm_")) %>%
  # Since dataframe for ggm is symmetric, i.e. after n/2 where n is
  # the number of ggm edges, starting n/2 + 1 we have the reverse
  # reactions. This is not really needed, since it adds redundancy
  # in the visualization, where a reversible node can be used instead.
  dplyr::slice(1:(n()/2)) %>%
  dplyr::transmute(rxn, reactant = from, product = to,
                   # since there are no cofactors annotated to ggm
                   # edges, just set to FALSE by default
                   connected = FALSE,
                   reversible = TRUE)

rxn_steps <-
  bind_rows(tmp_recon_forward,
            tmp_recon_backward,
            tmp_ggm) %>%
  distinct()