# This script sources all internal tracing codes and all Shiny module files,
# and then loads various pre-computed files needed for the app

# source code files ------------------------------------------------------------

work_dir <- getwd()

formals(source)$local <- TRUE
source(glue::glue(work_dir, "/src/tracing/tracing.R"), local = TRUE)
source(glue::glue(work_dir, "/src/gene_ranking/generanking.R"), local = TRUE)

source(glue::glue(work_dir, "/src/shiny_app/module_blacklist.R"))
source(glue::glue(work_dir, "/src/shiny_app/module_gene_ranking.R"))
source(glue::glue(work_dir, "/src/shiny_app/module_node_input.R"))
source(glue::glue(work_dir, "/src/shiny_app/module_results.R"))
source(glue::glue(work_dir, "/src/shiny_app/module_dl_generank.R"))
source(glue::glue(work_dir, "/src/shiny_app/busy_indicator.R"), local = TRUE)


# random seed settings
preset_seed <- 1
set.seed(preset_seed)




# create lookup tables -----------------------------------------------------

# Metabolites
metabolites <- met_network %>%
  select(from, to) %>%
  unlist() %>%
  unique()

# define recon metabolites
recon_metabolites <- recon_lookup$id
names(recon_metabolites) <- recon_lookup$name
recon_metabolites <- recon_metabolites[which(recon_metabolites %in% metabolites)]


# Genes
genes <<- ppitf_network %>%
  select(from, to) %>%
  unlist() %>%
  unique() %>%
  sort()
names(genes) <- genes


# Create metabolite database for the search bar
# Computation is rather lengthy, so create it once, save it and load on init

priority <- c("c", "m", "e", "r", "n", "g", "x", "l", "i")

# load IDs and metabolite lists
# id_map <- read_rds(glue::glue(work_dir, "/precalculated_data/id_map.rds")) # no need to preload, can be generated very fast... see further down in this file
full_metabolite_database <- read_rds(glue::glue(work_dir, "/precalculated_data/full_metabolite_database.rds"))
extracellular_mets <- read_rds(glue::glue(work_dir, "/precalculated_data/extracellular_mets.rds"))
cytoplasm_mets <- read_rds(glue::glue(work_dir, "/precalculated_data/cytoplasm_mets.rds"))



# load Omnipath files ------------------------------------------------------

omnipathdb_edges <-
  data.table::fread(glue::glue(work_dir, "/precalculated_data/omnipath_webservice_interactions__recent.tsv")) %>%
  # ncbi taxonomy id for humans: 9606
  filter(ncbi_tax_id_source == 9606 & ncbi_tax_id_target == 9606) %>%
  filter(consensus_direction == 1) %>%
  select("source_genesymbol", "target_genesymbol", "sources", "references") # %>%
# rowwise() %>%
# mutate(sources = str_split(sources, ";") %>% unlist() %>% list()) %>%
# unnest(sources) %>%
# rowwise() %>%
# mutate(references = str_split(references, ";") %>% unlist() %>% list()) %>%
# unnest(references) %>%
# filter(references != "") %>%
# separate(references, c("source", "source_id"), ":")
# .[, c("source_genesymbol", "target_genesymbol", "source", "source_id")]

omnipathdb_source_list <- readr::read_rds(glue::glue(work_dir, "/precalculated_data/omnipathdb_source_list.rds"))


# Function to write Excel file/sheet with openxlsx package, which works on all platforms
# openxlsx::writexlsx does not support adding sheets to existing files; this is implemented in this wrapper
#
# Last update: JK, 17.09.2017
#
writexlsx_append <- function(x, file, sheetName, append=T, rowNames=T, colNames=T) {
  
  # if no append => make sure file does not exist
  if (!append) unlink(file)
  # if file still exists => load it, otherwise, create one
  if (file.exists(file)) {
    wb <- openxlsx::loadWorkbook(file)
  } else {
    wb <- openxlsx::createWorkbook()
  }
  
  # write and save 
  ws=openxlsx::addWorksheet(wb,sheetName=sheetName)
  openxlsx::writeData(wb=wb, sheet=sheetName, x=x, colNames = colNames, rowNames = rowNames)
  openxlsx::saveWorkbook(wb, file, overwrite=T)
  
}


###########################################################################
# CODE ARCHIVE ----
###########################################################################

# translate_names <- TRUE

# toast_options <- list(
#   progressBar = FALSE,
#   positionClass = "toast-top-center",
#
#   # same as defaults
#   newestOnTop = TRUE,
#   preventDuplicates = FALSE,
#   showDuration = 500,
#   hideDuration = 1000,
#   extendedTimeOut = 1000,
#   showEasing = "linear",
#   hideEasing = "linear",
#   showMethod = "fadeIn",
#   hideMethod = "fadeOut"
# )


# search_bar_choices <- read_rds(glue::glue(work_dir, "/precalculated_data/search_bar_choices.rds"))


# generate ID map
id_map_withcomp <<- read_feather("precalculated_data/app_lookup.feather") %>%
  left_join(recon_lookup, by = c("recon_id" = "id")) %>%
  #filter(recon_id %in% metabolites) %>%
  # this part smartly figures out what is HMDB, PubChem, or KEGG ID
  dplyr::mutate(comp_id = str_match(recon_id, "\\[([a-z])\\]")[,2],
                id_type = case_when(
                  #is.idna(id) ~ "Recon ID",
                  str_detect(id, "HMDB") ~ "HMDB",
                  !str_detect(id, "[A-z]") ~ "PubChem",
                  TRUE ~ "KEGG"
                ))  %>%
  # add compartment to id column
  mutate(id=paste0(id, "[", comp_id, "]")) %>%
  # also add recon ID without compartment to table
  mutate(recon_id_nocomp = str_replace_all(recon_id, "\\[([a-z])\\]",""))

# infer ID map with no compartments
id_map_nocomp <<- 
  id_map_withcomp %>% 
  # delete compartment info from recon and id
  mutate(recon_id=str_replace_all(recon_id, "\\[([a-z])\\]","")) %>%
  mutate(id=str_replace_all(id, "\\[([a-z])\\]","")) %>%
  mutate(name=trimws(str_replace_all(name, "\\[([a-z])\\]",""))) %>%
  # drop unneeded columns
  dplyr::select(-comp_id) %>%
  # delete duplicate rows
  dplyr::distinct()



# the following steps have been removed ... they keep only one compartment version of each metabolites
# %>%
#   group_by(id) %>%
#   dplyr::arrange(match(comp_id, priority)) %>%
#   dplyr::arrange(id) %>%
#   dplyr::slice(1)
#
# id_map %>%
#   write_rds(glue::glue(work_dir, "/precalculated_data/id_map.rds"))
#





#
# missing_metabolites_df <- data.frame(
#   recon_id = recon_metabolites %>%
#     setdiff(id_map$recon_id)
# ) %>%
#   mutate(
#     id = NA,
#     name = names(recon_metabolites)[match(recon_id, recon_metabolites)],
#     comp_id = substring(name, nchar(name) - 1, nchar(name) - 1),
#     id_type = "recon"
#   )
#
# full_metabolite_database <<- id_map %>%
#   full_join(missing_metabolites_df) %>%
#   dplyr::mutate(name = name %>%
#                   # str_split(" \\[[a-z]\\]+") %>%
#                   # sapply("[[", 1, USE.NAMES = FALSE),
#                   str_remove(" \\[[a-z]\\]"),
#                 recon_id = recon_id %>%
#                   # str_split("\\[[a-z]\\]+") %>%
#                   # sapply("[[", 1, USE.NAMES = FALSE)
#                   str_remove(" \\[[a-z]\\]")
#   )
#
# full_metabolite_database %>%
#   write_rds(path = glue::glue(work_dir, "/precalculated_data/full_metabolite_database.rds"))


# Creates a named character vector containing all recon IDs of metabolites, their
# KEGG/HMDB/PubChem equivalent, and genes. It takes a while to process, so in
# order to reduce loading times, it's saved in the data directory and loaded
# for future use

# mets_with_other_ids <- id_map %>%
#   filter(!str_detect(id_type, "recon")) %>%
#   .$name %>%
#   sapply(function(met_name) {
#     suggestion_value <- id_map %>%
#       filter(name == met_name) %>%
#       .$recon_id %>%
#       .[1]
#
#     names(suggestion_value) <- id_map %>%
#       filter(name == met_name) %>%
#       arrange(match(id_type, c("KEGG", "HMDB", "PubChem"))) %>%
#       .$id %>%
#       paste(collapse = ", ") %>%
#       sprintf("%s (%s)",
#               str_remove(met_name, "\\s\\[[a-z]\\]"),
#               .)
#
#     suggestion_value
#   },
#   USE.NAMES = FALSE) %>%
#   # c(recon_metabolites[which(!recon_metabolites %in% id_map$recon_id)],
#   #   genes) %>%
#   .[!duplicated(.)]
#
# missing_metabolites <- full_metabolite_database %>%
#   filter(str_detect(id_type, "recon")) %>%
#   .$recon_id
#
# names(missing_metabolites) <- full_metabolite_database %>%
#   filter(str_detect(id_type, "recon")) %>%
#   .$name %>%
#   str_remove("\\s\\[[a-z]\\]")
#
# search_bar_choices <- c(
#   mets_with_other_ids,
#   missing_metabolites,
#   genes
# ) %>%
#   .[order(.)]
#
# search_bar_choices %>%
#   str_which("\\[e\\]") %>% # str_subset doesn't keep names
#   search_bar_choices[.] %>%
#   names() %>%
#   unique() %>%
#   write_rds(path = glue::glue(work_dir, "/precalculated_data/extracellular_mets.rds"))
# search_bar_choices %>%
#   str_which("\\[c\\]") %>%
#   search_bar_choices[.] %>%
#   names() %>%
#   unique() %>%
#   write_rds(path = glue::glue(work_dir, "/precalculated_data/cytoplasm_mets.rds"))
