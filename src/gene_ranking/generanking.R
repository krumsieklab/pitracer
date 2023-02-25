# Entry point function to initialize gene ranking machinery

# settings
options(dplyr.summarise.inform = FALSE)
if (!exists("work_dir")) work_dir <- getwd()

# source functions
source(glue::glue(work_dir, "/src/gene_ranking/generanking_functions.R"), local = TRUE)

# load precalculated files
load("precalculated_data/generank_permanent_objects.rds", envir = globalenv()) 

in_distances_noppi_met_undir <- glue::glue(work_dir, "/precalculated_data/in_distances_noppi_met_undir.rds") %>%
  readRDS()
no_ppi_met_undir_ins <- readRDS(glue::glue(work_dir, "/precalculated_data/noppi_met_undir_ins.rds"))

