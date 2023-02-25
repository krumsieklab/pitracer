# Run gene ranking from code, outside of Shiny app
# This script can take a long time (potentially several hours, depending on the number of metabolites used as input)

# Runtime with 1 core: 3614.054 sec elapsed

# Load requires functions and objects ----
work_dir <- getwd()
source(glue::glue(work_dir, "/src/global.R")) # also sources all generanking_ codes in this folder

# Order in which compartments are attempted to be assigned for metabolites that 
# do not have an explicit compartment specification [x]. If the first compartment
# is not found, the next one will be tried.
# c = cytosol, e=extracellular, m=mitochondrial, etc.
comp_priority <- c("c", "m", "e", "r", "n", "g", "x", "l", "i")

# Add path to file containing two columns: one column for metabolite names
# (recognized are Recon, KEGG, HMDB, PubChem IDs and metabolite names); and one
# column that reflects statistical significance (recognized are "t"/"true" for
# statistically significant, "f"/"false" for stat. insignificant,
# case-insensitive). Supports .csv, .xls and .xlsx formats.
# metabolites_input_path <- "www/example_gene_targets_short.xlsx
metabolites_input_path <- "www/example_gene_targets.xlsx"

# Note: This file contains the metabolites used in the paper. Results might 
# slightly differ numerically, which is due to the non-deterministic shortest 
# path algorithm used in the algorithm.

# verify that the file exists
if (!file.exists(metabolites_input_path)) stop(sprintf("Input file %s does not exist!", metabolites_input_path))

# If the provided metabolites file has a header
file_header <- TRUE




# =============================================================================
# From here on out the code can be run without making any changes
# =============================================================================

# Run calculation -------------------------------------------------------------

# read metabolite list
extension <- tools::file_ext(metabolites_input_path)
# when reading semicolon separated files,
# having a comma separator causes `read.csv` to error
metabolites_input <- tryCatch({
  
  if (extension == "csv") {
    metabolites_input <- data.table::fread(metabolites_input_path,
                                           header = file_header)
  } else if (str_detect(extension, "^xlsx$|^xls$")) {
    metabolites_input <- readxl::read_excel(metabolites_input_path,
                                            col_names = file_header)
  }
  
},
error = function(e) {
  stop(safeError(e))
})

# run everything, with timing
library(tictoc)
tictoc::tic()
# main functional call
results <- calculate_generanks(
  user_input = metabolites_input,
  comp_priority = comp_priority
)
tictoc::toc()



# Save results ----------------------------------------------------------------

# Save gene scores
writexlsx_append(
  results$main_results,
  file.path(getwd(), "generank_results.xlsx"),
  sheetName = "Gene scores",
  append = FALSE,
  rowNames = FALSE
)

# Save metabolite contributions
writexlsx_append(
  results$met_contributions,
  file.path(getwd(), "generank_results.xlsx"),
  sheetName = "Metabolite contributions",
  append = TRUE,
  rowNames = FALSE
)

# Save pathway scores
writexlsx_append(
  results$pathway_scores,
  file.path(getwd(), "generank_results.xlsx"),
  sheetName = "Pathway scores",
  append = TRUE,
  rowNames = FALSE
)

# Save gene scores
writexlsx_append(
  results$user_input_mapping,
  file.path(getwd(), "generank_results.xlsx"),
  sheetName = "Input comparison",
  append = TRUE,
  rowNames = FALSE
)






# View results ----------------------------------------------------------------

# View gene scores
View(results$main_results)
# View metabolite contributions
View(results$met_contributions)
# View pathway scores
View(results$pathway_scores)
# show top 3 pathways per gene
View(results$pathway_scores %>% top_n(3))
# View table comparing user input and what the function recognized
View(results$user_input_mapping)
