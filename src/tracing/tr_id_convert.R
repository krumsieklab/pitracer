# Functions ----------------------------------------------------------------

# detect compartment or return default compartment
tr_conv_comp_priority <- function(df_comp) {

  if (!any(str_detect(df_comp, "\\["))) return(df_comp)

  comp_priority <- c("c", "m", "e", "r", "n", "g", "x", "l", "i")

  for (comp in comp_priority) {

    met_comp <-
      str_extract(df_comp, comp)

    if (!is.na(met_comp)) return(met_comp)
  }

}

tr_conv_id_converter <- function(trace_datatable, lookup) {

  trace_datatable %>%
    left_join(lookup, by = c("from" = "hmdb_id")) %>%
    left_join(lookup, by = c("to" = "hmdb_id"), suffix = c("_from", "_to")) %>%
    dplyr::mutate(original_from = from,
                  original_to = to) %>%
    transmute(from = if_else(!is.na(recon_id_from), recon_id_from, from),
              to = if_else(!is.na(recon_id_to), recon_id_to, to),
              weight,
              original_from, original_to)

}


# Load and create lookup --------------------------------------------------------

work_dir <- getwd()
id_lookup <- read_excel(glue::glue(work_dir, "/precalculated_data/metabolon_lookup_table.xlsx"))

hmdb <-
  id_lookup %>%
  filter(externalID_type == "HMDB",
         !is.na(externalID),
         !is.na(chemical_id)) %>%
  dplyr::transmute(chemical_id, hmdb_id = externalID)

recon <-
  id_lookup %>%
  filter(externalID_type == "Recon",
         !is.na(externalID),
         !is.na(chemical_id)) %>%
  dplyr::transmute(chemical_id, recon_id = externalID)

hmdb_recon <-
  inner_join(hmdb, recon, by = "chemical_id") %>%
  distinct() %>%
  # rows [446, 447, 24159, 24160] don't have compartment info!
  separate(recon_id, c("recon_id", "compartment"), "--") %>%
  filter(!is.na(compartment)) %>%
  rowwise() %>%
  dplyr::mutate(compartment = tr_conv_comp_priority(compartment)) %>%
  ungroup() %>%
  transmute(hmdb_id,
            recon_id = str_c(recon_id, "[", compartment, "]")) %>%
  distinct()

# Select first match (since multiple recon IDs may be mapped to one HMDB)
