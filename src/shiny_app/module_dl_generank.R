# Drug prediction list download UI element and server codes


download_generank_UI <- function(id) {
  downloadButton(NS(id, "download_results"),
                 label = "Download results",
                 style = "width: 100%")
}

download_generank_server <- function(id, prediction_results, download_overlay) {
  moduleServer(id, function(input, output, session) {
    
    output$download_results <- downloadHandler(
      
      filename = sprintf("%s_target_prediction.zip", Sys.Date()),
      
      content = function(zip_fn) {
        
        download_overlay$show()
        fn_list <- c()
        
        tmpdir <- tempdir()
        setwd(tmpdir)
        
        generank_results <- prediction_results() %>% isolate()
        
        # prep main results
        main_results_fn <- "main_results.xlsx"
        fn_list <- c(fn_list, main_results_fn)
        generank_results$main_results %>%
          openxlsx::write.xlsx(main_results_fn, row.names = FALSE)
        
        # prep metabolite contributions
        met_contributions_fn <- "met_contributions.xlsx"
        fn_list <- c(fn_list, met_contributions_fn)
        generank_results$met_contributions %>%
          openxlsx::write.xlsx(met_contributions_fn, row.names = FALSE)
        
        # prep pathway scores
        pathway_scores_fn <- "pathway_scores.xlsx"
        fn_list <- c(fn_list, pathway_scores_fn)
        generank_results$pathway_scores %>%
          openxlsx::write.xlsx(pathway_scores_fn, row.names = FALSE)
        
        # prep input mapping
        user_input_mapping_fn <- "input_comparison.xlsx"
        fn_list <- c(fn_list, user_input_mapping_fn)
        generank_results$user_input_mapping %>%
          openxlsx::write.xlsx(user_input_mapping_fn, row.names = FALSE)
        
        download_overlay$hide()
        
        system2("zip", args = paste(zip_fn, fn_list, sep = " "))
        
      },
      
      contentType = "application/zip"
    )
  })
}