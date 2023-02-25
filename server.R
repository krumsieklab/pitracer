# Server part of piTracer Shiny app
Sys.setenv(TZ=Sys.timezone()) # otherwise, we get errors in the dockerized Shiny server

library(igraph)
library(ipc)
library(future)
library(promises)
plan(multicore)

formals(source)$local <- TRUE
work_dir <- getwd()

source(glue::glue(work_dir, "/src/global.R"), local = TRUE)


server <- function(input, output, session) {
  
  hideTab(inputId = "main_panel", target = "ID map")
  download_overlay <- Waiter$new(
    html = tagList(
      spin_loaders(12, color = "grey"),
      h4("Preparing download...", style = "color: grey")
    ),
    color = waiter::transparent(.7)
  )
  
  observe_helpers()
  future::plan(multicore)
  
  results <- reactiveVal(list())
  running <- reactiveVal(FALSE)
  id_mapping_df <- reactiveVal(data.frame())
  options(future.globals.maxSize = 7480000000)
  
  # Load modules ----
  trace_params <- node_input_server("global")
  
  results_params <- results_server("global", results, id_mapping_df)
  
  retraced_results <- blacklist_server("global",
                                       results,
                                       reactive(results_params$trace_selector),
                                       reactive(results_params$visualization_selected),
                                       reactive(trace_params))
  
  # Output: Trace button / Cancel button ----
  output$trace_or_cancel <- renderUI({
    if (running()) {
      loadingButton("cancel_trace",
                    "Cancel",
                    class = "btn btn-danger",
                    loadingSpinner = "spinner",
                    loadingLabel = "Cancelling...",
                    style = "width: 100%")
    } else {
      loadingButton("trace",
                    "Trace",
                    class = "btn btn-primary",
                    loadingLabel = "Tracing...",
                    style = "width: 100%")
    }
  })
  
  # Button click:: Trace ----
  observeEvent(input$trace, {
    
    tryCatch({
      
      running(TRUE)
      
      toggleButtons("tracing")
      
      Sys.sleep(1)
      
      # figure out priority order, move the one the user selected to the beginning
      def_comp <- trace_params$compartment_default %>% substr(start = 2, stop = 2)
      ind <- which(priority==def_comp)
      priority_resorted <- c(priority[ind],priority[-ind])
      
      
      # split up input by lines, delete empty entries
      lst_start <- trace_params$sn %>% str_split("\n") %>% unlist() %>% purrr::discard(~ .x == "")
      lst_end <- trace_params$en %>% str_split("\n") %>% unlist() %>% purrr::discard(~ .x == "")
      
      
      # write back to user interface, so the user sees the filtered list
      # this currently doesn't work
      # txt_start <- reactiveValues(txt=lst_start %>% paste0(collapse = "\n"))
      # txt_end   <- reactiveValues(txt=lst_end %>% paste0(collapse = "\n"))
      # updateTextAreaInput(session, "sn", value=txt_start$txt)
      # updateTextAreaInput(session, "en", value=txt_end$txt)
       
      
      # map input to an identifier that pitracer understands later
      # run on both start and end nodes as multi-line texts
      nodes_to_trace <- list(start=lst_start, end=lst_end) %>%
        # go through two lists
        lapply(function(ls) {
          ls %>%
            # fix HMDB5, throws meaningful errors to the user if bad HMDB entries are given
            hmdb5_fixer() %>%
            # map to recon IDs
            to_recon_mapper(comp_priority = priority_resorted)
        })
      
      
      
      # convert this mapped list to a direct "user_input" -> "translation" paired list, which can be given back to the user
      met_translation_df <- data.frame(
        # actual input
        user_input = list(start = lst_start,
                          end = lst_end) %>% unname() %>% unlist(),
        # translated version
        translation = nodes_to_trace  %>% unname() %>%  unlist()
      )
      
      # throw errors for all unmapped ones
      if (any(is.na(met_translation_df$translation))) {
        stop(sprintf("Could not find the following input identifiers in the database: %s",
                     paste0(met_translation_df$user_input[is.na(met_translation_df$translation)], collapse = ", ")))
      }

      
      # set in shiny app
      id_mapping_df(met_translation_df)
      

      # verify that the user inputs are viable, throw meaningful errors
      tr_gen_input_checker(start_nodes = nodes_to_trace$start,
                           end_nodes = nodes_to_trace$end,
                           all_vs_all = trace_params$all_vs_all)
      
      # set GUI parameters
      non_reactive_k <- trace_params$k
      non_reactive_cluster_paths <- trace_params$cluster_paths
      # run tracing
      trace_future <<- future({
        
        results <- net_tracer(start_nodes = nodes_to_trace$start,
                              end_nodes = nodes_to_trace$end,
                              all_vs_all = trace_params$all_vs_all,
                              net = trace_net,
                              k = non_reactive_k,
                              k_cluster_min = if_else(non_reactive_cluster_paths, 2, 0),
                              show_clusters = TRUE)
        
        results
        
      })
      
      # extract results, check if anything went wrong
      prom <- trace_future %...>% results()
      prom <- catch(trace_future,
                    function(err) {
                      print(err$message)
                      if (str_detect(err$message, "Failed to retrieve the result of MulticoreFuture")) {
                        alert_title <- ""
                        alert_message <- "Tracing cancelled by user"
                        alert_type = "info"
                      } else {
                        alert_title <- "Error"
                        alert_message <- err$message
                        alert_type = "error"
                      }
                      shinyalert(title = alert_title,
                                 text = alert_message,
                                 type = alert_type,
                                 confirmButtonCol = "#2196F3"
                      )
                    })
      
      trace_future <- finally(trace_future, function() {
        toggleButtons("traced")
        running(FALSE)
      })
      
      NULL
    },
    error = function(err) {
      running(FALSE)
      results() %>%
        isolate() %>%
        identical(list()) %>%
        if_else(toggleButtons("init"),
                toggleButtons("traced"))
      
      shinyalert(title = "Error",
                 text = err$message,
                 type = "error",
                 confirmButtonCol = "#2196F3"
      )
    })
  })
  
  # Reaction: Cancel trace ----
  observeEvent(input$cancel_trace, {
    stopMulticoreFuture(trace_future)
  })
  
  # Button: Download results ----
  output$download_results <- downloadHandler(
    
    filename = paste0(Sys.Date(), "_tracing_results.zip"),
    
    content = function(zip_fn) {
      
      download_overlay$show()
      fn_list <- c()
      
      tmpdir <- tempdir()
      setwd(tmpdir)
      
      # Create html with all visualizations
      fn <- "trace_visualizations.html"
      fn_list <- c(fn_list, fn)
      
      results <- results() %>% isolate()
      
      template_file <- "www/results_main.Rmd"
      tempReport <- file.path(tmpdir, str_remove(template_file, "www/"))
      file.copy(file.path(work_dir, template_file), tempReport, overwrite = TRUE)
      
      # markdown options
      params <- list(work_dir = work_dir,
                     trace_names = names(results$visualizations),
                     trace_visualizations = results$visualizations %>%
                       lapply(function(trace) {
                         trace %>%
                           visOptions(
                             width = "900px"
                           )
                       }))
      
      # build report
      rmarkdown::render(tempReport,
                        output_file = fn,
                        output_dir = tmpdir,
                        params = params,
                        envir = new.env(parent = globalenv()))
      
      
      
      # Create excel file with results table
      table_fn <- "results_table.xlsx"
      fn_list <- c(fn_list, table_fn)
      
      results$results_table %>%
        unnest(association_path) %>%
        rowwise() %>%
        dplyr::mutate(trace_length = length(association_path)) %>%
        ungroup() %>%
        dplyr::mutate(trace = map2(association_path, has_trace,
                                   function(ap, ht){
                                     if(ht){
                                       translate_metabolite_name_list(ap, recon_lookup) %>% str_c(collapse = "|")
                                     } else {
                                       NA
                                     }
                                   }) %>% unlist() %>% as.character(),
                      trace_reconIDs = ifelse(!has_trace, NA, association_path),
                      trace_reconIDs = map(trace_reconIDs, function(p) str_c(p, collapse = "|")) %>% as.character()) %>%
        group_by(from, to) %>%
        dplyr::mutate(trace_rank = dense_rank(trace_length),
                      trace_order = row_number()) %>%
        ungroup() %>%
        dplyr::mutate(trace_rank       = ifelse(!has_trace, NA, trace_rank),
                      trace_order      = ifelse(!has_trace, NA, trace_order),
                      trace_length     = ifelse(!has_trace, NA, trace_length))  %>%
        dplyr::select(from:has_trace, trace, trace_reconIDs,
                      everything(), -association_path) %>%
        distinct() %>%
        as.data.frame() %>%
        openxlsx::write.xlsx(table_fn, rowNames = FALSE)
      
      # Add translation table if nodes were entered as KEGG, HMDB or PubChem
      id_mapping_df <- id_mapping_df() %>% isolate()
      if (nrow(id_mapping_df) > 0) {
        id_map_fn <- "metabolite_translation.xlsx"
        fn_list <- c(fn_list, id_map_fn)
        id_mapping_df %>%
          openxlsx::write.xlsx2(id_map_fn, rowNames = FALSE)
      }
      
      # Add .graphml objects for each visualization
      for (trace in names(results$visualizations)) {
        graphml_fn <- trace %>%
          str_replace(" ---> ", " to ") %>%
          str_remove_all("\\s\\[[a-z]\\]") %>%
          str_remove("\\s\\([0-9]*\\)") %>%
          str_replace_all(" ", "_") %>%
          sprintf("%s.graphml", .)
        
        fn_list <- c(fn_list, graphml_fn)
        
        vis <- results$visualizations[[trace]]
        
        edges <- vis$x$edges %>%
          mutate(
            edge.lty = dashes %>%
              if_else(2, 1), # 2 is dashed in igraph, 1 is solid
            .keep = "unused"
          ) %>%
          rename(
            edge.width = width,
            edge.color = color
          )
        
        vertices <- vis$x$nodes %>%
          mutate(
            vertex.size = size/10,
            vertex.shape = case_when(
              shape == "diamond" ~ "Diamond",
              shape == "dot" ~ "Ellipse",
              shape == "star" ~ "Triangle"
            ),
            .keep = "unused"
          )
        
        graph_from_data_frame(edges, vertices = vertices, directed = TRUE) %>%
          write_graph(graphml_fn, format = "graphml")
      }
      
      download_overlay$hide()
      
      system2("zip", args = paste(zip_fn, fn_list, sep = " "))
      
    },
    
    contentType = "application/zip"
  )
  
  observeEvent(retraced_results$trigger, {
    results(retraced_results$results)
  })
  
  drug_effects_server("drugs")
  
  # Output: Guide ----
  output$guide <- renderUI({
    HTML(markdown::markdownToHTML(
      knit("www/guide.Rmd", quiet = TRUE)
    ))
  })
  
  toggleButtons("init")
}