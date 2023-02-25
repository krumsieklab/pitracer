# Blacklist definition-related UI elements and server codes

blacklist_UI <- function(id) {

  tagList(
    fluidRow(
      bsCollapse(id = "collapse_blacklist",
                 bsCollapsePanel(
                   # title = tags$div(
                     HTML("<i class='fa fa-caret-right pull-left' style='position: absolute; padding-top: 1px;'></i>
                         <i class='fa fa-caret-down pull-left' style='position: absolute;'></i><div style='position: relative; left: 15px;'>Blacklist</div>"),
                   # ),

                   style = "info",

                   selectizeInput(NS(id, "blacklist"),
                                  "Blacklist (metabolites to be excluded)",
                                  choices = NULL,
                                  multiple = TRUE,
                                  options = list(closeAfterSelect = TRUE)
                   ) %>%
                     shinyInput_label_embed(
                       icon("question-circle") %>%
                         bs_embed_tooltip("You can either search for molecules to blacklist, upload your own blacklist, or blacklist nodes directly from inside the visualization by selecting a node and pressing \"blacklist node\".
                                                        To apply the blacklist, press \"Retrace\".")
                     ),

                   fluidRow(
                     column(9,
                            actionButton(NS(id, "add_node"),
                                         "Blacklist selected node",
                                         class = "btn-primary",
                                         style = "width: 100%")
                     ),
                     column(3,
                            actionButton(NS(id, "clear_blacklist"),
                                         "Clear",
                                         style = "width: 100%")
                     )
                   ),

                   br(),
                   br(),

                   fluidRow(style = "vertical-align: top;",
                            column(6,
                                   radioButtons(NS(id, "blacklist_apply"),
                                                "Apply uploaded blacklist to traces:",
                                                choices = c("current", "all"),
                                                selected = character(0),
                                                inline = TRUE)
                            ),
                            column(6,
                                   fileInput(NS(id, "user_blacklist"),
                                             "Upload a blacklist",
                                             accept = c(".txt"),
                                             placeholder = "") %>%
                                     bs_embed_tooltip("First, please choose if your blacklist file should be applied on the currently selected trace or on all traces"),
                                   tagList(a("Download example file", href="./example_blacklist.txt"), " (File format documented in User Guide)")
                            )
                        
                   ),

                   br(),
                   br(),

                   fluidRow(
                     column(6, uiOutput(NS(id, "retrace_or_cancel"))),
                     column(6,
                            downloadButton(NS(id, "download_blacklist"),
                                           "Download blacklist",
                                           class = "btn-primary",
                                           style = "width: 100%"
                            )
                     )
                   )
                 )
      )
    )
  )
}

blacklist_server <- function(id, results, selected_trace, selected_node, trace_params) {
  moduleServer(id, function(input, output, session) {

    retrace_running <- reactiveVal(FALSE)
    return_list <- reactiveValues(results = NULL, trigger = 0)
    trace_selector_choices <- reactiveVal()
    blacklists <- reactiveValues()
    net_blacklist <- reactiveValues()

    # Reaction: Create empty blacklists ----
    observeEvent(results(), {
      results()$visualizations %>%
        names() %>%
        trace_selector_choices()

      results()$visualizations %>%
        names() %>%
        isolate() %>%
        sapply(function(trace) {
          blacklists[[trace]] <- NULL
          net_blacklist[[trace]] <- NULL
        })

    })

    # Reactive: Trace nodes for the blacklist input ----
    trace_nodes_list <- eventReactive(results(), {
      results()$visualizations %>%
        names() %>%
        sapply(function(trace) {
          list(
            metabolites = filter(results()$visualizations[[trace]]$x$nodes,
                                 node_type == "recon"
            ) %>%
              .$id,
            genes = filter(results()$visualizations[[trace]]$x$nodes,
                           node_type == "gene"
            ) %>%
              .$id,
            interactions = filter(results()$visualizations[[trace]]$x$nodes,
                                  node_type == "rxn"
            ) %>%
              .$id
          )
        },
        simplify = FALSE,
        USE.NAMES = TRUE) %>%
        lapply(function(lst) {
          Filter(length, lst)
        })
    })

    # Reaction: Use all nodes from current trace as options for the blacklist
    observeEvent(selected_trace(), {

      req(trace_nodes_list(), selected_trace())

      updateSelectizeInput(session, "blacklist",
                           choices = trace_nodes_list()[[selected_trace()]],
                           selected = blacklists[[selected_trace()]]
      )

    })

    # Reaction: Changes in blacklist ----
    observeEvent(input$blacklist, {

      req(selected_trace())

      blacklists[[selected_trace()]] <- input$blacklist

    })

    # Button: Add node to blacklist ----
    observeEvent(input$add_node, {

      req(selected_node())

      blacklists[[selected_trace()]] <- c(blacklists[[selected_trace()]], selected_node())

      net_blacklist[[selected_trace()]] <- results()$main_node_dicts[[selected_trace()]] %>%
        dplyr::filter(translated_name %in% blacklists[[selected_trace()]],
                      !endpoint_node) %>%
        .$original_name

      updateSelectizeInput(session, "blacklist",
                           choices = trace_nodes_list()[[selected_trace()]],
                           selected = blacklists[[selected_trace()]]
      )

    })

    # Reaction: En-/disable uploading a blacklist if input$blacklist_apply is not selected ----
    observeEvent(input$blacklist_apply, {

      if (length(input$blacklist_apply) > 0) {
        shinyjs::enable("user_blacklist")
      }

    })

    # Button: Clear blacklist ----
    observeEvent(input$clear_blacklist, {
      blacklists[[selected_trace()]] <- NULL
      net_blacklist[[selected_trace()]] <- NULL

      updateSelectizeInput(session, "blacklist",
                           choices = trace_nodes_list()[[selected_trace()]],
                           selected = blacklists[[selected_trace()]]
      )
    })

    # Output: Retrace button / Cancel button ----
    output$retrace_or_cancel <- renderUI({
      if (retrace_running()) {
        loadingButton(NS(id, "cancel_retrace"),
                      "Cancel",
                      class = "btn btn-danger",
                      loadingSpinner = "spinner",
                      loadingLabel = "Cancelling...",
                      style = "width: 100%")
      } else {
        loadingButton(NS(id, "retrace"),
                      "Re-Trace",
                      class = "btn btn-primary",
                      loadingLabel = "Re-Tracing...",
                      style = "width: 100%") %>%
          if (identical(results(), list())) shinyjs::disabled(.) else .
      }
    })

    # Button: Retrace ----
    observeEvent(input$retrace, {

      tryCatch({

        retrace_running(TRUE)

        toggleButtons("tracing")
        Sys.sleep(1)

        old_results <- results() %>% isolate()

        # Fetch start and end node IDs
        start_node <-
          old_results$main_node_dicts[[selected_trace()]] %>%
          dplyr::filter(
            translated_name == str_split(selected_trace(), " ---> ") %>%
              sapply("[[", 1),
            endpoint_node) %>%
          .$original_name

        translated_start_node <-
          old_results$main_node_dicts[[selected_trace()]] %>%
          dplyr::filter(
            translated_name == str_split(selected_trace(), " ---> ") %>%
              sapply("[[", 1),
            endpoint_node) %>%
          .$translated_name

        end_node <-
          old_results$main_node_dicts[[selected_trace()]] %>%
          dplyr::filter(
            translated_name == str_split(selected_trace(), " ---> ") %>%
              sapply("[[", 2) %>%
              str_split(" \\(") %>%
              sapply("[[", 1),
            endpoint_node) %>%
          .$original_name

        translated_end_node <-
          old_results$main_node_dicts[[selected_trace()]] %>%
          dplyr::filter(
            translated_name == str_split(selected_trace(), " ---> ") %>%
              sapply("[[", 2) %>%
              str_split(" \\(") %>%
              sapply("[[", 1),
            endpoint_node) %>%
          .$translated_name

        bl_interactions <- trace_nodes_list()[[selected_trace()]]$interactions %>%
          subset(. %in% input$blacklist) %>%
          str_split(" \\(#") %>%
          sapply("[[", 1) %>%
          unique()

        rxn_ids <- subset(rxn_readable_names, rxn_name %in% bl_interactions) %>%
          .$rxn

        rxn_reac_prod <- subset(rxn_steps, rxn %in% rxn_ids) %>%
          .[, c("reactant", "product")] %>%
          unlist() %>%
          unname() %>%
          unique()

        blacklisted_net <-
          trace_net %>%
          filter(!(from %in% net_blacklist[[selected_trace()]]),
                 !(to %in% net_blacklist[[selected_trace()]]),
                 !(from %in% rxn_reac_prod),
                 !(to %in% rxn_reac_prod))

        trace_net %>% nrow() %>% paste("Length of old trace net:", .) %>% print()
        blacklisted_net %>% nrow() %>% paste("Length of new trace net:", .) %>% print()

        non_reactive_all_vs_all <- trace_params()$all_vs_all
        non_reactive_k <- trace_params()$k
        non_reactive_cluster_paths <- trace_params()$cluster_paths

        retrace_future <- future({

          retrace_result <-
            net_tracer(start_nodes = start_node,
                       end_nodes = end_node,
                       all_vs_all = non_reactive_all_vs_all,
                       net = blacklisted_net,
                       k = non_reactive_k,
                       k_cluster_min = if_else(non_reactive_cluster_paths, 2, 0),
                       show_clusters = TRUE)

          results_index <- paste(translated_start_node, "--->", translated_end_node) %>%
            grep(names(old_results$visualizations), fixed = TRUE)

          retrace_index <- paste(translated_start_node, "--->", translated_end_node) %>%
            grep(names(retrace_result$visualizations), fixed = TRUE)

          old_results$visualizations[[results_index]] <- retrace_result$visualizations[[retrace_index]]
          old_results$main_node_dicts[[results_index]] <- retrace_result$main_node_dicts[[retrace_index]]
          old_results$results_table <- old_results$results_table %>%
            dplyr::mutate(has_trace = replace(has_trace,
                                              from == retrace_result$results_table$from &
                                                to == retrace_result$results_table$to,
                                              retrace_result$results_table$has_trace),

                          association_path = replace(association_path,
                                                     from == retrace_result$results_table$from &
                                                       to == retrace_result$results_table$to,
                                                     retrace_result$results_table$association_path))

          old_results
        })

        prom <- retrace_future %...>% {. -> return_list$results}
        prom <- catch(retrace_future,
                      function(err) {
                        print(err$message)
                        shinyalert(title = "",
                                   text = "Tracing cancelled by user.",
                                   type = "info",
                                   confirmButtonCol = "#2196F3"
                        )
                      })

        retrace_future <- finally(retrace_future, function() {
          toggleButtons("traced")
          retrace_running(FALSE)
          return_list$trigger <- return_list$trigger + 1
        })

        NULL
      },

      error = function(err) {
        retrace_running(FALSE)
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

    # Reaction: Cancel retrace ----
    observeEvent(input$cancel_retrace, {
      stopMulticoreFuture(retrace_future)
    })

    # Button: Download blacklist ----
    output$download_blacklist <- downloadHandler(

      filename = function() {
        paste0(selected_trace(), Sys.Date(), ".txt")
      },

      content = function(file) {

        blacklist_nodes <- input$blacklist %>%
          unlist()

        # Filter which nodes to translate
        metabolite_ids <- blacklist_nodes[blacklist_nodes %in% recon_lookup$id]

        for (metabolite in metabolite_ids) {
          blacklist_nodes <- replace(blacklist_nodes,
                                     which(blacklist_nodes == metabolite),
                                     recon_lookup$name[which(recon_lookup$id == metabolite)])
        }

        blacklist_nodes %>%
          cat(file = file, sep = "\n")
      },

      contentType = "text/csv"
    )

    # Reaction: User uploaded blacklist ----
    observeEvent(input$user_blacklist, {

      req(selected_trace())

      df <- input$user_blacklist %>%
        .$datapath %>%
        read.delim(header = FALSE)

      print(na.omit(df[[1]]))

      if (!is.na(selected_trace())) {
        if (input$blacklist_apply == "current") {
          blacklists[[selected_trace()]] <<- na.omit(df[[1]])
        } else if (input$blacklist_apply == "all") {
          for (trace in trace_selector_choices()) {
            blacklists[[trace]] <<- na.omit(df[[1]])
          }
        }
      }

      blacklists[[selected_trace()]] <- na.omit(df[[1]])

    })

    return(return_list)
  })
}