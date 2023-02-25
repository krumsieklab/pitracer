# Main piTracer network result UI elements and server codes

results_ui <- function(id) {
  tagList(
    mainPanel(
      tabsetPanel(type = "tabs",
                  id = NS(id, "main_panel"),
                  tabPanel("Trace selection",
                           fluidRow(
                             column(9,
                                    style = "display: inline-block; vertical-align: center; padding: 10px;",
                                    selectizeInput(NS(id, "trace_selector"),
                                                   "Select trace:",
                                                   choices = "",
                                                   width = "100%",
                                                   options = list(
                                                     placeholder = "",
                                                     onInitialize = I('function() { this.setValue(""); }'))
                                    )
                             ),
                             column(3,
                                    style = "display: inline-block; vertical-align: center; padding-top: 14px;",
                                    selectizeInput(NS(id, "sort_options"),
                                                   " ",
                                                   choices = c("By: start/end nodes",
                                                               "By: path length (ascending)",
                                                               "By: path length (descending)"),
                                                   selected = NULL,
                                                   multiple = TRUE,
                                                   width = "100%",
                                                   options = list(placeholder = "Sort by: start/end nodes",
                                                                  maxItems = 1))
                             )
                           ),

                           visNetworkOutput(NS(id, "visualization"),
                                            width = "100%",
                                            height = "60vh"),

                           shinyjs::hidden(imageOutput(NS(id, "visLegend")))
                  ),
                  tabPanel("Trace list",
                           DT::dataTableOutput(NS(id, "results_table"))
                  ),
                  tabPanel("ID mapping",
                           DT::dataTableOutput(NS(id, "id_mapping"))
                  )
      )
    )
  )

}

results_server <- function(id, results, id_mapping_df) {
  moduleServer(id, function(input, output, session) {

    zoom_overlay <- Waiter$new(id = NS(id, "visualization"),
                               html = h4("Scroll to zoom", style = "color: grey"),
                               color = waiter::transparent(.5)
    )
    zoom_shown <<- FALSE

    trace_selector_choices <- reactiveVal()

    # Reactive: Trace selector choices ----
    observeEvent(results(), {
      results()$visualizations %>%
        names() %>%
        trace_selector_choices()
    })

    # Reaction: Update trace selector when traces change ----
    observeEvent(trace_selector_choices(), {
      req(trace_selector_choices())

      updateSelectizeInput(session, "trace_selector",
                           choices = trace_selector_choices(),
                           selected = trace_selector_choices()[1]
      )

    })

    # Reaction: Trace selector ----
    observeEvent(input$trace_selector, {

      req(input$trace_selector)

      output$visualization <-
        renderVisNetwork({
          results()$visualizations[[input$trace_selector]]
        })

      shinyjs::show("visLegend")

      if (!zoom_shown) {
        shinyjs::delay(10, {
          zoom_overlay$show()
          Sys.sleep(3)
          zoom_overlay$hide()
        })
        zoom_shown <<- TRUE
      }

    })

    # Output: Results table ----
    output$results_table <- DT::renderDataTable(
      {
        if (!identical(results(), list())) {
          results()$results_table %>%
            ungroup() %>%
            rowwise() %>%
            transmute(from,
                      to,
                      has_trace = ifelse(has_trace, "True", "False"),
                      number_of_paths = length(association_path)
            )
        } else {
          data.table()
        }
      },
      options = list(scroller = TRUE,
                     scrollY = "400px",
                     paging = FALSE,
                     lengthChange = FALSE),
      rownames = FALSE,
      colnames = c("From", "To", "Has Trace", "Number of Paths")
    )

    # Reaction: Show/hide ID translation tab ----
    observeEvent(id_mapping_df(), {
      req(id_mapping_df())

      if (nrow(id_mapping_df()) > 0) {
        showTab(inputId = "main_panel", target = "ID mapping")
      } else {
        hideTab(inputId = "main_panel", target = "ID mapping")
      }
    })

    # Output: ID translation ----
    output$id_mapping <- DT::renderDataTable(
      id_mapping_df(),
      options = list(scroller = TRUE,
                     scrollY = "400px",
                     paging = FALSE,
                     lengthChange = FALSE),
      rownames = FALSE,
      colnames = c("User submitted ID", "Internal database ID")
    )

    # Reaction: Sort options ----
    observeEvent(input$sort_options, {

      req(trace_selector_choices())

      switch(input$sort_options,
             "By: start/end nodes" = trace_selector_choices() %>%
               .[order(.)],

             "By: path length (ascending)" = trace_selector_choices() %>%
               sapply(str_split, pattern = "\\(") %>%
               sapply(tail, 1) %>%
               sapply(str_split, pattern = "\\)") %>%
               sapply(head, 1) %>%
               sapply(as.numeric) %>%
               .[order(.)] %>%
               names(),

             "By: path length (descending)" = trace_selector_choices() %>%
               sapply(str_split, pattern = "\\(") %>%
               sapply(tail, 1) %>%
               sapply(str_split, pattern = "\\)") %>%
               sapply(head, 1) %>%
               sapply(as.numeric) %>%
               .[order(., decreasing = TRUE)] %>%
               names()
      ) %>%
        trace_selector_choices()

    })

    # Reaction: Change legend image ----
    output$visLegend <- renderImage({

      image_file <- "./www/visLegend_wo_ggm.png"

      list(src = image_file,
           height = "125px",
           style = "display: block; margin-left: auto; margin-right: auto;")
    }, deleteFile = FALSE)


    return(input)
  })

}