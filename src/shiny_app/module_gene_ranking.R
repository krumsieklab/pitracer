# Drug effect prediction main interface, UI and server codes

library(DT)
library(future)
library(promises)
plan(multicore)

drug_effects_ui <- function(id) {
  tagList(

    tags$head(tags$style(HTML("
      a[aria-expanded=true] .fa-caret-right {
        display: none;
      }
      a[aria-expanded=false] .fa-caret-down {
        display: none;
      }
      "))),

    sidebarLayout(
      sidebarPanel(

        div(style = "height: auto; max-height: 500px; overflow-x: hidden; overflow-y: scroll",
            DTOutput(NS(id, "mets_table"),
                     height = "100%")
        ),

        br(),

        fluidRow(
          column(4, style = "padding-right: 0px;",
                 actionButton(NS(id, "add_row"),
                              "Add row",
                              class = "btn btn-primary",
                              width = "100%")
          ),
          column(4, style = "padding: 0px 4px;",
                 actionButton(NS(id, "delete_rows"),
                              "Delete rows",
                              class = "btn btn-primary",
                              width = "100%")
          ),
          column(4, style = "padding-left: 0px;",
                 actionButton(NS(id, "clear_table"),
                              "Clear",
                              class = "btn btn-secondary",
                              width = "100%")
          )
        ),

        br(),

        fluidRow(
          bsCollapse(id = "collapse_group",
                     bsCollapsePanel(
                       # title = tags$div(
                         HTML("<i class='fa fa-caret-right pull-left' style='position: absolute; padding-top: 1px;'></i>
                           <i class='fa fa-caret-down pull-left' style='position: absolute;'></i><div style='position: relative; left: 15px;'>File upload options</div>"),
                       # ),
                       "Input file (CSV or Excel)",
                       style = "info",

                       checkboxInput(NS(id, "header"), "Header", TRUE),

                       fileInput(NS(id, "user_file"),
                                 label = "Choose file",
                                 multiple = TRUE,
                                 accept = c(".xlsx", ".xls", ".csv")
                       ),
                       tagList(a("Download example file", href="./example_gene_targets.xlsx"), " (File format documented in User Guide)")
                       
                     )
          )
        ),
        
        # compartment panel
        fluidRow(
          column(5, h6("Default compartment:")),
          column(1,
                 div(
                   icon("question-circle") %>%
                     bs_embed_tooltip("Compartment of metabolites where no compartment has been specified. If the metabolite does not exist in this compartment, a standard priority order will be followed.",
                                      placement = "right"),
                   style = "display: inline-block;
                                        vertical-align: center;
                                        float: right;
                                        padding-top: 9px;"
                 )
          ),
          column(6, style = "padding-top: 4px;",
                 selectInput(NS(id, "drugrank_compartment"),
                             label = NULL,
                             choices = c("[c] - cytoplasm",
                                         "[e] - extracellular",
                                         "[m] - mitochondrion",
                                         "[r] - endoplasmic reticulum",
                                         "[n] - nucleus",
                                         "[g] - Golgi",
                                         "[x] - peroxisome",
                                         "[l] - lysosome")
                 )
          )
        ),
        
        # buttons
        fluidRow(
          column(6,
                 uiOutput(NS(id, "predict_or_cancel"))
          ),
          column(6,
                 # downloadButton(NS(id, "download_results"),
                 #                label = "Download results",
                 #                style = "width: 100%")
                 download_generank_UI(NS(id, "download_generank"))
          )
        ),
        fluidRow(
          column(12, p("This might take a long time, up to several hours"))
        )
      ),
      mainPanel(
        tabsetPanel(type = "tabs",
                    id = NS(id, "main_panel"),
                    tabPanel("Main results",
                             DT::DTOutput(NS(id, "main_results"))
                    ),
                    tabPanel("Metabolite contributions",
                             DT::DTOutput(NS(id, "met_contributions"))
                    ), 
                    tabPanel("Pathway scores",
                                DT::DTOutput(NS(id, "pathway_scores"))
                    ),
                    tabPanel("User input",
                             DT::DTOutput(NS(id, "user_input_mapping"))
                    )
        )
      )
    )
  )
}

drug_effects_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    future::plan(multicore)

    download_overlay <- Waiter$new(
      html = tagList(
        spin_loaders(12, color = "grey"),
        h4("Preparing download...", style = "color: grey")
      ),
      color = waiter::transparent(.7)
    )

    rv <- reactiveValues(
      input_dt = data.table(
        metabolite_name = c("enter a metabolite"),
        is_significant = c("")
      )
    )

    prediction_results <- reactiveVal(list())
    prediction_running <- reactiveVal(FALSE)

    output$mets_table <- renderDT(
      rv$input_dt,
      server = FALSE,
      editable = "cell",
      rownames = FALSE,
      colnames = c("Metabolite name", "significant"),
      options = list(
        dom = "t",
        scrollY = "300px",
        paging = FALSE,
        lengthChange = FALSE
      )
    )

    # Reaction: Store edits made to input table ----
    observeEvent(input$mets_table_cell_edit, {

      row <- input$mets_table_cell_edit$row
      col <- input$mets_table_cell_edit$col+1

      rv$input_dt[row, col] <- input$mets_table_cell_edit$value

    })

    # Button: Add row ----
    observeEvent(input$add_row, {

      rv$input_dt <- rv$input_dt %>%
        rbind(
          data.table(
            metabolite_name = "",
            is_significant = "false"
          )
        )

    })

    # Button: Delete rows ----
    observeEvent(input$delete_rows, {

      req(input$mets_table_rows_selected)

      rv$input_dt <- rv$input_dt[-input$mets_table_rows_selected, ]

    })

    # Button: Clear table ----
    observeEvent(input$clear_table, {

      rv$input_dt <- data.table(
        metabolite_name = "",
        is_significant = "false"
      )
    })

    # Button: User file ----
    observeEvent(input$user_file, {
      req(input$user_file)

      extension <- tools::file_ext(input$user_file$datapath)
      # when reading semicolon separated files,
      # having a comma separator causes `read.csv` to error
      tryCatch(
        {

          if (extension == "csv") {
            df <- data.table::fread(input$user_file$datapath,
                                    header = input$header)
          } else if (str_detect(extension, "^xlsx$|^xls$")) {
            df <- readxl::read_excel(input$user_file$datapath,
                                     col_names = input$header)
          }

        },
        error = function(e) {
          # return a safeError if a parsing error occurs
          stop(safeError(e))
        }
      )

      validate(
        need(ncol(df) == 2,
             message = "Please provide a file with two columns: one for metabolites, one for whether they are significant or not"
        )
      )

      validate(
        need(is.character(df[[1]]), "The first column needs to be text"),
        need(is.character(df[[2]]) | is.logical(df[[2]]), "The second column needs to be either logical or text")
      )

      names(df) <- c("metabolite_name", "is_significant")
      rv$input_dt <- df

    })

    # Output: Predict button / Cancel button ----
    output$predict_or_cancel <- renderUI({
      if (prediction_running()) {
        loadingButton(NS(id, "cancel_prediction"),
                      "Cancel",
                      class = "btn btn-danger",
                      loadingSpinner = "spinner",
                      loadingLabel = "Cancelling...",
                      style = "width: 100%")
      } else {
        loadingButton(NS(id, "calculate"),
                      "Predict targets",
                      class = "btn btn-primary",
                      loadingLabel = "Calculating...",
                      style = "width: 100%")
      }
    })

    # New Button: RUN!  Calculate gene scores ----
    observeEvent(input$calculate, {

      tryCatch({
        
        # verify user input
        validate(need(nrow(rv$input_dt) > 0, message = "User input empty"))

        # initialize UI
        prediction_running(TRUE)
        generank_progress <- ipc::AsyncProgress$new(
          session,
          message = "Running identification process... (this could take a long time)",
          min = 0,
          max = 100
        )

        # get user input
        user_input <- rv$input_dt %>% isolate()
        
        # figure out priority order, move the one the user selected to the beginning
        def_comp <- input$drugrank_compartment %>% substr(start = 2, stop = 2)
        ind <- which(priority==def_comp)
        priority_resorted <- c(priority[ind],priority[-ind])

        # run gene ranking
        prediction_future <<- future({
          # execute actual function
          calculate_generanks(user_input, comp_priority = priority_resorted, shinyProgress = generank_progress)
        })

        # UI controls
        prom <- prediction_future %...>% prediction_results()
        prom <- catch(prediction_future,
                      function(err) {
                        print(err$message)
                        if (str_detect(err$message, "Failed to retrieve the result of MulticoreFuture")) {
                          alert_title <- ""
                          alert_message <- "Process cancelled by user"
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

        prediction_future <- finally(prediction_future, function() {
          prediction_running(FALSE)
          generank_progress$close()
        })

        NULL
      },
      error = function(err) {
        prediction_running(FALSE)
        generank_progress$close()

        shinyalert(title = "Error",
                   text = err$message,
                   type = "error",
                   confirmButtonCol = "#2196F3"
        )
      })
    })

    # Reaction: Cancel prediction ----
    observeEvent(input$cancel_prediction, {
      stopMulticoreFuture(prediction_future)
    })

  
    download_generank_server("download_generank", prediction_results, download_overlay)

    # Reaction: Hide user input tab when no results exist ----
    observeEvent(prediction_results(), {

      if (!identical(prediction_results(), list())) {
        showTab(inputId = "main_panel", target = "User input")
      } else {
        hideTab(inputId = "main_panel", target = "User input")
      }

    })

    # Output: Main results ----
    output$main_results <- DT::renderDataTable(
      {
        if (!identical(prediction_results(), list())) {
          prediction_results()$main_results %>%
            mutate(
              score = score %>% round(2),
              .keep = "unused",
              .after = gene
            )
        } else {
          data.frame()
        }
      },
      rownames = FALSE,
      colnames = c("Gene", "Score", "Gene rank"),
      options = list(
        scrollY = "500px",
        pageLength = 25
      )
    )

    # Output: Metabolite contributions ----
    output$met_contributions <- DT::renderDataTable(
      {
        if (!identical(prediction_results(), list())) {
          prediction_results()$met_contributions %>%
            mutate(
              score = score %>% round(2),
              weighted_score = weighted_score %>% round(2),
              .keep = "unused",
              .after = met_name
            )
        } else {
          data.frame()
        }
      },
      rownames = FALSE,
      colnames = c("Gene", "Recon ID", "Metabolite name", "Score", "Weighted score", "Pathways"),
      options = list(
        scrollY = "500px",
        pageLength = 25
      )
    )
    
    # Output: Pathway scores ----
    output$pathway_scores <- DT::renderDataTable(
      {
        if (!identical(prediction_results(), list())) {
          prediction_results()$pathway_scores 
        } else {
          data.frame()
        }
      },
      rownames = FALSE,
      colnames = c("Gene", "Pathway", "Score"),
      options = list(
        scrollY = "500px",
        pageLength = 25
      )
    )

    # Output: User input comparison ----
    output$user_input_mapping <- DT::renderDataTable(
      {
        if (!identical(prediction_results(), list())) {
          prediction_results()$user_input_mapping
        } else {
          data.frame()
        }
      },
      rownames = FALSE,
      colnames = c("Input name", "Recognized as", "Stat. significance"),
      options = list(
        scrollY = "500px",
        pageLength = 25
      )
    )

  })

}