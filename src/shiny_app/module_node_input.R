# Main piTracer node input UI elements and server codes

node_input_UI <- function(id) {

  tagList(
    tags$head(tags$style(HTML(
      "#global-user_file {
        margin-left: 10px;
      }
      "
    ))),

    fluidRow(
      column(6, h4("Search for nodes")),
      column(6, h5("Add as"))
    ),

    fluidRow(
      column(6, style = 'padding-top: 6px;',
             selectizeInput(NS(id, "search_bar"),
                            label = NULL,
                            choices = NULL,
                            multiple = TRUE,
                            options = list(
                              closeAfterSelect = TRUE,
                              placeholder = "Type to search for nodes"
                            ),
                            width = "100%"
             )
      ),
      column(3, style = 'padding-top: 6px; padding-right: 2px;',
             actionButton(NS(id, "add_to_sn"),
                          label = "Start node",
                          class = "btn-primary",
                          style = "width: 100%; font-size: 85%;")
      ),
      column(3, style = 'padding-top: 6px; padding-left: 2px;',
             actionButton(NS(id, "add_to_en"),
                          label = "End node",
                          class = "btn-primary",
                          style = "width: 100%; font-size: 85%")
      )
    ),
    
    fluidRow(
      column(5, h6("Compartment of added node:")),
      column(1,
             div(
               icon("question-circle") %>%
                 bs_embed_tooltip("The metabolic model background contains metabolites in different compartments. This option defines the compartment in which metabolites will be searched in the field above, and the compartment specification that will be added to the search list below.",
                                  placement = "right"),
               style = "display: inline-block;
                                        vertical-align: center;
                                        float: right;
                                        padding-top: 9px;"
             )
      ),
      column(6, style = "padding-top: 4px;",
             selectInput(NS(id, "compartment"),
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
    selectizeInput(
      NS(id, "examples"),
      label = "Examples",
      choices = c(
        "",
        "Citrate → Oxaloacetate",
        "ALDH18A1 → L-Citrulline",
        # `"MCCC1 → 3-Hydroxy-Isovaleryl Carnitine",
        "DBH → 3-Methoxy-4-Hydroxymandelate",
        # "LDHA/GTF2H1/SAA1 → 2-Hydroxy-Isovalerate",
        "Signalling Examples"
      ),
      options = list(
        placeholder = "Choose an example",
        onInitialize = I('function() { this.setValue(""); }')
      )
    ),

    br(),

    fluidRow(
      column(6,
             textAreaInput(NS(id, "sn"), h4("Start nodes"),
                           height = "200px",
                           value = "",
                           resize = "vertical",
                           placeholder = "Enter gene symbols (HGNC), KEGG/HMDB/PubChem IDs or Recon names/IDs for metabolites"
             )
      ),
      column(6,
             textAreaInput(NS(id, "en"), h4("End nodes"),
                           height = "200px",
                           value = "",
                           resize = "vertical"
             )
      )
    ),

    fluidRow(
      column(6),
      column(6,
             actionButton(NS(id, "clear_nodes"), "Clear nodes",
                          class = "btn-primary",
                          style = "width: 100%")
      )
    ),

    br(),
    
    fluidRow(
      bsCollapse(id = "collapse_group",
                 bsCollapsePanel(
                   # title = tags$div(
                   HTML("<i class='fa fa-caret-right pull-left' style='position: absolute; padding-top: 1px;'></i>
                           <i class='fa fa-caret-down pull-left' style='position: absolute;'></i><div style='position: relative; left: 15px;'>File upload</div>"),
                   # ),
                   "Input file:",
                   style = "info",
                   
                   checkboxInput(NS(id, "header"), "Header", TRUE),
                   
                   fileInput(NS(id, "user_file"),
                             label = "Choose file",
                             multiple = TRUE,
                             accept = c(".xlsx", ".xls", ".csv")
                   ),
                   tagList(a("Download example file", href="./example_shin_associations.xlsx"), " (File format documented in User Guide)")
                 )
      )
    ),
    
    blacklist_UI("global"),
    
    fluidRow(
      bsCollapse(id = "collapse_trace_options",
                 
                 bsCollapsePanel(
                   # title = tags$div(
                   HTML("<i class='fa fa-caret-right pull-left' style='position: absolute; padding-top: 1px;'></i>
                           <i class='fa fa-caret-down pull-left' style='position: absolute;'></i>
                           <div style='position: relative; left: 15px;'>Tracing options</div>"),
                   # ),
                   style = "info",
                   
                   fluidRow(
                     column(6,
                            checkboxInput(NS(id, "all_vs_all"),
                                          label = "All vs. all",
                                          value = TRUE) %>%
                              shinyInput_label_embed(
                                icon("question-circle") %>%
                                  bs_embed_tooltip("Find traces between all combinations of start & end nodes, else pairwise",
                                                   placement = "right")
                              )
                     ),
                     column(6,
                            checkboxInput(NS(id, "cluster_paths"),
                                          label = "Cluster paths",
                                          value = TRUE) %>%
                              shinyInput_label_embed(
                                icon("question-circle") %>%
                                  bs_embed_tooltip("Checking this option will cluster paths into biologically similar paths. Clusters represent groups of paths that may share biological characteristics, similar to biological pathways.",
                                                   placement = "right")
                              )
                     )
                   ),
                   
                   fluidRow(
                     column(1, style = "display: inline-block; vertical-align: center;",
                            div(h4("k:"))
                     ),
                     column(4, style = "padding-top: 7px;",
                            tags$style("#k {background-color: white;}"),
                            numericInput(NS(id, "k"),
                                         label = NULL,
                                         value = 10,
                                         min = 1)
                     ),
                     column(1,
                            div(
                              icon("question-circle") %>%
                                bs_embed_tooltip("Maximum distance between start & end node",
                                                 placement = "right"),
                              style = "display: inline-block;
                                        vertical-align: center;
                                        float: right;
                                        padding-top: 13px;"
                            )
                     )
                   ),
                   
                   
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
                            selectInput(NS(id, "compartment_default"),
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
                   )
                   
                   
                 )
      )
    )

    

  )
}

node_input_server <- function(id) {
  moduleServer(id, function(input, output, session) {

    # Button: Add nodes to start nodes ----
    observeEvent(input$add_to_sn, {

      req(input$search_bar)

      added_nodes <- input$search_bar %>%
        sapply(function(node) {
          print(node)
          node_wo_other_ids <- node %>%
            str_split(" \\(\\s*(?=[^(]+$)") %>%
            sapply("[[", 1)

          if (node_wo_other_ids %in% genes) {
            return(node_wo_other_ids)
          } else {
            sprintf("%s [%s]",
                    node_wo_other_ids,
                    input$compartment %>%
                      substr(start = 2, stop = 2)) %>%
              return()
          }
        }) %>%
        paste(collapse = "\n")

      if (input$sn != "") {
        new_start_nodes <- paste(input$sn, added_nodes, sep = "\n")
      } else {
        new_start_nodes <- added_nodes
      }

      updateTextAreaInput(session, "sn",
                          value = new_start_nodes)

      updateSelectizeInput(session, "search_bar",
                           choices = c(met_choices, genes),
                           selected = NULL,
                           server = TRUE)

    })

    # Button: Add nodes to end nodes ----
    observeEvent(input$add_to_en, {

      req(input$search_bar)

      added_nodes <- input$search_bar %>%
        sapply(function(node) {
          node_wo_other_ids <- node %>%
            str_split(" \\(\\s*(?=[^(]+$)") %>%
            sapply("[[", 1)

          if (node_wo_other_ids %in% genes) {
            return(node_wo_other_ids)
          } else {
            sprintf("%s [%s]",
                    node_wo_other_ids,
                    input$compartment %>%
                      substr(start = 2, stop = 2)) %>%
              return()
          }
        }) %>%
        paste(collapse = "\n")

      if (input$en != "") {
        new_end_nodes <- paste(input$en, added_nodes, sep = "\n")
      } else {
        new_end_nodes <- added_nodes
      }

      updateTextAreaInput(session, "en",
                          value = new_end_nodes)

      updateSelectizeInput(session, "search_bar",
                           choices = c(met_choices, genes),
                           selected = NULL,
                           server = TRUE)

    })

    observeEvent(input$examples, {

      if (input$examples == "Citrate → Oxaloacetate") {

        updateTextAreaInput(session, "sn",
                            value = "Citrate [e]")
        updateTextAreaInput(session, "en",
                            value = "Oxaloacetate [e]")
        updateNumericInput(session, "k",
                           value = 10)

      } else if (input$examples == "ALDH18A1 → L-Citrulline") {

        updateTextAreaInput(session, "sn",
                            value = "ALDH18A1")
        updateTextAreaInput(session, "en",
                            value = "L-Citrulline [e]")
        updateNumericInput(session, "k",
                           value = 10)

      } else if (input$examples == "MCCC1 → 3-Hydroxy-Isovaleryl Carnitine") {

        updateTextAreaInput(session, "sn",
                            value = "MCCC1")
        updateTextAreaInput(session, "en",
                            value = "3-Hydroxy-Isovaleryl Carnitine [e]")
        updateNumericInput(session, "k",
                           value = 10)

      } else if (input$examples == "DBH → 3-Methoxy-4-Hydroxymandelate") {

        updateTextAreaInput(session, "sn",
                            value = "DBH")
        updateTextAreaInput(session, "en",
                            value = "3-Methoxy-4-Hydroxymandelate [e]")

      } else if (input$examples == "LDHA/GTF2H1/SAA1 → 2-Hydroxy-Isovalerate") {

        updateTextAreaInput(session, "sn",
                            value = "LDHA\nGTF2H1\nSAA1")
        updateTextAreaInput(session, "en",
                            value = "2-Hydroxy-Isovalerate [e]")
        updateCheckboxInput(session, "all_vs_all",
                            value = TRUE)
        updateNumericInput(session, "k",
                           value = 10)

      } else if (input$examples == "Signalling Examples") {

        updateTextAreaInput(session, "sn",
                            value = "NOD1\nNODAL\nSTK3")
        updateTextAreaInput(session, "en",
                            value = "NFKB1\nLEFTY1\nTEAD1")
        updateCheckboxInput(session, "all_vs_all",
                            value = FALSE)
        updateNumericInput(session, "k",
                           value = 10)

      }
    })

    # Button: Clear nodes ----
    observeEvent(input$clear_nodes, {

      updateTextAreaInput(session, "sn",
                          value = "")

      updateTextAreaInput(session, "en",
                          value = "")
      
      updateSelectizeInput(session, "examples", selected = '')
      
    })

    # Dropdown: Compartment location ----
    observeEvent(input$compartment, {

      compartment_location <- input$compartment %>%
        substr(start = 2, stop = 2)

      if (compartment_location == "e") {
        met_choices <<- extracellular_mets
      } else if (compartment_location == "c") {
        met_choices <<- cytoplasm_mets
      } else {
        comp_metabolites <- full_metabolite_database %>%
          filter(str_detect(comp_id, compartment_location)) %>%
          arrange(name)

        met_choices <<- comp_metabolites$name %>%
          unique() %>%
          sapply(function(met_name) {

            non_recon_ids <- comp_metabolites %>%
              filter(name == met_name, !is.na(id))

            if (nrow(non_recon_ids) > 0) {
              non_recon_ids %>%
                filter(name == met_name) %>%
                arrange(match(id_type, c("KEGG", "HMDB", "PubChem"))) %>%
                .$id %>%
                paste(collapse = ", ") %>%
                sprintf("%s (%s)", met_name, .) %>%
                return()
            } else {
              return(met_name)
            }


          }, USE.NAMES = FALSE)
      }

      updateSelectizeInput(session, "search_bar",
                           choices = c(met_choices, genes),
                           server = TRUE)
    })

    # Reaction: k slider ----
    observeEvent(input$k, {

      req(input$k)

      if (input$k < 1) {
        showFeedbackWarning(
          inputId = "k",
          text = "Please input a k with at least a value of 1"
        )
      } else {
        hideFeedback("k")
      }

    })

    # Reaction: User uploaded file ----
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
             message = "Please provide a file with two columns: one for start nodes, one for end nodes"
        )
      )
      print(df)

      updateCheckboxInput(session, "all_vs_all", value = FALSE)

      updateTextAreaInput(session, "sn",
                          value = str_c(na.omit(df[[1]]), collapse = "\n"))

      updateTextAreaInput(session, "en",
                          value = str_c(na.omit(df[[2]]), collapse = "\n"))

    })

    return(input)
  })

}