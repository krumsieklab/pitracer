# UI part of piTracer Shiny app
Sys.setenv(TZ=Sys.timezone())

library(shinyBS)
library(bsplus)
library(shinyFeedback)
library(shinythemes)
library(shinyhelper)
library(shinybusy)
library(visNetwork)
library(waiter)
library(shinyjs)
library(knitr)
library(shinyalert)

formals(source)$local <- TRUE
work_dir <- getwd()

source(glue::glue(work_dir, "/src/global.R"), local = TRUE)


ui <- tagList(
  tags$style(type = "text/css",
             "body {padding-top: 70px;}"),
  tags$head(
    tags$link(
      rel = "stylesheet", type = "text/css", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css"
    ),
    tags$script(
      "function openGitHub()
      {
        window.open('https://github.com/krumsieklab/pitracer', '_blank');
      }"
    )
  ),

  useShinyFeedback(), # include shinyFeedback
  useShinyalert(),
  useShinyjs(),
  use_waiter(),
  # waiter_show_on_load(spin_loaders(12)),

  navbarPage("piTracer",

             theme = shinytheme("paper"),

             tabPanel("Cascade reconstruction",

                      use_bs_tooltip(),

                      tags$script(
                        HTML("var header = $('.navbar> .container-fluid');
                       header.append('<div style=\"float:right; margin:13px auto;\" type=\"button\" class=\"btn btn-default action-button\" onclick=\"openGitHub();\"><i class=\"fa-brands fa-github\"></i> piTracer GitHub</div>');
                       console.log(header)"
                        )
                      ),

                      # CSS ----
                      tags$head(tags$style(HTML('
                            ::-webkit-input-placeholder {
                              color: grey;
                            }
                            ::-webkit-placeholder {
                              color: grey;
                            }
                            ::-placeholder {
                              color: grey;
                            }
                            ::-moz-placeholder { /* Firefox 19+ */
                              color: grey;
                            }
                            :-ms-input-placeholder { /* IE 10+ */
                              color: grey;
                            }
                            textArea {
                              background-color: #ffffff !important;
                              border: none;
                            }
                            input,
                            .actual-edit-overlay {
                              -webkit-text-fill-color: #000 !important;
                            }
                            hr {border-top: 1px solid #bbbbbb;}
                            .tooltip {
                              transform: translateZ(0);
                            }
                            .tooltip-inner {
                              min-width: 180px;
                              white-space: pre-line;
                              text-align: justify;
                            }
                            .tooltip {font-size: 14px;}
                            .radio-inline {
                              margin-left: 5px;
                            }
                            a[aria-expanded=true] .fa-caret-right {
                               display: none;
                            }
                            a[aria-expanded=false] .fa-caret-down {
                               display: none;
                            }
                            body {
                                margin-bottom: 25px;
                            }
                            #footer {
                                position: fixed;
                                background-color: white;
                                height: 20px;
                                bottom: 0px;
                                left: 0px;
                                right: 0px;
                                margin-bottom: 0px;
                                margin-left: 4px;
                            }
                            .shiny-notification {
                              position: fixed;
                              top: calc(50%);
                              left: calc(50%);
                              width: 400px;
                            }
                            '))
                      ),

                      # sidebarlayout ----
                      sidebarLayout(
                        # Sidebar panel ----
                        sidebarPanel(

                          node_input_UI("global"),

                          #blacklist_UI("global"), # has been moved inside of node_input_UI

                          fluidRow(

                            column(6,
                                   uiOutput("trace_or_cancel")
                            ),

                            column(6,
                                   shinyjs::disabled(
                                     downloadButton(
                                       "download_results",
                                       "Download Results",
                                       class = "btn-primary",
                                       style = "width: 100%"
                                     )
                                   )
                            )
                          )
                        ),

                        # Main panel ----
                        results_ui("global")

                      )

             ),
             tabPanel("Gene target identification",
                      drug_effects_ui("drugs")
             ),
             tabPanel("User Guide",
                      tags$iframe(src = './guide.html',
                                  width = '100%', height = '900px',
                                  frameborder = 0, scrolling = 'auto'
                      )
             ),

             footer = tagList(
               tags$section("© KrumsiekLab", id = "footer")
             ),

             position = "fixed-top"
  )
)