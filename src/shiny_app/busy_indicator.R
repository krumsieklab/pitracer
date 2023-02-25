# Busy indicators for Shiny app
# source: https://github.com/daattali/advanced-shiny/tree/master/busy-indicator


withBusyIndicatorCSS <- "
.btn-loading-container {
margin-left: 10px;
font-size: 1.2em;
}
.btn-done-indicator {
color: green;
}
.btn-err {
margin-top: 10px;
color: red;
}
"

withBusyIndicatorUI <- function(button) {
  id <- button[['attribs']][['id']]
  div(
    shinyjs::useShinyjs(),
    singleton(tags$head(
      tags$style(withBusyIndicatorCSS)
    )),
    `data-for-btn` = id,
    button,
    span(
      class = "btn-loading-container",
      shinyjs::hidden(
        icon("spinner", class = "btn-loading-indicator fa-spin"),
        icon("check", class = "btn-done-indicator")
      )
    ),
    shinyjs::hidden(
      div(class = "btn-err",
          div(icon("exclamation-circle"),
              tags$b("Error: "),
              span(class = "btn-err-msg")
          )
      )
    )
  )
}

# Call this function from the server with the button id that is clicked and the
# expression to run when the button is clicked
withBusyIndicatorServer <- function(buttonId, expr) {
  # UX stuff: show the "busy" message, hide the other messages, disable the button
  loadingEl <- sprintf("[data-for-btn=%s] .btn-loading-indicator", buttonId)
  doneEl <- sprintf("[data-for-btn=%s] .btn-done-indicator", buttonId)
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  shinyjs::disable(buttonId)
  shinyjs::show(selector = loadingEl)
  shinyjs::hide(selector = doneEl)
  shinyjs::hide(selector = errEl)
  on.exit({
    shinyjs::enable(buttonId)
    shinyjs::hide(selector = loadingEl)
  })

  # Try to run the code when the button is clicked and show an error message if
  # an error occurs or a success message if it completes
  tryCatch({
    value <- expr
    shinyjs::show(selector = doneEl)
    shinyjs::delay(2000, shinyjs::hide(selector = doneEl, anim = TRUE, animType = "fade",
                                       time = 0.5))
    value
  }, error = function(err) { errorFunc(err, buttonId) })
}

# When an error happens after a button click, show the error
errorFunc <- function(err, buttonId) {
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  errElMsg <- sprintf("[data-for-btn=%s] .btn-err-msg", buttonId)
  errMessage <- gsub("^ddpcr: (.*)", "\\1", err$message)
  shinyjs::html(html = errMessage, selector = errElMsg)
  shinyjs::show(selector = errEl, anim = TRUE, animType = "fade")
}

toggleButtons <- function(mode = c("init", "tracing", "traced")) {
  if (mode == "init") {
      shinyjs::enable(NS("global", "search_bar"))
      shinyjs::enable(NS("global", "add_to_sn"))
      shinyjs::enable(NS("global", "add_to_en"))
      shinyjs::enable(NS("global", "ex_met_met"))
      shinyjs::enable(NS("global", "ex_gene_met"))
      shinyjs::enable(NS("global", "ex_gene_gene"))
      shinyjs::enable(NS("global", "sn"))
      shinyjs::enable(NS("global", "en"))
      shinyjs::enable(NS("global", "clear_nodes"))
      shinyjs::enable(NS("global", "compartment"))
      shinyjs::enable(NS("global", "all_vs_all"))
      shinyjs::enable(NS("global", "cluster_paths"))
      shinyjs::enable(NS("global", "k"))
      shinyjs::enable(NS("global", "header"))
      shinyjs::enable(NS("global", "user_file"))
      shinyjs::disable(NS("global", "blacklist"))
      shinyjs::disable(NS("global", "add_node"))
      shinyjs::disable(NS("global", "clear_blacklist"))
      shinyjs::disable(NS("global", "blacklist_apply"))
      shinyjs::disable(NS("global", "retrace"))
      shinyjs::disable(NS("global", "user_blacklist"))
      shinyjs::disable(NS("global", "download_blacklist"))
      shinyjs::enable(NS("global", "trace"))
      shinyjs::disable("download_results")
  } else if (mode == "tracing") {
      shinyjs::disable(NS("global", "search_bar"))
      shinyjs::disable(NS("global", "add_to_sn"))
      shinyjs::disable(NS("global", "add_to_en"))
      shinyjs::disable(NS("global", "ex_met_met"))
      shinyjs::disable(NS("global", "ex_gene_met"))
      shinyjs::disable(NS("global", "ex_gene_gene"))
      shinyjs::disable(NS("global", "sn"))
      shinyjs::disable(NS("global", "en"))
      shinyjs::disable(NS("global", "clear_nodes"))
      shinyjs::disable(NS("global", "compartment"))
      shinyjs::disable(NS("global", "all_vs_all"))
      shinyjs::disable(NS("global", "cluster_paths"))
      shinyjs::disable(NS("global", "k"))
      shinyjs::disable(NS("global", "header"))
      shinyjs::disable(NS("global", "user_file"))
      shinyjs::disable(NS("global", "blacklist"))
      shinyjs::disable(NS("global", "add_node"))
      shinyjs::disable(NS("global", "clear_blacklist"))
      shinyjs::disable(NS("global", "blacklist_apply"))
      shinyjs::disable(NS("global", "retrace"))
      shinyjs::disable(NS("global", "user_blacklist"))
      shinyjs::disable(NS("global", "download_blacklist"))
      shinyjs::disable(NS("global", "trace"))
      shinyjs::disable("download_results")
    } else if (mode == "traced") {
      shinyjs::enable(NS("global", "search_bar"))
      shinyjs::enable(NS("global", "add_to_sn"))
      shinyjs::enable(NS("global", "add_to_en"))
      shinyjs::enable(NS("global", "ex_met_met"))
      shinyjs::enable(NS("global", "ex_gene_met"))
      shinyjs::enable(NS("global", "ex_gene_gene"))
      shinyjs::enable(NS("global", "sn"))
      shinyjs::enable(NS("global", "en"))
      shinyjs::enable(NS("global", "clear_nodes"))
      shinyjs::enable(NS("global", "compartment"))
      shinyjs::enable(NS("global", "all_vs_all"))
      shinyjs::enable(NS("global", "cluster_paths"))
      shinyjs::enable(NS("global", "k"))
      shinyjs::enable(NS("global", "header"))
      shinyjs::enable(NS("global", "user_file"))
      shinyjs::enable(NS("global", "blacklist"))
      shinyjs::enable(NS("global", "add_node"))
      shinyjs::enable(NS("global", "clear_blacklist"))
      shinyjs::enable(NS("global", "retrace"))
      shinyjs::enable(NS("global", "blacklist_apply"))
      shinyjs::enable(NS("global", "user_blacklist"))
      shinyjs::enable(NS("global", "download_blacklist"))
      shinyjs::enable(NS("global", "trace"))
      shinyjs::enable("download_results")
    }
}
