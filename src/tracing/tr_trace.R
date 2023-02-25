library(yenpathy)
library(igraph)
library(data.table)
library(tidyverse)

work_dir <- getwd()
source(glue::glue(work_dir, "/src/tracing/tr_load_trace_helper.R"))

# Calculate shortest paths
tr_trace_path <- function(trace_datatable, net,
                          k = 1) {

  # pairs of nodes should be in the trace, else it will be set to FALSE
  net_nodes <- base::union(net$from, net$to)

  trace_datatable <-
    trace_datatable %>%
    dplyr::mutate(from_in_net = from %in% net_nodes,
                  to_in_net = to %in% net_nodes,
                  usable = if_else(!from_in_net | !to_in_net,
                                   FALSE, TRUE)) %>%
    as.data.table()

  if (all(!trace_datatable$usable)) {
    stop("There are NO traces between any from-to pairs.")
  }

  met_regex <- c("\\[" = "\\\\[", "\\]" = "\\\\]")
  trace_datatable <-
    trace_datatable %>%
    rowwise() %>%
    # NOTE: is relying on reachable column enough to
    # assure that a path will be found?
    dplyr::mutate(reachable =
                    ifelse(usable,
                           tr_get_subcomponent(in_distances, to) %>%
                             str_detect(
                               sprintf("^%s$",
                                       if_else(str_detect(from, "[\\[\\]]"),
                                               str_replace_all(from, met_regex),
                                               from))) %>%
                             any(),
                           FALSE))


  if (all(!trace_datatable$reachable)) {
    stop("There are NO traces between any from-to pairs.")
  }

  trace_datatable <-
    trace_datatable %>%
    # invoking arrange to enable NA assignment in the next step
    dplyr::arrange(desc(reachable)) %>%
    dplyr::mutate(
      association_path =
        ifelse(reachable,
               map2(from, to,
                    function(f, t) {
                      print(glue::glue("tracing {f} -> {t} for k = {k}"))
                      k_shortest_paths(graph_df = tr_prune_network(net, f, t),
                                       from = f,
                                       to = t,
                                       k = k)
                    }),
               c(NA))) %>%

    # NOTE: unnesting transforms NAs to NULLs, which could be a source
    # of error later?
    unnest(association_path)

  trace_datatable %>%
    dplyr::mutate(gene = from,
                  metabolite = to,
                  met_ids = to,
                  association_path) %>%
    dplyr::select(gene,
                  metabolite,
                  met_ids,
                  association_path,
                  reachable) %>%
    ungroup()

}

# Define custom piTracer function for speed --------------------------------

# This bypasses clustering and a slower step, visualization

fast_trace <- function(trace_datatable, # Either supply a data table
                       all_vs_all = TRUE,
                       start_nodes,     # or a pair of start and end node lists
                       end_nodes,
                       net,
                       k) {             # Network to use for tracing

  trace_table <-
    tr_gen_input_checker(trace_datatable = trace_datatable,
                         start_nodes = start_nodes,
                         end_nodes = end_nodes,
                         all_vs_all = all_vs_all) %>%
    tr_trace_path(net = net, k = k)

  trace_stats <-
    trace_table %>%
    dplyr::rename(from = gene,
                  to = metabolite,
                  has_trace = reachable) %>%
    group_by(from, to, has_trace) %>%
    summarise(association_path = list(association_path)) %>%
    arrange(desc(has_trace))

  trace_stats$from <- trace_stats$from %>%
    translate_metabolite_name_list(recon_lookup)

  trace_stats$to <- trace_stats$to %>%
    translate_metabolite_name_list(recon_lookup)

  trace_stats

}

fast_trace_loop <- function(start_nodes, end_nodes, network, k) {

  start_nodes %>%
    sapply(function(metabo) {
        tryCatch(
          {
            fast_trace(start_nodes = metabo,
                       end_nodes = setdiff(end_nodes, metabo),
                       all_vs_all = TRUE,
                       net = network,
                       k = k)

          },
          
          error = function(cond) {
            # warning(cond) # this used to be message(), which cause Shiny to crash
            return(NULL)
          }
        )
    }, simplify = FALSE, USE.NAMES = TRUE)

}
