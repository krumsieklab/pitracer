# This is the central code file that initializes the tracing machinery.
# It loads and preprocessed various variables and tables needed for tracing, and loads all helper scripts.
# This script also defines the central "net_tracer" function.

library(tidyverse)
library(readxl)
library(feather)



# Load networks ------------------------------------------------------------

work_dir <- getwd()

print("loading networks...")

if (!exists("ppitf_network")) ppitf_network <- readRDS(glue::glue(work_dir, "/precalculated_data/ppi_tf_edges.rds"))
if (!exists("met_network")) met_network <- readRDS(glue::glue(work_dir, "/precalculated_data/recon_edges.rds"))

# bind into one network
trace_net <-
  bind_rows(ppitf_network, met_network) %>%
  dplyr::select(-moa) %>%
  as.data.frame()


print("loading pre-calculated distance matrices...")

if (!exists("in_distances")) in_distances <- read_rds(glue::glue(work_dir, "/precalculated_data/distance_matrix_in.rds") )
if (!exists("out_distances")) out_distances <- read_rds(glue::glue(work_dir, "/precalculated_data/distance_matrix_out.rds") )
if (!exists("cluster_distance_matrix")) cluster_distance_matrix <- read_rds(glue::glue(work_dir, "/precalculated_data/distance_matrix_all.rds") )

print("loading finished.")


# Source all scripts ------------------------------------------------------

source(glue::glue(work_dir, "/src/tracing/visualization_params.R"))
source(glue::glue(work_dir, "/src/tracing/tr_cluster.R"))
source(glue::glue(work_dir, "/src/tracing/tr_general_helpers.R"))
source(glue::glue(work_dir, "/src/tracing/tr_id_convert.R"))
source(glue::glue(work_dir, "/src/tracing/tr_trace.R"))
source(glue::glue(work_dir, "/src/tracing/tr_visualize.R"))



# Define main function ----------------------------------------------------

net_tracer <- function(trace_datatable, # Either supply a data table
                       start_nodes,     # or a pair of start and end node lists
                       end_nodes,
                       all_vs_all = TRUE, # Whether to trace between row/list pairs
                       net,             # Network to use for tracing
                       blacklist,       # A blacklist of node to apply to net
                       k,               # number of paths to trace per node pair
                       k_cluster_min = 2, # minimum number of paths per node pair
                       # before clustering paths
                       show_clusters = TRUE) {


  if(!missing(blacklist)) net <- tr_gen_blacklist(net, blacklist)

  traces <-
    tr_gen_input_checker(trace_datatable = trace_datatable,
                         start_nodes = start_nodes,
                         end_nodes = end_nodes,
                         all_vs_all = all_vs_all) %>%

    tr_trace_path(net = net, k = k) %>%

    {. ->> trace_table} %>%

    filter(reachable) %>%
    group_by(gene, metabolite) %>%
    group_split() %>%

    # NOTE: all_distance just means it is a distance matrix calculated
    # on the undirected network
    sapply(tr_add_cluster, cluster_distance_matrix,
           k_cluster_min = k_cluster_min,
           simplify = FALSE, USE.NAMES = TRUE) %>%

    sapply(tr_visualize_trace, show_clusters = show_clusters,
           simplify = FALSE, USE.NAMES = TRUE) %>%

    setNames(lapply(., function(tr) tr$trace_name))

  trace_visualizaitons <-
    traces %>%
    sapply(function(tr) tr$visualization,
           simplify = FALSE, USE.NAMES = TRUE)

  trace_main_node_dicts <-
    traces %>%
    sapply(function(tr) tr$main_nodes_dict,
           simplify = FALSE, USE.NAMES = TRUE)


  trace_stats <-
    trace_table %>%
    dplyr::rename(from = gene,
                  to = metabolite,
                  has_trace = reachable) %>%
    group_by(from, to, has_trace) %>%
    summarise(association_path = list(association_path)) %>%
    dplyr::arrange(desc(has_trace))

  trace_stats$from <- trace_stats$from %>%
    translate_metabolite_name_list(recon_lookup)

  trace_stats$to <- trace_stats$to %>%
    translate_metabolite_name_list(recon_lookup)

  names(trace_stats$association_path) <- names(trace_visualizaitons)

  list(visualizations = trace_visualizaitons,
       results_table = trace_stats,
       main_node_dicts = trace_main_node_dicts)

}
