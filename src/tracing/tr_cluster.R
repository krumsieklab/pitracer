library(dynamicTreeCut)
library(tidyverse)
library(data.table)

work_dir <- getwd()
source(glue::glue(work_dir, "/src/tracing/tr_general_helpers.R"))

tr_get_distances <- function(k_paths, dist_matrix) {

  # NOTE: path_nodes were used as main_nodes, this will cause ERRORS
  path_nodes <- k_paths %>% unlist() %>% unique()

  node_distances <-
    dist_matrix[path_nodes, path_nodes]

  # Calculate path distances
  path_dists <- c()
  for (path_i in 1:length(k_paths)) {

    path_dists_internal <- c()
    for(path_j in path_i:length(k_paths)) {

      dist_matrix <- node_distances[k_paths[[path_i]], k_paths[[path_j]]]

      row_size <- nrow(dist_matrix)
      col_size <- ncol(dist_matrix)

      row_idx <- 1
      row_idx_prev <- 0

      dist_score <- 0
      for (col_idx in 1:col_size) {

        row_idx <- row_idx_prev + (dist_matrix[row_idx:row_size, col_idx] %>%
                                     which.min() %>%
                                     as.numeric())
        row_idx_prev <-  row_idx - 1

        dist_score <- dist_score + dist_matrix[row_idx, col_idx]

      }

      # normalize by longer path length
      path_dists_internal[[path_j]] <- dist_score / (max(row_size, col_size) + 1)

    }

    path_dists[[path_i]] <- path_dists_internal

  }

  do.call(cbind, path_dists)

}


tr_cluster_paths <- function(path_distances) {

  minClusterSize = 2 # 5

  # this step is necessary, since "NULL" causes as.dist
  # to throw an error.
  path_distances[path_distances == "NULL"] <- NA

  # create a distance object variable
  dist_obj <- as.dist(path_distances, diag = TRUE)

  # Perform hierarchical clustering
  hclust_obj <- hclust(dist_obj, method = "complete")

  # Use dynamic tree cutting to get optimal clusters
  cutreeHybrid(
    # Input data: basic tree cutiing
    hclust_obj,

    dist_obj %>% as.matrix(),
    # Branch cut criteria and options
    cutHeight = NULL, minClusterSize = minClusterSize, deepSplit = TRUE,

    # External (user-supplied) measure of branch split
    externalBranchSplitFnc = NULL, minExternalSplit = NULL,
    externalSplitOptions = list(),
    externalSplitFncNeedsDistance = NULL,
    assumeSimpleExternalSpecification = TRUE,
    # PAM stage options
    pamStage = TRUE, pamRespectsDendro = TRUE,
    useMedoids = FALSE,
    maxPamDist = NULL,
    respectSmallClusters = TRUE,
    # Various options
    verbose = 1, indent = 0
  ) %>%
    .$labels %>%
    str_pad(width = 2, side = "left", pad = "0")

}



tr_add_cluster <- function(trace_df, dist_matrix, k_cluster_min = 2) {


  if (k_cluster_min == 0) {

    sprintf("assigning one cluster per path for trace: %s",
            tr_gen_get_endpoints(trace_df)) %>%
      print()

    clustered_paths <-
      trace_df %>%
      dplyr::mutate(cluster_num = row_number(),
                    cluster_num = str_pad(
                      cluster_num,
                      width =
                        max(cluster_num) %>%
                        nchar(),
                      side = "left",
                      pad = 0))

  } else {

    if(nrow(trace_df) < k_cluster_min) {

      sprintf("setting no clusters for trace: %s",
              tr_gen_get_endpoints(trace_df)) %>%
        print()

      clustered_paths <-
        trace_df %>%
        cbind(cluster_num = "01")

    } else {

      sprintf("assigning clusters to trace: %s",
              tr_gen_get_endpoints(trace_df)) %>%
        print()

      k_paths <- trace_df$association_path

      path_dist_mat <- tr_get_distances(k_paths, dist_matrix)
      path_clusters <- tr_cluster_paths(path_dist_mat)

      clustered_paths <-
        trace_df %>%
        cbind(cluster_num = path_clusters)

    }

  }

  clustered_paths

}