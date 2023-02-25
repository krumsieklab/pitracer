tr_get_subcomponent <- function(distance_matrix, node, radius_bound = Inf) {

  dists <- tryCatch({ distance_matrix[node,] },
                    error = function(e) {
                      glue::glue("Node {node} does not exist in the network. Check your distance matrix for reference.") %>%
                        stop()
                    })

  dists[dists < radius_bound] %>% names()

}

tr_get_intersection <- function(node1, node2, radius1, radius2) {

  intersect(tr_get_subcomponent(out_distances, node1, radius1),
            tr_get_subcomponent(in_distances, node2, radius2))

}

tr_radius_search <- function(node1, node2, min_overlap_size = 5) {

  i <- 0
  j <- 0

  inc_i <- TRUE
  inc_j <- FALSE

  overlap_size <- 0
  # throw error if no nodes can be found in common after a max distance
  while (overlap_size < min_overlap_size) {

    if (inc_i) {

      i <- i + 1
      inc_i <- FALSE

    } else {

      j <- j + 1
      inc_i <- TRUE

    }

    overlap_size <-
      tr_get_intersection(node1, node2, i, j) %>% length()

  }

  list(node1_radius = i, node2_radius = j)

}

# This function basically calculates all the shortest distance lengths between
# a starting/ending node and interface genes going through the intermediate
# nodes
tr_get_path_lengths <- function(
  node_reachables,
  node,
  boundary__genes,
  direction,
  # n_paths_min ensures that each boundary gene
  # selects a path_length_min that has at least
  # n_paths_min paths
  n_paths_min = 5
) {

  if (direction == "from") {

    dist_mat <- out_distances[node, node_reachables, drop = FALSE]
    dist_mat_complement <- in_distances[boundary__genes, node_reachables, drop = FALSE]

  } else if (direction == "to") {

    dist_mat <- in_distances[node, node_reachables, drop = FALSE]
    dist_mat_complement <- out_distances[boundary__genes, node_reachables, drop = FALSE]

  } else {
    stop('input either "from" or "to" for direction parameter')
  }


  lapply(boundary__genes, function(bg) {

    bg_reach <- dist_mat_complement[bg, node_reachables, drop = FALSE]
    bg_reach <- bg_reach[, bg_reach < Inf, drop = FALSE]

    path_lengths <- bg_reach + dist_mat[, colnames(bg_reach), drop = FALSE]
    # added plus 1 to go above the shortest length, which usually includes the node
    # itself, the start/end node, and a very few other nodes with the same length

    n_paths <- 0
    adder <- 0
    while (n_paths <= n_paths_min) {

      # ensure that you get back results if there really are
      # only ncol(path_lengths) <= n_paths_min
      if (ncol(path_lengths) <= n_paths_min) return(path_lengths %>% colnames())

      path_length_min <- min( path_lengths[, path_lengths > 0]) + adder

      res <- path_lengths[, path_lengths <= path_length_min, drop = FALSE] %>% colnames()
      if (path_length_min == Inf) break

      n_paths <- res %>% length()

      adder <- adder + 1
    }

    res

  }) %>%
    unlist() %>%
    unique()

}

tr_get_pruned_nodes <- function(boundary__genes,
                                node,
                                node_reachables,
                                node_radius,
                                direction,
                                dist_matrix) {

  # drop ensures that row name is retained for subsetting one row
  dist_subset <- dist_matrix[boundary__genes, node_reachables, drop = FALSE]

  node_bg_reachables <-
    boundary__genes %>%
    lapply(function(bg) tr_get_subcomponent(dist_subset, bg, node_radius)) %>%
    unlist() %>%
    unique()

  tr_get_path_lengths(node_reachables,
                      node,
                      boundary__genes,
                      direction = direction)

}

tr_prune_network_simply <- function(net, node_1, node_2) {

  from_reachables <- tr_get_subcomponent(out_distances, node_1)
  to_reachables <- tr_get_subcomponent(in_distances, node_2)

  common_nodes <- intersect(from_reachables, to_reachables)

  net %>%
    filter(from %in% common_nodes,
           to   %in% common_nodes)

}


tr_prune_network_extensively <- function(net, node_1, node_2) {

  node_radii <- tr_radius_search(node_1, node_2)

  node1_radius <- node_radii$node1_radius
  node2_radius <- node_radii$node2_radius

  node1_reachables <- tr_get_subcomponent(out_distances, node_1, node1_radius)
  node2_reachables <- tr_get_subcomponent(in_distances, node_2, node2_radius)

  boundary_genes <- intersect(node1_reachables, node2_reachables)

  print(glue::glue("prunning network extensively for trace {node_1} --> {node_2}"))
  node1_pruned <- tr_get_pruned_nodes(boundary_genes,
                                      node_1,
                                      node1_reachables,
                                      node1_radius,
                                      direction = "from",
                                      in_distances)


  node2_pruned <- tr_get_pruned_nodes(boundary_genes,
                                      node_2,
                                      node2_reachables,
                                      node2_radius,
                                      direction = "to",
                                      out_distances)

  commons <- union(node1_pruned, node2_pruned)

  pruned_net <-
    net %>%
    filter(from %in% commons,
           to   %in% commons)

  print(glue::glue("pruned network down to {nrow(pruned_net)} edges"))

  pruned_net

}

tr_prune_network <- function(net,
                             node_1, node_2,
                             # when network has > n_edge_threshold
                             # edges, then prune network extensively
                             n_edge_threshold = 20000
                             ) {

  pruned_net <- tr_prune_network_simply(net, node_1, node_2)

  if (nrow(pruned_net) > n_edge_threshold) {

    print(glue::glue("network currently has {nrow(pruned_net)} edges"))
    pruned_net <- tr_prune_network_extensively(net, node_1, node_2)

  }

  pruned_net

}