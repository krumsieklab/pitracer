# Defines various hard-coded visualization parameters

library(tidyverse)
library(data.table)


# Which nodes to show? ----------------------------------------

translate_names <- TRUE
show_neighboring_rxns <- FALSE

# remove only metaabolite or genes option as well
remove_minor_met_nodes <- TRUE
remove_main_secondary_met_nodes <- FALSE

remove_minor_gene_nodes <- TRUE
remove_main_secondary_gene_nodes <- FALSE



# Use short reaction names? -----------------------------------------------

use_short_rxn_names <- TRUE


# Shape and size parameters -----------------------------------------------

## EDGE WIDTH SIZES
edge_main_width <- 1.5
edge_main_secodary_width <- 0.5
edge_minor_width <- 0.5

## NODE SHAPES
node_rxn_shape <- "diamond"
node_endpoint_shape <- "star"
node_all_others_shape <- "dot"

## NODE SHAPE SIZES
node_endpoint_size <- 70
node_main_size <- 20
node_main_secondary_size <- 10
node_rxn_size <- 15
node_minor_size <- 10

## FONT SIZES
node_endpoint_font_size <- 50
node_main_font_size <- 10
node_main_secondary_font_size <- 10
node_minor_font_size <- 10




# Color parameters ------------------------------------------------

## EDGE COLORS
edge_main_color <- "black"
edge_main_secondary_color <- "#9E9E9E"
edge_minor_color <- "#ECECEC" # "#D4D4D4"# "#BEBEBE"


## NODE COLORS
node_source_color <- "#9911AA"
node_target_color <- "#F34F41"
node_action_color <- "#FFCF35"

node_minor_metabolite_color <- "#ECECEC" # "#FFDFBE"
node_minor_gene_color <- "#D4D4D4" # "#E0C6E2"
node_minor_action_color <- "#FFF979"
node_minor_reversible_border_color <- "#F7E5EA"

node_metabolite_ggm_color <- "#FFB347"
node_metabolite_recon_color <- "#FFB1B0"
node_gene_color <- "#CC99FF"

node_reversible_border_color <- "#FFB347"
node_crowded_border_color <- "#FF6961"

## NODE FONT COLORS
node_main_font_color <- "black"
node_main_secondary_font_color <- "gray"
node_minor_font_color <- "#DEDEDE"



# Visualization legends ---------------------------------------------------

# nodes
lnodes <-
  rbind(
    c("Recon\nMetabolite",            "dot",     node_metabolite_recon_color, 0,   "#973365", 10, 20),
    c("GGM\nMetabolite",              "dot",     node_metabolite_ggm_color,   0,   "#973365", 10, 20),
    c("Gene\n",                       "dot",     node_gene_color,             0,   "#973365", 10, 20),
    c("In >5\nReactions",             "dot",     node_metabolite_recon_color, 5,   node_crowded_border_color, 10, 20),
    c("Interaction\nNode",            "diamond", node_action_color,           0,   "#973365", 10, 20),
    c("Reversible\nInteraction\nNode","diamond", node_action_color,           5,   node_reversible_border_color, 10, 20),
    c("Minor\nInteraction\nNode",     "diamond", node_minor_action_color,     0,   "#973365", 10, 20),
    c("Start",                        "star",    node_source_color,           0,   "#973365", 10, 20),
    c("Target",                       "star",    node_target_color,           0,   "#973365", 10, 20)
  ) %>%
  as.data.table() %>%
  set_names(
    c("label", "shape", "color.background", "borderWidth", "color.border", "font.size", "size")
  ) %>%
  dplyr::mutate(label = as.factor(label),
                shape = as.factor(shape),
                color.background = as.factor(color.background),
                borderWidth = as.numeric(borderWidth),
                color.border = as.factor(color.border),
                font.size = as.numeric(font.size),
                size = as.numeric(size))

# edges
ledges <-
  data.frame(color = c("black", edge_minor_color, "black"),
             label = c("Main\nPath", "Minor\nPath", "GGM\nEdges"),
             dashes = c(FALSE, FALSE, TRUE),
             arrows = c("to", "to", "to"),
             width = c(3, 0.5, 1),
             font.size = 10)