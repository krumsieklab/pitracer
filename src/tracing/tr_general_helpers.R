library(tidyverse)


# This function will find all strings that start with "HMDB" and make sure they
# come back in the old 5 digit, form.
#
# - Ignores any string not starting with HMBD
# - Throws an error if the string is not 9 or 11 characters long ("HMDB" plus 5 or 7 digits)
# - Throws an error if it's a 7-digit version, and the first two digits are not zero
# - Does some format checking whether characters instead of numbers are used in the 7 digit version
# - Converts 7-digit version to 5-digit version
#
# - This function ignores [x] compartment information at the end of the string... it'll simply add them back
hmdb5_fixer <- function(v) {
  
  # find all strings starting with HMDB
  ind <- which(startsWith(v, "HMDB"))
  
  # work through all of them
  v[ind] <- sapply(v[ind], function(s_full) {
    
    # compartment info?
    has_comp <- !is.na(str_match(s_full, "\\[([a-z])\\]")[,2])
    if (has_comp) {
      splitup <- strsplit(s_full, "\\[")
      s <- splitup[[1]][1]
      s_comp <- paste0("[",splitup[[1]][2])
      
    } else {
      s <- s_full
      s_comp <- ""
    }
    
    # must be 9 or 11 characters long (5 or 7 digits)
    if (nchar(s)!=9 && nchar(s)!=11) {
      stop(sprintf("Invalid HMDB identifier: %s... must have 5 or 7 digits", s))
    }
    # if it's a 7 digit version, then the first two ones cannot be used
    if (nchar(s) == 11) {
      # check that there are numbers at all
      if(any(is.na(suppressWarnings(as.numeric(substr(s,5,11)))))) {
        stop(sprintf("Invalid HMDB identifier: %s", s))
      }
      # check that numbers are zero
      if (as.numeric(substr(s, 5,6)) > 0) {
        stop(sprintf("Unsupported HMDB identifier: %s. For 7-digit identifiers, the first two digits have to be zero (HMDB00xxxxx). Note that the database does not contain any metabolites from the 7-digit HMDB space.", s))     
      }
    }
    # if it's a valid 7 digit version, remove the first two zeros
    if (nchar(s) == 11) {
      s <- sprintf("HMDB%s", substr(s,7,11))
    }
    # 5 digits ones are untouched
    
    # add back compartment info?
    if (has_comp) {
      s <- paste0(s, s_comp)
    }
    
    # return string
    s
    
    
    
  }) 
  
  # return vector
  v
}

# use global id_map data frame to map any string to recon
# will take care of compartment prioritization, i.e. picking them in the order of the vector that comes in
to_recon_mapper <- function(lst, comp_priority) {
 
  lst %>% sapply(function(node) {
    
    if (node %in% genes) {
      # if it's in the gene list, leave as is
      node
    } else {
      # not a gene, must be a metabolite
      
      # determine if the string already has an [x] compartment identifier in it
      has_comp <- !is.na(str_match(node, "\\[([a-z])\\]")[,2])
      
      if (has_comp) {
        # it has a compartment annotation in it, look up in respective id_map
        
        # go through cases
        case_when(
          # if it's in the compartmentalized recon names, leave as is
          node %in% id_map_withcomp$recon_id ~ node,
          # if it's in the recon names list, map it over to recon id
          node %in% id_map_withcomp$name ~ id_map_withcomp$recon_id[match(node, id_map_withcomp$name)],
          # if it's found in id (containing HMDB, KEGG, and PubChem), map it over from there
          node %in% id_map_withcomp$id ~ id_map_withcomp$recon_id[match(node, id_map_withcomp$id)]
        )
        # note: if no case is fulfilled, NA will be returned
        
      } else {
        # does not have compartment information look up in respective id_map

        # go through casess
        this_recon_id <- case_when(
          # if it's in the recon names list, leave as is
          node %in% id_map_nocomp$recon_id ~ node,
          # if it's in the recon names list, map it over to recon id
          node %in% id_map_nocomp$name ~ id_map_nocomp$recon_id[match(node, id_map_nocomp$name)],
          # if it's found in id (containing HMDB, KEGG, and PubChem), map it over from there
          node %in% id_map_nocomp$id ~ id_map_nocomp$recon_id[match(node, id_map_nocomp$id)]
        )
        
        # now go through priorities
        if (!is.na(this_recon_id)) {
          # metabolite was mapped

          # find all entries with compartments
          entries <- id_map_withcomp %>% 
            dplyr::filter(recon_id_nocomp == this_recon_id) %>%
            # clean up to have unique entries
            dplyr::select(-id, -id_type) %>%
            distinct()
          # go through priorities list and use first one
          m <- comp_priority %in% entries$comp_id
          # first TRUE in priority order will be used
          use_comp <- comp_priority[min(which(m))]
          # add to recon_id and return
          paste0(this_recon_id, "[", use_comp, "]")

        } else {
          # just return the NA
          NA
        }

      }
    }
  })

}

# test code
# "glc_D[c]" %>% hmdb5_fixer() %>% to_recon_mapper()
# "HMDB0000660[c]" %>% hmdb5_fixer() %>% to_recon_mapper()
# "glc_D" %>% hmdb5_fixer() %>% to_recon_mapper(comp_priority = priority)
# "HMDB0000660" %>% hmdb5_fixer() %>% to_recon_mapper(comp_priority = priority)
# "D-Fructose [c]" %>% hmdb5_fixer() %>% to_recon_mapper()
# "D-Fructose" %>% hmdb5_fixer() %>% to_recon_mapper(comp_priority = priority)


# verify that the user inputs are viable
# thow errors otherwise
tr_gen_input_checker <- function(trace_datatable,    # input comes either as this OR as the other three parameters
                                 start_nodes,
                                 end_nodes,
                                 all_vs_all = TRUE) {
  
  # can only input either a list of nodes or a dataframe
  if (!missing("trace_datatable")) {
    
    if (!missing("start_nodes") | !missing("end_nodes")) {
      stop("Either supply trace_datatable or start_nodes and end_nodes
           to use for the tracing, NOT both")
    }
    
    trace_datatable_nodes <- base::union(trace_datatable$from,
                                         trace_datatable$to)
    
    if (all_vs_all) {
      
      # Remove NAs if there are any
      start_nodes <- trace_datatable %>% filter(!is.na(from)) %>% .$from
      end_nodes <- trace_datatable %>% filter(!is.na(to)) %>% .$to
      
      trace_datatable <-
        expand.grid(from = start_nodes, to = end_nodes, weight = 1,
                    stringsAsFactors = FALSE) %>%
        as.data.table()
      
    } else {
      
      if(any(is.na(trace_datatable_nodes))) {
        stop("Expecting no NAs in the trace_datatable when all_vs_all == TRUE.
             Remove the NAs from trace_datatable to proceed")
      }
      
      # ensure that we have a data.table object. This removes factors, which
      # tibble usually assigns to strings
      trace_datatable <-
        trace_datatable %>%
        # currently set weights to 1
        dplyr::mutate(weight = 1) %>%
        as.data.table()
      
    }
    
    
  } else {
    
    # if not all-vs-all, there has to be an equal number of start and end nodes.
    if (!all_vs_all & length(start_nodes) != length(end_nodes)) {
      stop("There is an unequal number of start and end nodes. Check 'Tracing options' > 'All vs. all'")
    }
    
    if (!all_vs_all) {
      
      trace_datatable <- data.table(from = start_nodes,
                                    to = end_nodes,
                                    weight = 1)
      
    } else {
      
      trace_datatable <-
        expand.grid(from = start_nodes, to = end_nodes, weight = 1,
                    stringsAsFactors = FALSE) %>%
        as.data.table()
      
    }
    
    
  }
  
  # Make sure start and end nodes are not the same
  trace_datatable %>%
    tr_conv_id_converter(hmdb_recon) %>%
    filter(from != to) %>%
    distinct()
  
}


tr_gen_get_endpoints <- function(trace_df) {
  
  start_node <- trace_df$gene %>% unique()
  end_node <- trace_df$metabolite %>% unique()
  
  sprintf("%s ---> %s", start_node, end_node)
  
}


tr_gen_blacklist <- function(net, blacklist) {
  
  if (blacklist == "") return(net)
  
  net %>%
    filter(!str_detect(from, str_c(blacklist, collapse = "|")),
           !str_detect(to,   str_c(blacklist, collapse = "|")))
  
}
