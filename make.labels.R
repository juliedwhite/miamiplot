# Make probe/snp/gene labels for adding to miami plot.
make.labels <- function(data, p.name, hits.label, top.n.hits = NULL){
  # If the user has requested labeling the top n hits: 
  # Check that top.n.hits is a positive natural number and
  # Check that the column name provided for labels is in the df.
  if(testCount(top.n.hits) & 
     testNames(hits.label, type = "named", subset.of = colnames(data))){
    # If the user has given one column name in "hits.label"
    if(length(hits.label) == 1){
      # Filter the data to find top n values and label with user requested cols.
      label.df <- data %>%
        arrange(desc(-log10(!!sym(p.name)))) %>%
        mutate(label = !!sym(hits.label)) %>%
        select(rel.pos, !!sym(p.name), label) %>%
        slice(1:top.n.hits)
      # Return label df
      return(label.df)
      
      # Or, if the user has given two column names in "hits.label"
    } else if(length(hits.label) == 2){
      # Filter the data to find top n values and label with user requested cols.
      label.df <- data %>%
        arrange(desc(-log10(!!sym(p.name)))) %>%
        mutate(label = paste0(!!sym(hits.label[1]), "\n", !!sym(hits.label[2]))) %>%
        select(rel.pos, !!sym(p.name), label) %>%
        slice(1:top.n.hits)
      # Return label.df
      return(label.df)
    } 
    
  # If instead hte user has given a character vector of which SNPs/probes/genes to label
  } else if(is.null(top.n.hits) & !is.null(hits.label)) {
    # Find which column contains the hits.label items
    hits.col.indx <- unname(which(sapply(data, function(x) any(x %in% hits.label))))
    # Filter the data to get the position of these labels
    label.df <- data %>%
      filter(!!sym(colnames(data)[hits.col.indx]) %in% hits.label) %>% 
      mutate(label = !!sym(colnames(data)[hits.col.indx])) %>%
      select(rel.pos, !!sym(p.name), label)
    # Return label.df 
    return(label.df)
  }
}

