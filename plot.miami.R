# MiamiPlot.R
# Function to generate ggplot2 based miami plot
# Julie D. White

plot.miami <- function(data, chr=NULL, pos=NULL) {
  # Necessary packages
  require(rlang)
  require(dplyr)
  require(ggplot2)
  require(patchwork)
  require(ggrepel)
  
  # Get the dataframe 
  df <- data
  
  # Identify the column containing chromosome information.
  if(!is.null(chr)){
    # Chromosome column specified. Find column index matching this value.
    chr.name <- chr
    chr.indx <- grep(pattern = chr.name, x = names(df), ignore.case = T)
  }  else {
    # Chromosome column not specified. Find column index matching "chr" while ignoring case.
    chr.indx <- grep(pattern = "chr", x = names(df), ignore.case = T)
    chr.name <- names(df)[chr.indx]
    if(is_empty(chr.indx)){
      stop("I cannot find a chromosome column named 'chr'. Please specify using chr argument.")
    }
  }
  
  # Figure out if the chromosome column is numeric or not.
  if(!is.numeric(df[,chr.indx])) {
    stop("Please make sure your chromosome column is numeric.")
  }
  
  # Idenitfy the column containing position information. 
  if(!is.null(pos)){
    # Position column specified. Find column index matching this value.
    pos.name <- pos
    pos.indx <- grep(pattern = pos.name, x = names(df), ignore.case = T)
  }  else {
    # Position column not specified. Find column index matching "pos" while ignoring case.
    pos.indx <- grep(pattern = "pos", x = names(df), ignore.case = T)
    pos.name <- names(df)[pos.indx]
    if(is_empty(pos.indx)){
      stop("I cannot find a position column named 'pos'. Please specify using pos argument.")
    }
  }
  
  # Figure out if the position column is numeric or not.
  if(!is.numeric(df[,pos.indx])) {
    stop("Please make sure your position column is numeric.")
  }
  
  # To make the plot, we need to know the cumulative position of each probe/snp across the whole genome.
  df <- df %>% 
    group_by(!!sym(chr.name)) %>% 
    # Compute chromosome size
    summarise(chrlength = max(!!sym(pos.name))) %>%  
    # Calculate cumulative position of each chromosome
    mutate(cumulativechrlength = cumsum(as.numeric(chrlength))-chrlength) %>% 
    # Remove chr length column
    select(-chrlength) %>%
    # Temporarily add the cumulative length of each chromosome to the initial dataset 
    left_join(df, ., by=chr.name) %>%
    # Sort by chr then position 
    arrange(!!sym(chr.name), !!sym(pos.name)) %>%
    # Add the position to the cumulative chromosome length to get the position of this probe relative to all other probes
    mutate(rel.pos = !!sym(pos.name) + cumulativechrlength) %>%
    # Remove cumulative chr length column
    select(-cumulativechrlength)
  
}


