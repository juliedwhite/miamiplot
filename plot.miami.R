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
    if(is_empty(chr)){
      stop("I cannot find a chromosome column named 'chr'. Please specify using chr argument.")
    }
  }
  
  # Figure out if the chromosome column is numeric or not.
  if(!is.numeric(df[,chr.indx])) {
    stop("Please make sure your chromosome column is numeric.")
  }
  
  # Idenitfy the column containing position information. 
  if(!is.null(pos)){
    # Chromosome column specified. Find column index matching this value.
    pos <- grep(pattern = pos, x = names(df), ignore.case = T)
  }  else {
    # Chromosome column not specified. Find column index matching "pos" while ignoring case.
    pos <- grep(pattern = "pos", x = names(df), ignore.case = T)
    if(is_empty(pos)){
      stop("I cannot find a position column named 'pos'. Please specify using pos argument.")
    }
  }
  
  # Figure out if the position column is numeric or not.
  if(!is.numeric(df[,pos])) {
    stop("Please make sure your position column is numeric.")
  }
  
  # To make the plot, we need to know the cumulative position of each probe/snp across the whole genome.
  df <- df %>% 
    group_by(chr) %>% 
    # Compute chromosome size
    summarise(chrlength = max(pos)) %>%  
    # Calculate cumulative position of each chromosome
    mutate(cumulativechrlength = cumsum(as.numeric(chrlength))-chrlength) %>% 
    # Remove chr length column
    select(-chrlength) %>%
    # Temporarily add the cumulative length of each chromosome to the initial dataset 
    left_join(df, ., by=c(chr=chr)) %>%
    # Sort by chr then position 
    arrange(chr, pos) %>%
    # Add the position to the cumulative chromosome length to get the position of this probe relative to all other probes
    mutate(rel.pos = pos + cumulativechrlength) %>%
    # Remove cumulative chr length column
    select(-cumulativechrlength)
  
  
}


