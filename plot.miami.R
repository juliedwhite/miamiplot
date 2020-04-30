# MiamiPlot.R
# Function to generate ggplot2 based miami plot
# Julie D. White

plot.miami <- function(data, split.by, split.at, chr="chr", pos="pos") {
  # Necessary packages
  require(rlang)
  require(checkmate)
  require(dplyr)
  require(ggplot2)
  require(patchwork)
  require(ggrepel)
  
  # Rename data
  df <- data
  
  # Establish a new argument check object
  check <- makeAssertCollection()
  
  # Check the input 
  check.miami.input(data = data, split.by = split.by, split.at = split.at)
  
  # Identify column containing chromosome information.
  if(testNames(chr, type = "named", subset.of = colnames(df))) {
    # Chromosome column specified and is subset of colnames. Find column index matching this value.
    chr.name <- chr
    chr.indx <- which(colnames(df)==chr.name)
  } else {
    # Chromosome column name is not in the right format or is not a subset of colnames. 
    assertNames(chr, type = "named", subset.of = colnames(df), add = check)
  }
  
  if(testNames(pos, type = "named", subset.of = colnames(df))) {
    # Position column specified is a subset of colnames. Find column index matching this value.
    pos.name <- pos
    pos.indx <- which(colnames(df)==pos.name)
  } else {
    # Position column name is not in the right format or is not a subset of colnames. 
    assertNames(pos, type = "named", subset.of = colnames(df), add = check)
  }
  
  # At this point, check if there are any errors. 
  if(!check$isEmpty()) {
    reportAssertions(check)
  }
  
  # Test if chromosome column is numeric. 
  if(!testMultiClass(df[,chr.indx], classes = c("numeric", "integer"))){
    check$push(paste0("Your chromosome column is of class ",
                      class(df[,chr.indx]),
                      ". Please make sure your chromosome column is numeric or integer."))
  }

  # Test if position column is numeric. 
  if(!testMultiClass(df[,pos.indx], classes = c("numeric", "integer"))){
    check$push(paste0("Your position column is of class ",
                      class(df[,pos.indx]),
                      ". Please make sure your position column is numeric or integer."))
  }

  # At this point, check if there are any errors. 
  if(!check$isEmpty()) {
    reportAssertions(check)
  }
  
  # To make the plot, we need to know the cumulative position of each probe/snp across the whole genome.
  df <- df %>% 
    #Group by chromosome
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
  
  # However, we don't want to print the cumulative position on the x-axis, so we need to provide the chromosome position labels relative to the entire genome.
  axis.df = df %>% 
    #Group by chromosome
    group_by(!!sym(chr.name)) %>% 
    #Find the center of the chromosome
    summarize(chr.center=(max(rel.pos) + min(rel.pos))/2)
  
  # Depending on what the user has input for split.by and split.at, make top and bottom df.
  if(is.numeric(split.at)){
    # Create top df using numeric info.
    top.df <- df %>%
      filter(!!sym(split.by) >= split.at)
    # Create bottom df using numeric info
    bot.df <- df %>%
      filter(!!sym(split.by) < split.at)
  } else if(is.character(split.at)){
    # Create top df using character specified by user
    top.df <- df %>%
      filter(!!sym(split.by) == split.at)
    # Create bottom df using character not specified by user
    bot.df <- df %>%
      filter(!!sym(split.by) != split.at)
  }
  

}


