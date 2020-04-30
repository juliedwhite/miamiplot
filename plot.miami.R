# MiamiPlot.R
# Function to generate ggplot2 based miami plot
# Julie D. White

plot.miami <- function(data, split.by, split.at, chr="chr", pos="pos", p = "P") {
  # Necessary packages
  require(rlang)
  require(checkmate)
  require(dplyr)
  require(ggplot2)
  require(patchwork)
  require(ggrepel)
  
  # Establish a new argument check object
  check <- makeAssertCollection()
  
  # Check the input 
  check.miami.input(data = data, split.by = split.by, split.at = split.at)
  
  # Identify column index for chromosome information. 
  chr.indx <- find.col.indx(data = data, x = chr)
  chr.name <- chr
  
  # Identify column index for position information
  pos.indx <- find.col.indx(data = data, x = pos)
  pos.name <- pos
  
  # Identify column index for p-value information
  p.indx <- find.col.indx(data = data, x = p)
  p.name <- p
  
  # Test if chromosome, position, pvalue columns are numeric or integer. 
  check.miami.numeric(data = data, col.indx = chr.indx)
  check.miami.numeric(data = data, col.indx = pos.indx)
  check.miami.numeric(data = data, col.indx = p.indx)
  
  # To make the plot, we need to know the cumulative position of each probe/snp across the whole genome.
  data <- data %>% 
    #Group by chromosome
    group_by(!!sym(chr.name)) %>% 
    # Compute chromosome size
    summarise(chrlength = max(!!sym(pos.name))) %>%  
    # Calculate cumulative position of each chromosome
    mutate(cumulativechrlength = cumsum(as.numeric(chrlength))-chrlength) %>% 
    # Remove chr length column
    select(-chrlength) %>%
    # Temporarily add the cumulative length of each chromosome to the initial dataset 
    left_join(data, ., by=chr.name) %>%
    # Sort by chr then position 
    arrange(!!sym(chr.name), !!sym(pos.name)) %>%
    # Add the position to the cumulative chromosome length to get the position of this probe relative to all other probes
    mutate(rel.pos = !!sym(pos.name) + cumulativechrlength) %>%
    # Remove cumulative chr length column
    select(-cumulativechrlength)
  
  # However, we don't want to print the cumulative position on the x-axis, so we need to provide the chromosome position labels relative to the entire genome.
  axis.data = data %>% 
    #Group by chromosome
    group_by(!!sym(chr.name)) %>% 
    #Find the center of the chromosome
    summarize(chr.center=(max(rel.pos) + min(rel.pos))/2)
  
  # Depending on what the user has input for split.by and split.at, make top and bottom data.
  if(is.numeric(split.at)){
    # Create top data using numeric info.
    top.data <- data %>%
      filter(!!sym(split.by) >= split.at)
    # Create bottom data using numeric info
    bot.data <- data %>%
      filter(!!sym(split.by) < split.at)
  } else if(is.character(split.at)){
    # Create top data using character specified by user
    top.data <- data %>%
      filter(!!sym(split.by) == split.at)
    # Create bottom data using character not specified by user
    bot.data <- data %>%
      filter(!!sym(split.by) != split.at)
  }
  
  
  
  
}


