# MiamiPlot.R
# Function to generate ggplot2 based miami plot
# Julie D. White

plot.miami <- function(data, split.by, split.at, chr=NULL, pos=NULL) {
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
  
  # Check that df is a dataframe
  assertDataFrame(df, col.names = "named", add = check)
  
  # The user should supply the column name by which they want to split their data using split.by.
    # For example: beta values, study name, genotype status. 
  assertNames(split.by, type = "named", subset.of = colnames(df),  add = check)
  
  # The user should supply the point at which they want to split their data using split.at, as an integer or character.
    # For example: if they want positive beta values on top, then split.at = 0
    # For example, if they want "genotyped" SNPs on top, then split.at = "genotyped"
  # If the user has supplied a single character for split.at
  if(!testCharacter(split.at, len = 1)){
    if(!testNumeric(split.at, len = 1)){
      check$push("split.at must be either character or numeric, with length = 1")
    }
  } 
  
  # Identify the column containing chromosome information.
  if(testNull(chr)){
    # Chromosome column not specified. Find column index matching "chr" while ignoring case.
    if(testNames("chr", subset.of = tolower(colnames(df)))){
      # If this is true, then there is a column named "chr" or some permutation of it
      chr.indx <- grep(pattern = "chr", x = colnames(df), ignore.case = T)
      chr.name <- colnames(df)[chr.indx]
    } else {
      # I didn't find a column named chr. Ask the user to specify this.
      check$push("I cannot find a chromosome column named 'chr'. Please specify using chr argument.")
    }
  } else if(testNames(chr, type = "named", subset.of = colnames(df))) {
    # Chromosome column specified and is subset of colnames. Find column index matching this value.
    chr.name <- chr
    chr.indx <- which(colnames(df)==chr.name)
  } else {
    # Chromosome column specified, but not in the right format or not a subset of colnames.
    assertNames(chr, type = "named", subset.of = colnames(df), add = check)
  }
  
  # Identify the column containing position information.
  if(testNull(pos)){
    # Position column not specified. Find column index matching "pos" while ignoring case.
    if(testNames("pos", subset.of = tolower(colnames(df)))){
      # If this is true, then there is a column named "pos" or some permutation of it
      pos.indx <- grep(pattern = "pos", x = colnames(df), ignore.case = T)
      pos.name <- colnames(df)[pos.indx]
    } else {
      # I didn't find a column named pos. Ask the user to specify this.
      check$push("I cannot find a position column named 'pos'. Please specify using pos argument.")
    }
  } else if(testNames(pos, type = "named", subset.of = colnames(df))) {
    # Position column specified and is subset of colnames. Find column index matching this value.
    pos.name <- pos
    pos.indx <- which(colnames(df)==pos.name)
  } else {
    # Position column specified, but not in the right format or not a subset of colnames.
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


