# MiamiPlot.R
# Function to generate ggplot2 based miami plot
# Julie D. White

plot.miami <- function(data, chr=NULL) {
  # Necessary packages
  require(rlang)
  require(dplyr)
  require(ggplot2)
  require(patchwork)
  require(ggrepel)
  
  # Identify the column containing chromosome information.
  if(!is.null(chr)){
    # Chromosome column specified. Find column index matching this value.
    chr <- grep(pattern = chr, x = names(data), ignore.case = T)
  }  else {
    # Chromosome column not specified. Find column index matching "chr" while ignoring case
    chr <- grep(pattern = "chr", x = names(data), ignore.case = T)
    if(is_empty(chr)){
      stop("I cannot find a chromosome column named 'chr'. Please specify using chr argument.")
    }
  }
  
}
