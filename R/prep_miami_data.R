#' Prepare data for miami plot
#'
#' @param data A data.frame object. Required.
#' @param split.by A character vector. The name of the column to use for
#'   splitting into top and bottom sections of the Miami plot. Required.
#' @param split.at A character or numeric vector. If numeric, the top plot will
#'   contain your results where the values in the \code{split.by} column are
#'   >= \code{split.at}. If character, top plot will contain your results where
#'   the values in the \code{split.by} column are equal to \code{split.at}. Required.
#' @param chr The name of the column containing your chromosome information.
#'   Defaults to "chr"
#' @param pos The name of the column containing your position information.
#'   Defaults to "pos"
#' @param p The name of the column containing your p-value information.
#'   Defaults to "P"
#' @param genomewideline Should we draw a genome-wide significance line? Set to
#'   FALSE if you do not want this line. Defaults to 5e-8, or supply your own number.
#' @param suggestiveline Should we draw a suggestive significance line? Set to
#'   FALSE if you do not want this line. Defaults to 1e-5, or supply your own number.
#' @param hits.label Either the name of the column(s), max. 2, to use for
#'   automatically labeling n hits, determined using \code{top.n.hits}, or a
#'   character vector of probes/genes/SNPs to label. If supplying a list of
#'   genes, it is helpful to also supply \code{top.n.hits} to limit the labels
#'   to the top N results, though this isn't necessary.
#'   Defaults to FALSE: labels aren't displayed.
#' @param top.n.hits How many of the top hits do you want to label? Defaults to
#'   FALSE: the number of labels aren't filtered.
#' @export
#' @return A list containing the data needed for the top plot, the bottom plot,
#'   the axes, and the maximum p-value (for sizing the plot)
#' @author Julie White
#' @references \url{https://github.com/juliedwhite/miamiplot}
#'
#' @importFrom magrittr %>%
#'

prep_miami_data <- function(
  data,
  split.by,
  split.at,
  chr="chr",
  pos="pos",
  p = "P",
  ...) {

  # Check the required input
  check_miami_input(data = data, split.by = split.by, split.at = split.at)

  # Identify column index for chromosome information and check class
  chr.indx <- find_col_info(data = data, col = chr)
  chr.name <- chr

  # Identify column index for position information and check class
  pos.indx <- find_col_info(data = data, col = pos)
  pos.name <- pos

  # Identify column index for p-value information and check class
  p.indx <- find_col_info(data = data, col = p)
  p.name <- p

  # To make the plot, we need to know the cumulative position of each probe/snp
  # across the whole genome.
  data <- data %>%
    #Group by chromosome
    dplyr::group_by(!!rlang::sym(chr.name)) %>%
    # Compute chromosome size
    dplyr::summarise(chrlength = max(!!rlang::sym(pos.name))) %>%
    # Calculate cumulative position of each chromosome
    dplyr::mutate(cumulativechrlength = cumsum(as.numeric(chrlength))-chrlength) %>%
    # Remove chr length column
    dplyr::select(-chrlength) %>%
    # Temporarily add the cumulative length of each chromosome to the initial
    # dataset
    dplyr::left_join(data, ., by=chr.name) %>%
    # Sort by chr then position
    dplyr::arrange(!!rlang::sym(chr.name), !!rlang::sym(pos.name)) %>%
    # Add the position to the cumulative chromosome length to get the position
    # of this probe relative to all other probes
    dplyr::mutate(rel.pos = !!rlang::sym(pos.name) + cumulativechrlength) %>%
    # Remove cumulative chr length column
    dplyr::select(-cumulativechrlength)

  # However, we don't want to print the cumulative position on the x-axis, so we
  # need to provide the chromosome position labels relative to the entire genome.
  axis.data = data %>%
    #Group by chromosome
    dplyr::group_by(!!rlang::sym(chr.name)) %>%
    #Find the center of the chromosome
    dplyr::summarize(chr.center=(max(rel.pos) + min(rel.pos))/2)

  # To create a symmetric looking plot, calculate the maximum p-value
  maxp <- ceiling(max(-log10(data[,p.indx])))

  # Depending on what the user has input for split.by and split.at, make top
  # and bottom data.
  if(is.numeric(split.at)){
    # Create top data using numeric info.
    top.data <- data %>%
      dplyr::filter(!!rlang::sym(split.by) >= split.at) %>%
      dplyr::mutate(loggedp = -log10(.[,p.indx])) %>%
      dplyr::rename(chr = as.name(chr.name)) %>%
      dplyr::select(chr, rel.pos, loggedp)
    # Create bottom data using numeric info
    bot.data <- data %>%
      dplyr::filter(!!rlang::sym(split.by) < split.at) %>%
      dplyr::mutate(loggedp = -log10(.[,p.indx])) %>%
      dplyr::rename(chr = as.name(chr.name)) %>%
      dplyr::select(chr, rel.pos, loggedp)

  } else if(is.character(split.at)){
    # Create top data using character specified by user
    top.data <- data %>%
      dplyr::filter(!!rlang::sym(split.by) == split.at) %>%
      dplyr::mutate(loggedp = -log10(.[,p.indx])) %>%
      dplyr::rename(chr = as.name(chr.name)) %>%
      dplyr::select(chr, rel.pos, loggedp)

    # Create bottom data using character not specified by user
    bot.data <- data %>%
      dplyr::filter(!!rlang::sym(split.by) != split.at) %>%
      dplyr::mutate(loggedp = -log10(.[,p.indx])) %>%
      dplyr::rename(chr = as.name(chr.name)) %>%
      dplyr::select(chr, rel.pos, loggedp)
  }

  miami_data_list <- list("top" = top.data, "bottom" = bot.data,
                          "axis" = axis.data, "maxp" = maxp)
  return(miami_data_list)
}
