#' Prepare data for miami plot
#'
#' @param data A data.frame object. Required.
#' @param split_by A character vector. The name of the column to use for
#'   splitting into upper and lower sections of the Miami plot. Required.
#' @param split_at A character or numeric vector. If numeric, the upper plot
#'   will contain your results where the values in the \code{split_by} column
#'   are >= \code{split_at}. If character, upper plot will contain your results
#'   where the values in the \code{split_by} column are equal to
#'   \code{split_at}. Required.
#' @param chr The name of the column containing your chromosome information.
#'   Defaults to "chr"
#' @param pos The name of the column containing your position information.
#'   Defaults to "pos"
#' @param p The name of the column containing your p-value information.
#'   Defaults to "p"
#' @examples
#'   To create plot data where results are split with positive beta values in
#'   the upper plto and negative beta values in the lower plot:
#'     plot_data <- prep_miami_data(data = df, split_by = "Beta", split_at = 0)
#' @export
#' @return A list containing the data needed for the upper and lower plots,
#'   axes, and the maximum p-value (for sizing the plot)
#' @author Julie White
#' @references \url{https://github.com/juliedwhite/miamiplot}
#'
#' @importFrom magrittr %>%
#'

prep_miami_data <- function(
  data,
  split_by,
  split_at,
  chr = "chr",
  pos = "pos",
  p  =  "p") {

  # Check the required input
  check_miami_input(data = data, split_by = split_by, split_at = split_at,
                    chr = chr, pos = pos, p = p)

  # To make the plot, we need to know the cumulative position of each probe/snp
  # across the whole genome.
  data <- data %>%
    #Group by chromosome
    dplyr::group_by(!!rlang::sym(chr)) %>%
    # Compute chromosome size
    dplyr::summarise(chrlength = max(!!rlang::sym(pos))) %>%
    # Calculate cumulative position of each chromosome
    dplyr::mutate(cumulativechrlength = cumsum(as.numeric(chrlength)) -
                    chrlength) %>%
    # Remove chr length column
    dplyr::select(-chrlength) %>%
    # Temporarily add the cumulative length of each chromosome to the initial
    # dataset
    dplyr::left_join(data, ., by = chr) %>%
    # Sort by chr then position
    dplyr::arrange(!!rlang::sym(chr), !!rlang::sym(pos)) %>%
    # Add the position to the cumulative chromosome length to get the position
    # of this probe relative to all other probes
    dplyr::mutate(rel_pos = !!rlang::sym(pos) + cumulativechrlength) %>%
    # Remove cumulative chr length column
    dplyr::select(-cumulativechrlength)

  # However, we don't want to print the cumulative position on the x-axis, so we
  # need to provide the chromosome position labels relative to the entire
  # genome.
  axis_data <- data %>%
    #Group by chromosome
    dplyr::group_by(!!rlang::sym(chr)) %>%
    #Find the center of the chromosome
    dplyr::summarize(chr_center = (max(rel_pos) + min(rel_pos)) / 2)

  # To make it easier for ourselves later, create function-named chromosome and
  # logged p-value columns
  data <- data %>%
    dplyr::mutate(logged_p = -log10(!!rlang::sym(p))) %>%
    dplyr::rename(chr = as.name(chr))

  # To create a symmetric looking plot, calculate the maximum p-value
  maxp <- ceiling(max(data$logged_p, na.rm = TRUE))

  # Depending on what the user has input for split_by and split_at, make upper
  # and lower plot data.
  if (is.numeric(split_at)) {
    # Create upper plot data using numeric info.
    upper_data <- data %>%
      dplyr::filter(!!rlang::sym(split_by) >= split_at)

    # Create lower plot data using numeric info
    lower_data <- data %>%
      dplyr::filter(!!rlang::sym(split_by) < split_at)

  } else if (is.character(split_at)) {
    # Create upper plot data using character specified by user
    upper_data <- data %>%
      dplyr::filter(!!rlang::sym(split_by) == split_at)

    # Create lower plot data using whatever is left in that column
    lower_data <- data %>%
      dplyr::filter(!!rlang::sym(split_by) != split_at)

  }

  miami_data_list <- list("upper" = upper_data, "lower" = lower_data,
                          "axis" = axis_data, "maxp" = maxp)
  return(miami_data_list)
}
