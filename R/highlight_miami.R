#' Make highlight data.frame for miami plot
#'
#' @param data A data.frame object. Required.
#' @param highlight A vector of SNPs or gene names you would like to highlight.
#'   Defaults to NULL. If you specify this, you must also specify highlight_col
#'   so that we know where to find the items.
#' @param highlight_col The column where the values you provide in
#'   \code{highlight} can be found.
#' @param highlight_color The color that you would like the highlighted SNPs or
#'   genes to be. Defaults to "green."
#' @param logged_p The name of the column containing your logged p-value
#'   information. Defaults to "logged_p," assuming data preparation with
#'   \code{prep_miami_data}
#' @param rel_pos The name of the column containing the positon of each
#'   SNP/probe relative to all other SNPs/probes in the genome. Defaults to
#'   "rel_pos," assuming data preparation with \code{prep_miami_data}.
#' @examples
#'   upper_highlight <- highlight_miami(data = plot_data$upper,
#'                                      highlight = "snps_of_interest"
#'                                      highlight_col = "rsid")
#'
#'  When highlighting, the items in 'highlight' must all come from a single
#'  column. Specifying multiple columns in highlight_col will return an error.
#'
#'  \dontrun{
#'  highlight_miami(data = df, split_by = "study", split_at = "A",
#'                  highlight = "snps_of_interest",
#'                  highlight_col = c("rsid", "gene"))
#'  }
#'
#' @export
#' @return A data.frame with three columns: relative position, logged p-value,
#'   and highlight color.
#' @author Julie White
#' @references \url{https://github.com/juliedwhite/miamiplot}
#'
#' @importFrom magrittr %>%
#'

highlight_miami <- function(data, highlight, highlight_col,
                            highlight_color = "green", logged_p = "logged_p",
                            rel_pos = "rel_pos") {

  # Check that the column given in highlight_col is named in the dataframe
  checkmate::assertNames(highlight_col, type = "named",
                         subset.of = colnames(data))

  # and has length = 1, because we can only match one column.
  if (length(highlight_col) > 1) {
    stop("You have given two column names in highlight_col. This function
         cannot currently match items in more one column, so please only
         provide one column where the items in highlight can be found.")
  }

  # Check if the column contains special characters, which makes it difficult
  # to identify and sort the labels.
  if (!rlang::is_empty(grep("[[:punct:]]", data[, highlight_col]))) {
    warning("The column {", highlight_col, "} has special characters in it.
            This function only respects exact matches, so your highlights
            may be incorrect. Please consider removing the special
            characters from this column or providing exact matches.")
  }

  # Filter the data to get the position of these labels
  highlight_df <- data %>%
    dplyr::filter(!!rlang::sym(highlight_col) %in% highlight) %>%
    dplyr::mutate(color = highlight_color) %>%
    dplyr::rename(rel_pos = !!rlang::sym(rel_pos)) %>%
    dplyr::rename(logged_p = !!rlang::sym(logged_p)) %>%
    dplyr::select(rel_pos, logged_p, color)

  # If label_df is empty, it means we did not find any matches so throw a
  # warning.
  if (nrow(highlight_df) == 0) {
    warning("I could not find any matches for {",
            head(paste(highlight, collapse = ", ")), "} in {", highlight_col,
            "}. Perhaps there was a typo?")
  }

  return(highlight_df)
}

