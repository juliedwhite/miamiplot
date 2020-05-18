#' Make label data.frame for miami plot
#'
#' @param data A data.frame object. Required.
#' @param hits_label_col Either the name of the column(s), max = 2, to use for
#'   automatically labeling n hits, determined using \code{top_n_hits}, or the
#'   column where the values you provide in \code{hits_label} can be found.
#'   Required.
#' @param hits_label A user-specified character vector of probes/genes/SNPs
#'   to label. Defaults to NULL.
#' @param top_n_hits How many of the top hits do you want to label? Defaults to
#'   5. Set to NULL to turn off this filtering, but this is not recommended
#'   because the plot can get cluttered easily.
#' @param loggedp The name of the column containing your logged p-value
#'   information. Defaults to "loggedp," assuming data preparation with
#'   \code{prep_miami_data}
#' @param rel_pos The name of the column containing the positon of each
#'   SNP/probe relative to all other SNPs/probes in the genome. Defaults to
#'   "rel_pos," assuming data preparation with \code{prep_miami_data}.
#' @export
#' @return A data.frame with three columns: relative position, logged p-value,
#'   and the label to add to the miami plot. Designed to be used with ggrepel.
#' @author Julie White
#' @references \url{https://github.com/juliedwhite/miamiplot}
#'
#' @importFrom magrittr %>%
#'

make_miami_labels <- function(data, hits_label_col, hits_label = NULL,
                              top_n_hits = 5, loggedp = "loggedp",
                              rel_pos = "rel_pos") {
  # If the user has not requested specific labels
  if (is.null(hits_label)) {
    # Check that the columns given in hits_label_col are named in the dataframe.
    checkmate::assertNames(hits_label_col, type = "named",
                          subset.of = colnames(data))

    # If the user has given one column name in "hits_label"
    if (length(hits_label_col) == 1) {
      label_df <- data %>%
        dplyr::mutate(label = !!rlang::sym(hits_label_col)) %>%
        dplyr::select(!!rlang::sym(rel_pos), !!rlang::sym(loggedp), label)

      # Or, if the user has given two column names in "hits_label"
    } else if (length(hits_label_col) == 2) {
      label_df <- data %>%
        dplyr::mutate(label = paste0(!!rlang::sym(hits_label_col[1]), "\n",
                                     !!rlang::sym(hits_label_col[2]))) %>%
        dplyr::select(!!rlang::sym(rel_pos), !!rlang::sym(loggedp), label)
    }

  # If instead the user has given a character vector of which SNPs/probes/genes
  # to label
  } else if (!is.null(hits_label)) {
    # Check that the column given in hits_label_col is named in the dataframe
    checkmate::assertNames(hits_label_col, type = "named",
                           subset.of = colnames(data))

    # and has length = 1, because we can only match one column.
    if (length(hits_label_col) > 1) {
      stop("You have given two column names in addition to specifying ",
           "hits_label. This function cannot currently match items in more ",
           "one column, so please only provide one column where the items in ",
           "hits_label can be found.")
    }

    # Check if the column contains special characters, which makes it difficult
    # to identify and sort the labels.
    if (!rlang::is_empty(grep("[[:punct:]]", data[, hits_label_col]))) {
      warning("The column {", hits_label_col, "} has special characters in it.",
              " This function only respects exact matches, so your labels may ",
              "be incorrect. Please consider removing the special characters ",
              "from this column or providing exact matches.")
    }

    # Filter the data to get the position of these labels
    label_df <- data %>%
      dplyr::filter(!!rlang::sym(hits_label_col) %in% hits_label) %>%
      dplyr::mutate(label = !!rlang::sym(hits_label_col)) %>%
      dplyr::select(!!rlang::sym(rel_pos), !!rlang::sym(loggedp), label)

    # If label_df is empty, it means we did not find any matches so throw a
    # warning.
    if (nrow(label_df) == 0) {
      warning("I could not find any matches for {",
              paste(hits_label, collapse = ", "), "} in {", hits_label_col,
              "}. Perhaps there was a typo?")
    }
  }

  # If top_n_hits is a number, re-arrange the dataframe and select the top n
  if (checkmate::testCount(top_n_hits)) {
    label_df <- label_df %>%
      dplyr::arrange(desc(!!rlang::sym(loggedp))) %>%
      dplyr::slice(1:top_n_hits)

    # Return the df
    return(label_df)

  } else if (is.null(top_n_hits)) {
    # Top_n_hits is null, so just return full df
    return(label_df)
  }
}
