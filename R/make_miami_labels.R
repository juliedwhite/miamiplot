#' Make label data.frame for miami plot
#'
#' @param data A data.frame object. Required.
#' @param hits.label.col Either the name of the column(s), max. 2, to use for
#'   automatically labeling n hits, determined using \code{top.n.hits}, or the
#'   column where the values you provide in \code{hits.label} can be found. Required.
#' @param hits.label A user-specified character vector of probes/genes/SNPs
#'   to label. Defaults to NULL.
#' @param top.n.hits How many of the top hits do you want to label? Defaults to
#'   5. Set to NULL to turn off this filtering, but this is not recommended
#'   because the plot can get cluttered easily.
#' @param loggedp The name of the column containing your logged p-value information.
#'   Defaults to "loggedp," assuming data preparation with `prep_miami_data`
#' @param rel.pos The name of the column containing the positon of each
#'   SNP/probe relative to all other SNPs/probes in the genome. Defaults to
#'   "rel.pos," assuming data preparation with `prep_miami_data`
#' @export
#' @return A data.frame with three columns: relative position, logged p-value,
#'   and the label to add to the miami plot.
#' @author Julie White
#' @references \url{https://github.com/juliedwhite/miamiplot}
#'
#' @importFrom magrittr %>%
#'

make_miami_labels <- function(data, hits.label.col, hits.label = NULL,
                              top.n.hits = 5, loggedp = "loggedp",
                              rel.pos = "rel.pos"){
  # If the user has not requested specific labels
  if(is.null(hits.label)){
    # Check that the columns given in hits.label.col are named in the dataframe.
    checkmate::assertNames(hits.label.col, type = "named",
                          subset.of = colnames(data))

    # If the user has given one column name in "hits.label"
    if(length(hits.label.col) == 1){
      label.df <- data %>%
        dplyr::mutate(label = !!rlang::sym(hits.label.col)) %>%
        dplyr::select(!!rlang::sym(rel.pos), !!rlang::sym(loggedp), label)

      # Or, if the user has given two column names in "hits.label"
    } else if(length(hits.label.col) == 2){
      label.df <- data %>%
        dplyr::mutate(label = paste0(!!rlang::sym(hits.label.col[1]), "\n",
                                     !!rlang::sym(hits.label.col[2]))) %>%
        dplyr::select(!!rlang::sym(rel.pos), !!rlang::sym(loggedp), label)
    }

  # If instead the user has given a character vector of which SNPs/probes/genes
  # to label
  } else if(!is.null(hits.label)){
    # Check that the column given in hits.label.col is named in the dataframe
    checkmate::assertNames(hits.label.col, type = "named",
                           subset.of = colnames(data))

    # and has length = 1, because we can only match one column.
    if(length(hits.label.col) > 1){
      stop("You have given two column names in addition to specifying ",
           "hits.label. This function cannot currently match items in more ",
           "one column, so please only provide one column where the items in ",
           "hits.label can be found.")
    }

    # Check if the column contains special characters, which makes it difficult
    # to identify and sort the labels.
    if(!rlang::is_empty(grep('[[:punct:]]', data[,hits.label.col]))){
      warning("The column {", hits.label.col, "} has special characters in it.",
              " This function only respects exact matches, so your labels may ",
              "be incorrect. Please consider removing the special characters ",
              "from this column or providing exact matches.")
    }

    # Filter the data to get the position of these labels
    label.df <- data %>%
      dplyr::filter(!!rlang::sym(hits.label.col) %in% hits.label) %>%
      dplyr::mutate(label = !!rlang::sym(hits.label.col)) %>%
      dplyr::select(!!rlang::sym(rel.pos), !!rlang::sym(loggedp), label)

    # If label.df is empty, it means we did not find any matches so throw a
    # warning.
    if(nrow(label.df) == 0){
      warning("I could not find any matches for {",
              paste(hits.label, collapse = ", "), "} in {", hits.label.col,
              "}. Perhaps there was a typo?")
    }
  }

  # If top.n.hits is a number, re-arrange the dataframe and select the top n
  if(checkmate::testCount(top.n.hits)){
    label.df <- label.df %>%
      dplyr::arrange(desc(!!rlang::sym(loggedp))) %>%
      dplyr::slice(1:top.n.hits)

    # Return the df
    return(label.df)

  } else if(is.null(top.n.hits)) {
    # Top.n.hits is null, so just return full df
    return(label.df)
  }
}

