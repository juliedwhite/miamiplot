#' Make label data.frame for miami plot
#'
#' @param data A data.frame object. Required.
#' @param p The name of the column containing your logged p-value information.
#'   Defaults to "loggedp," assuming data preparation with `prep_miami_data`
#' @param hits.label Either the name of the column(s), max. 2, to use for
#'   automatically labeling n hits, determined using \code{top.n.hits}, or a
#'   character vector of probes/genes/SNPs to label. If supplying a list of
#'   genes, it is helpful to also supply \code{top.n.hits} to limit the labels
#'   to the top N results, though this isn't necessary.
#' @param top.n.hits How many of the top hits do you want to label? Defaults to
#'   5. Set to NULL to turn off this filtering.
#' @export
#' @return A data.frame with
#' @author Julie White
#' @references \url{https://github.com/juliedwhite/miamiplot}
#'
#' @importFrom magrittr %>%
#'


make_miami_labels <- function(data, p = "loggedp", hits.label, top.n.hits = 5){
  # Check that top.n.hits is a positive natural number and
    # if hits.label is a column in the df.
  if(checkmate::testCount(top.n.hits) &
     checkmate::testNames(hits.label, type = "named",
                          subset.of = colnames(data))){

    # If the user has given one column name in "hits.label"
    if(length(hits.label) == 1){
      # Filter the data to find top n values and label with user requested cols.
      label.df <- data %>%
        dplyr::mutate(label = !!rlang::sym(hits.label)) %>%
        dplyr::select(rel.pos, !!rlang::sym(p), label) %>%
        dplyr::arrange(desc(!!rlang::sym(p))) %>%
        dplyr::slice(1:top.n.hits)
      # Return label df
      return(label.df)

      # Or, if the user has given two column names in "hits.label"
    } else if(length(hits.label) == 2){
      # Filter the data to find top n values and label with user requested cols.
      label.df <- data %>%
        dplyr::mutate(label = paste0(!!rlang::sym(hits.label[1]), "\n",
                                     !!rlang::sym(hits.label[2]))) %>%
        dplyr::select(rel.pos, !!rlang::sym(p), label) %>%
        dplyr::arrange(desc(!!rlang::sym(p))) %>%
        dplyr::slice(1:top.n.hits)
      # Return label.df
      return(label.df)
    }

  # If instead the user has given a character vector of which SNPs/probes/genes
    # to label, this should not be a subset of colnames.
  } else if(!checkmate::testNames(hits.label, type = "named",
                                  subset.of = colnames(data))) {
    # If they don't want us to also filter by top n, then top.n.hits should be null
    if(is.null(top.n.hits)){
      hits.col.indx <- unname(which(sapply(data, function(x) any(x %in% hits.label))))
      # Filter the data to get the position of these labels
      label.df <- data %>%
        dplyr::filter(!!rlang::sym(colnames(data)[hits.col.indx]) %in% hits.label) %>%
        dplyr::mutate(label = !!rlang::sym(colnames(data)[hits.col.indx])) %>%
        dplyr::select(rel.pos, !!rlang::sym(p), label)
      # Return label.df
      return(label.df)

      # If they also give a positive natural number in top.n.hits, then we will
        # only find the top n values matching the hits.label.
    } else if(checkmate::testCount(top.n.hits)){
      hits.col.indx <- unname(which(sapply(data, function(x) any(x %in% hits.label))))
      # Filter the data to get the position of these labels
      label.df <- data %>%
        dplyr::filter(!!rlang::sym(colnames(data)[hits.col.indx]) %in% hits.label) %>%
        dplyr::mutate(label = !!rlang::sym(colnames(data)[hits.col.indx])) %>%
        dplyr::select(rel.pos, !!rlang::sym(p), label) %>%
        dplyr::arrange(desc(!!rlang::sym(p))) %>%
        slice(1:top.n.hits)
      # Return label.df
      return(label.df)
    }
  }
}

