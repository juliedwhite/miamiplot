#' Make label data.frame for miami plot
#'
#' @param data A data.frame object. Required.
#' @param loggedp The name of the column containing your logged p-value information.
#'   Defaults to "loggedp," assuming data preparation with `prep_miami_data`
#' @param hits.label Either the name of the column(s), max. 2, to use for
#'   automatically labeling n hits, determined using \code{top.n.hits}, or a
#'   character vector of probes/genes/SNPs to label. If supplying a list of
#'   genes, it is helpful to also supply \code{top.n.hits} to limit the labels
#'   to the top N results, though this isn't necessary.
#' @param top.n.hits How many of the top hits do you want to label? Defaults to
#'   5. Set to NULL to turn off this filtering, but this is not recommended.
#' @export
#' @return A data.frame with
#' @author Julie White
#' @references \url{https://github.com/juliedwhite/miamiplot}
#'
#' @importFrom magrittr %>%
#'

make_miami_labels <- function(data, loggedp = "loggedp", hits.label, top.n.hits = 5){
  # If hits.label is a column in the df.
  if(checkmate::testNames(hits.label, type = "named",
                          subset.of = colnames(data))){

    # If the user has given one column name in "hits.label"
    if(length(hits.label) == 1){
      label.df <- data %>%
        dplyr::mutate(label = !!rlang::sym(hits.label)) %>%
        dplyr::select(rel.pos, !!rlang::sym(p), label)

      # Or, if the user has given two column names in "hits.label"
    } else if(length(hits.label) == 2){
      label.df <- data %>%
        dplyr::mutate(label = paste0(!!rlang::sym(hits.label[1]), "\n",
                                     !!rlang::sym(hits.label[2]))) %>%
        dplyr::select(rel.pos, !!rlang::sym(loggedp), label)
    }

  # If instead the user has given a character vector of which SNPs/probes/genes
    # to label, this should not be a subset of colnames.
  } else if(!checkmate::testNames(hits.label, type = "named",
                                  subset.of = colnames(data))) {

    # Find the column where the labels they have given match the data.
    hits.col.indx <- unname(which(sapply(data, function(x) any(x %in% hits.label))))
    if(rlang::is_empty(hits.col.indx)){
      stop(paste0("I cannot find ", hits.label, " in your dataframe."))
    }

    # Filter the data to get the position of these labels
    label.df <- data %>%
      dplyr::filter(!!rlang::sym(colnames(data)[hits.col.indx]) %in% hits.label) %>%
      dplyr::mutate(label = !!rlang::sym(colnames(data)[hits.col.indx])) %>%
      dplyr::select(rel.pos, !!rlang::sym(loggedp), label)
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

