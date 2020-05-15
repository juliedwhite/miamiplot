#' Check input to miami plot
#'
#' @param data A data.frame object. Required.
#' @param split.by A character vector. The name of the column to use for
#'   splitting into top and bottom sections of the Miami plot. Required.
#' @param split.at A character or numeric vector. If numeric, the top plot will
#'   contain your results where the values in the \code{split.by} column are
#'   >= \code{split.at}. If character, top plot will contain your results where
#'   the values in the \code{split.by} column are equal to \code{split.at}. Required.
#'
#' @return Will return an error message if there are problems.
#' @examples
#' # If you would like to put positive beta values in the top of your plot, and
#' # negative values on the bottom.
#' check_miami_input(data = df, split.by = "Beta", split.at = 0)
#'
#' # If you want genotyped SNPs on top, and imputed SNPs on the bottom plot.
#' check_miami_input(data = df, split.by = "GenoStatus", split.at = "Genotyped")
#'
#' # If you want to compare results from two studies, with study A on top.
#' check_miami_input(data = df, split.by = "Study", split.at = "A")
#'
#' @author Julie White
#' @references \url{https://github.com/juliedwhite/miamiplot}
#'

check_miami_input <- function(data, split.by, split.at) {

  # Establish a new argument check object so that we can nicely report errors.
  check <- checkmate::makeAssertCollection()

  # Check that the data is a data.frame, so that it plays nicely with the tidyverse.
  checkmate::assertDataFrame(data, col.names = "named", add = check)

  # The split.by column should exist in the data
  checkmate::assertNames(split.by, type = "named", subset.of = colnames(data),
                         add = check)

  # For now, only character and numeric options are supported, so that people
  # can split on things like genotyping status or numeric p-values.
  if(!checkmate::testCharacter(split.at, len = 1)){
    if(!checkmate::testNumeric(split.at, len = 1)){
      check$push("split.at must be either character or numeric, with length = 1")
    }
  }
  checkmate::reportAssertions(check)
}

