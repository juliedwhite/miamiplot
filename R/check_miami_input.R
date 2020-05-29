#' Check input to miami plot
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
#' @examples check_miami_input(data = gwas_results, split_by = "beta",
#'                             split_at = 0, p = "pval")
#' @export
#' @return If no problems, executes silently. Will return an error message if
#'   there are problems.
#' @author Julie White
#' @references \url{https://github.com/juliedwhite/miamiplot}
#'

check_miami_input <- function(data, split_by, split_at, chr = "chr",
                              pos = "pos", p = "p") {

  # Establish a new argument check object so that we can nicely report errors.
  check <- checkmate::makeAssertCollection()

  # Check that the data is a data.frame, so that it plays nicely with the
  # tidyverse.
  checkmate::assertDataFrame(data, col.names = "named", add = check)

  # The split.by column should exist in the data
  checkmate::assertNames(split_by, type = "named", subset.of = colnames(data),
                         add = check)

  # For now, only character and numeric options are supported, so that people
  # can split on things like genotyping status or numeric p-values.
  if (!checkmate::testCharacter(split_at, len = 1)) {
    if (!checkmate::testNumeric(split_at, len = 1)) {
      check$push("split_at must be either character or numeric,
                 with length = 1")
    }
  }

  # Check the chromosome, position, and p-value column.
  for (i in c(chr, pos, p)) {
    # Column should exist in the data
    if (!checkmate::testNames(i, type = "named", subset.of = colnames(data))) {
      stop(paste0("The column ", i, " must be a subset of set {",
                  paste(colnames(data), collapse = ","), "}."))
    }

    if (!checkmate::testMultiClass(data[, i],
                                   classes = c("numeric", "integer"))) {
      stop(paste0("Your ", i, " column is of class ", class(data[, i]),
                  ". Please make sure this column is numeric or integer"))
    }
  }

  checkmate::reportAssertions(check)
}
