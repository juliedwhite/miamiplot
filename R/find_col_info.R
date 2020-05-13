#' Find column index matching a given string and verify numeric or integer class.
#'
#' @param data A data.frame object. Required.
#' @param col An object name or character vector. The name of the column to
#'   search the data.frame for and check class. Required.
#' @export
#' @return Will return an error message if there are problems. Otherwise returns
#'   the index of the named column in the dataframe.
#' @examples
#' If you are searching for the "chr" column.
#' find_col_info(data = df, col = "chr")
#'
#' @author Julie White
#' @references \url{https://github.com/juliedwhite/miamiplot}
#'

# Find column information and check class
find_col_info <- function(data, col) {

  # Identify column containing the x value supplied
  if (checkmate::testNames(col, type = "named", subset.of = colnames(data))) {
    # Column specified and is subset of colnames.
    # Find column index matching this value.
    col.indx <- which(colnames(data)==col)

    # Check that column is either numeric or integer
    if (!testMultiClass(data[,col.indx], classes = c("numeric", "integer"))) {
      stop(paste0("Your ", colnames(data)[col.indx], " column is of class ",
                  class(data[,col.indx]), ". Please make sure this column is ",
                  "numeric or integer."))
    }

    # Return column index
    return(col.indx)

  } else {
    # Column name is not in the right format or is not a subset of colnames.
    checkmate::assertNames(col, type = "named", subset.of = colnames(data))
  }
}




