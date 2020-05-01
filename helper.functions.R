# Helper functions

# Check that required args exist and are of appropriate data types.
check.miami.input <- function(data, split.by, split.at) {
  require(checkmate)

  # Establish a new argument check object
  check <- makeAssertCollection()
  
  # Check that data is a dataframe
  assertDataFrame(data, col.names = "named", add = check)
  
  # The user should supply the column name by which they want to split their data using split.by.
  # For example: beta values, study name, genotype status. 
  assertNames(split.by, type = "named", subset.of = colnames(data),  add = check)
  
  # The user should supply the point at which they want to split their data using split.at, as an integer or character.
  # For example: if they want positive beta values on top, then split.at = 0
  # For example, if they want "genotyped" SNPs on top, then split.at = "genotyped"
  # If the user has supplied a single character for split.at
  if(!testCharacter(split.at, len = 1)){
    if(!testNumeric(split.at, len = 1)){
      check$push("split.at must be either character or numeric, with length = 1")
    }
  }
  reportAssertions(check)
}


# Find column information 
find.col.indx <- function(data, col) {
  require(checkmate)
  
  # Identify column containing the x value supplied
  if(testNames(col, type = "named", subset.of = colnames(data))) {
    # Column specified and is subset of colnames. Find column index matching this value.
    col.indx <- which(colnames(data)==col)
    # Return column index
    return(col.indx)
  } else {
    # Column name is not in the right format or is not a subset of colnames. 
    assertNames(col, type = "named", subset.of = colnames(data))
  }
}


# Check that columns are either numeric or integer.
check.miami.numeric <- function(data, col.indx){
  require(checkmate)
  
  # Test if given column is numeric. 
  if(!testMultiClass(data[,col.indx], classes = c("numeric", "integer"))){
    stop(paste0("Your ", colnames(data)[col.indx], " column is of class ", 
                class(data[,col.indx]), ". Please make sure this column is numeric or integer."))
  }
}

