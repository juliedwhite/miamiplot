# Checking functions

check.miami.input <- function(data, split.by, split.at) {
  require(checkmate)
  
  # Rename data
  df <- data
  
  # Establish a new argument check object
  check <- makeAssertCollection()
  
  # Check that df is a dataframe
  assertDataFrame(df, col.names = "named", add = check)
  
  # The user should supply the column name by which they want to split their data using split.by.
  # For example: beta values, study name, genotype status. 
  assertNames(split.by, type = "named", subset.of = colnames(df),  add = check)
  
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

