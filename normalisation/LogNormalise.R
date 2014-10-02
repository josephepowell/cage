#-------------------------------------------------------------------------------
# File  : LogNormalise.R
# Author: Alex Holloway
# Date  : 29/09/2014
#-------------------------------------------------------------------------------
LogNormalise <- function(x, neg.rm = TRUE, row.names = TRUE) {
  # Log transform numerical data, ommitting NA values.
  #
  # Args:
  #   x: data to normalise.
  #   neg.rm: logical, should non-positive values be set to NA.
  #   row.names: logical, does x contain a column of character values.
  #
  # Returns:
  #   data.frame of log transformed numeric values.
  if (row.names) {
    rownames(x) <- x[, 1]  # move column of character entries to row.names
    x           <- x[, -1] # strip character column
  }

  if (neg.rm) {
    x[x <= 0] <- NA       # remove negative and zero values
  } else {
    x[x <= 0] <- 10^(-10) # replace negative and zero values with small value
  }

  x[!is.na(x)] <- log2(x[!is.na(x)]) # ommit NA values
  x[x < 0]     <- 0                  # zero negative values
  return(data.frame(x))
}
