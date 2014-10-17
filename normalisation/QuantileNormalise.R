#-------------------------------------------------------------------------------
# File  : QuantileNormalise.R
# Author: Alex Holloway
# Date  : 30/09/2014
#-------------------------------------------------------------------------------
QuantileNormalise <- function(x, row.names = FALSE) {
  # Perform quantile normalisation on numerical data using the Bioconductor
  # package `limma`.
  #
  # Args:
  #   x: data to normalise.
  #   row.names: logical, does x contain a column of character values.
  #
  # Returns:
  #   data.frame of numeric values, normalised across chips.
  stopifnot(require("limma"))  # hard check for requisite package

  if (row.names) {
    rownames(x) <- x[, 1]  # move column of character entries to row.names
    x           <- x[, -1]  # strip character column
  }

  x <- as.matrix(x)  # convert data to matrix
  if (!all(x[which(!is.na(x))] >= 0)) {
    stop("\n  x contains negative values. Perform log normalisation.")
  }

  x <- normalizeBetweenArrays(x)
  return(data.frame(x))
}
