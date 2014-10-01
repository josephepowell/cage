QuantileNormalise <- function(x) {
  # Perform quantile normalisation on numerical data using the Bioconductor
  # package `limma`.
  #
  # Args:
  #   x: data to normalise.
  #
  # Returns:
  #   data.frame of numeric values, normalised across chips.
  stopifnot(require("limma")) # hard check for requisite package
  x <- as.matrix(x)           # convert data to matrix
  if (!all(x[which(!is.na(x))] >= 0)) {
    stop("\n  x contains negative values. Perform log normalisation.")
  }
  x <- normalizeBetweenArrays(x)
  return(data.frame(x))
}
