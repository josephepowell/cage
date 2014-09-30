QuantileNormalise <- function(x) {
  # Perform quantile normalisation on numerical data using the Bioconductor
  # package `limma`.
  #   x: data to normalise.
  # Returns: data.frame of numeric values, normalised across chips.
  stopifnot(require("limma")) # hard check for requisite package
  x <- as.matrix(x)           # convert data to numerical matrix
  if (!all(x[which(!is.na(x))] >= 0)) {
    stop("\n  x contains negative values. Perform log normalisation.")
  }
  x <- data.frame(normalizeBetweenArrays(x))
  return(x)
}
