QuantileNormalise <- function(x) {
  # Perform quantile normalisation on numerical data using the Bioconductor
  # package `limma`.
  #   x: data to normalise.
  # Returns: data.frame of numeric values, normalised across chips.
  stopifnot(require("limma"))
  x <- as.matrix(x)
  if (!all(x >= 0)) {
    stop("\n  x contains negative values. Perform log normalisation.")
  }
  x <- normalizeBetweenArrays(x)
  x <- data.frame(x)
  return(x)
}
