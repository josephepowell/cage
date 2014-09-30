QuantileNormalise <- function(x) {
  # Perform quantile normalisation on numerical data using the Bioconductor
  # package `limma`.
  #   x: data to normalise.
  # Returns: matrix of numeric values, normalised across chips.
  x <- as.matrix(x)
  if (!all(x >= 0)) {
    stop("\n  x contains negative values. Perform log normalisation.")
  }
  x <- limma::normalizeBetweenArrays(x)
  return(x)
}
