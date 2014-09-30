QuantileNormalise <- function(x) {
  # Perform quantile normalisation on numerical data using the Bioconductor
  # package `limma`.
  #   x: data to normalise.
  # Returns: matrix of numeric values, normalised across chips.
  stopifnot(all(x >= 0))
  # TODO: exit status for data that has not been log normalised
  x <- as.matrix(x)
  x <- normalizeBetweenArrays(x))
}
