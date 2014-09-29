LogNormalise <- function(x) {
  # Log transform numerical data, ommitting NA values
  #   x: data to normalise.
  # Returns: Log transformed data.
  if (class(x[ ,1]) == "character") {
    x <- x[, -1]
  }
  x[x <= 0]    <- 10^(-10)
  x[!is.na(x)] <- log2(x[!is.na(x)])
  x[x < 0]     <- 0
  return(x)
}
