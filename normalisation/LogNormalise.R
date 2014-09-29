LogNormalise <- function(x) {
  # Log transform numerical data, ommitting NA values
  #   x: data to normalise.
  # Returns: Log transformed data.
  if (class(x[ ,1]) == "character") {
    rownames(x) <- x[, 1]
    x <- x[, -1]
  }
  x[x <= 0]    <- 10^(-10)           # manage negative and zero values
  x[!is.na(x)] <- log2(x[!is.na(x)]) # ommit NA values
  x[x < 0]     <- 0                  # zero negative values
  return(x)
}
