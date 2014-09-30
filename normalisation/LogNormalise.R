LogNormalise <- function(x, neg.rm = TRUE) {
  # Log transform numerical data, ommitting NA values
  #   x: data to normalise.
  #   neg.rm: logical, should non-positive values be set to NA.
  # Returns: Log transformed data.
  if (class(x[, 1]) == "character") {
    rownames(x) <- x[, 1]
    x <- x[, -1]
  }
  if (neg.rm) {
    x[x <= 0] <- NA       # remove negative and zero values
  } else {
    x[x <= 0] <- 10^(-10) # replace negative and zero values with small value
  }
  x[!is.na(x)] <- log2(x[!is.na(x)]) # ommit NA values
  x[x < 0]     <- 0                  # zero negative values
  return(x)
}
