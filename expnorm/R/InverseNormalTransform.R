#-------------------------------------------------------------------------------
# File  : InverseNormal.R
# Author: Alex Holloway
# Date  : 09/10/2014
#-------------------------------------------------------------------------------
InverseNormalTransform <- function(x) {
  # Perform rank-based inverse normal transformation on a vector of values.
  #
  # Args:
  #   x: vector of numeric values.
  #
  # Returns:
  #   vector of numeric values, transformed according to Blom (1958).
  c  <- 3 / 8  # arbitrary constant
  i  <- !is.na(x)
  n  <- length(x)
  zn <- (x[i] - mean(x[i])) / sd(x[i])
  if (all(zn[!is.na(zn)] == 0)) {  # check for vectors containing only 0 or NA
    x[i] <- zn
    x[is.nan(x)] <- 0
    return(x)
  } 
  r  <- rank(zn[!is.na(zn)], ties.method = "average")
  nrm  <- qnorm((r - c) / (n - 2 * c + 1))
  x[i] <- nrm
  return(x)
}
