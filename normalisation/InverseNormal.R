#-------------------------------------------------------------------------------
# File  : InverseNormal.R
# Author: Alex Holloway
# Date  : 09/10/2014
#-------------------------------------------------------------------------------
InverseNormal <- function(x) {
  # Transform the values in a vector to an inverse-normal distribution.
  # i.e. scaled to have values between 0 and 1.
  #
  # Args:
  #   x: vector of numeric values.
  #
  # Returns:
  #   vector of numeric values, with range in 0 to 1.
  x <- (x - min(x)) / (max(x) - min(x))
  return(x)
}
