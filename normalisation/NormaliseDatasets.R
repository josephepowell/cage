#-------------------------------------------------------------------------------
# File  : NormaliseDatasets.R
# Author: Alex Holloway
# Date  : 09/10/2014
# 
# Functions to perform end-to-end normalisation of each individual dataset in
# the CAGE consortium. Studies that contain data for multiple tissue types are 
# processed simultaneously, but returned as independent datasets.
#-------------------------------------------------------------------------------
stopifnot(require("ppca"))  # hard check for requisite package
data("cage.expression")  # load CAGE gene expression (included in package)
#-------------------------------------------------------------------------------
NormaliseBsgsMain <- function() {
  # Convenience function for performing end-to-end normalisation of expression 
  # data from the main BSGS dataset.
  #
  # Args:
  #   null.
  #
  # Returns:
  #  data.frame of normalised gene expression values.
  exp <- exp.bsgs
  exp <- LogNormalise(exp)
  exp <- QuantileNormalise(exp)
  exp <- exp[!apply(exp, 1, function(x) all(is.na(x))), ]
  pc  <- ppca(t(exp), nPcs = 25)
  exp <- CorrectByPca(exp, pc)
  return(exp)
}

# TODO: NormaliseBsgsPilot

# TODO: NormaliseCad

NormaliseChdwb <- function() {
  # Convenience function for performing end-to-end normalisation of expression 
  # data from the CHDWB dataset.
  #
  # Args:
  #   null.
  #
  # Returns:
  #   data.frame of normalised gene expression values.
  exp <- exp.chdwb
  exp <- QuantileNormalise(exp)
  exp <- exp[!apply(exp, 1, function(x) all(is.na(x))), ]
  pc  <- ppca(t(exp), nPcs = 25)
  exp <- CorrectByPca(exp, pc)
  return(exp)
}

# TODO: NormaliseEgcut

# TODO: NormaliseMorocco

# TODO: NormaliseMuTHER
