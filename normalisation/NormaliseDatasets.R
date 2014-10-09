#-------------------------------------------------------------------------------
# File  : NormaliseDatasets.R
# Author: Alex Holloway
# Date  : 09/10/2014
# 
# Functions to perform end-to-end normalisation of each individual dataset in
# the CAGE consortium. Studies that contain data for multiple tissue types are 
# processed simultaneously, but returned as independent datasets.
#-------------------------------------------------------------------------------
stopifnot(require("ppca") & require("limma"))  # hard check for requisite packages
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
  #   data.frame of normalised gene expression values.
  exp <- exp.bsgs
  exp <- LogNormalise(exp)
  exp <- QuantileNormalise(exp)
  exp <- exp[!apply(exp, 1, function(x) all(is.na(x))), ]
  pc  <- ppca(t(exp), nPcs = 25)
  exp <- CorrectByPca(exp, pc)
  return(exp)
}

NormaliseBsgsPilot <- function() {
  # Convenience function for performing end-to-end normalisation of expression 
  # data from the main BSGS dataset.
  #
  # Args:
  #   null.
  #
  # Returns:
  #   list of data.frames, containing normalised gene expression values.
  exp.1 <- exp.bsgs.lcl
  exp.2 <- exp.bsgs.pbmc
  exp.1 <- LogNormalise(exp.1)
  exp.2 <- LogNormalise(exp.2)
  exp.1 <- QuantileNormalise(exp.1)
  exp.2 <- QuantileNormalise(exp.2)
  exp.1 <- exp.1[!apply(exp, 1, function(x) all(is.na(x))), ]
  exp.2 <- exp.2[!apply(exp, 1, function(x) all(is.na(x))), ]
  pc.1  <- ppca(t(exp.1), nPcs = 25)
  pc.2  <- ppca(t(exp.2), nPcs = 25)
  exp.1 <- CorrectByPca(exp.1, pc.1)
  exp.2 <- CorrectByPca(exp.2, pc.2)
  result <- list (exp.1, exp.2)
  return(result)
}

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
