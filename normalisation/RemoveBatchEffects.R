#-------------------------------------------------------------------------------
# File  : RemoveBatchEffects.R
# Author: Alex Holloway
# Date  : 01/10/2014
#-------------------------------------------------------------------------------
RemoveBatchEffects <- function(exp, process.info, sample.info) {
  # Normalise gene expression data by accounting for various batch effects.
  #
  # Args:
  #   exp: expression data.
  #   process.info: microarray batch process data (i.e. process and extraction
  #                 dates).
  #   sample.info: individual sample data (i.e. Sentrix IDs and positions).
  #
  # Returns:
  #   
  # TODO: add arguments for individual vectors of data
  # TODO: check dimensions of all inputs
  ex.date <- as.factor(process.info$RNA_EXTRACT_DATE)
  sen.id  <- as.factor(sample.info$SENTRIX_ID)
  sen.pos <- as.factor(sample.info$SENTRIX_POSITION)

  exp.nrm <- matrix(nrow = nrow(exp),
                    ncol = ncol(exp),
                dimnames = dimnames(exp))

  for (i in 1:nrow(exp)) {
    non.na <- array(!is.na(exp[i, ]))
    if (length(which(non.na)) < 16) {  # arbitrary threshold for NAs, prevents
      exp.nrm[i, ] <- NA               # error during model construction
      next
    }

    mod.1  <- lm(as.numeric(exp[i, ]) ~ sen.id + sen.pos + ex.date)
    fit.1  <- summary(mod.1)
    exp.nrm[i, non.na] <- fit.1$residuals
  }

  exp.nrm <- data.frame(exp.nrm)
  return(exp.nrm)
}
