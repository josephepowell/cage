#-------------------------------------------------------------------------------
# File  : CorrectByPca.R
# Author: Alex Holloway
# Date  : 08/10/2014
#-------------------------------------------------------------------------------
CorrectByPca <- function(exp, n.pcs) {
  # Attempt to remove batch processing effects using principal components.
  #
  # Args:
  #   exp: expression data.
  #   n.pcs: number of principal components to perform analysis with.
  #
  # Returns:
  #   data.frame of expression levels, adjusted for principal component scores.
  stopifnot(require("pcaMethods"))  # hard check for requisite package
  exp.nrm <- matrix(nrow = nrow(exp),
                    ncol = ncol(exp),
                dimnames = dimnames(exp))
  pca <- ppca(t(exp), n.pcs)
  for (i in 1:nrow(exp)) {
    non.na <- array(!is.na(exp[i, ]))
    if (length(which(non.na)) < 16) {  # arbitrary threshold for NAs, prevents
      exp.nrm[i, ] <- NA               # error during model construction
      next
    }
    mod <- lm(as.numeric(exp[i, ]) ~ scores(pca)[, 1:24])
    fit <- summary(mod)
    exp.nrm[i, non.na] <- fit$residuals
  }
  exp.nrm <- data.frame(exp.nrm)
}
