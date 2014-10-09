#-------------------------------------------------------------------------------
# File  : Normalise.R
# Author: Alex Holloway
# Date  : 09/10/2014
#-------------------------------------------------------------------------------
Normalise <- function(x,
                      log = TRUE,
                      quantile = TRUE,
                      pca = TRUE,
                      inv.norm = TRUE,
                      n.pcs = 25,
                      neg.rm = TRUE,
                      row.names = FALSE) {
  # Perform end-to-end normalisation of gene expression data, through log2 
  # transformation, quantile normalisation, inverse normal transformation, 
  # and principal components analysis. By default, all normalisations are 
  # performed.
  #
  # Args:
  #   x: matrix or data.frame of expression data to be normalised.
  #   log: logical, should log2 transformation be performed.
  #   quantile: logical, should quantile normalisation be performed.
  #   inv.norm: logical, should the data be transformed have range in 0 to 1.
  #   pca: logical, should batch effects be removed by principal components analysis.
  #   n.pcs: numeric, quantity of principal components to fit.
  #   neg.rm: logical, should negative values be removed during log2 transformation.
  #   row.names: logical, does the matrix contain a character column of probe IDs.
  #
  #  Returns:
  #    data.frame of gene expression data, normalised according to selected arguments.
  if (log) {
    x <- LogNormalise(x, neg.rm, row.names)
  }
  if (quantile) {
    x <- QuantileNormalise(x, row.names)
  }
  if (pca) {
    x  <- x[apply(x, 1, function(x) all(is.na(x))), ]
    pc <- ppca(t(x), n.pcs)
    x  <- CorrectByPca(x, pc)
  }
  if (inv.norm) {
    x <- apply(x, 1, InverseNormal)
  }
  return(x)
}
