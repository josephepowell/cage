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
  exp     <- as.matrix(exp)
  prc.inf <- as.matrix(process.info)
  smp.inf <- as.matrix(sample.info)
  
  ex.date <- as.factor(prc.inf$RNA_EXTRACT_DATE)
  sen.id  <- as.factor(smp.inf$SENTRIX_ID)
  sen.pos <- as.factor(smp.inf$SENTRIX_POSITION)

  mod <- lm(as.numeric(x[1, ]) ~ sen.id + sen.pos + ex.date)
  fit <- summary(mod)
  # TODO: extract summary stats
}
