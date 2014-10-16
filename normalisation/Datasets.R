#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------
stopifnot(require("limma") & require("pcaMethods"))
source("/clusterdata/uqahollo/scripts/cage/normalisation/LogTransform.R")
source("/clusterdata/uqahollo/scripts/cage/normalisation/QuantileNormalise.R")
source("/clusterdata/uqahollo/scripts/cage/normalisation/CorrectByPca.R")
source("/clusterdata/uqahollo/scripts/cage/normalisation/InverseNormalTransform.R")
source("/clusterdata/uqahollo/scripts/cage/normalisation/RemoveBatchEffects.R")
source("/clusterdata/uqahollo/scripts/cage/normalisation/Normalise.R")
inpath <- "/ibscratch/wrayvisscher/xander/CAGE/data/"
#-------------------------------------------------------------------------------
# BSGS main
#-------------------------------------------------------------------------------
exp.bsgs <- read.table(file.path(inpath, "BSGSmain/clean/BSGSmain_exp.txt"), header = TRUE, sep = "\t")
bsgs.prc <- read.table(file.path(inpath, "BSGSmain/clean/BSGSmain_process_info.txt"), header = TRUE, sep = "\t")
bsgs.smp <- read.table(file.path(inpath, "BSGSmain/clean/BSGSmain_cov.txt"), header = TRUE, sep = "\t")
exp.bsgs <- LogTransform(exp.bsgs, neg.rm = TRUE, row.names = TRUE)
exp.bsgs <- QuantileNormalise(exp.bsgs)
exp.bsgs <- RemoveBatchEffects(exp.bsgs, bsgs.prc, bsgs.smp)[[1]]
exp.bsgs <- t(apply(exp.bsgs, 1, InverseNormalTransform))
#-------------------------------------------------------------------------------
# BSGS pilot
#-------------------------------------------------------------------------------
exp.bsgs.lcl  <- read.table(file.path(inpath, "BSGSpilot/clean/BSGSpilot_exp_LCL.txt"), header = TRUE, sep = "\t")
exp.bsgs.lcl  <- LogTransform(exp.bsgs.lcl, neg.rm = TRUE, row.names = TRUE)
exp.bsgs.lcl  <- QuantileNormalise(exp.bsgs.lcl)
exp.bsgs.lcl  <- CorrectByPca(exp.bsgs.lcl, n.pcs = 25)
exp.bsgs.lcl  <- t(apply(exp.bsgs.lcl, 1, InverseNormalTransform))
#-------------------------------------------------------------------------------
exp.bsgs.pbmc <- read.table(file.path(inpath, "BSGSpilot/clean/BSGSpilot_exp_PBMC.txt"), header = TRUE, sep = "\t")
exp.bsgs.pbmc <- LogTransform(exp.bsgs.pbmc, neg.rm = TRUE, row.names = TRUE)
exp.bsgs.pbmc <- QuantileNormalise(exp.bsgs.pbmc)
exp.bsgs.pbmc <- CorrectByPca(exp.bsgs.pbmc, n.pcs = 25)
exp.bsgs.pbmc  <- t(apply(exp.bsgs.pbmc, 1, InverseNormalTransform))
#-------------------------------------------------------------------------------
# CAD
#-------------------------------------------------------------------------------
exp.cad.1 <- read.table(file.path(inpath, "CAD/clean/CAD_exp_phase1.txt"), header = TRUE, sep = "\t")
exp.cad.1 <- LogTransform(exp.cad.1, neg.rm = TRUE, row.names = TRUE)
exp.cad.1 <- QuantileNormalise(exp.cad.1)
exp.cad.1 <- CorrectByPca(exp.cad.1, n.pcs = 25)
exp.cad.1 <- t(apply(exp.cad.1, 1, InverseNormalTransform))
#-------------------------------------------------------------------------------
exp.cad.2 <- read.table(file.path(inpath, "CAD/clean/CAD_exp_phase2.txt"), header = TRUE, sep = "\t")
exp.cad.2 <- LogTransform(exp.cad.2, neg.rm = TRUE, row.names = TRUE)
exp.cad.2 <- QuantileNormalise(exp.cad.2)
exp.cad.2 <- CorrectByPca(exp.cad.2, n.pcs = 25)
exp.cad.2 <- t(apply(exp.cad.2, 1, InverseNormalTransform))
#-------------------------------------------------------------------------------
# CHDWB
#-------------------------------------------------------------------------------
exp.chdwb <- read.table(file.path(inpath, "CHDWB/clean/CHDWB_exp_log2.txt"), header = TRUE, sep = "\t")
exp.chdwb <- QuantileNormalise(exp.chdwb, row.names = TRUE)
exp.chdwb <- CorrectByPca(exp.chdwb, n.pcs = 25)
exp.chdwb <- t(apply(exp.chdwb, 1, InverseNormalTransform))
#-------------------------------------------------------------------------------
# EGCUT
#-------------------------------------------------------------------------------
exp.egcut <- read.table(file.path(inpath, "EGCUT/clean/EGCUT_exp.txt"), header = TRUE, sep = "\t")
exp.egcut <- LogTransform(exp.egcut, neg.rm = TRUE, row.names = TRUE)
exp.egcut <- QuantileNormalise(exp.egcut)
exp.egcut <- CorrectByPca(exp.egcut, n.pcs = 25)
exp.egcut <- t(apply(exp.egcut, 1, InverseNormalTransform))
#-------------------------------------------------------------------------------
# Morocco
#-------------------------------------------------------------------------------
exp.moroc <- read.table(file.path(inpath, "Morocco/clean/Morocco_exp.txt"), header = TRUE, sep = "\t")
exp.moroc <- LogTransform(exp.moroc, neg.rm = TRUE, row.names = TRUE)
exp.moroc <- QuantileNormalise(exp.moroc)
exp.moroc <- CorrectByPca(exp.moroc, n.pcs = 25)
exp.moroc <- t(apply(exp.moroc, 1, InverseNormalTransform))
#-------------------------------------------------------------------------------
# MuTHER
#-------------------------------------------------------------------------------
exp.mut.fat <- read.table(file.path(inpath, "MuTHER/clean/MuTHER_exp_fat.txt"), header = TRUE, sep = "\t")
exp.mut.fat <- LogTransform(exp.mut.fat, neg.rm = TRUE, row.names = TRUE)
exp.mut.fat <- QuantileNormalise(exp.mut.fat)
exp.mut.fat <- CorrectByPca(exp.mut.fat, n.pcs = 25)
exp.mut.fat <- t(apply(exp.mut.fat, 1, InverseNormalTransform))
#-------------------------------------------------------------------------------
exp.mut.skin <- read.table(file.path(inpath, "MuTHER/clean/MuTHER_exp_skin.txt"), header = TRUE, sep = "\t")
exp.mut.skin <- LogTransform(exp.mut.skin, neg.rm = TRUE, row.names = TRUE)
exp.mut.skin <- QuantileNormalise(exp.mut.skin)
exp.mut.skin <- CorrectByPca(exp.mut.skin, n.pcs = 25)
exp.mut.skin <- t(apply(exp.mut.skin, 1, InverseNormalTransform))
#-------------------------------------------------------------------------------
exp.mut.lcl <- read.table(file.path(inpath, "MuTHER/clean/MuTHER_exp_LCL.txt"), header = TRUE, sep = "\t")
exp.mut.lcl <- LogTransform(exp.mut.lcl, neg.rm = TRUE, row.names = TRUE)
exp.mut.lcl <- QuantileNormalise(exp.mut.lcl)
exp.mut.lcl <- CorrectByPca(exp.mut.lcl, n.pcs = 25)
exp.mut.lcl <- t(apply(exp.mut.lcl, 1, InverseNormalTransform))
