setwd("/ibscratch/wrayvisscher/xander/CAGE/data")
source("/clusterdata/uqahollo/scripts/cage/normalisation/LogNormalise.R")
source("/clusterdata/uqahollo/scripts/cage/normalisation/QuantileNormalise.R")
library("pcaMethods")

exp.bsgs <- read.table("BSGSmain/clean/BSGSmain_exp.txt", header = TRUE, sep = "\t")
exp.bsgs <- QuantileNormalise(LogNormalise(exp.bsgs))
exp.bsgs <- exp.bsgs[!apply(exp.bsgs, 1, function(x) all(is.na(x))), ]
pc.bsgs  <- ppca(t(exp.bsgs), nPcs = 25)

exp.bsgs.lcl <- read.table("BSGSpilot/clean/BSGSpilot_exp_LCL.txt", header = TRUE, sep = "\t")
exp.bsgs.lcl <- QuantileNormalise(LogNormalise(exp.bsgs.lcl))
exp.bsgs.lcl <- exp.bsgs.lcl[!apply(exp.bsgs.lcl, 1, function(x) all(is.na(x))), ]
pc.bsgs.lcl  <- ppca(t(exp.bsgs.lcl), nPcs = 25)

exp.bsgs.pbmc <- read.table("BSGSpilot/clean/BSGSpilot_exp_PBMC.txt", header = TRUE, sep = "\t")
exp.bsgs.pbmc <- QuantileNormalise(LogNormalise(exp.bsgs.pbmc))
exp.bsgs.pbmc <- exp.bsgs.pbmc[!apply(exp.bsgs.pbmc, 1, function(x) all(is.na(x))), ]
pc.bsgs.pbmc  <- ppca(t(exp.bsgs.pbmc), nPcs = 25)

exp.egcut <- read.table("EGCUT/clean/EGCUT_exp.txt", header = TRUE, sep = "\t")
exp.egcut <- QuantileNormalise(LogNormalise(exp.egcut))
exp.egcut <- exp.egcut[!apply(exp.egcut, 1, function(x) all(is.na(x))), ]
pc.egcut  <- ppca(t(exp.egcut), nPcs = 25)

exp.chdwb <- read.table("CHDWB/clean/CHDWB_exp-log2.txt", header = TRUE, sep = "\t")
exp.chdwb <- QuantileNormalise(exp.chdwb, row.names = TRUE)  # already log2 normalised
exp.chdwb <- exp.chdwb[!apply(exp.chdwb, 1, function(x) all(is.na(x))), ]
pc.chdwb  <- ppca(t(exp.chdwb), nPcs = 25)

exp.cad.1 <- read.table("CAD/clean/CAD_exp_phase1.txt", header = TRUE, sep = "\t")
exp.cad.1 <- QuantileNormalise(LogNormalise(exp.cad.1))
exp.cad.1 <- exp.cad.1[!apply(exp.cad.1, 1, function(x) all(is.na(x))), ]
pc.cad.1  <- ppca(t(exp.cad.1), nPcs = 25)

exp.cad.2 <- read.table("CAD/clean/CAD_exp_phase2.txt", header = TRUE, sep = "\t")
exp.cad.2 <- QuantileNormalise(LogNormalise(exp.cad.2))
exp.cad.2 <- exp.cad.2[!apply(exp.cad.2, 1, function(x) all(is.na(x))), ]
pc.cad.2  <- ppca(t(exp.cad.2), nPcs = 25)

exp.moroc <- read.table("Morocco/clean/Morocco_exp.txt", header = TRUE, sep = "\t")
exp.moroc <- QuantileNormalise(LogNormalise(exp.moroc))
exp.moroc <- exp.moroc[!apply(exp.moroc, 1, function(x) all(is.na(x))), ]
pc.moroc  <- ppca(t(exp.moroc), nPcs = 25)

exp.mut.fat <- read.table("MuTHER/clean/MuTHER_exp_fat.txt", header = TRUE, sep = "\t")
exp.mut.fat <- QuantileNormalise(LogNormalise(exp.mut.fat))
exp.mut.fat <- exp.mut.fat[!apply(exp.mut.fat, 1, function(x) all(is.na(x))), ]
pc.mut.fat  <- ppca(t(exp.mut.fat), nPcs = 25)

exp.mut.lcl <- read.table("MuTHER/clean/MuTHER_exp_LCL.txt", header = TRUE, sep = "\t")
exp.mut.lcl <- QuantileNormalise(LogNormalise(exp.mut.lcl))
exp.mut.lcl <- exp.mut.lcl[!apply(exp.mut.lcl, 1, function(x) all(is.na(x))), ]
pc.mut.lcl  <- ppca(t(exp.mut.lcl), nPcs = 25)

exp.mut.skin <- read.table("MuTHER/clean/MuTHER_exp_skin.txt", header = TRUE, sep = "\t")
exp.mut.skin <- QuantileNormalise(LogNormalise(exp.mut.skin))
exp.mut.skin <- exp.mut.skin[!apply(exp.mut.skin, 1, function(x) all(is.na(x))), ]
pc.mut.skin  <- ppca(t(exp.mut.skin), nPcs = 25)
