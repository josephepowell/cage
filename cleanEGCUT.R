                             ####################
                             # Clean EGCUT data #
                             ####################

#--------------------------------### Setup ###----------------------------------
inpath  <- "/ibscratch/wrayvisscher/xander/CAGE/data/EGCUT/raw"
outpath <- "/ibscratch/wrayvisscher/xander/CAGE/data/EGCUT/clean"
require(snpStats)
source("/clusterdata/uqahollo/scripts/Write.R")
#---------------------------### Expression data ###-----------------------------
exp.raw <- read.table(file = file.path(inpath, "ExpressionData.txt"),
                       sep = "\t",
                    header = TRUE,
               check.names = FALSE)
info    <- read.table(file = file.path(inpath, "EGC_Proteomics_Vcodes_Pcodes_180711.txt"),
                       sep = "\t",
                    header = TRUE)
PROBE_ID <- exp.raw[, 1]
exp   <- exp.raw[, -1]
# get Vcode from probe ID
names <- sapply(colnames(exp), function(x) info[which(info[, 1] == x),2] )
colnames(exp) <- names
exp   <- exp[, order(colnames(exp))]
exp   <- cbind(PROBE_ID, exp)
#---------------------------### Genotype data ###--------------------------------
plink.1 <- read.plink(file.path(inpath, "EGCUT_370CNV_ExpSample_290814"))
# Remove CNV genotypes that do not overlap with expression data
gen.1   <- plink.1$genotypes@.Data
gen.1   <- gen.1[which(rownames(gen.1) %in% names), ]
fam.1   <- plink.1$fam
fam.1   <- fam.1[rownames(gen.1), ]
# All Omni genotypes overlap with expression data
plink.2 <- read.plink(file.path(inpath, "EGCUT_OmniX_ExpSample_290814"))
gen.2   <- plink.2$genotypes@.Data
fam.2   <- plink.2$fam
# Remove expression data for samples with no genotypes
names <- unique( c(rownames(plink.1$genotypes), rownames(plink.2$genotypes)) )
names <- names[which(names %in% colnames(exp))]
exp   <- exp[, names]
exp   <- cbind(PROBE_ID, exp)
# TODO: Map probe IDs to ILMN_ID
ill.path <- "/ibscratch/wrayvisscher/xander/CAGE/annotation/humanht-12_v3_0_r3_11283641_a_txt.zip"
ill.anno <- read.delim(pipe(paste0("gunzip -c ", annotation)),
                       skip = 8, header = TRUE)
ill.ids  <- ill.anno[which(ill$Array_Address_Id %in% exp$PROBE_ID),"Probe_Id"]
exp$PROBE_ID <- ill.ids
#-------------------------------------------------------------------------------
# Perform outer join on genotype matrices, favouring Omni genotypes
# names.1 <- rownames(gen.1)
# names.2 <- rownames(gen.2)
# gen.1 <- apply(gen.1, 2, as.numeric)
# gen.2 <- apply(gen.2, 2, as.numeric)
# gen.1 <- data.frame(ID = names.1, gen.1)
# gen.2 <- data.frame(ID = names.2, gen.2)
# gen <- plyr::join(gen.2, gen.1, type = "full", by = "ID")
# rownames(gen) <- gen$ID
# Strip ID row and column
# gen <- gen[!rownames(gen) %in% "ID",!colnames(gen) %in% "ID"]
# gen <- gen[order(rownames(gen)), ] # Sorting is undesirable because it breaks
# gen <- gen[, order(colnames(gen))] # up the contiguous blocks of missing values
# TODO: Merge fam tables
# TODO: Merge map tables
#-------------------------------------------------------------------------------
# Create `SnpMatrix` objects
gen.1 <- new("SnpMatrix", gen.1)
gen.2 <- new("SnpMatrix", gen.2)
#----------------------------### Write to file ###------------------------------
Write(exp, "EGCUT_expression_signals.txt")
write.plink(file.base = file.path(outpath, "EGCUT_CNV"),
            snp.major = TRUE,
                 snps = gen.1,
         subject.data = fam.1,
             snp.data = plink.1$map)

write.plink(file.base = file.path(outpath, "EGCUT_Omni"),
            snp.major = TRUE,
                 snps = gen.2,
         subject.data = fam.2,
             snp.data = plink.2$map)
