                              ##################
                              # Clean CAD data #
                              ##################

#--------------------------------### Setup ###----------------------------------
inpath  <- "/ibscratch/wrayvisscher/xander/CAGE/data/CAD/raw"
outpath <- "/ibscratch/wrayvisscher/xander/CAGE/data/CAD/clean"
require(snpStats)
source("/clusterdata/uqahollo/scripts/Write.R")

PadWithZeroes <- function(id) {
  # Pad the numeric portion of sample IDs with zeroes, to five characters.
  #   id: sample ID number to pad
  # Returns: sample ID, with numeric portion padded to length 5.
  sep <- strsplit( gsub("([A-Z]{3})", "\\1~", id), "~")
  pre <- sep[[1]][1]
  suf <- sep[[1]][2]
  suf <- formatC(as.numeric(suf), width = 5, format = "d", flag = "0")
  id  <- paste0(pre, suf)
}

#---------------------------### Expression data ###-----------------------------
exp.raw.1 <- read.table(file = file.path(inpath, "Cardiology_phase1_raw.txt"),
                      header = TRUE,
                         sep = "\t")
# phase 2 expression data is labelled with IDs in the form GG.0### for ### in (372:534)
# these IDs match to the ColumnName in the covariate file
exp.raw.2 <- read.table(file = file.path(inpath, "Cardiology_phase2_raw.txt"),
                      header = TRUE,
                         sep = "\t")
# Extract probe info from phase 1 matrix
PROBE_ID.1 <- exp.raw.1[, 1]
exp.1      <- exp.raw.1[, order(colnames(exp.raw.1))]
info.probe <- exp.1[, 1:10]
exp.1      <- exp.1[, -c(1:11)] # strip probe info
# Extract STUDY_ID from phase 1 expression data
exp.id <- colnames(exp.1)
exp.id <- gsub("X[0-9]*_|\\.[A-z]*", "", exp.id)
exp.id <- sapply(exp.id, PadWithZeroes)
colnames(exp.1) <- exp.id
#---------------------------### Genotype data ###-------------------------------
gen.plink <- read.plink("/ibscratch/wrayvisscher/xander/CAGE/data/CAD/raw/genotypes/CAD")
gen <- t(gen.plink$genotypes@.Data)
# Extract STUDY_ID from genotype data
gen.id <- gsub("X[0-9]*_|_C", "", colnames(gen))
gen.id <- sapply(gen.id, PadWithZeroes)
colnames(gen) <- gen.id
# re-label FAM matrix
fam <- gen.plink$fam
fam <- fam[order(rownames(fam)), ]
rownames(fam) <- gen.id
fam$pedigree  <- gen.id
fam$member    <- gen.id
# remove SNPs that lack allele info
map <- gen.plink$map
map <- map[which(complete.cases(map[, 5:6])), ]
gen <- gen[rownames(map), ]
#------------------------------### Covariates ###-------------------------------
cov <- read.table(file = file.path(inpath, "Cardiology_exptdes_bothphases.txt"),
                   sep = "\t",
                header = TRUE)
PROBE_ID.2 <- exp.raw.2[, 1]
exp.2 <- exp.raw.2[, -1]
colnames(exp.2) <- gsub("\\.[A-Z]*_[A-z]*", "", colnames(exp.2))
colnames(exp.2) <- gsub("\\.", "_", colnames(exp.2))
colnames(exp.2) <- cov[which(cov$ColumnName %in% colnames(exp.2)),"STUDY_ID"]
cov <- cov[which(cov$STUDY_ID %in% c(colnames(exp.1), colnames(exp.2))), ]
# TODO: remove unnecessary columns
cov <- cov[, -c(1,3)]
# TODO: re-label columns for consistency across datasets
colnames(cov)[1] <- "SAMPLE_ID"
#------------------------------### Cleanup ###----------------------------------
names <- gen.id[which(gen.id %in% exp.id)]
exp.1 <- exp.1[, names]           # only keep samples found in exp and gen data
exp.1 <- exp.1[, order(colnames(exp.1))]
exp.1 <- cbind(PROBE_ID.1, exp.1) # append probe IDs
exp.2 <- exp.2[, order(colnames(exp.2))]
exp.2 <- cbind(PROBE_ID.2, exp.2)
colnames(exp.1)[1] <- "PROBE_ID"
colnames(exp.2)[1] <- "PROBE_ID"
gen   <- gen[, names]
gen   <- gen[, order(colnames(gen))]
gen   <- gen[order(rownames(gen)), ]
gen   <- new("SnpMatrix", t(gen)) # create SnpMatrix object for output
#----------------------------### Write to file ###------------------------------
Write(exp.1, "CAD_exp_phase1.txt")
Write(exp.2, "CAD_exp_phase2.txt")
Write(cov, "CAD_cov.txt")
Write(info.probe, "CAD_probe_info.txt")
write.plink(file.base = file.path(outpath, "genotypes/CAD_phase1"),
            snp.major = TRUE,
                 snps = gen,
         subject.data = fam,
             snp.data = map)
