                           #########################
                           # Clean BSGS pilot data #
                           #########################

#--------------------------------### Setup ###----------------------------------
require(snpStats)
inpath  <- "/ibscratch/wrayvisscher/xander/CAGE/data/BSGSpilot/raw"
outpath <- "/ibscratch/wrayvisscher/xander/CAGE/data/BSGSpilot/clean"
source("/clusterdata/uqahollo/scripts/Write.R")
#---------------------------### Expression data ###-----------------------------
exp.1   <- read.csv(file.path(inpath, "Probe_signal_unnormalised_LCL_only.csv"), header = TRUE)
exp.2   <- read.csv(file.path(inpath, "Probe_signal_unnormalised_PBMC_only.csv"), header = TRUE)
prb.inf <- read.csv(file.path(inpath, "Probe_information.csv"), header = TRUE)
# Sample information / covariates
cov.1 <- read.csv(file.path(inpath, "Sample_info_LCL_only.csv"), header = TRUE)
cov.2 <- read.csv(file.path(inpath, "Sample_info_PBMC_only.csv"), header = TRUE)
# Get project-specific sample IDs
colnames(exp.1) <- cov.1[, 2]
colnames(exp.2) <- cov.2[, 2]
# Sort by sample IDs
exp.1 <- exp.1[, order(colnames(exp.1))]
exp.2 <- exp.2[, order(colnames(exp.2))]
# Extract Illumina probe IDs
PROBE_ID <- prb.inf[, 1]
exp.1    <- cbind(PROBE_ID, exp.1)
exp.2    <- cbind(PROBE_ID, exp.2)

# Remove redundant columns
cov.1 <- cov.1[, -c(1,3)]
cov.2 <- cov.2[, -c(1,3)]
# Alter column names for consistency across files
names    <- sapply(colnames(cov.1), toupper)
names[1] <- "SAMPLE_ID"
colnames(cov.1) <- names
colnames(cov.2) <- names
# Re-order rows by sample ID
cov.1 <- cov.1[order(cov.1[, 1]), ]
cov.2 <- cov.2[order(cov.2[, 1]), ]
#----------------------------### Genotype data ###------------------------------
labels  <- c("AA", "AB", "BB")
alleles <- c("A", "C", "G", "T")
#---------------------------------## LCL ##-------------------------------------
# Read raw genotype data
plink.1 <- read.plink(file.path(inpath, "cleaned_geno_lcl.bed"))
gen.1 <- t(plink.1$genotypes@.Data)
# Recode SNPs as AA, AB or BB 
gen.1 <- apply(gen.1, 2, as.numeric)
gen.1[gen.1 == 0] <- NA
RS_ID <- colnames(plink.1$genotypes)
SAMPLE_ID <- plink.1$fam[, 2]
# 
map.1 <- plink.1$map
rownames(map.1) <- c(1:nrow(map.1))
map.1 <- cbind(RS_ID, map.1[, -c(2, 3)])
# Resolve allele IDs to nucleotide bases
alleles    <- c("A", "C", "G", "T")
map.1[, 4] <- alleles[map.1[, 4]]
map.1[, 5] <- alleles[map.1[, 5]]
# Enumerate actual allele combinations
alleles <- matrix(nrow = nrow(map.1),
                  ncol = 3)
for (i in 1:nrow(map.1)) {
  alleles[i, ] <- c(paste0(map.1[i, 4], map.1[i, 4]),
                    paste0(map.1[i, 4], map.1[i, 5]),
                    paste0(map.1[i, 5], map.1[i, 5]))
}
# Encode genotypes using allele combinations
gen.res <- matrix(nrow = nrow(gen.1),
                  ncol = ncol(gen.1))
for (i in 1:nrow(gen.1)) {
  gen.res[i, ] <- alleles[i, gen.1[i, ]]
}
# Set genotype row and column names
colnames(gen.res) <- SAMPLE_ID
gen.res <- cbind(RS_ID, gen.res)
gen.1   <- gen.res
# Correct column labelling for consistency between datasets
names   <- sapply(colnames(map.1), toupper)
names   <- gsub("\\.", "_", names)
colnames(map.1) <- names
#---------------------------------## PBMC ##------------------------------------
# Read raw genotype data
plink.2 <- read.plink(file.path(inpath, "cleaned_geno_pbmc.bed"))
gen.2 <- t(plink.2$genotypes@.Data)
# Recode SNPs as 0, 1, or 2 for AA, AB, or BB, respectively
gen.2 <- apply(gen.2, 2, as.numeric)
gen.2[gen.2 == 0] <- NA
RS_ID <- colnames(plink.2$genotypes)
SAMPLE_ID  <- plink.2$fam[, 2]
# Metadata
map.2 <- plink.2$map
rownames(map.2) <- c(1:nrow(map.2))
map.2 <- cbind(RS_ID, map.2[, -c(2, 3)])
# Resolve allele IDs to nucleotide bases
alleles <- c("A", "C", "G", "T")
map.2[, 4] <- alleles[map.2[, 4]]
map.2[, 5] <- alleles[map.2[, 5]]
# Resolve allele labels to nucleotide bases
# Enumerate actual allele combinations
alleles <- matrix(nrow = nrow(map.2),
                  ncol = 3)
for (i in 1:nrow(map.2)) {
  alleles[i, ] <- c(paste0(map.2[i, 4], map.2[i, 4]),
                    paste0(map.2[i, 4], map.2[i, 5]),
                    paste0(map.2[i, 5], map.2[i, 5]))
}
# Encode genotypes using allele combinations
gen.res <- matrix(nrow = nrow(gen.2),
                  ncol = ncol(gen.2))
for (i in 1:nrow(gen.2)) {
  gen.res[i, ] <- alleles[i, gen.2[i, ]]
}
# Set genotype row and column names
colnames(gen.res) <- SAMPLE_ID
gen.res <- cbind(RS_ID, gen.res)
gen.2   <- gen.res
# Correct column labelling for consistency between datasets
colnames(map.2) <- names
#---------------------### Make data entries consistent ###----------------------
#---------------------------------## LCL ##-------------------------------------
cols  <- which(!colnames(exp.1) %in% colnames(gen.1))
cols  <- cols[-1]
exp.1 <- exp.1[, -cols]
cols  <- cols - 1
cov.1 <- cov.1[-cols, ]
#---------------------------------## PBMC ##------------------------------------
cols  <- which(!colnames(exp.2) %in% colnames(gen.2))
cols  <- cols[-1]
exp.2 <- exp.2[, -cols]
cols  <- cols - 1
cov.2 <- cov.2[-cols, ]
#----------------------------### Write to file ###------------------------------
Write(exp.1, "BSGSpilot_exp_LCL.txt")
Write(exp.2, "BSGSpilot_exp_PBMC.txt")
Write(prb.inf, "BSGSpilot_probe_info.txt")
Write(cov.1, "BSGSpilot_sample_info_LCL.txt")
Write(cov.2, "BSGSpilot_sample_info_PBMC.txt")
Write(gen.1, "genotypes/BSGSpilot_gen_LCL.txt")
Write(gen.2, "genotypes/BSGSpilot_gen_PBMC.txt")
Write(map.1, "genotypes/BSGSpilot_gen_map_LCL.txt")
Write(map.2, "genotypes/BSGSpilot_gen_map_PBMC.txt")
# TODO: copy and rename PLINK files in clean directory
#system(paste0("bash for file in ", inpath, "/cleaned_geno_lcl.*
#do
#  cp '$file' 'test/BSGSpilot_LCL.${file#*.}'\n
#done"))
#------------------------------### Clean up ###---------------------------------
rm(inpath, outpath, PROBE_ID, names, plink.1, plink.2, gen.res, alleles,
   labels, RS_ID, SAMPLE_ID, cols, i, Write)
