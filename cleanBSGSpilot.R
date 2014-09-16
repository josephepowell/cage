                           #########################
                           # Clean BSGS pilot data #
                           #########################

#--------------------------------### Setup ###----------------------------------
require(snpStats)
inpath  <- "/ibscratch/wrayvisscher/xander/CAGE/data/BSGSpilot/raw"
outpath <- "/ibscratch/wrayvisscher/xander/CAGE/data/BSGSpilot/clean"
source("/clusterdata/uqahollo/scripts/Write.R")
#---------------------------### Expression data ###-----------------------------
# Read raw expression data
exp.1    <- read.csv(file.path(inpath, "Probe_signal_unnormalised_LCL_only.csv"), header = TRUE)
exp.2    <- read.csv(file.path(inpath, "Probe_signal_unnormalised_PBMC_only.csv"), header = TRUE)
exp.meta <- read.csv(file.path(inpath, "Probe_information.csv"), header = TRUE)
#---------------------------### Phenotype data ###------------------------------
# Read raw phenotype data
phen.1 <- read.csv(file.path(inpath, "Sample_info_LCL_only.csv"), header = TRUE)
phen.2 <- read.csv(file.path(inpath, "Sample_info_PBMC_only.csv"), header = TRUE)
# Get project-specific sample IDs
colnames(exp.1) <- phen.1[ ,2]
colnames(exp.2) <- phen.2[ ,2]
# Sort by sample IDs
exp.1 <- exp.1[ ,order(colnames(exp.1))]
exp.2 <- exp.2[ ,order(colnames(exp.2))]
# Extract Illumina probe IDs
PROBE_ID <- exp.meta[ ,1]
exp.1    <- cbind(PROBE_ID, exp.1)
exp.2    <- cbind(PROBE_ID, exp.2)

# Remove redundant columns
phen.1 <- phen.1[ ,-c(1,3)]
phen.2 <- phen.2[ ,-c(1,3)]
# Alter column names for consistency across files
names    <- sapply(colnames(phen.1), toupper)
names[1] <- "SAMPLE_ID"
colnames(phen.1) <- names
colnames(phen.2) <- names
# Re-order rows by sample ID
phen.1 <- phen.1[order(phen.1[ ,1]), ]
phen.2 <- phen.2[order(phen.2[ ,1]), ]
#----------------------------### Genotype data ###------------------------------
labels <- c("AA", "AB", "BB")
alleles <- c("A", "C", "G", "T")
#---------------------------------## LCL ##-------------------------------------
# Read raw genotype data
gen.raw.1 <- read.plink(file.path(inpath, "cleaned_geno_lcl.bed"))
gen.1 <- gen.raw.1$genotypes
# Recode SNPs as AA, AB or BB 
gen.1 <- apply(gen.1@.Data, 2, as.numeric)
gen.1[gen.1 == 0] <- NA
# Transpose genotype matrix for consistency between datasets
gen.1 <- t(gen.1)
RS_ID <- colnames(gen.raw.1$genotypes)
SAMPLE_ID <- gen.raw.1$fam[ ,2]
# Metadata
gen.meta.1 <- gen.raw.1$map
rownames(gen.meta.1) <- c(1:nrow(gen.meta.1))
gen.meta.1 <- cbind(RS_ID, gen.meta.1[ ,-c(2,3)])
# Resolve allele IDs to nucleotide bases
alleles <- c("A", "C", "G", "T")
gen.meta.1[ ,4] <- alleles[gen.meta.1[ ,4]]
gen.meta.1[ ,5] <- alleles[gen.meta.1[ ,5]]
# Enumerate actual allele combinations
alleles <- matrix(nrow=nrow(gen.meta.1), ncol=3)
for (i in 1:nrow(gen.meta.1)) {
    alleles[i, ] <- c(paste(gen.meta.1[i,4], gen.meta.1[i,4], sep=""),
                      paste(gen.meta.1[i,4], gen.meta.1[i,5], sep=""),
                      paste(gen.meta.1[i,5], gen.meta.1[i,5], sep=""))
}

# Encode genotypes using allele combinations
gen.res <- matrix(nrow=nrow(gen.1), ncol=ncol(gen.1))
for (i in 1:nrow(gen.1)) {
    gen.res[i, ] <- alleles[i,gen.1[i, ]]
}
# Set genotype row and column names
colnames(gen.res) <- SAMPLE_ID
gen.res <- cbind(RS_ID, gen.res)
gen.1   <- gen.res
# Correct column labelling for consistency between datasets
names   <- sapply(colnames(gen.meta.1), toupper)
names   <- gsub("\\.", "_", names)
colnames(gen.meta.1) <- names
#---------------------------------## PBMC ##------------------------------------
# Read raw genotype data
gen.raw.2 <- read.plink(file.path(inpath, "cleaned_geno_pbmc.bed"))
gen.2 <- gen.raw.2$genotypes
# Recode SNPs as 0, 1, or 2 for AA, AB, or BB, respectively
gen.2 <- apply(gen.2@.Data, 2, as.numeric)
gen.2[gen.2 == 0] <- NA
# Transpose genotype matrix for consistency between datasets
gen.2 <- t(gen.2)
RS_ID <- colnames(gen.raw.2$genotypes)
SAMPLE_ID  <- gen.raw.2$fam[ ,2]
# Metadata
gen.meta.2 <- gen.raw.2$map
rownames(gen.meta.2) <- c(1:nrow(gen.meta.2))
gen.meta.2 <- cbind(RS_ID, gen.meta.2[ ,-c(2,3)])
# Resolve allele IDs to nucleotide bases
alleles <- c("A", "C", "G", "T")
gen.meta.2[ ,4] <- alleles[gen.meta.2[ ,4]]
gen.meta.2[ ,5] <- alleles[gen.meta.2[ ,5]]
# Resolve allele labels to nucleotide bases
# Enumerate actual allele combinations
alleles <- matrix(nrow=nrow(gen.meta.2), ncol=3)
for (i in 1:nrow(gen.meta.2)) {
    alleles[i, ] <- c(paste(gen.meta.2[i,4], gen.meta.2[i,4], sep=""),
                      paste(gen.meta.2[i,4], gen.meta.2[i,5], sep=""),
                      paste(gen.meta.2[i,5], gen.meta.2[i,5], sep=""))
}

# Encode genotypes using allele combinations
gen.res <- matrix(nrow=nrow(gen.2), ncol=ncol(gen.2))
for (i in 1:nrow(gen.2)) {
    gen.res[i, ] <- alleles[i,gen.2[i, ]]
}
# Set genotype row and column names
colnames(gen.res) <- SAMPLE_ID
gen.res <- cbind(RS_ID, gen.res)
gen.2   <- gen.res
# Correct column labelling for consistency between datasets
colnames(gen.meta.2) <- names
#---------------------### Make data entries consistent ###----------------------
#---------------------------------## LCL ##-------------------------------------
cols   <- which(!colnames(exp.1) %in% colnames(gen.1))
cols   <- cols[-1]
exp.1  <- exp.1[ ,-cols]
cols   <- cols - 1
phen.1 <- phen.1[-cols, ]
#---------------------------------## PBMC ##------------------------------------
cols   <- which(!colnames(exp.2) %in% colnames(gen.2))
cols   <- cols[-1]
exp.2  <- exp.2[ ,-cols]
cols   <- cols - 1
phen.2 <- phen.2[-cols, ]
#----------------------------### Write to file ###------------------------------
if (all(colnames(exp.1)[-1] == colnames(gen.1)[-1]) &&
  all(colnames(exp.2)[-1] == colnames(gen.2)[-1]) &&
  all(colnames(exp.1)[-1] == phen.1[ ,1])         &&
  all(colnames(exp.2)[-1] == phen.2[ ,1])) {
  Write(exp.1, "BSGSpilot_expression_signals-LCL.txt")
  Write(exp.2, "BSGSpilot_expression_signals-PBMC.txt")
  Write(phen.1, "BSGSpilot_phenotypes-LCL.txt")
  Write(phen.2, "BSGSpilot_phenotypes-PBMC.txt")
  Write(gen.1, "BSGSpilot_genotypes-LCL.txt")
  Write(gen.2, "BSGSpilot_genotypes-PBMC.txt")
  Write(gen.meta.1, "BSGSpilot_genotypes-meta-LCL.txt")
  Write(gen.meta.2, "BSGSpilot_genotypes-meta-PBMC.txt")
} else {
  print("Error: Sample labelling is inconsistent between matrices.")
}
#------------------------------### Clean up ###---------------------------------
rm(inpath, outpath, PROBE_ID, names, gen.raw.1, gen.raw.2, gen.res, alleles,
   labels, RS_ID, SAMPLE_ID, cols, i, write)
