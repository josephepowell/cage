                           ########################
                           # Clean BSGS main data #
                           ########################

#--------------------------------### Setup ###----------------------------------
require(snpStats)
inpath  <- "/ibscratch/wrayvisscher/xander/CAGE/data/BSGSmain/raw"
outpath <- "/ibscratch/wrayvisscher/xander/CAGE/data/BSGSmain/clean"
# Read "raw" data
load(file.path(inpath, "Pre_normalisation_data.Rdata"))
source("/clusterdata/uqahollo/scripts/Write.R")
#---------------------------### Expression data ###-----------------------------
# Rename and sort existing data structures
exp   <- probe_signal
names <- gsub("(X|\\..*)", "", colnames(exp))
colnames(exp) <- names
exp   <- exp[, order(colnames(exp))]
PROBE_ID <- probe_info[, 3]
exp   <- cbind(PROBE_ID, exp)
# Metadata
probe.info <- cbind(probe_info[, c(3,2,4:length(probe_info))])
colnames(probe.info)[6] <- "RS_ID"
colnames(probe.info) <- sapply(colnames(probe.info), toupper)
# p-values
exp.p <- probe_pval
names <- gsub("(X|\\..*)", "", colnames(exp.p))
colnames(exp.p) <- names
exp.p <- exp.p[, order(colnames(exp.p))]
exp.p <- cbind(PROBE_ID, exp.p)
#---------------------------### Phenotype data ###------------------------------
cov <- sample_info
colnames(cov) <- gsub("\\.", "_", colnames(cov))
colnames(cov) <- sapply(colnames(cov), toupper)
cov <- cov[order(cov[, 1]), ]
#----------------------------### Process info ###-------------------------------
info <- process_info
colnames(info)[c(3,4)] <- c("PROCESS_DATE", "RNA_EXTRACT_DATE")
#----------------------------### Genotype data ###------------------------------
# Read raw genotype data
gen.raw <- read.plink(file.path(inpath, "clean_geno_final.bed"))
gen     <- gen.raw$genotypes
# Encode genotypes as a number between 1 and 3 and remove null values
gen <- apply(gen@.Data, 2, as.numeric)
gen[gen == 0] <- NA
# Transpose genotype matrix for consistency between datasets
gen <- t(gen)
# Extract SNP and sample IDs
RS_ID     <- colnames(gen.raw$genotypes)
SAMPLE_ID <- gen.raw$fam[, 2]
# Metadata
gen.map <- gen.raw$map
rownames(gen.map) <- c(1:nrow(gen.map))
gen.map <- cbind(RS_ID, gen.map[, -c(2,3)])
# Resolve allele IDs to nucleotide bases
labels       <- c("A", "C", "G", "T")
gen.map[, 4] <- labels[gen.map[, 4]]
gen.map[, 5] <- labels[gen.map[, 5]]
# Enumerate actual allele combinations
alleles <- matrix(nrow = nrow(gen.map),
                  ncol = 3)
for (i in 1:nrow(gen.map)) {
  alleles[i, ] <- c(paste0(gen.map[i, 4], gen.map[i, 4]),
                    paste0(gen.map[i, 4], gen.map[i, 5]),
                    paste0(gen.map[i, 5], gen.map[i, 5]))
}

# Encode genotypes using allele combinations
gen.res <- matrix(nrow = nrow(gen),
                  ncol = ncol(gen))
for (i in 1:nrow(gen)) {
  gen.res[i, ] <- alleles[i, gen[i, ]]
}
# Set genotype row and column names
colnames(gen.res) <- SAMPLE_ID
gen.res <- cbind(RS_ID, gen.res)
gen     <- gen.res
# Correct column labelling for consistency between datasets
colnames(gen.map) <- sapply(colnames(gen.map), toupper)
colnames(gen.map) <- gsub("\\.", "_", colnames(gen.map))
#---------------------### Make data entries consistent ###----------------------
cols  <- which(!colnames(exp) %in% colnames(gen))
cols  <- cols[-1]
exp   <- exp[, -cols]
exp.p <-  exp.p[, -cols]
cols  <- cols - 1
cov   <- cov[-cols, ]
#----------------------------### Write to file ###------------------------------
Write(exp, "BSGSmain_exp.txt")
Write(exp.p, "BSGSmain_exp_pval.txt")
Write(probe.info, "BSGSmain_probe_info.txt")
Write(cov, "BSGSmain_cov.txt")
Write(info, "BSGSmain_process_info.txt")
Write(gen, "genotypes/BSGSmain_gen.txt")
Write(gen.map, "genotypes/BSGSmain_gen_map.txt")
# Copy PLINK files to clean directory
system(paste0("cp ", file.path(inpath, "*.bed "), file.path(outpath, "genotypes/BSGSmain.bed")))
system(paste0("cp ", file.path(inpath, "*.bim "), file.path(outpath, "genotypes/BSGSmain.bim")))
#------------------------------### Clean up ###---------------------------------
rm(inpath, outpath, PROBE_ID, RS_ID, SAMPLE_ID, cols, names, probe_info, probe_pval,
   probe_signal, process_info, sample_info, gen.res, alleles, labels, i, write)
