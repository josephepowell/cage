                           ########################
                           # Clean BSGS main data #
                           ########################

#--------------------------------### Setup ###----------------------------------
require(snpStats)
inpath  <- "/ibscratch/wrayvisscher/xander/CAGE/data/BSGSmain/raw"
outpath <- "/ibscratch/wrayvisscher/xander/CAGE/data/BSGSmain/clean"
load(file.path(inpath, "Pre_normalisation_data.Rdata")) # read "raw" data
source("/clusterdata/uqahollo/scripts/Write.R")
source("/clusterdata/uqahollo/scripts/ConvertToDosage.R")
#---------------------------### Expression data ###-----------------------------
# Rename and sort existing data structures
exp   <- probe_signal
names <- gsub("(X|\\..*)", "", colnames(exp))
colnames(exp) <- names
exp   <- exp[, order(colnames(exp))]
PROBE_ID <- probe_info[, 3]
exp   <- cbind(PROBE_ID, exp)
# Metadata
probe.info <- cbind(probe_info[, c(3, 2, 4:length(probe_info))])
colnames(probe.info)[6] <- "RS_ID"
colnames(probe.info) <- sapply(colnames(probe.info), toupper)
# p-values
exp.p <- probe_pval
names <- gsub("(X|\\..*)", "", colnames(exp.p))
colnames(exp.p) <- names
exp.p <- exp.p[, order(colnames(exp.p))]
exp.p <- cbind(PROBE_ID, exp.p)
#-----------------------------### Covariates ###--------------------------------
cov  <- sample_info
colnames(cov) <- gsub("\\.", "_", colnames(cov))
colnames(cov) <- sapply(colnames(cov), toupper)
cov  <- cov[order(cov[, 1]), ]
info <- process_info
colnames(info)[c(1, 3, 4)] <- c("SAMPLE_ID", "PROCESS_DATE", "RNA_EXTRACT_DATE")
#----------------------------### Genotype data ###------------------------------
gen.plink <- read.plink(file.path(inpath, "clean_geno_final.bed"))
gen       <- ConvertToDosage(gen.plink)
gen       <- gen + 1 # add 1 to index alleles correctly
# Extract SNP and sample IDs
RS_ID     <- colnames(gen.plink$genotypes)
SAMPLE_ID <- gen.plink$fam[, 2]
# Metadata
gen.map <- gen.plink$map
rownames(gen.map) <- c(1:nrow(gen.map))
gen.map <- cbind(RS_ID, gen.map[, -c(2, 3)])
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
info  <- info[which(info$SAMPLE_ID %in% colnames(exp)), ]
#----------------------------### Write to file ###------------------------------
Write(exp, "BSGSmain_exp.txt")
Write(exp.p, "BSGSmain_exp_pval.txt")
Write(probe.info, "BSGSmain_probe_info.txt")
Write(cov, "BSGSmain_cov.txt")
Write(info, "BSGSmain_process_info.txt")
Write(gen, "genotypes/BSGSmain_gen.txt")
Write(gen.map, "genotypes/BSGSmain_gen_map.txt")
# Copy PLINK files to clean directory
system(paste0("cp ", file.path(inpath, "*.bed "), file.path(outpath, "genotypes/plink/BSGSmain.bed")))
system(paste0("cp ", file.path(inpath, "*.bim "), file.path(outpath, "genotypes/plink/BSGSmain.bim")))
#------------------------------### Clean up ###---------------------------------
rm(inpath, outpath, PROBE_ID, RS_ID, SAMPLE_ID, cols, names, probe_info, probe_pval,
   probe_signal, process_info, sample_info, gen.res, alleles, labels, i)
