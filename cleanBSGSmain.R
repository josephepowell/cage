                           ########################
                           # Clean BSGS main data #
                           ########################

#--------------------------------### Setup ###----------------------------------
require(snpStats)
inpath  <- "/ibscratch/wrayvisscher/xander/CAGE/data/BSGSmain/raw"
outpath <- "/ibscratch/wrayvisscher/xander/CAGE/data/BSGSmain/clean"
# Read "raw" data
load(file.path(inpath, "Pre_normalisation_data.Rdata"))

#---------------------------### Expression data ###-----------------------------
# Rename and sort existing data structures
exp   <- probe_signal
names <- gsub("(X|\\..*)", "", colnames(exp))
colnames(exp) <- names
exp   <- exp[ ,order(colnames(exp))]
PROBE_ID <- probe_info[ ,3]
exp   <- cbind(PROBE_ID, exp)
# Metadata
exp.meta <- cbind(probe_info[ ,c(3,2,4:length(probe_info))])
colnames(exp.meta)[6] <- "RS_ID"
colnames(exp.meta) <- sapply(colnames(exp.meta), toupper)
# p-values
exp.p <- probe_pval
names <- gsub("(X|\\..*)", "", colnames(exp.p))
colnames(exp.p) <- names
exp.p <- exp.p[ ,order(colnames(exp.p))]
exp.p <- cbind(PROBE_ID, exp.p)

#---------------------------### Phenotype data ###------------------------------
phen <- sample_info
colnames(phen) <- gsub("\\.", "_", colnames(phen))
colnames(phen) <- sapply(colnames(phen), toupper)
phen <- phen[order(phen[ ,1]), ]

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
SAMPLE_ID <- gen.raw$fam[ ,2]
# Metadata
gen.meta <- gen.raw$map
rownames(gen.meta) <- c(1:nrow(gen.meta))
gen.meta <- cbind(RS_ID, gen.meta[ ,-c(2,3)])
# Resolve allele IDs to nucleotide bases
labels        <- c("A", "C", "G", "T")
gen.meta[ ,4] <- labels[gen.meta[ ,4]]
gen.meta[ ,5] <- labels[gen.meta[ ,5]]
# Enumerate actual allele combinations
alleles <- matrix(nrow = nrow(gen.meta),
                  ncol = 3)
for (i in 1:nrow(gen.meta)) {
    alleles[i, ] <- c(paste(gen.meta[i,4], gen.meta[i,4], sep=""),
                      paste(gen.meta[i,4], gen.meta[i,5], sep=""),
                      paste(gen.meta[i,5], gen.meta[i,5], sep=""))
}

# Encode genotypes using allele combinations
gen.res <- matrix(nrow = nrow(gen),
                  ncol = ncol(gen))
for (i in 1:nrow(gen)) {
    gen.res[i, ] <- alleles[i,gen[i, ]]
}
# Set genotype row and column names
colnames(gen.res) <- SAMPLE_ID
gen.res <- cbind(RS_ID, gen.res)
gen     <- gen.res
# Correct column labelling for consistency between datasets
colnames(gen.meta) <- sapply(colnames(gen.meta), toupper)
colnames(gen.meta) <- gsub("\\.", "_", colnames(gen.meta))

#---------------------### Make data entries consistent ###----------------------
cols  <- which(!colnames(exp) %in% colnames(gen))
cols  <- cols[-1]
exp   <- exp[ ,-cols]
exp.p <-  exp.p[ ,-cols]
cols  <- cols - 1
phen  <- phen[-cols, ]

#----------------------------### Write to file ###------------------------------
# Convenience function for easily changing output file format
write <- function(table, filename) write.table(table, file.path(outpath, filename),
                                               sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE)
if (all(colnames(exp)[-1] == colnames(gen)[-1]) &&
    all(colnames(exp)[-1] == phen[ ,1])) {
    write(exp, "BSGSmain_expression_signals.txt")
    write(exp.meta, "BSGSmain_expression_signals-meta.txt")
    write(exp.p, "BSGSmain_expression_signals-pval.txt")
    write(phen, "BSGSmain_phenotypes.txt")
    write(info, "BSGSmain_process_info.txt")
    write(gen, "genotypes/BSGSmain_genotypes.txt")
    write(gen.meta, "genotypes/BSGSmain_genotypes-meta.txt")
} else {
    print("Error: Sample labelling is inconsistent between matrices.")
}

# Copy PLINK files to clean directory
system(paste("cp ", file.path(inpath, "*.bed "), file.path(outpath, "genotypes/BSGSmain.bed"), sep = ""))
system(paste("cp ", file.path(inpath, "*.bim "), file.path(outpath, "genotypes/BSGSmain.bim"), sep = ""))

#------------------------------### Clean up ###---------------------------------
rm(inpath, outpath, PROBE_ID, RS_ID, SAMPLE_ID, cols, names, probe_info, probe_pval,
   probe_signal, process_info, sample_info, gen.res, alleles, labels, i, write)
