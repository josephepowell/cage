                             #####################
                             # Clean MuTHER data #
                             #####################

#--------------------------------### Setup ###----------------------------------
inpath  <- "/ibscratch/wrayvisscher/xander/CAGE/data/MuTHER/raw"
outpath <- "/ibscratch/wrayvisscher/xander/CAGE/data/MuTHER/clean"

#---------------------------### Expression data ###-----------------------------
# Read raw expression data
exp.raw.1 <- read.delim(file.path(inpath, "MuTHER_Fat_main_raw_data_newIds.txt"), sep = "\t", fill = TRUE)
exp.raw.2 <- read.delim(file.path(inpath, "MuTHER_LCL_main_raw_data_newIds.txt"), sep = "\t", fill = TRUE)
exp.raw.3 <- read.delim(file.path(inpath, "MuTHER_SKIN_main_raw_data_newIds.txt"), sep = "\t", fill = TRUE)
# Replace column names
colnames(exp.raw.1) <- gsub("_.*", "", colnames(exp.raw.1))
colnames(exp.raw.2) <- gsub("_.*", "", colnames(exp.raw.2))
colnames(exp.raw.3) <- gsub("_.*", "", colnames(exp.raw.3))
# Extract probe labels
PROBE_ID.1 <- exp.raw.1[ ,1]
PROBE_ID.2 <- exp.raw.2[ ,1]
PROBE_ID.3 <- exp.raw.3[ ,1]
# Separate expression signals from bead counts
exp.1  <- exp.raw.1[ ,seq(2, length(exp.raw.1), by = 2)]
bead.1 <- exp.raw.1[ ,seq(3, length(exp.raw.1), by = 2)]
exp.2  <- exp.raw.2[ ,seq(2, length(exp.raw.2), by = 2)]
bead.2 <- exp.raw.2[ ,seq(3, length(exp.raw.2), by = 2)]
exp.3  <- exp.raw.3[ ,seq(2, length(exp.raw.3), by = 2)]
bead.3 <- exp.raw.3[ ,seq(3, length(exp.raw.3), by = 2)]
# Re-order samples
exp.1 <- exp.1[ ,order(colnames(exp.1))]
exp.1 <- exp.1[order(exp.1[ ,1]), ]
exp.2 <- exp.2[ ,order(colnames(exp.2))]
exp.3 <- exp.3[ ,order(colnames(exp.3))]
# Split replicates
exp.1 <- exp.1[ ,grep(".*\\.+.", colnames(exp.1), invert = TRUE)]
exp.2 <- exp.2[ ,grep(".*\\.+.", colnames(exp.2), invert = TRUE)]
exp.3 <- exp.3[ ,grep(".*\\.+.", colnames(exp.3), invert = TRUE)]
# Add probe IDs
exp.1 <- cbind(PROBE_ID.1, exp.1)
colnames(exp.1)[1] <- "PROBE_ID"
rownames(bead.1)   <- PROBE_ID.1
exp.2 <- cbind(PROBE_ID.2, exp.2)
colnames(exp.2)[1] <- "PROBE_ID"
rownames(bead.2)   <- PROBE_ID.2
exp.3 <- cbind(PROBE_ID.3, exp.3)
colnames(exp.3)[1] <- "PROBE_ID"
rownames(bead.3)   <- PROBE_ID.3

#----------------------------### Genotype data ###------------------------------
# Genotype data is available for 2,601 sampels of various tissue types, in the 
# form of posterior probabilities. These dosages are not labelled by sample 
# ID in the raw form, but the order of their occurrence is described in the 
# associated .sample file. There are three columns per sample (one for each 
# combination of alleles).

#--------------------------## Setup for iteration ##----------------------------
# Read sample info, applicable to all genotypes
cov <- read.delim(file.path(inpath, "genotypes/MuTHER_Covariates_File_Public.txt"), sep = "\t", header = TRUE)
info <- read.delim(file.path(inpath, "genotypes/MuTHER_release_100912.sample"), sep = " ", header = TRUE)
info <- info[-1,1]
# Allele labels
labels <- c("AA", "AB", "BB")
# Convenience function for easily changing output file format
write <- function(table, filename) write.table(table, file.path(outpath, filename),
                                               sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE)
#-------------------## Iterate over all chromosome files ##---------------------
gen.files <- list.files(file.path(inpath, "genotypes"))
sequence <-  seq(1, length(gen.files), by = 2)
for (i in 1:ceiling((length(gen.files) - 3) / 2)) {
    # Extract name of chromosome
    chrom <- strsplit(gen.files[sequence[i]], "\\_")[[1]][1]
    # Read raw data
    gen.raw <- read.delim(file.path(inpath,paste("genotypes/", chrom, "_MuTHER_release_100912_impute.gen.gz", sep="")),
                          sep=" ", header=FALSE)
    gen.meta <- read.delim(file.path(inpath,paste("genotypes/", chrom, "_MuTHER_release_100912_impute.info.gz", sep="")),
                          sep=" ", header=FALSE)
    # Create matrix to store genotype posterior probabilities
    gen.post <- data.matrix( gen.raw[ ,6:(length(gen.raw))] )
    gen.post <- gen.post[ ,-ncol(gen.post)]
    temp <- matrix(nrow = 1,
                   ncol = ncol(gen.post))
    colnames(temp) <- colnames(gen.post)
    gen.post <- rbind(temp, gen.post)
    # Label posterior probabilities by alleles
    gen.post[1, ] <- labels
    times <- rep(3, length(info))
    names <- rep(info, times)
    colnames(gen.post) <- names
    gen.post <- cbind(c("", gen.raw[ ,2]), gen.post)
    colnames(gen.post)[1] <- "RS_ID"
    # Create genotype matrix
    gen.clean <- gen.post[-1,-1]
    # Flatten matrix, and transpose for column-wise processing
    gen.flat <- as.vector(t(gen.clean))
    gen.flat <- vapply(split(gen.flat, ceiling(1:length(gen.flat) / 3)),
                       which.max,
                       -1,
                       USE.NAMES=FALSE)
    # Convert maximum indices back to SNP labels
    gen.res <- matrix(gen.flat,
                      byrow = TRUE,
                      ncol  = ncol(gen.clean) / 3)
    # Enumerate actual allele combinations
    alleles <- matrix(nrow = nrow(gen.raw),
                      ncol = 3)
    # TODO: Vectorise for-loop
    for (i in 1:nrow(gen.raw)) {
        alleles[i, ] <- c(paste(gen.raw[i,4], gen.raw[i,4], sep=""),
                          paste(gen.raw[i,4], gen.raw[i,5], sep=""),
                          paste(gen.raw[i,5], gen.raw[i,5], sep=""))
    }
    # Encode genotypes using allele combinations
    gen <- matrix(nrow = nrow(gen.res),
                  ncol = ncol(gen.res))
    #TODO: Vectorise for-loop
    for (i in 1:nrow(gen.res)) {
	    gen[i, ] <- alleles[i,gen.res[i, ]]
    }
    # Set column and row names
    colnames(gen) <- info
    gen <- cbind(gen.raw[ ,2], gen)
    colnames(gen)[1] <- "RS_ID"
    gen <- gen[order(gen[ ,1]), ]
    #-----------------------### Format PLINK files ###--------------------
    ped <- matrix(nrow = ncol(gen) - 1, ncol = 6 + nrow(gen))
    for (i in 1:ncol(gen) - 1) {
        ped[i,1:6] <- c("", colnames(gen)[i + 1], "0", "0", "0", "0")
        ped[i,7:ncol(ped)] <- sapply(gen[ ,i + 1],
                                     function(x) paste(substring(x, 1, 1), substring(x, 2)))
    }
    map        <- matrix(nrow = nrow(gen), ncol = 4)
    map[ ,1:4] <- c(rep(chrom, nrow(map)),
                    gen[ ,1],
                    rep("0", nrow(map)),
                    rep("0", nrow(map)))

    #--------------------------### Write to file ###----------------------
    write(gen, paste("genotypes/MuTHER_genotypes-", chrom, ".txt", sep = ""))
    write(gen.post, paste("genotypes/MuTHER_genotypes-", chrom, "-posteriors.txt", sep = ""))
    write.table(ped, paste(outpath, "/genotypes/plink/MuTHER-", chrom, ".ped", sep = ""),
                sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(map, paste(outpath, "/genotypes/plink/MuTHER-", chrom, ".map", sep = ""),
                sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# TODO: Merge all PLINK files before write
#system("bash ~/scripts/cage/plink/mergeMuTHER.sh")

#----------------------------### Write to file ###------------------------------
write(exp.1, "MuTHER_expression_signals-fat.txt")
write(exp.2, "MuTHER_expression_signals-LCL.txt")
write(exp.3, "MuTHER_expression_signals-skin.txt")
write(bead.1, "MuTHER_bead_counts-fat.txt")
write(bead.2, "MuTHER_bead_counts-LCL.txt")
write(bead.3, "MuTHER_bead_counts-skin.txt")

#------------------------------### Clean up ###---------------------------------
rm(inpath, outpath, ids, names, order, temp, labels, gen.files, sequence,
   i, j, chrom, times, alleles, gen.raw, gen.clean, gen.flat, gen.res, write)
