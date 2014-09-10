                            ######################
                            # Clean Morocco data #
                            ######################

#--------------------------------### Setup ###----------------------------------
inpath  <- "/ibscratch/wrayvisscher/xander/CAGE/data/Morocco/raw"
outpath <- "/ibscratch/wrayvisscher/xander/CAGE/data/Morocco/clean"

#---------------------------### Expression data ###-----------------------------
# Read raw expression data
exp.raw  <- read.delim(file.path(inpath, "Group_Probe_Profile_Gibson_1-3_11-12-08.txt"), sep = "\t", fill = TRUE)
# Extract probe labels and metadata from expression matrix
PROBE_ID <- exp.raw[ ,3]
exp      <- exp.raw[ ,4:length(exp.raw)]
exp.info <- exp.raw[ ,1:2]
# Replace default column names with project-specific sample IDs
colnames(exp) <- gsub("\\..*", "", colnames(exp))
# Sort samples by project-specific sample ID
exp      <- exp[ ,order(colnames(exp))]
# Insert platform-specific probe IDs
exp.info <- cbind(PROBE_ID, exp.info)
# Change column names for consistency between datasets
colnames(exp.info)[2:3] <- c("TARGET_ID", "PROBE_NUMBER")

#----------------------------### Genotype data ###------------------------------
# Read raw genotype data
gen.raw <- read.delim(file.path(inpath, "Morocco_194_genotypes_598k_tall.txt"), sep = "\t", header = TRUE)
# Extract SNP labels
RS_ID   <- gen.raw[ ,1]
gen     <- gen.raw[ ,-1]
# Sort samples by project-specific sample ID
gen     <- gen[ ,order(colnames(gen))]

#---------------------### Make data entries consistent ###----------------------
# Remove samples which have only one of expression or genotype data entries
cols <- which(!colnames(exp) %in% colnames(gen))
exp  <- exp[ ,-(cols)]
cols <- which(!colnames(gen) %in% colnames(exp))
gen  <- gen[ ,-(cols)]
# Insert platform-specific probe IDs
exp  <- cbind(PROBE_ID, exp)
gen  <- cbind(RS_ID, gen)

#-----------------------------### PLINK files ###-------------------------------
gen <- apply(gen, 2, function(x) ifelse(x == "--", "00", x))
ped <- matrix(nrow = ncol(gen) - 1, ncol = 6 + nrow(gen))
for (i in 1:nrow(ped)) {
  ped[i,1:6] <- c("", colnames(gen)[i + 1], "0", "0", "0", "0")
  ped[i,7:ncol(ped)] <- sapply(gen[ ,i + 1],
                               function(x) paste(substring(x, 1, 1), substring(x, 2)))
}
# Map file
map         <- matrix(nrow = nrow(gen), ncol = 4)
#alleles     <- data.frame(matrix(nrow = nrow(gen), ncol = 2))
#get.alleles <- function(x) { 
#                   temp <- unique(unlist(x)[-1])
#                   temp <- temp[order(temp)]
#                   temp <- strsplit(temp[2], "")
#                   temp
#}
#for (i in 1:nrow(map)) {
#    alleles[i, ] <- get.alleles(gen[i, ])[[1]]
#}

# TODO: fix NA values in .map file
map <- cbind(rep("0", nrow(gen)), # missing chromosomes
             gen[ ,1],
             rep("0", nrow(gen)),
             rep("0", nrow(gen)))
#              alleles)

#----------------------------### Write to file ###------------------------------
# Convenience function for easily changing output file format
write <- function(table, filename) write.table(table, file.path(outpath, filename),
                                               sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE)

if (all(colnames(exp)[-1] == colnames(gen)[-1])) {
  write(exp, "Morocco_expression_signals.txt")
  write(exp.info, "Morocco_probe_info.txt")
  write(gen, "genotypes/Morocco_genotypes.txt")
  write.table(ped, paste(outpath, "/genotypes/Morocco.ped", sep = ""), sep = "\t", eol = "\n",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(map, paste(outpath, "/genotypes/Morocco.map", sep = ""), sep = "\t", eol = "\n",
              quote = FALSE, row.names = FALSE, col.names = FALSE)
} else {
  print("Error: Sample labelling is not consistent between expression and genotype data")
}

#------------------------------### Clean up ###---------------------------------
#rm(inpath, outpath, PROBE_ID, RS_ID, names, cols, write, alleles, getAlleles)
