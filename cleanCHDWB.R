                             ####################
                             # Clean CHDWB data #
                             ####################

#--------------------------------### Setup ###----------------------------------
inpath  <- "/ibscratch/wrayvisscher/xander/CAGE/data/CHDWB/raw_2014-09-01"
outpath <- "/ibscratch/wrayvisscher/xander/CAGE/data/CHDWB/clean"
require(snpStats)
source("/clusterdata/uqahollo/scripts/Write.R")

FixId <- function(id) {
  # Correct formatting of GG_#### sample IDs found in the batch 2 and 3 ID map
  # files to match GG2_####, as found in the full study design map. GG1_######
  # format IDs are left intact.
  #   id: the ID number to process.
  # Returns: ID number either in GG2_####, or GG1_###### (unchanged) format.
  if (grepl("GG_", id)[1]) {
    id <- gsub("GG_", "GG2_", id)
  }
  id
}

FixIds <- function(map) {
  # Correct formatting of "GG_####" sample IDs found in the ID map files to 
  # match "GG1_######" sample IDs found in expression data.
  #   map: ID map file, read from .csv file
  # Returns: ID map data.frame, with correctly formatted IDs in appended column.
  names <- map$ID_3
  names <- gsub("GG.?_", "", names)
  names <- formatC(as.numeric(names), width = 6, format = "d", flag = "0")
  names <- paste0("GG1_", names)
  map$ID_FIX <- names
  map
}

RemoveIncomplete <- function(map) {
  # Convenience function for removing inomplete cases in ID map data.frames.
  # Used to strip samples that do not have a valid mapping from genotype to 
  # expression ID.
  #   map: ID map file, read from .csv file
  # Returns: ID map data.frame, with incomplete rows removed.
  map[map == ""] <- NA
  map[complete.cases(map), ]
}
#---------------------------### Expression data ###-----------------------------
exp.raw <- read.table(file = file.path(inpath, "chdwb_final_log2.txt"),
                       sep = "\t",
                    header = TRUE)
# Read ID map files
map.1 <- read.csv(file.path(inpath, "genotypes/batch1_IDs.csv"), header = TRUE)
map.1 <- RemoveIncomplete(map.1)
map.2 <- read.csv(file.path(inpath, "genotypes/batch2_IDs.csv"), header = TRUE)
map.2 <- RemoveIncomplete(map.2)
map.3 <- read.csv(file.path(inpath, "genotypes/batch3_IDs.csv"), header = TRUE)
map.3 <- RemoveIncomplete(map.3)
# Complete ID map resolves all ID associations
map.new <- read.csv(file.path(inpath, "chdwb_all_design.csv"), header = TRUE)
map.new$Sample  <- gsub("GG2-", "GG2_", map.new$Sample) # subtle difference in ID
map.new$CHD_ID  <- toupper(map.new$CHD_ID)
# Correct formatting of IDs to match expression data
map.2 <- FixIds(map.2)
map.3 <- FixIds(map.3)
# Get all unique sample names in "GG1_######" format
names <- unique(c(map.1$ID_3, map.2$ID_FIX, map.3$ID_FIX))
# Correct expression IDs
PROBE_ID <- exp.raw[, 2]
exp      <- exp.raw[, -c(1:2)]
exp      <- exp[, which(colnames(exp) %in% map.new$Sample)]
exp.col  <- map.new[which(map.new$Sample %in% colnames(exp)),1]
colnames(exp) <- exp.col
exp      <- cbind(PROBE_ID, exp)
#------------------------------### Metadata ###-------------------------------
# sample info
info    <- read.table(file = file.path(inpath, "chdwb_final_exptdes.txt"),
                       sep = "\t",
                    header = TRUE)
info$ColumnName <- map.new$CHD_ID
colnames(info)  <- c("SAMPLE_ID", toupper(colnames(info)[-1]))
# covariates
cov <- read.csv(file = file.path(inpath, "sample_info_CHDWB.csv"),
              header = TRUE)
col <- grep("CHDWB", names(cov)) # move sample ID to first column
cov <- cov[, c(col, (1:ncol(cov))[-col])]
# TODO: split column names with underscores
colnames(cov) <- c("SAMPLE_ID", toupper(colnames(cov)[-1]))
cov$SAMPLE_ID <- gsub("WB", "", cov$SAMPLE_ID)
cov <- cov[which(cov$SAMPLE_ID %in% colnames(exp)), ]
# probe info
probe <- read.table(file = file.path(inpath, "CHDWB_probe_info.txt"),
                     sep = "\t",
                  header = TRUE)
probe <- probe[which(probe$PROBE_ID %in% exp$PROBE_ID), ] # strip to match exp
#---------------------------### Genotype data ###-----------------------------
plink.2 <- read.plink(file.path(inpath, "genotypes/batch2"))
plink.3 <- read.plink(file.path(inpath, "genotypes/batch3"))
gen.2   <- plink.2$genotypes@.Data
gen.2   <- gen.2[which(rownames(gen.2) %in% map.2$sample.ID_1), ]
gen.3   <- plink.3$genotypes@.Data
# Map genotype sample IDs to `map.new` IDs
names.2 <- map.2[which(map.2$sample.ID_1 %in% rownames(gen.2)),"ID_3"]
names.3 <- map.3[which(map.3$sample.ID_1 %in% rownames(gen.3)),"ID_3"]
names.2 <- sapply(names.2, FixId)
names.3 <- sapply(names.3, FixId)
names.2 <- map.new[which(map.new$Sample %in% names.2),"CHD_ID"]
names.3 <- map.new[which(map.new$Sample %in% names.3),"CHD_ID"]
# Set row names with Sample IDs
rownames(gen.2) <- names.2
rownames(gen.3) <- names.3
# Create SnpMartix objects
gen.2 <- new("SnpMatrix", gen.2)
gen.3 <- new("SnpMatrix", gen.3)
# Label batch 2 sample data
fam.2 <- plink.2$fam
fam.2 <- fam.2[which(rownames(fam.2) %in% map.2$sample.ID_1), ]
fam.2$pedigree  <- names.2
fam.2$member    <- names.2
rownames(fam.2) <- names.2
# Label batch 3 sample data
fam.3 <- plink.3$fam
fam.3$pedigree  <- names.3
fam.3$member    <- names.3
rownames(fam.3) <- names.3
#----------------------------### Write to file ###------------------------------
Write(exp, "CHDWB_exp_log2.txt")
Write(cov, "CHDWB_cov.txt")
Write(info, "CHDWB_sample_info.txt")
Write(probe, "CHDWB_probe_info.txt")
write.plink(file.base = file.path(outpath, "genotypes/CHDWB_batch2"),
            snp.major = TRUE,
                 snps = gen.2,
         subject.data = fam.2,
             snp.data = plink.2$map)
write.plink(file.base = file.path(outpath, "genotypes/CHDWB_batch3"),
            snp.major = TRUE,
                 snps = gen.3,
         subject.data = fam.3,
             snp.data = plink.3$map)
