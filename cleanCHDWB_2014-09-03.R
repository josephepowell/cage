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
  id <- map.2[map.2 == id,3]
  if (grepl("GG_", id)) {
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
info    <- read.table(file = file.path(inpath, "chdwb_final_exptdes.txt"),
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
#---------------------------### Genotype data ###-----------------------------
plink.2 <- read.plink(file.path(inpath, "genotypes/batch2"))
plink.3 <- read.plink(file.path(inpath, "genotypes/batch3"))
gen.2   <- plink.2$genotypes@.Data
gen.3   <- plink.3$genotypes@.Data
# TODO: Rectify sample retention according to `map.new` IDs
# TODO: Retain both GG1_###### and GG2_#### IDs for `map.new`
gen.2   <- gen.2[which(rownames(gen.2) %in% map.2$sample.ID_1), ]
rownames(gen.2) <- map.2[which(map.2$sample.ID_1 %in% rownames(gen.2)),"ID_FIX"]
gen.2   <- gen.2[which(rownames(gen.2) %in% map.new$Sample), ]
rownames(gen.2) <- map.new[which(map.new$Sample %in% rownames(gen.2)),"CHD_ID"]





# Store original rownames for sample ID mapping
names.2 <- rownames(gen.2)
names.3 <- rownames(gen.3)
# TODO: Mapping of sample IDs
names.3 <- gsub(".*_", "", names.3)
# rownames(gen.2) <- map.2["sample.ID_1" %in% names.2, "ID_2"]
#----------------------------### Write to file ###------------------------------
# TODO: Write to file
