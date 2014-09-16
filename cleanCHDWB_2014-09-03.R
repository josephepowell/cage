                             ####################
                             # Clean CHDWB data #
                             ####################

#--------------------------------### Setup ###----------------------------------
inpath  <- "/ibscratch/wrayvisscher/xander/CAGE/data/CHDWB/raw_2014-09-01"
outpath <- "/ibscratch/wrayvisscher/xander/CAGE/data/CHDWB/clean"
source("/clusterdata/uqahollo/scripts/Write.R")

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
# Correct formatting of IDs to match expression data
map.3 <- FixIds(map.3)
map.2 <- FixIds(map.2)
# Get all unique sample names in "GG1_######" format
names <- unique(c(map.1$ID_3, map.2$ID_FIX, map.3$ID_FIX))

batch.2 <- snpStats::read.plink(file.path(inpath, "genotypes/batch2"))
batch.3 <- snpStats::read.plink(file.path(inpath, "genotypes/batch3"))
gen.2   <- batch.2$genotypes@.Data
gen.3   <- batch.3$genotypes@.Data




names.2 <- rownames(gen.2)

names.3 <- rownames(gen.3)
names.3 <- gsub(".*_", "", names.3)


rownames(gen.2) <- map.2["sample.ID_1" %in% names.2, "ID_2"]
