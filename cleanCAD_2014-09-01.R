                              ##################
                              # Clean CAD data #
                              ##################

#--------------------------------### Setup ###----------------------------------
inpath  <- "/ibscratch/wrayvisscher/xander/CAGE/data/CAD/raw"
outpath <- "/ibscratch/wrayvisscher/xander/CAGE/data/CAD/clean"

PadWithZeroes <- function(id) {
  # Pad the numeric portion of sample IDs with zeroes, to five characters.
  #   id: sample ID number to pad
  # Returns: sample ID, with numeric portion padded to length 5.
  sep <- strsplit( gsub("([A-Z]{3})", "\\1~", id), "~")
  pre <- sep[[1]][1]
  suf <- sep[[1]][2]
  suf <- formatC(as.numeric(suf), width = 5, format = "d", flag = "0")
  id  <- paste0(pre, suf)
}

#---------------------------### Expression data ###-----------------------------
exp.raw.1 <- read.table(file = file.path(inpath, "Cardiology_phase1_raw.txt"),
                      header = TRUE,
                         sep = "\t")
exp.raw.2 <- read.table(file = file.path(inpath, "Cardiology_phase2_raw.txt"),
                      header = TRUE,
                         sep = "\t")
# Extract probe info from phase 1 matrix
PROBE_ID   <- exp.raw.1[, 1]
exp.1      <- exp.raw.1[, order(colnames(exp.raw.1))]
exp.info.1 <- exp.1[, 1:10]
exp.1      <- exp.1[, -c(1:11)] # strip probe info

# Extract STUDY_ID from phase 1 expression data
exp.id <- colnames(exp.1)
exp.id <- gsub("X[0-9]*_|\\.[A-z]*", "", exp.id)
exp.id <- sapply(exp.id, PadWithZeroes)
colnames(exp.1) <- exp.id

phen <- read.table(file.path(inpath, "Cardiology_exptdes_bothphases.txt"), sep = "\t", header = TRUE)

# Extract STUDY_ID from genotype data
gen.plink <- snpStats::read.plink("/ibscratch/wrayvisscher/xander/CAGE/data/CAD/clean/CAD")
gen <- t(gen.plink$genotypes@.Data)
gen.id <- gsub("X[0-9]*_|_C", "", colnames(gen))
gen.id <- sapply(gen.id, PadWithZeroes)
colnames(gen) <- gen.id

# re-label FAM matrix
fam <- gen.plink$fam
fam <- fam[order(rownames(fam)), ]
rownames(fam) <- gen.id
fam$pedigree  <- gen.id
fam$member    <- gen.id
# remove SNPs that lack allele info
map <- gen.plink$map
map <- map[which(complete.cases(map[, 5:6])), ]
gen <- gen[rownames(map), ]

names <- gen.id[which(gen.id %in% exp.id)]
exp.1 <- exp.1[, names]
gen   <- gen[, names]
# Append probe IDs
exp.1 <- cbind(PROBE_ID, exp.1)

#----------------------------### Write to file ###------------------------------
write.table(exp.1, file.path(outpath, "CAD_exp.txt"), sep = "\t", eol = "\n", quote = FALSE, row.names = FALSE)

write.plink(file.base = file.path(outpath, "CAD"),
            snp.major = TRUE,
                 snps = gen,
         subject.data = fam,
             snp.data = map)
