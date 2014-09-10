                           #########################
                           #   Combine CAGE data   #
                           # into unified matrices #
                           #########################

#--------------------------------### Setup ###----------------------------------
inpath <- "/ibscratch/wrayvisscher/xander/CAGE/data"
# Convenience function for reading matrices
read <- function(filename) read.table(file.path(inpath, filename), sep = "\t", header = TRUE)

#---------------------------### Expression data ###-----------------------------
exp.bsgsmain       <- read("BSGSmain/clean/BSGSmain_expression_signals.txt")
exp.bsgspilot.lcl  <- read("BSGSpilot/clean/BSGSpilot_expression_signals-LCL.txt")
exp.bsgspilot.pbmc <- read("BSGSpilot/clean/BSGSpilot_expression_signals-PBMC.txt")
exp.chdwb          <- read("CHDWB/clean/CHDWB_expression_signals-log2.txt")
exp.morocco        <- read("Morocco/clean/Morocco_expression_signals.txt")
exp.muther.lcl     <- read("MuTHER/clean/MuTHER_expression_signals-LCL.txt")
exp.muther.fat     <- read("MuTHER/clean/MuTHER_expression_signals-fat.txt")
exp.muther.skin    <- read("MuTHER/clean/MuTHER_expression_signals-skin.txt")
#exp.cad            <- read("CAD/clean/CAD_expression_signals.txt")
#exp.egcut          <-

# CHDWB names contain . instead of - for unknown reason
colnames(exp.chdwb) <- gsub("\\.", "-", colnames(exp.chdwb))

# Get a list of unique probes across all matrices
dfs    <- sapply( ls(pattern = "exp"), get)
probes <- sapply(dfs, function(x) c(x[, 1]))
probes <- unique( do.call(c, list(unlist(probes))) )
probes <- probes[order(probes)]
# Add empty rows for missing probes
dfs <- llply(dfs,
             function(x) {
                labels <- probes[which(!probes %in% x[, 1])]
                mat <- matrix(nrow = length(labels), ncol = ncol(x))
                mat[, 1] <- labels
                colnames(mat) <- colnames(x)
                x <- rbind(x, mat)
                x <- x[order(x[, 1]), ]
             })
# Merge dataframes
e.merged      <- data.frame(matrix(nrow = nrow(dfs[[1]])))
e.merged[, 1] <- probes
e.merged[, 1] <- "PROBE_ID"
for (i in 1:length(dfs)) {
    e.merged <- cbind(e.merged, dfs[[i]][, -1])
}
e.merged      <- e.merged[, -1]

#----------------------------### Genotype data ###------------------------------
gen.bsgsmain       <- read("BSGSmain/clean/genotypes/BSGSmain_genotypes.txt")
gen.bsgspilot.lcl  <- read("BSGSpilot/clean/genotypes/BSGSpilot_genotypes-LCL.txt")
gen.bsgspilot.pbmc <- read("BSGSpilot/clean/genotypes/BSGSpilot_genotypes-PBMC.txt")
gen.chdwb          <- read("CHDWB/clean/CHDWB_genotypes.txt")
gen.morocco        <- read("Morocco/clean/genotypes/Morocco_genotypes.txt")

# MuTHER genotypes are across all chromosomes and must be handled individually
gen.muther   <- list()
files    <- list.files(file.path(inpath, "MuTHER/clean/genotypes"))
sequence <- seq(2, length(files), by = 2)
for (i in 1:22) {
   # Extract name of chromosome
    chrom <- strsplit(files[sequence[i]], "-|\\.")[[1]][2]
    temp  <- read(paste("MuTHER/clean/genotypes/", files[sequence[i]], sep = ""))
    gen.muther[[length(gen.muther) + 1]] <- temp
}

# Same applies for CAD
#gen.cad <- list()
#files <- list.files(file.path(inpath, "CAD/clean/genotypes/"))
#for (i in 1:25) {
#    chrom <- strsplit(files[i], "-|\\.")[[1]][2]
#    temp <- read(paste("CAD/clean/genotypes/", files[i], sep = ""))
#    gen.cad[[length(gen.cad) + 1]] <- temp
#}

# Get all genotype matrices in a single list
dfs <- sapply( ls(pattern = "gen"), get )

#-----------------------------## Set column names ##----------------------------
# Store dataset name and sample type in vector
#codes <- list("BSm", "BSp.LCL", "BSp.PBMC", "CAD", "CHD", "Mor", "MuT.fat", "MuT.LCL", "MuT.PBMC")
codes <- list("BSm", "BSp.LCL", "BSp.PBMC", "CHD", "Mor", "MuT.fat", "MuT.LCL", "MuT.PBMC")
e.names  <- unique(gsub("\\..*", "", colnames(e.merged)))
e.names.1 <- colnames(e.merged)
cage.names <- paste("CAGE", formatC(seq(1:length(e.names)), width = 6, format = "d", flag = "0"), sep = "")

map <- cbind(e.names, cage.names)
rows <- lapply(gsub("\\..*", "", colnames(e.merged)), function(x) which(map[, 1] == x))
names <- map[unlist(rows),2]

# Append dataset code to CAGE IDs
from <- 1
sets  <- ls(pattern = "exp")
setnames <- list()
for (i in 1:length(sets)) {
    dataset <- mget(sets[i])
    to <- from + length(dataset[[1]]) - 2 # subtract 2 to account for base 1 indexing and probe ID column
    unlabelled <- names[from:to]
    labelled <- paste(unlabelled, codes[i], sep = "_")
    colnames(e.merged)[from:to] <- labelled
    setnames <- c(setnames, rep(codes[i], (to - from)))
    from <- to + 1 # increment start index for next iteration
}

# TODO: Replace plaintext genotype matrix labels with unique CAGE IDs
sets <- ls(pattern = "gen")
dfs.g <- list()
for (i in 1:length(sets)) {
    dataset <- mget(sets[[i]])
    # CAD and MuTHER contain lists of matrices by chromosome
    if (class(dataset[[1]]) == "list") {
#        next
        rows <- lapply(colnames(dataset[[1]]), function(x) which(map[, 1] == x))
        names <- c("RS_ID", map[unlist(rows),2])
        for (i in 1:length(dataset[[1]])) {
            colnames(dataset[[1]][[i]]) <- names
        }
    # All other datasets are represented by a single matrix
    } else {
        rows <- lapply(colnames(dataset[[1]]), function(x) which(map[, 1] == x))
        names <- c("RS_ID", map[unlist(rows),2])
        colnames(dataset[[1]]) <- names
    }
    dfs.g <- c(dfs.g, dataset)
    rm(dataset, rows, names)
}

# TODO: Replace PLINK genotype labels with unique CAGE IDs

# TODO: Write to file
e.merged <- cbind(probes, e.merged)
colnames(e.merged)[1] <- "PROBE_ID"
write.table(e.merged, "/ibscratch/wrayvisscher/xander/CAGE/data/CAGE_gene_expressions.txt",
            sep = "\t", eol = "\n", row.names = FALSE)
