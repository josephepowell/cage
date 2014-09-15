
#--------------------------------### Setup ###----------------------------------
inpath  <- "/ibscratch/wrayvisscher/xander/CAGE/data/EGCUT/raw"
outpath <- "/ibscratch/wrayvisscher/xander/CAGE/data/EGCUT/clean"

#---------------------------### Expression data ###-----------------------------
exp.raw <- read.table(file = file.path(inpath, "ExpressionData.txt"),
                       sep = "\t",
                    header = TRUE,
               check.names = FALSE)
info    <- read.table(file = file.path(inpath, "EGC_Proteomics_Vcodes_Pcodes_180711.txt"),
                       sep = "\t",
                    header = TRUE)
PROBE_ID <- exp.raw[, 1]
exp   <- exp.raw[, -1]
names <- sapply(colnames(exp), function(x) info[which(info[, 1] == x),2] )
colnames(exp) <- names
exp   <- exp[, order(colnames(exp))]
exp   <- cbind(PROBE_ID, exp)

plink.1 <- snpStats::read.plink(file.path(inpath, "EGCUT_370CNV_ExpSample_290814"))
# Remove CNV genotypes that do not overlap with expression data
gen.1 <- plink.1$genotypes@.Data
gen.1 <- gen.1[which(rownames(gen.1) %in% names), ]
fam.1 <- plink.1$fam
fam.1 <- fam.1[rownames(gen.1), ]


plink.1$genotypes <- plink.1$genotypes[which(rownames(plink.1$genotypes) %in% names), ]
plink.1$fam       <- plink.1$fam[rownames(plink.1$genotypes), ]
# All Omni genotypes overlap with expression data
plink.2 <- snpStats::read.plink(file.path(inpath, "EGCUT_OmniX_ExpSample_290814"))

# Remove expression data for samples with no genotypes
names <- unique( c(rownames(plink.1$genotypes), rownames(plink.2$genotypes)) )
exp   <- exp[, names]
exp   <- cbind(PROBE_ID, exp)

# TODO: Handle duplication of SNPs between genotype runs
# gen.1$ID <- rownames(gen.1)
# gen.2$ID <- rownames(gen.2)
# gen <- plyr::join(gen.2, gen.1, type = "full", by = "ID")
# TODO: Handle missingness of differing SNPs between genotype runs

#----------------------------### Write to file ###------------------------------
write <- function(x, filename) {
  # Convenience function for easily changing output file format.
  #   x: data.frame or matrix object to write.
  #   file: file outpath.
  # Returns: void.
  write.table(x = x,
           file = file.path(outpath, filename),
            sep = "\t",
            eol = "\n",
          quote = FALSE,
      row.names = FALSE)
}

write(exp, "EGCUT_expression_signals.txt")
write.plink(file.base = file.path(outpath, "EGCUT-CNV"),
            snp.major = TRUE,
                 snps = plink.1$genotypes,
         subject.data = plink.1$fam,
             snp.data = plink.1$map)

write.plink(file.base = file.path(outpath, "EGCUT-Omni"),
            snp.major = TRUE,
                 snps = plink.2$genotypes,
         subject.data = plink.2$fam,
             snp.data = plink.2$map)


