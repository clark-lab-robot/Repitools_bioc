library(Repitools)
library(aroma.affymetrix)
library(GenomicRanges)

# Get annotations of genes.
hg18info <- makeTranscriptDbFromUCSC(genome = "hg18", tablename = "refGene")
genes <- transcripts(hg18info)

# Normalise arrays.
cdfFile <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02")
cdfFileU <- getUniqueCdf(cdfFile)
celSet <- AffymetrixCelSet$byName("Methylation", cdf = cdfFile)
MATNormalise <- MatNormalization(celSet)
celSetNorm <- process(MATNormalise)
celSetNormUniq <- convertToUnique(celSetNorm)

# Create a design matrix

diff.design <- matrix(c(-1, -1, 1, 1),
                       dimnames = list(getNames(celSetNormalisedU), c("CancerVsNormal")))


# Get statistics results for differences around -2500 to 500 of TSS.

results <- blocksStats(celSetNormUniq, genes, 2500, 500, design = diff.design)
