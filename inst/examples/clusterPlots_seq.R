library(Repitools)
library(GenomicRanges)
library(GenomicFeatures)

# Load a GRangesList of different epigenetic marks in the same cell type.
load("cancerReads.Rdata")

# > names(readsList)
# [1] "H3K27me3" "H3K27me3" "H3K27me3" "H3K27me3" "MBD2IP"   "MBD2IP"   "H3K4me3" 
# [8] "H3K36me3"

# Merge the biological replicates.
merged <- mergeReplicates(readsList, names(readsList))

# > names(merged)
# [1] "H3K27me3" "H3K36me3" "H3K4me3"  "MBD2IP"

# Get annotations of genes.
hg18info <- makeTranscriptDbFromUCSC(genome = "hg18", tablename = "refGene")
genes <- transcripts(hg18info)

# Get the coverage of each mark type at regular intervals 5000 bases either
# side of each TSS.
covs <- featureScores(merged, genes, up = 5000, down = 5000, freq = 1000, s.width = 1000)

# Load pre-computed expression values for the genes.
load("cancerExpr.Rdata")

# Now look at side-by-side heatmaps of all of the marks.
pdf("marksCombo.pdf", height = 10, width = 28)
clusterPlots(covs, function(x) sqrt(x), expr = cancerExpr, plot.type = "heatmap",
             t.name = "Four Marks In Cancer")
dev.off()
