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
