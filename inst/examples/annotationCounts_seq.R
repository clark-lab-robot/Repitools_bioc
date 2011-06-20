library(Repitools)
library(GenomicRanges)
library(GenomicFeatures)

# Load a saved GRangesList 'grl' into memory.
load("/home/data/Public/Ruike_MeDIPSeq/grl_subset.Rdata")

# Get annotation of genes.
hg18info <- makeTranscriptDbFromUCSC(genome = "hg18", tablename = "refGene")
genes <- transcripts(hg18info)

# Make counts 2500 bases upstream and 500 bases downstream of each gene.
# Use extended short reads, that are extended to 300 bases.
counts <- annotationCounts(grl, genes, 2500, 500, 300)

# Make counts within each gene.
# Use extended short reads, that are extended to 300 bases.
counts <- annotationBlocksCounts(grl, genes, 300)
