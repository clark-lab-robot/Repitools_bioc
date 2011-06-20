library(Repitools)
library(BSgenome.Hsapiens.UCSC.hg18)

# load a saved GRangesList
load("/home/data/Public/Ruike_MeDIPSeq/grl_subset.Rdata")

pdf("enrich.pdf",10,10)
enrichmentPlot(grl, seq.len=300, organism=Hsapiens, total.lib.size=FALSE)
dev.off()

