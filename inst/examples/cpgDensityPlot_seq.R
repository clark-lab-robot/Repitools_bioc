library(Repitools)
library(BSgenome.Hsapiens.UCSC.hg18)

# load a saved GRangesList
load("/home/data/Public/Ruike_MeDIPSeq/grl_subset.Rdata")

pdf("CDP.pdf",10,10)
cpgDensityPlot(grl, cols = c("black", "red", "green","blue","orange"), 
               xlim = c(0,30), wFunction = "none", organism = Hsapiens,
               seqLen = 300, lwd = 4, verbose=TRUE)
dev.off()

