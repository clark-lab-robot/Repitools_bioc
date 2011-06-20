library(Repitools)
library(GenomicRanges)
library(GenomicFeatures)

# load saved GRangesList
load("/home/data/Public/Ruike_MeDIPSeq/grl_subset.Rdata")

#> class(grl)
#[1] "GRangesList"
#attr(,"package")
#[1] "GenomicRanges"
#> names(grl)
#[1] "HMEC_1"        "HMEC_2"        "HMEC_I"        "MCF7_1_single"
#[5] "MCF7_2"        "MCF7_4"        "MCF7_EMT_1"    "MCF7_EMT_2"   
#[9] "MCF7_I"      

# grab annotation
refseq = makeTranscriptDbFromUCSC('hg18', tablename='refGene')
refseqGR = transcripts(refseq)
keep <- seqnames(refseqGR) %in% paste("chr", c(1:22,"X","Y"),sep="")
refseqGR <- refseqGR[keep]
refseqDF <- as.data.frame(refseqGR)
w <- !duplicated(refseqDF[,c(1:3,5)])
refseqDF <- refseqDF[w,]  # remove duplicates
colnames(refseqDF)[ colnames(refseqDF)=="seqnames" ] <- "chr"

# expression data vector/matrix that matches refseqGR (replace with real data)
expr <- matrix(rnorm( nrow(refseqDF)  ), ncol=1)

# create design matrix 
mm <- matrix(0,nrow=length(grl),ncol=3)
rownames(mm) <- names(grl)
colnames(mm) <- c("HMEC","MCF7","MCF7-HMEC")
mm[,1] <- c(1,1,0, 0, 0, 1,0,0,0)
mm[,2] <- c(0,0,0, 1, 1, 1,0,0,0)
mm[,3] <- mm[,1]-mm[,2]

#binPlots(grl, refseqGR, design=mm[,1,drop=FALSE], verbose=TRUE, seqLen=300, nbins=50,
#           ordering=expr, ordLabel="expression",plotType="heatmap")

binPlots(grl, refseqDF, design=mm[,1,drop=FALSE], verbose=TRUE, seqLen=300, nbins=10,
           ordering=expr, ordLabel="expression",plotType="line")

