library(Repitools)
library(aroma.affymetrix)

# set to local directory where affy files are pre-preprocessed
setwd("~/projects/microarray/mbd2")

# assumes appropriate files are at annotationData/chipTypes/Hs_PromPR_v02/
cdf <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02",verbose=-20)
cdfU <- getUniqueCdf(cdf,verbose=-20)

# assumes appropriate files are at rawData/mbd2/Hs_PromPR_v02/
cs <- AffymetrixCelSet$byName("mbd2",cdf=cdf,verbose=-20)
mn <- MatNormalization(cs)
csMN <- process(mn,verbose=-50)
csMNU <- convertToUnique(csMN,verbose=-20)

#> getNames(cs)
# [1] "elution_5_LNCaP"     "elution_5_PrEC"      "elution_6_LNCaP"
# [4] "elution_6_PrEC"      "LNCaP_input"         "Lncap1_MeDNA_Input1"
# [7] "Lncap1_MeDNA_IP1"    "Lncap1_MeDNA_IP2"    "Lncap2_MeDNA_Input1"
#[10] "Lncap2_MeDNA_IP1"    "Lncap2_MeDNA_IP2"    "PrEC_input"
#[13] "Prec1_MeDNA_Input1"  "Prec1_MeDNA_IP1"     "Prec1_MeDNA_IP2"
#[16] "Prec2_MeDNA_Input1"  "Prec2_MeDNA_IP1"     "Prec2_MeDNA_IP2"
#[19] "prot_k_LNCaP"

des <- matrix( c(0,0,1,-1,rep(0,length(cs)-4)), ncol=1, dimnames=list(getNames(cs),"elut5_L-P") )

annoFile <- system.file("data","chr21genes.csv", package="Repitools")
annoDF <- read.csv(annoFile)

ms <- MatSmoothing(csMNU, design = des, probeWindow = 300, 
                   tag = "300bp_smoothing", nProbes = 10)
csTS <- process(ms, units = NULL, verbose=TRUE)

binPlots(csTS, coordinatesTable=annoDF, verbose=TRUE, nbins=10,
           ordering=expr, ordLabel="expression",plotType="line")

