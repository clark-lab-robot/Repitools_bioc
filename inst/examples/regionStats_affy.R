library(Repitools)
library(aroma.affymetrix)

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

design <- matrix( c(1,-1,rep(0,length(cs)-2)), ncol=1, dimnames=list(getNames(cs),"elut5_L-P") )

# just get indices of chr7 here
ind <- getCellIndices(cdfU, unit = indexOf(cdfU, "chr7F"), unlist = TRUE, useNames = FALSE)

regs <- regionStats(csMNU, design, ind = ind, probeWindow = 500, verbose = TRUE)
