library(Repitools)
library(GenomicFeatures)
library(aroma.affymetrix)

# Load a set of 1.0R arrays and do typical MAT normalisation steps on them.

cdf <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02")
cdfU <- getUniqueCdf(cdf)
celSet <- AffymetrixCelSet$byName("el6", cdf = cdf)
MATNormalise <- MatNormalization(celSet)
cs.norm <- process(MATNormalise)
cs.norm.uniq <- convertToUnique(cs.norm)

# > cs.norm.uniq
# AffymetrixCelSet:
# Name: el6
# Tags: MN,lm,UNQ
# Path: probeData/el6,MN,lm,UNQ/Hs_PromPR_v02
# Platform: Affymetrix
# Chip type: Hs_PromPR_v02,unique
# Number of arrays: 4
# Names: LNCaP_MBD2IP_Elution_6, LNCaP_MBD2IP_Input, PrEC_MBD2IP_Elution_6, PrEC_MBD2IP_Input
# Time period: 2009-11-24 10:23:57 -- 2009-11-24 10:24:14
# Total file size: 166.31MB
# RAM: 0.01MB

# Do smoothing of the array signal.
design <- matrix(c(0, 0, 1, 0, 1, 0, 0, 0), nrow = 4, ncol = 2)
colnames(design) <- c("PrEC Elution 6", "LNCaP Elution 6")
smoothingProcess <- MatSmoothing(cs.norm.uniq, design = design, probeWindow = 300,
                                 tag = "el6example", nProbes = 10)
smoothedScores <- process(smoothingProcess, verbose = -10)

# > smoothedScores
# AffymetrixCelSet:
# Name: el6
# Tags: MN,lm,UNQ,el6example
# Path: probeData/el6,MN,lm,UNQ,el6example/Hs_PromPR_v02
# Platform: Affymetrix
# Chip type: Hs_PromPR_v02,unique
# Number of arrays: 2
# Names: PrEC Elution 6, LNCaP Elution 6
# Time period: 2009-11-24 10:23:57 -- 2009-11-24 10:23:57
# Total file size: 83.16MB
# RAM: 0.01MB

# Get the annotation.
hg18info <- loadFeatures("hg18.txdb")
genes <- as.data.frame(transcripts(hg18info))[, -4]
colnames(genes)[1] <- "chr"

# Get equally spaced intensities around TSSs.
el6scores <- featureScores(smoothedScores, anno = genes, up = 5000, down = 2000, freq = 500,
                           chrs = paste("chr", c(1:22, 'X', 'Y', 'M'), sep = ''))

# > el6scores
# An object of class 'ScoresList'.
# Tables: PrEC Elution 6, LNCaP Elution 6.
# Features:
# GRanges with 37844 ranges and 2 elementMetadata values
#         seqnames               ranges strand   |     tx_id      tx_name
#            <Rle>            <IRanges>  <Rle>   | <integer>     <factor>
#     [1]     chr1     [ 58954,  59871]      +   |      1050 NM_001005484
#     [2]     chr1     [313755, 318444]      +   |      1051    NR_028322
#     [3]     chr1     [313755, 318444]      +   |      1052    NR_028325
#     [4]     chr1     [313755, 318444]      +   |      1056    NR_028327
#     [5]     chr1     [357522, 358460]      +   |      1053 NM_001005221
#     [6]     chr1     [357522, 358460]      +   |      1054 NM_001005224
#     [7]     chr1     [357522, 358460]      +   |      1055 NM_001005277
#     [8]     chr1     [752927, 779603]      +   |      1062    NR_015368
#     [9]     chr1     [850984, 869824]      +   |      1065    NM_152486
#     ...      ...                  ...    ... ...       ...          ...
# [37836]     chrY [24765502, 24770366]      -   |     19536    NR_002195
# [37837]     chrY [25318604, 25369027]      -   |     19546    NM_020420
# [37838]     chrY [25318604, 25369027]      -   |     19547    NM_020364
# [37839]     chrY [25586438, 25607639]      -   |     19548    NM_004678
# [37840]     chrY [25586438, 25607639]      -   |     19549 NM_001002761
# [37841]     chrY [25586438, 25607639]      -   |     19550 NM_001002760
# [37842]     chrY [25739178, 25740308]      -   |     19551    NR_001526
# [37843]     chrY [25739178, 25740308]      -   |     19552    NR_002179
# [37844]     chrY [25739178, 25740308]      -   |     19553    NR_002180
# 
# seqlengths
#            chr1   chr1_random         chr10 ...   chrX_random          chrY
#             NA            NA            NA ...            NA            NA
# Region: 5000 bases up to 2000 bases down.
# Window Width : 500 bases.

######################################
# Example of featureScores on a matrix
######################################

intens.matrix <- extract(smoothedScores, 1:2)
el6scores <- featureScores(intens.matrix, anno = genes, up = 5000, down = 2000, freq = 500,
                           chrs = paste("chr", c(1:22, 'X', 'Y', 'M'), sep = ''))

