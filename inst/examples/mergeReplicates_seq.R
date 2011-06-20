library(Repitools)
library(GenomicRanges)

# Load a saved GRangesList 'grl' into memory.
load("/home/data/Public/Ruike_MeDIPSeq/grl_subset.Rdata")

# > names(grl)
# [1] "HMEC_1"        "HMEC_2"        "HMEC_I"        "MCF7_1_single"
# [5] "MCF7_2"        "MCF7_4"        "MCF7_EMT_1"    "MCF7_EMT_2"   
# [9] "MCF7_I"

grl.types <- c("HMEC", "HMEC", "HMECinput", "MCF7", "MCF7", "MCF7", "MCF7EMT", "MCF7EMT", "MCF7input")

merged <- mergeReplicates(grl, grl.types)

# > names(merged)
# [1] "HMEC"      "HMECinput" "MCF7"      "MCF7EMT"   "MCF7input"
