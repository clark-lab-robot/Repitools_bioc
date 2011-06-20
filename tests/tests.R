require("Repitools")
require("BSgenome.Hsapiens.UCSC.hg18")
options(warn = -1)

probes <- data.frame(chr = c("chr1", "chr9", "chrY", "chr1", "chr21", "chr6", "chr6", "chr2", "chrX", "chr11"), position = c(10000, 5500, 100000, 11000, 20000000, 500100, 499900, 700000, 9900, 90000), strand = c('+', '+', '+', '+', '-', '-', '-', '+', '-', '+'), stringsAsFactors = FALSE)
genes <- data.frame(chr = c("chr1", "chr9", "chr11", "chr1", "chr11", "chr6", "chr6", "chr22", "chrY", "chr21"), start = c(10000, 7900, 950000, 10500, 74000000, 450000, 5000000, 44000000, 1500, 9800000), end = c(12500, 9500, 1000000, 14500, 75000000, 500000, 9000000, 45000000, 3000, 10000000), strand = c('+', '-', '-', '+', '-', '-', '-', '+', '+', '-'), name = paste("Gene", 1:10), stringsAsFactors = FALSE)

crossMatch <- annotationLookup(probes, genes, 5000, 5000)
correctCrossMatch <- list(indexes = list(`Gene 1` = as.integer(c(1, 4)), `Gene 2` = as.integer(2), `Gene 3` = integer(), `Gene 4` = as.integer(c(1, 4)), `Gene 5` = integer(), `Gene 6` = as.integer(c(7, 6)), `Gene 7` = integer(), `Gene 8` = integer(), `Gene 9` = integer(), `Gene 10` = integer()), offsets = list(`Gene 1` = as.integer(c(0, 1000)), `Gene 2` = as.integer(4000), `Gene 3` = numeric(), `Gene 4` = as.integer(c(-500, 500)), `Gene 5` = numeric(), `Gene 6` = as.integer(c(100, -100)), `Gene 7` = numeric(), `Gene 8` = numeric(), `Gene 9` = numeric(), `Gene 10` = numeric()))
names(correctCrossMatch$offsets$`Gene 1`) <- c(1, 4)
names(correctCrossMatch$offsets$`Gene 2`) <- c(2)
names(correctCrossMatch$offsets$`Gene 4`) <- c(1, 4)
names(correctCrossMatch$offsets$`Gene 6`) <- c(7, 6)
names(correctCrossMatch$offsets$`Gene 8`) <- character()
names(correctCrossMatch$offsets$`Gene 9`) <- character()

if(!isTRUE(all.equal(crossMatch, correctCrossMatch))) 
	stop("Error in annotationLookup function.")
cat("anontationLookup tested fine.\n")

lookupTable <- makeWindowLookupTable(crossMatch$indexes, crossMatch$offsets, starts = seq(-5000, 4900, 100), ends = seq(-4900, 5000, 100))
correctLookupTable <- matrix(NA, nrow = 10, ncol = 100, dimnames = list(genes$names, seq(-4950, 4950, 100)))
correctLookupTable[1, c(50, 51)] <- 1
correctLookupTable[1, c(60, 61)] <- 4
correctLookupTable[2, c(90, 91)] <- 2
correctLookupTable[4, c(45, 46)] <- 1
correctLookupTable[4, c(55, 56)] <- 4
correctLookupTable[6, c(49, 50)] <- 6
correctLookupTable[6, c(51, 52)] <- 7

if(!all(lookupTable == correctLookupTable, na.rm = TRUE))
    stop("Error in makeWindowLookupTable function")
cat("makeWindowLookupTable tested fine.\n")

cpgDensity <- cpgDensityCalc(genes, organism = Hsapiens, window = 500, w.function="linear")
if(!isTRUE(all.equal(cpgDensity, c(5.784, 7.620, 5.828, 2.928, 2.080, 1.252, 0.000, 7.404, 3.928, 0.000))))
    stop("cpgDensityCalc not working for window = 500, scaling = linear")
cpgDensity <- cpgDensityCalc(genes, window = 100, w.function = "log", organism = Hsapiens)
if(!isTRUE(all.equal(round(cpgDensity, 3), c(2.424, 1.882, 1.436, 0.084, 0.379, 0.000, 0.000, 0.263, 1.392, 0.000))))
    stop("cpgDensityCalc not working for window = 100, scaling = log")
cpgDensity <- cpgDensityCalc(genes, window = 1000, w.function = "exp", organism = Hsapiens)
if(!isTRUE(all.equal(round(cpgDensity, 3), c(4.874, 5.835, 4.999, 2.239, 1.567, 0.851, 0.054, 5.589, 3.229,0.062))))
    stop("cpgDensityCalc not working for window = 1000, scaling = exp")
cpgDensity <- cpgDensityCalc(genes, window = 500, w.function = "none", organism = Hsapiens)
if(!isTRUE(all.equal(cpgDensity, c(11, 14, 16, 6, 4, 2, 0, 15, 9, 0))))
    stop("cpgDensityCalc not working for window = 500, scaling = none")
cat("cpgDensityCalc tested fine.\n")

GCpercent <- gcContentCalc(genes, Hsapiens, 500)
if(!isTRUE(all.equal(GCpercent, c(0.504, 0.586, 0.560, 0.470, 0.540, 0.304, 0.356, 0.638, 0.444, 0.388))))
    stop("Error in gcContentCalc function")
cat("gcContentCalc tested fine.\n")

findsCount <- sequenceCalc(genes, Hsapiens, 500, pattern = "AATT")
if(!isTRUE(all.equal(findsCount, c(1, 1, 0, 2, 1, 8, 2, 0, 4, 10))))
    stop("Error in sequenceCalc function counting task")

findsPlaces <- sequenceCalc(genes, Hsapiens, 500, pattern = "AATT", positions = TRUE)
correctPlaces <- list(-62, 181, NULL, c(-140, -98), 231, c(-219, -146, -88, -12, 12, 182, 209, 214), c(-61, 60), NULL, c(-115, -30, 11, 80), c(-238, -228, -202, -189, -177, -106, -21, 148, 158, 238))
if(!isTRUE(all.equal(findsPlaces, correctPlaces)))
    stop("Error in sequenceCalc function positions task")
cat("sequenceCalc tested fine.\n")
cat("All tests passed.\n")
