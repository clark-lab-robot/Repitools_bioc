# Get a table of genes with a t-statistic column.
stats <- read.csv("K9AcDiff.csv", stringsAsFactors = FALSE)

# Ensure the genes are already sorted in positional order.
pos <- ifelse(stats$strand == '+', stats$start, stats$end)
pos.order <- order(stats$chr, stats$pos)
stats <- stats[pos.order, ]

# T t-statistic is in column 7 of the data frame. Find all gene clusters that
# have a candidate gene within a window of 5 genes,  with at least 2 candidate
# genes in a window, and a minimum window size of 3. Try cutoffs
# 1, 2, ..., 9, 10, and use the lowest cutoff that gives an FDR < 0.05.
# Look for genes that are consitently positive in the t-statistic.

findClusters(stats, 7, 5, 2, 3, seq(1, 10, 1), trend = "up", n.perm = 2)
