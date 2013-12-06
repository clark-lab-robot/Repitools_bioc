empBayes <- function(x, ngroups=100, ncomp=1, maxBins=50000, method="beta", controlMethod=list(mode="full", weights=c(0.1, 0.8, 0.1), param=c(1,1)), ncpu=NULL, verbose=FALSE){

    if(!class(x) == "BayMethList"){
        stop("\n\tObject must be of class BayMethList.\n\n")
    }
    if(!(ncomp %in% c(1, 2, 3))){
        stop("\n\tThe number of mixture components ",  
            "should be 1, 2, or 3.\n\n")
    }
    if(!(method %in% c("DBD", "beta"))){
        stop("\n\tMethod must be DBD or beta.\n\n")
    }
    if(method == "DBD"){
        if(is.null(controlMethod$mode)){
            controlMethod$mode <- "full"
        } else {
            if(!is.element(controlMethod$mode, c("full", "fixedWeights", "fixedBeta"))){
                stop("controlMethod$mode must be either \"full\", \"fixedWeights\" or \"fixedBeta\".")
            }
        }
        if(is.null(controlMethod$weights)){
            controlMethod$weights <- c(0.1, 0.8, 0.1)
        } else {
            if(sum(.belongsTo(controlMethod$weights, 0, 1)) != 3){
                    stop("The weights of the mixture must be between 0 and 1")
            }
            if(sum(controlMethod$weights) != 1){
                    stop("The weights of the mixture must sum up to 1")
        }
        }
        if(is.null(controlMethod$param)){
            controlMethod$param <- c(1,1)
        } else {
            if(sum(.belongsTo(controlMethod$param, 1e-20, Inf))!=2){
                stop("The parameters for the beta component must be positive")
            }
        }
    }

    f <- fOffset(x)
    cpgdens <- cpgDens(x)

    ## for the empirical Bayes use only the bins that are not masked out
    cpgdens_filtered <- cpgdens[!maskEmpBayes(x)]
    sI_filtered <- as.matrix(sampleInterest(x)[!maskEmpBayes(x),])
    if(nrow(f) > 1){
        f_filtered <- f[!maskEmpBayes(x), , drop=FALSE]
    } else {
        ## here we have only one offset independent of the bin
        f_filtered <- f
    }

    ## find the CpG groups
    lastgroup <- quantile(cpgdens_filtered, prob=0.9999)
    ## these are our intervals limits
    cu <- c(seq(min(cpgdens), lastgroup, length.out=ngroups), max(cpgdens))

    ## split the bins according to the CpG density groups
    cpggroup <- cut(cpgdens_filtered, cu, include.lowest=TRUE)
    ## get a list of length ngroups, containing for each group the corresponding bin indices
    inds <- split(1:length(cpggroup), cpggroup)

    ## prepare for each CpG group the corresponding bin indices
    red_idx <- list()
    for( i in 1:ngroups){
        ## just use the indices for CpG group i
        idx <-  inds[[i]]
        ## in particular for the low CpG density groups there are many genomic
        ## regions which we do not all need, so use an upper limit maxBins.
        len <- length(idx)
        if(len > maxBins){
            sidx <- sample(1:len, maxBins, replace=TRUE)
        } else {
            sidx <- 1:len
        }
        red_idx[[i]] <- idx[sidx]
    }

    if(is.null(ncpu)){
        ncpu <- .getCpu(maxCPU=FALSE)
    }
    paramTab <- list()
    ## save for all bins the information to which CpG group it belongs
    ## although we don't use it here, but later for deriving the methylation estimates
    paramTab[["CpG groups"]] <- cut(cpgdens, cu, include.lowest=TRUE)
    paramTab[["ncomp"]] <- ncomp
    paramTab[["method"]] <- method

    controlMat <- as.matrix(control(x))
    ## check whether we have control information
    if(nrow(controlMat) != 1){
        co_filtered <- as.matrix(controlMat[!maskEmpBayes(x),])
        co_tmp <- co_filtered[,1]
        # co_list contains for first sample the control values splitted by CpG density group
        co_list <- sapply(1:ngroups, function(u){co_tmp[red_idx[[u]]]})
        len_control <- ncontrol(x)
    } else {
        if(verbose){
            message("\nNOTE: There is no control information that can be taken into account\n")
        }
        len_control <- 1
        co_list <- lapply(1:ngroups, function(u){0})
    }
        
    for(j in 1:nsampleInterest(x)){
        if(verbose){
            message(paste0("Sample ",j, ":\tPrior parameters are determined\t"))
        }
        f_list <- list()
        sI_tmp <- sI_filtered[,j]
        sI_list <-  sapply(1:ngroups, function(u){sI_tmp[red_idx[[u]]]})

        if(length(f_filtered[,j]) == 1){
            f_list <- sapply(1:ngroups, function(u){rep(f_filtered[,j], length(red_idx[[u]]))})
        } else {
            f_tmp <- f_filtered[,j]
            f_list <- sapply(1:ngroups, function(u){f_tmp[red_idx[[u]]]})
        }

        if(len_control > 1){
            co_tmp <- co_filtered[,j]
            co_list <- sapply(1:ngroups, function(u){co_tmp[red_idx[[u]]]})
        }
        if(method=="beta"){
            tmp <- mclapply(1:ngroups, .myoptimize, 
                sI_list, co_list, f_list, ncomp, mc.cores=ncpu)
            paramTab[[j+3]] <- t(do.call(rbind, tmp))
        } else {
            tmp <- mclapply(1:ngroups, .myoptimizeDirac, 
                sI_list, co_list, f_list, controlMethod, mc.cores=ncpu)
            paramTab[[j+3]] <- t(do.call(rbind, tmp))
        }

        gc()

        if(controlMethod$mode=="fixedBeta"){
            paramTab[[j+3]] <- rbind(paramTab[[j+3]][1:2,], controlMethod$param[1], controlMethod$param[2], paramTab[[j+3]][3:4,])
        }
        if(controlMethod$mode=="fixedWeights"){
            paramTab[[j+3]] <- rbind(paramTab[[j+3]][1:4,], controlMethod$weights[2], controlMethod$weights[3])
        }
    }
    names(paramTab) <- c("CpG groups", "ncomp", "method", colnames(sampleInterest(x)))

    priorTab(x) <- paramTab
    return(x)
}
