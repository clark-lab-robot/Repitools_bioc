empBayes <- function(x, ngroups=100, ncomp=1, maxBins=50000, ncpu=NULL){

    if(!class(x) == "BayMethList"){
        stop("Object must be of class BayMethList.")
    }
    if(ncomp != 1){
        stop("At the moment we do not support ncomp > 1, but it is in progress")
    }

    f <- fOffset(x)
    cpgdens <- cpgDens(x)
    cpgdens_filtered <- cpgdens[!maskEmpBayes(x)]
    len_control <- ncontrol(x)

    ## find the CpG groups
    lastgroup <- quantile(cpgdens_filtered, prob=0.9999)
    ## these are our intervals limits
    cu <- c(seq(min(cpgdens), lastgroup, length.out=ngroups), max(cpgdens))

    ## for the empirical Bayes use only the bins that are not masked out
    sI_filtered <- as.matrix(sampleInterest(x)[!maskEmpBayes(x),])
    co_filtered <- as.matrix(control(x)[!maskEmpBayes(x),])
    if(nrow(f) > 1){
        f_filtered <- f[!maskEmpBayes(x)]
    } else {
        f_filtered <- f
    }

    ## split the bins according to the CpG density groups
    cpggroup <- cut(cpgdens_filtered, cu, include.lowest=TRUE)
    ## these are the group indices
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
        ncpu <- Repitools:::.getCpu(maxCPU=FALSE)
    }
    paramTab <- list()
    paramTab[["CpG groups"]] <- cut(cpgdens, cu, include.lowest=TRUE)
    paramTab[["ncomp"]] <- ncomp

    co_tmp <- co_filtered[,1]
    co_list <- sapply(1:length(red_idx), function(u){co_tmp[red_idx[[u]]]})
    for(j in 1:nsampleInterest(x)){
        
        message("Prior parameters are determined for sample of interest: ", j, "\n")
 
        f_list <- list()
        sI_tmp <- sI_filtered[,j]
        sI_list <-  sapply(1:length(red_idx), function(u){sI_tmp[red_idx[[u]]]})

        if(length(f_filtered[,j]) == 1){
            f_list <- sapply(1:length(red_idx), function(u){rep(f_filtered[,j], length(red_idx[[u]]))})
        } else {
            f_tmp <- f_filtered[,j]
            f_list <- sapply(1:length(red_idx), function(u){f_tmp[red_idx[[u]]]})
        }

        if(len_control != 1){
            co_tmp <- co_filtered[,j]
            co_list <- sapply(1:length(red_idx), function(u){co_tmp[red_idx[[u]]]})
        }
        sfInit(parallel=TRUE,cpus=ncpu)
        sfLibrary("Rsolnp", character.only=TRUE )
        sfLibrary("Repitools", character.only=TRUE )
        sfExport(".marg", namespace="Repitools")
        sfExport(".myoptimize", namespace="Repitools")
        #sfExport("x")
        sfExport("co_list")
        sfExport("sI_list")
        sfExport("f_list")
        sfExport("maxBins")
        sfExport("ngroups")
        sfExport("ncomp")
        sfExport("inds")    
        sfExport("len_control")    

            if(len_control == 1){
                paramTab[[j+2]] <- snowfall:::sfSapply(1:ngroups,
                Repitools:::.myoptimize, 
                sI_list, co_list, f_list, ncomp, maxBins)
            } else {
                paramTab[[j+2]] <- snowfall:::sfSapply(1:ngroups,
                Repitools:::.myoptimize, 
                sI_list, co_list, f_list, ncomp, maxBins)
            }
        sfStop()
        gc()
    }

    priorTab(x) <- paramTab
    return(x)
}
