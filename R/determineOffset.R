determineOffset <- function(x,  quantile=0.998, 
    controlPlot=list(show=FALSE, nsamp=50000, mfrow=c(1,1), ask=FALSE)){
  
    if(class(x) != "BayMethList"){
        stop("x must be a BayMethList object")
    }
    if(nrow(control(x)) == 1){
        stop("There is no control information available. Please specify the
            normalizing offset manually using ``fOffset <- '' ")
    }
    if(is.null(controlPlot$show)){
        controlPlot$show <- FALSE
    }
    if(is.null(controlPlot$nsamp)){
        controlPlot$nsamp <- 50000
    }
    if(is.null(controlPlot$mfrow)){
        controlPlot$mfrow <- c(1,1)
    }
    if(is.null(controlPlot$ask)){
        controlPlot$ask <- FALSE
    }

    ## get the number of samples
    nc <- ncontrol(x)
    ns <- nsampleInterest(x)

    cnames <- colnames(sampleInterest(x))

    message("\n\n\tYou provided ", nc, " SssI samples and ", ns, 
        " samples of interest\n")

    f <- matrix(NA, ncol=ns, nrow=1)
    colnames(f) <- colnames(sampleInterest(x))

    if(controlPlot$show){
        par(mfrow=controlPlot$mfrow, ask=controlPlot$ask)
    }

    ## compute the normalizing offset for each sample of interest
    tmpc <- control(x)[,1]
    for(i in 1:ns){
        if(nc != 1){
        tmpc <- control(x)[,i]
        }
        message("Offset is determined for sample of interest: ", i, "\n")
        tmps <- sampleInterest(x)[,i]
        A <- (log2(tmps) + log2(tmpc))/2
        M <- log2(tmps/tmpc)
        th <- quantile(A, probs=quantile, na.rm=TRUE)
        h <- median(M[A > th], na.rm=TRUE)
        f[1,i] <- 2^h

        if(controlPlot$show){
            ms <- sample(length(tmpc), controlPlot$nsamp)
            tmpcc <- tmpc[ms]
            tmpss <- tmps[ms]
        
            test <- cbind(tmpcc, tmpss)
            rs <- rowSums(test)

            mp <- edgeR:::maPlot(tmpcc[!is.na(rs)], tmpss[!is.na(rs)],
                normalize=FALSE, xlab="A=log2(sampleInterest*SssI)/2",
                ylab="M=log2(sampleInterest/SssI)", pch=19, cex=.4, main=cnames[i]) 

            grid(col="blue")
            abline(v=th, col=2, lwd=2, lty=2)
            abline(h=h, col=2, lwd=2)
            text(0.9*max(mp$A), h-0.3, paste("Offset =",round(f[1,i],3)), 
                col=2, lwd=4, cex=1.2)
        }
    }
    fOffset(x) <- f
    return(x)
}
