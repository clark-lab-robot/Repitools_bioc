

QdnaData <- function(counts,regions,design,cnv.offsets=NULL,neutral=NULL) {
	require(edgeR)
# do a slew of checks
	stopifnot( !is.null(counts) | !is.null(regions) | !is.null(design) )
	stopifnot( nrow(design)==ncol(counts) )
    stopifnot( is(regions,"GRanges") )
	stopifnot( length(regions)==nrow(counts) )
	
# populate offsets/neutral if not already
	if( is.null(cnv.offsets) )
	  cnv.offsets <- matrix(1,nrow=nrow(counts),ncol=ncol(counts))
	if( is.null(neutral) )
		neutral <- rep(TRUE,nrow(counts))
					   
# do more checks
	stopifnot( all(dim(counts)==dim(cnv.offsets)) )
    stopifnot( length(neutral)==nrow(counts) )
	
	suppressMessages(d <- DGEList(counts=counts, genes=as.data.frame(regions)[,-c(4:5)]))
			  
# create object
	new("QdnaData",list(DGEList=d,design=as.matrix(design),
		           cnv.offsets=as.matrix(cnv.offsets),neutral=as.logical(neutral),
				   sample.specific.calculated=FALSE))
}


.getQuantile <- function(u,min,quan) {
	top <- max(min,round((1-quan)*sum(!u$w)))
	o <- order(-u$A)[1:top]
	max(min(u$A[o]),max(u$A[u$w]))	
}

getSampleOffsets <- function(obj, ref=1, quantile=0.99, min.n=100, plot.it=FALSE, force=FALSE, ...) {
	if( !is(obj,"QdnaData") )
	    stop("Only operates on 'QdnaData' objects.")
	if( obj$sample.specific.calculated & !force )
	    stop("Sample-specific factors have already been calculated. Set force=TRUE to recalculate.")
	cnt <- obj$DGEList$counts
	nc <- ncol(cnt)
	nf <- rep(0,nc)
	idx <- setdiff(1:nc,ref)
	rs <- obj$DGEList$samples$lib.size
	neu <- obj$neutral
	for(i in idx) {
		map <- maPlot(cnt[neu,ref], cnt[neu,i], normalize=TRUE, plot.it=plot.it, ...)
		q <- .getQuantile(map,min.n,quantile)
		if(plot.it) {
			grid()
			abline(v=q,col="blue")
		}
		nf[i] <- median( map$M[map$A > q] )
		if(plot.it)
			abline(h=nf[i],col="red",lwd=4)
	}
	obj$DGEList$samples$norm.factors <- exp(nf-mean(nf))
	obj$sample.specific.calculated <- TRUE
    obj
}



plotQdnaByCN <- function(obj, cnv.group, idx.ref=1, idx.sam=2, min.n=100, quantile=0.99, ylim=c(-5,5), ...) {
    if( !is(obj,"QdnaData") )
        stop("Only operates on 'QdnaData' objects.")
	
	# effective library size
    els <- obj$DGEList$samples$lib.size*obj$DGEList$samples$norm.factors 
	cnt <- obj$DGEList$counts
	
	tc <- table(cnv.group)
	n <- length(tc)
	f <- rep(NA,n)
	
	layout(matrix(c(1:n,rep(n+1,n)),2,n,byrow=TRUE))


    for(i in 1:n) {
		cn <- names(tc)[i]
	    kk <- cnv.group==cn
	    map <- maPlot(cnt[kk,idx.ref]/els[idx.ref], 
					  cnt[kk,idx.sam]/els[idx.sam], 
		    		  normalize=FALSE, main=cn, ylim=ylim, ...); grid();
		q <- .getQuantile(map,min.n,quantile)
	    abline(v=q,col="blue"); 
	    f[i] <- median(map$M[map$A>q])
	    abline(h=f[i],col="red",lwd=4); 
	}

    f <- exp(f)
    plot(1:n,f,pch=19,ylab="Relative Depth by CNV",type="b",
             xlab="",xaxs="i",xlim=c(1-.5,n+.5))

}


abcdDNA <- function(obj, coef=ncol(obj$design), dispersion=NULL) {
	d <- obj$DGEList
	if( !obj$sample.specific.calculated )
	    message("Sample-specific factors have not been calculated.  First call getSampleOffsets().")
	stopifnot( !is.null(dispersion) )

	o <- outer( rep(1,nrow(d)), getOffset(d)) + log(cn)	
	fit <- glmFit(d, obj$design, offset=o, dispersion=dispersion)
	glmLRT(fit, coef=coef)
}


setCNVOffsets <- function(obj, cnv.offsets) {	
	stopifnot( all(dim(obj$DGEList$counts)==dim(cnv.offsets)) )
	obj$cnv.offsets <- cnv.offsets
	obj
}

