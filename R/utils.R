# Common short functions used by multiple function in Repitools.

setGeneric(".validate", signature = c("anno"), function(anno, up, down)
                                        {standardGeneric(".validate")})

setMethod(".validate", "GRanges", function(anno, up, down)
{
    if(is.null(up))
        stop("Mandatory argument 'up' not given.")
    if(is.null(down))
        stop("Mandatory argument 'down' not given.")

    str <- strand(anno)

    if(-up > down)
	stop("'up' and 'down' window boundaries cross over each other.\n")
	
    if(any(str == '*') && any(str %in% c('+', '-')))
	stop("Annotation contains mixed feature types.")

    # For unstranded features.
    if(any(str %in% '*') && up != down)
	stop("Different upstream and downstream window distances don't make
              sense for unstranded features.\n")
    NULL
})

.getNames <- function(annoGR)
{
    if("name" %in% names(elementMetadata(annoGR)))
        return(elementMetadata(annoGR)[, "name"])
    else
        return(1:length(annoGR))
}

.getCpu <- function(maxCPU=FALSE){

    numCores <- detectCores()
    if(maxCPU){
        message("\n\tInformation: The program will take advantage of ", 
            numCores, " CPUs\n\tIf you would like to change this ",
            "please cancel and set explicitly the argument 'ncpu'\n\n")
        return(numCores)
    } else {
        if(numCores <= 4){
            numCoresRed <- numCores/2
            message("\n\tInformation: The program will take advantage of ", 
                numCoresRed, " CPUs from total ", numCores, "\n\tIf you would",
                " like to change this please cancel and set explicitly the ",
                "argument 'ncpus'\n\n")
        } else {
            numCoresRed <- floor(numCores/3*2)
            message("\n\tInformation: The program will take advantage of ", 
                numCoresRed, " CPUs from total ", numCores, "\n\tIf you would",
                " like to change this please cancel and set explicitly the ",
                "argument 'ncpus'\n\n")
        }
        return(numCoresRed)
    }
}

## here we assume to get ONE sample of interest, 
## ONE control, and ONE value or vector for f.
.myoptimize <- function(i, sample, control, f, ncomp, maxBins){    

    ## lower bound for parameter values (zero not possible)
    eps <- 10^(-3)

    sum.tmp <- sample[[i]]+control[[i]]
    len <- length(sample[[i]])
    if(sum(is.na(sum.tmp))){
        w.na <- which(is.na(sum.tmp))
        message("\nCpG group", i, ":\n\tWe remove ", length(w.na), " bins out of ", len,  " in the empirical Bayes\n")
        message("\tThe reason is that there are NA's in one of the sample reads")
        print(head(cbind("sampleInterest"=sample[[i]][w.na], "control"=control[[i]][w.na])))
        cat("\n")
        f[[i]] <- f[[i]][-w.na]
        sample[[i]] <- sample[[i]][-w.na]
        if(length(control[[i]]) != 1){
            control[[i]] <- control[[i]][-w.na]
        }
    }


    ## solnp is in particular better for mixtures of beta priors as you 
    ## can give equality constraints so that the weights sum to one
    ## (not relevant now) 
    paramVec <- Rsolnp:::solnp(c(2.14,1.34), fun=Repitools:::.marg, 
        eqfun=NULL, eqB=NULL, LB=c(eps, eps), 
        UB=c(Inf, Inf), cons=f[[i]], 
        y1=sample[[i]], y2=control[[i]], ncomp=ncomp)$pars

# system.time(Rsolnp:::solnp(c(2.14,1.34), fun=Repitools:::.marg, 
#         eqfun=NULL, eqB=NULL, LB=c(eps, eps), 
#         UB=c(Inf, Inf), cons=f[[i]], 
#         y1=sample[[i]], y2=control[[i]], ncomp=ncomp)$pars)
# 
# system.time(
# optim(c(2.14,1.34), fn=Repitools:::.marg, method="L-BFGS-B",
#         lower=c(eps, eps), 
#         upper=c(Inf, Inf), cons=f[[i]], 
#         y1=sample[[i]], y2=control[[i]], ncomp=ncomp)$par)
# system.time(
# nlminb(start=c(2.14,1.34), objective=Repitools:::.marg, cons=f[[i]], 
#         y1=sample[[i]], y2=control[[i]], ncomp=ncomp,
#         lower=c(eps, eps), 
#         upper=c(Inf, Inf))$par)


    return(paramVec)
}

.marg <- function(arg, cons, y1, y2, ncomp){

    if(!(ncomp %in% c(1))){
        stop("\n\tThe proper support of mixtures of beta ",  
            "distributions is in progress.\n\n")
    }
    no_control <- length(y2) == 1


    ## parameters for lambda
    al <- arg[1]
    bl <- arg[2]

    ## uniform prior (alternatives in progress)
    w <- 1
    a <- 1
    b <- 1

    len <- length(y1)
    if(no_control){
        denom <- bl + cons
    } else {
        denom <- bl + 1 + cons
    }
    ## z-parameter
    argZ <- cons/denom
    ## a-parameter
    argA <- y1 + y2 + al

    ## evertyhing before the mixture is part 1 (do it on log-scale)
    part1 <- lgamma(argA) - (lgamma(al) + lgamma(y1 + 1) + lgamma(y2 + 1)) + 
        al*log(bl/denom) + y1*log(argZ) + y2*log(1/denom)
    ## replicate the z parameter to have one for each observation
    if(length(argZ) != len){
        argZ <- rep(argZ, len)
    }
    part2 <- 0
    ## go over the number of mixture components
    for(i in 1:length(a)){
        argB <- rep(b[i], len)
        argC <- y1 + a[i] + b[i]
        part2 <- part2 + exp(log(w[i]) + 
            (lgamma(a[i] + b[i]) + lgamma(y1 + a[i]))-
            (lgamma(a[i]) + lgamma(y1 + a[i] + b[i])) +
            log(hyperg2F1_vec(a=argA, b=argB, c=argC, z=argZ)))
    }
    s_tmp <- part1 + log(part2)

    ## check for problematic bins (usually with extreme high control
    ## values and low to zero sample of interest counts.
    w_na <- which(is.na(s_tmp))
    w_inf <- which(is.infinite(s_tmp))
    w_rm <- c(w_na, w_inf)
    
    if(length(w_rm) > 0){
        message("We remove ", length(w_rm), " bins out of ", len,  " in the empirical Bayes\n")
        message("The reason might be extreme read density in one of the samples\n")
        print(head(cbind("sampleInterest"=y1[w_rm], "control"=y2[w_rm])))
        s_tmp <- s_tmp[-w_rm]
    }

    res <- sum(s_tmp)
    
    return(-res)
}


## logit function
.logit <- function(x){
    return(log(x/{1-x}))
}

## first derivative of the logit function
.dlogit <- function(x){
    return(1/{x*{1-x}})
}

## inverse logit function
.invlogit <- function(x){
    return(1/{1+exp(-x)})
}

## get quantile or HPD based intervals
.getcredible <-  function(u, method, level, nmarg, y1, y2, al, bl, W, cons, my_sd){

    if(!is.element(method, c("quantile", "HPD"))){
        stop("\n\t Argument 'method' must be either equal to 'quantile' or 'HPD'!\n")
    }
    names <- c(rev(paste("lb_", level, sep="")), paste("ub_", level, sep=""))

    # find the mode
    my_max <- optimize(Repitools:::.mydmarginal, interval=c(0,1), y1=y1[u], y2=y2[u], al=al[u], bl=bl[u], W=W[u], cons=cons[u], maximum=T)$maximum
    # find support points which represent the density 
    # (take nmarg points: nmarg/3 equally spaced between 0 and 1
    #  nmarg/3*2 equally spaced in the higher probability mass, i.e
    # (max(0, mode - 2*sd), min(mode+2*sd, 1))
    n1 <- ceiling(nmarg/3)
    eps <- 1e-8
    x <- sort(c(seq(eps, 1-eps, length.out=n1), seq(max(my_max - 2*my_sd,eps),min(1-eps,my_max + 2*my_sd), 
            length.out=nmarg-n1)))
    # get the density
    y <- Repitools:::.mydmarginal(x, y1=y1[u], y2=y2[u], al=al[u], bl=bl[u], W=W[u], cons=cons[u])

    if(method=="quantile"){
        ll <- (1-level)/2
        q <- c(rev(ll), 1-ll)
        my.ci <- Repitools:::.myquantile( q=q, cbind(x,y))
        names(my.ci) <- names
        return(my.ci)
    }

    if(method=="HPD"){
        if(my_max > 1e-4 & my_max < 0.9999){
            my.hpd <- Repitools:::.myhpd(level, cbind(x,y))
            names(my.hpd) <- names
            if(sum(my.hpd < 0 | my.hpd > 1) > 0){
                cat("\n\n\t WARNING:  Something is wrong with the HPD intervals!\n\tThey are outside (0,1)!\n")
            }
            return(my.hpd)
        }
        
        if(my_max <= 1e-4){
            #cat("\n\tThe lower bound of the HPD interval is zero => Take the", level*100, "% quantile as upper bound\n\n")
            qn <- Repitools:::.myquantile(cbind(x,y), q=level)
            my.hpd <- c(rep(0, length(level)), qn)
            names(my.hpd) <- names
            return(my.hpd)
        } 
        if(my_max >= 0.9999){
            #cat("\n\tThe upper bound of the HPD interval is one => Take the", (1-level)*100, "% quantile as lower bound\n\n")
            qn <- rev(Repitools:::.myquantile(cbind(x,y), q=1-level))
            my.hpd <- c(qn, rep(1, length(level)))
            names(my.hpd) <- names
            return(my.hpd) 
        }
    }
}


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
.myquantile <- function( q, den){

    x <- den[,1]
    y <- den[,2]

    # approximate the cumulative distribution using the trapezoidal rule
    x.tmp <- x[-1] - x[-length(x)]
    y.tmp <- (y[-1] + y[-length(y)])/2

    dn <- cumsum(x.tmp*y.tmp)
    dn <- dn/dn[length(dn)]

    qn <- c()
    for(i in 1:length(q)){
        qn[i] <- x[which(dn >= q[i])[1]]
    }
    return(qn)
}



# for uniform prior for methylation level
.mydmarginal <- function(mu, y1, y2, al, bl, W, cons, log=F){

   if(log){
        if(mu < 0 || mu > 1){
            return(-Inf)
        }
        return(log( mu^{y1}/W* (1- (cons*(1-mu))/(bl + 1 +cons))^{-(al + y1 + y2)}))
    } else {
        if(mu < 0 || mu > 1){
            return(0)
        }
        return(mu^{y1}/W* (1- (cons*(1-mu))/(bl + 1 +cons))^{-(al + y1 + y2)})
    }
}

.mysmarginal <- function (marginal, log = FALSE, extrapolate = 0, keep.type = FALSE) 
{
    is.mat = is.matrix(marginal)
    m = Repitools:::.mymarginalfix(marginal)
    r = range(m$x)
    r = r[2] - r[1]
    ans = spline(m$x, log(m$y), xmin = min(m$x) - extrapolate * 
        r, xmax = max(m$x) + extrapolate * r, n = 10 * length(m$x))
    if (!log) {
        ans$y = exp(ans$y)
    }
    if (is.mat && keep.type) {
        return(cbind(ans$x, ans$y))
    }
    else {
        return(ans)
    }
}

.mysfmarginal <- function (marginal) 
{
    m = Repitools:::.mymarginalfix(marginal)
    r = range(m$x)
    return(list(range = r, fun = splinefun(m$x, log(m$y))))
}


.mymarginalfix <- function (marginal) 
{
    if (is.matrix(marginal)) {
        i = (marginal[, 2] > 0) & (abs(marginal[, 2]/max(marginal[,2])) > sqrt(.Machine$double.eps))
        m = list(x = marginal[i, 1], y = marginal[i, 2])
    }
    else {
        i = (marginal$y > 0) & (abs(marginal$y/max(marginal$y)) > sqrt(.Machine$double.eps))
        m = list(x = marginal$x[i], y = marginal$y[i])
    }
    return(m)
}

.myhpd <- function (p, marginal, len = 1024) 
{
    f = Repitools:::.mysfmarginal(Repitools:::.mysmarginal(marginal))
    xx = seq(f$range[1], f$range[2], length = len)
    d = cumsum(exp(f$fun(xx)))
    d = d/d[length(d)]
    eps = .Machine$double.eps * 1000
    for (val in c(0, 1)) {
        is.val = which(abs(d - val) <= eps)
        if (length(is.val) > 1) {
            is.val = is.val[-1]
            d = d[-is.val]
            xx = xx[-is.val]
        }
    }
    fq = splinefun(d, xx, method = "monoH.FC")
    np = length(p)
    pp = 1 - pmin(pmax(p, rep(0, np)), rep(1, np))
    f = function(x, posterior.icdf, conf) {
        return(posterior.icdf(1 - conf + x) - posterior.icdf(x))
    }
    tol = sqrt(.Machine$double.eps)
    result = matrix(NA, np, 2)
    for (i in 1:np) {
        out = optimize(f, c(0, pp[i]), posterior.icdf = fq, conf = pp[i], 
            tol = tol)
        result[i, ] = c(fq(out$minimum), fq(1 - pp[i] + out$minimum))
    }
    result <- c(rev(result[,1]), result[,2])
    return(result)
}
