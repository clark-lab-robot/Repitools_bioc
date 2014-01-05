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

    if (.Platform$OS.type == "windows"){
        message("\n\tCAUTION: Parallel computing using mclapply is not possible on Windows, we therefore set ncpu=1\n\n!")
        return(1)
    }

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
            numCoresRed <- floor(numCores/3)
            message("\n\tInformation: The program will take advantage of ", 
                numCoresRed, " CPUs from total ", numCores, "\n\tIf you would",
                " like to change this please cancel and set explicitly the ",
                "argument 'ncpus'\n\n")
        }
        return(numCoresRed)
    }
}



.myoptimizeDirac <- function(i, sample, control, f, controlMethod){

    ## lower bound for parameter values (zero not possible)
    eps <- 10^(-3)

    ## doesn't matter if control[[i]] = 0 (R does automatically enlarge it to a vector)
    sum.tmp <- sample[[i]]+control[[i]]
    ## get the number of bins included in this group
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

    if(controlMethod$mode == "full"){
        lb <- c(eps, eps, eps, eps, eps, eps)
        ub <- c(Inf, Inf, Inf, Inf, 1-eps, 1-eps)

        arg <- c(2.14, 1.3, 0.7,0.7, 0.1, 0.8)
        paramVec <- solnp(arg, fun=.margDirac, 
                    ineqfun=.ineqn, ineqLB=eps, ineqUB=1-eps, 
                    LB=lb, UB=ub,
                    y1=sample[[i]], y2=control[[i]], 
                    cons=f[[i]], control=list(trace=F))$pars
    }
    if(controlMethod$mode == "fixedWeights"){
        lb <- c(eps, eps, eps, eps)
        ub <- c(Inf, Inf, Inf, Inf)

        arg <- c(2.14, 1.3, 0.7,0.7)
        paramVec <- solnp(arg, fun=.margDirac_fW, 
                    LB=lb, UB=ub,
                    y1=sample[[i]], y2=control[[i]], 
                    cons=f[[i]], weights=controlMethod$weights, control=list(trace=F))$pars
    }
    if(controlMethod$mode == "fixedBeta"){
        lb <- c(eps, eps,  eps, eps)
        ub <- c(Inf, Inf, 1-eps, 1-eps)

        arg <- c(2.14, 1.3, 0.1, 0.8)
        paramVec <- solnp(arg, fun=.margDirac_fB, 
                    ineqfun=.ineqn2, ineqLB=eps, ineqUB=1-eps, 
                    LB=lb, UB=ub, 
                    y1=sample[[i]], y2=control[[i]], 
                    cons=f[[i]], param=controlMethod$param, control=list(trace=F))$pars
    }
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

.ineqn <- function(arg, y1, y2, cons){
    
    z <- arg[5] + arg[6] 

    return(z)
}


.ineqn2 <- function(arg, y1, y2, cons, param){
    
    z <- arg[3] + arg[4] 

    return(z)
}

.margDirac <- function(arg, y1, y2, cons){


    no_control <- length(y2) == 1

    ## parameters for lambda
    al <- arg[1]
    bl <- arg[2]
    ## parameters for beta
    a <- arg[3]
    b <- arg[4]
    ## weights are fixed (for beta component and one point-mass)
    w2 <- arg[5]
    w3 <- arg[6]

    len <- length(y1)
    if(no_control){
        denom <- bl + cons
    } else {
        denom <- bl + 1 + cons
    }

    ## z-parameter
    argZ <- cons/denom
    ## replicate the z parameter to have one for each observation
    if(length(argZ) != len){
        argZ <- rep(argZ, len)
    }
    ## a-parameter
    argA <- y1 + y2 + al

    ## evertyhing before the mixture is part 1 (do it on log-scale)
    part1 <- lgamma(argA) - (lgamma(al) + lgamma(y1 + 1) + lgamma(y2 + 1)) + 
        al*log(bl/denom) + y1*log(argZ) + y2*log(1/denom)

    argB <- rep(b, len)
    argC <- y1 + a + b

    part2 <- w3 + exp(log(w2) + 
        (lgamma(a + b) + lgamma(y1 + a))-
        (lgamma(a) + lgamma(y1 + a + b)) +
        log(hyperg2F1_vec(a=argA, b=argB, c=argC, z=argZ)))
    s_tmp <- part1 + log(part2)

    ## check for problematic bins (usually with extreme high control
    ## values and low to zero sample of interest counts.
    w_na <- which(is.na(s_tmp))
    w_inf <- which(is.infinite(s_tmp))
    w_rm <- c(w_na, w_inf)
    
    if(length(w_rm) > 0){
#         message("We remove ", length(w_rm), " bins out of ", len,  " in the empirical Bayes\n")
#         message("The reason might be extreme read density in one of the samples\n")
#         print(head(cbind("sampleInterest"=y1[w_rm], "control"=y2[w_rm])))
        s_tmp <- s_tmp[-w_rm]
    }

    res <- sum(s_tmp)
    
    return(-res)
}


.margDirac_fW <- function(arg, y1, y2, cons, weights){


    no_control <- length(y2) == 1

    ## parameters for lambda
    al <- arg[1]
    bl <- arg[2]
    ## parameters for beta
    a <- arg[3]
    b <- arg[4]
    ## weights are fixed (for beta component and one point-mass)
    w2 <- weights[2]
    w3 <- weights[3]

    len <- length(y1)
    if(no_control){
        denom <- bl + cons
    } else {
        denom <- bl + 1 + cons
    }

    ## z-parameter
    argZ <- cons/denom
    ## replicate the z parameter to have one for each observation
    if(length(argZ) != len){
        argZ <- rep(argZ, len)
    }
    ## a-parameter
    argA <- y1 + y2 + al

    ## evertyhing before the mixture is part 1 (do it on log-scale)
    part1 <- lgamma(argA) - (lgamma(al) + lgamma(y1 + 1) + lgamma(y2 + 1)) + 
        al*log(bl/denom) + y1*log(argZ) + y2*log(1/denom)

    argB <- rep(b, len)
    argC <- y1 + a + b

    part2 <- w3 + exp(log(w2) + 
        (lgamma(a + b) + lgamma(y1 + a))-
        (lgamma(a) + lgamma(y1 + a + b)) +
        log(hyperg2F1_vec(a=argA, b=argB, c=argC, z=argZ)))
    s_tmp <- part1 + log(part2)

    ## check for problematic bins (usually with extreme high control
    ## values and low to zero sample of interest counts.
    w_na <- which(is.na(s_tmp))
    w_inf <- which(is.infinite(s_tmp))
    w_rm <- c(w_na, w_inf)
    
    if(length(w_rm) > 0){
#         message("We remove ", length(w_rm), " bins out of ", len,  " in the empirical Bayes\n")
#         message("The reason might be extreme read density in one of the samples\n")
#         print(head(cbind("sampleInterest"=y1[w_rm], "control"=y2[w_rm])))
        s_tmp <- s_tmp[-w_rm]
    }

    res <- sum(s_tmp)
    
    return(-res)
}



.margDirac_fB <- function(arg, y1, y2, cons, param){


    no_control <- length(y2) == 1

    ## parameters for lambda
    al <- arg[1]
    bl <- arg[2]
    ## parameters for beta
    a <- param[1]
    b <- param[2]
    ## weights are fixed (for beta component and one point-mass)
    w2 <- arg[3]
    w3 <- arg[4]

    len <- length(y1)
    if(no_control){
        denom <- bl + cons
    } else {
        denom <- bl + 1 + cons
    }

    ## z-parameter
    argZ <- cons/denom
    ## replicate the z parameter to have one for each observation
    if(length(argZ) != len){
        argZ <- rep(argZ, len)
    }
    ## a-parameter
    argA <- y1 + y2 + al

    ## evertyhing before the mixture is part 1 (do it on log-scale)
    part1 <- lgamma(argA) - (lgamma(al) + lgamma(y1 + 1) + lgamma(y2 + 1)) + 
        al*log(bl/denom) + y1*log(argZ) + y2*log(1/denom)

    argB <- rep(b, len)
    argC <- y1 + a + b

    part2 <- w3 + exp(log(w2) + 
        (lgamma(a + b) + lgamma(y1 + a))-
        (lgamma(a) + lgamma(y1 + a + b)) +
        log(hyperg2F1_vec(a=argA, b=argB, c=argC, z=argZ)))
    s_tmp <- part1 + log(part2)

    ## check for problematic bins (usually with extreme high control
    ## values and low to zero sample of interest counts.
    w_na <- which(is.na(s_tmp))
    w_inf <- which(is.infinite(s_tmp))
    w_rm <- c(w_na, w_inf)
    
    if(length(w_rm) > 0){
#         message("We remove ", length(w_rm), " bins out of ", len,  " in the empirical Bayes\n")
#         message("The reason might be extreme read density in one of the samples\n")
#         print(head(cbind("sampleInterest"=y1[w_rm], "control"=y2[w_rm])))
        s_tmp <- s_tmp[-w_rm]
    }

    res <- sum(s_tmp)
    
    return(-res)
}



## here we assume to get ONE sample of interest, 
## ONE control, and ONE value or vector for f.
.myoptimize <- function(i, sample, control, f, ncomp){    

    ## lower bound for parameter values (zero not possible)
    eps <- 10^(-3)

    ## doesn't matter if control[[i]] = 0 (R does automatically enlarge it to a vector)
    sum.tmp <- sample[[i]]+control[[i]]
    ## get the number of bins included in this group
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

    if(ncomp > 1){
        if(ncomp==2){ 
            # a_gamma, b_gamma, w_1, w_2, a_1, b_1, a_2, b_2
            lb <-  c(eps, eps, eps,eps, rep(eps, 4))
            ub <- c(Inf, Inf, 1-eps, 1-eps, rep(Inf,4))
            arg <- c(2.14, 1.34, 0.5, 0.5, 1, 1, 1,1)
            paramVec <- solnp(arg, fun=.marg, 
                eqfun=.eqn2, eqB=1, LB=lb, UB=ub, cons=f[[i]], 
                y1=sample[[i]], y2=control[[i]], ncomp=ncomp, control=list(trace=F))$pars
        }
        if(ncomp==3){
            lb <- c(eps, eps, eps, eps, eps, rep(eps, 3))
            ub <- c(Inf, Inf, 1-eps, 1-eps, 1-eps, rep(Inf,3))
            arg <- c(2.14, 1.34, 0.2, 0.2, 0.6, 1, 1, 1)
            paramVec <- solnp(arg, fun=.marg, 
                eqfun=.eqn3, eqB=1, LB=lb, 
                UB=ub, cons=f[[i]], 
                y1=sample[[i]], y2=control[[i]], ncomp=ncomp, control=list(trace=F))$pars            
        }
    } else {
        ## solnp is in particular better for mixtures of beta priors as you 
        ## can give equality constraints so that the weights sum to one
        paramVec <- solnp(c(2.14,1.34), fun=.marg, 
            eqfun=NULL, eqB=NULL, LB=c(eps, eps), 
            UB=c(Inf, Inf), cons=f[[i]], 
            y1=sample[[i]], y2=control[[i]], ncomp=ncomp, control=list(trace=F))$pars
    }
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

    if(!(ncomp %in% c(1, 2, 3))){
        stop("\n\tThe number of mixture components ",  
            "should be 1, 2, or 3.\n\n")
    }
    no_control <- length(y2) == 1


    ## parameters for lambda
    al <- arg[1]
    bl <- arg[2]

    ## weights
    if(ncomp==3){
        w <- c(arg[3], arg[4], arg[5])
        a <- c(arg[6], arg[8], arg[7])
        b <- c(arg[7], arg[8], arg[6])
    } else {
        if(ncomp==2){
            w <- c(arg[3], arg[4])
            a <- c(arg[5], arg[7])
            b <- c(arg[6], arg[8])
        } else { # assume a uniform prior
            w <- 1
            a <- 1
            b <- 1
        }
    }

    len <- length(y1)
    if(no_control){
        denom <- bl + cons
    } else {
        denom <- bl + 1 + cons
    }
    ## z-parameter
    argZ <- cons/denom
    ## replicate the z parameter to have one for each observation
    if(length(argZ) != len){
        argZ <- rep(argZ, len)
    }
    ## a-parameter
    argA <- y1 + y2 + al

    ## evertyhing before the mixture is part 1 (do it on log-scale)
    part1 <- lgamma(argA) - (lgamma(al) + lgamma(y1 + 1) + lgamma(y2 + 1)) + 
        al*log(bl/denom) + y1*log(argZ) + y2*log(1/denom)
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
#         message("We remove ", length(w_rm), " bins out of ", len,  " in the empirical Bayes\n")
#         message("The reason might be extreme read density in one of the samples\n")
#         print(head(cbind("sampleInterest"=y1[w_rm], "control"=y2[w_rm])))
        s_tmp <- s_tmp[-w_rm]
    }

    res <- sum(s_tmp)
    
    return(-res)
}



## eqn3(...)
## Equality constraint function used in the optimization of the empirical Bayes step.
## For a mixture of three beta distributions it should gurantee that their
## weights sum up to one
##
## Arguments:
##############
## arg - parameter vector
## cons - multiplicative offset
## y1 - vector containing the read counts for the sample of interest
## y2 - vector containing the read counts for the SssI control
## ncomp - number of beta mixture components (using a uniform prior ncomp=1)
##
## Value:
##############
##
## The function returns the sum of the three mixture weights
.eqn3 <- function(arg, cons, y1, y2, ncomp){
    
    z <- arg[3] + arg[4] + arg[5]

    return(z)
}

## eqn2(...)
## Equality constraint function used in the optimization of the empirical Bayes step.
## For a mixture of two beta distributions it should guarantee that their
## weights sum up to one
##
## Arguments:
##############
## arg - parameter vector
## cons - multiplicative offset
## y1 - vector containing the read counts for the sample of interest
## y2 - vector containing the read counts for the SssI control
## ncomp - number of beta mixture components (using a uniform prior ncomp=1)
##
## Value:
##############
##
## The function returns the sum of the two mixture weights
.eqn2 <- function(arg, cons, y1, y2, ncomp){
    
    z <- arg[3] + arg[4]
    return(z)
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
.getcredible <-  function(u, method, level, nmarg, y1, y2, al, bl, W, control.available, w, a, b, cons){

    if(is.na(W[u]) || is.infinite(W[u])){
        return(c(NA, NA))
    }

    method <- match.arg(method, c("quantile", "HPD"))

    if(ncol(w) != 1){
        w_tmp <- w[,u]
        a_tmp <- a[,u]
        b_tmp <- b[,u]
    } else {
        w_tmp <- w[,1]
        a_tmp <- a[,1]
        b_tmp <- b[,1]
    }

    if(length(cons) != 1){
        cons <- cons[u]
    }

    if(control.available){
        y2 <- y2[u]
    } 

    eps <- 1e-6

    ## get idea where the probability mass is 
    x <- seq(eps, 1 - eps, length.out = 100)
    y <- .mydmarginal(x, y1 = y1[u], y2 = y2, al = al[u], 
        bl = bl[u], W = W[u], control.available = control.available, 
        w = w_tmp, a = a_tmp, b = b_tmp, cons = cons)

    lim <- 2*(1-level)
    lim2 <- 1- lim

    mq <- .myquantile(c(0.01, lim, lim2, 0.99),x,y)$x

    x <- c( seq(eps, mq[1], length.out=100), 
        seq(mq[1]+eps, mq[2], length.out=nmarg/2), 
        seq(mq[2]+eps, mq[3], length.out=100), 
        seq(mq[3]+eps, mq[4], length.out=nmarg/2), 
        seq(mq[4]+eps, 1-eps, length.out=100))

    y <- .mydmarginal(x, y1 = y1[u], y2 = y2, al = al[u], 
        bl = bl[u], W = W[u], control.available = control.available, 
        w = w_tmp, a = a_tmp, b = b_tmp, cons = cons)  

    if(method=="quantile"){
        ll <- (1-level)/2
        q <- c(rev(ll), 1-ll)
        return(.myquantile(q=q, x,y)$x)
    }

    if(method=="HPD"){
        return(.myhpd(level, x,y))
    }
}


## get quantile  based intervals
.getcredibleDBD <-  function(u, method, level, nmarg, y1, y2, al, bl, W, control.available, w, a, b, cons){

    method <- match.arg(method, c("quantile"))
    if(length(cons) != 1){
        cons <- cons[u]
    }

    if(control.available){
        y2 <- y2[u]
    } 
    eps <- 1e-6
#     x <- seq(eps, 1-eps, length.out=nmarg)
# 
#     y <- Repitools:::.mydmarginalDBD(x, y1=y1[u], y2=y2[u], al=al[u], 
#         bl=bl[u], W=W[u], control.available=control.available, 
#         w=w[,u], a=a[u], b=b[u], cons=cons[u])


    ## get idea where the probability mass is 
    x <- seq(eps, 1 - eps, length.out = 100)
    y <- .mydmarginalDBD(x, y1=y1[u], y2=y2, al=al[u], 
        bl=bl[u], W=W[u], control.available=control.available, 
        w=w[,u], a=a[u], b=b[u], cons=cons)

    lim <- 2*(1-level)
    lim2 <- 1- lim

    mq <- .myquantile(c(0.01, lim, lim2, 0.99),x,y)$x

    x <- c( seq(eps, mq[1], length.out=100), 
        seq(mq[1]+eps, mq[2], length.out=nmarg/2), 
        seq(mq[2]+eps, mq[3], length.out=100), 
        seq(mq[3]+eps, mq[4], length.out=nmarg/2), 
        seq(mq[4]+eps, 1-eps, length.out=100))

    y <- .mydmarginalDBD(x, y1=y1[u], y2=y2, al=al[u], 
        bl=bl[u], W=W[u], control.available=control.available, 
        w=w[,u], a=a[u], b=b[u], cons=cons)

    ll <- (1-level)/2
    q <- c(rev(ll), 1-ll)
    return(.myquantile( q=q, x,y)$x)
}


## function to numerially derive the HPD interval based on a grid of density values
.myhpd = function(level, x, y, delta=0.001, maxiter=15){

  delta <- delta^2

  # the lower HPD limit must be between zero
  lower1 <- x[1]
  # and the 1-level quantile
  tail <- 1-level
  tail_low <- .myquantile(tail, x, y)
  # if the y-value at the 1-level quantile is smaller than
  # the right most denisty value, the density is not well-behaved
  # and we return the 1-level quantile and 1 as result.
  if( tail_low$y < y[length(y)]){
    return(c(tail_low$x, 1))
  }
  # find the level qunantile and the corresponding y-value.
  # if the y-value at the level quantile is smaller than
  # the left most denisty value, the density is not well-behaved
  # and we return 0 and the level quantile as result.  
  tail_high <- .myquantile(level,x,y)
  if( tail_high$y < y[1]){
    return(c(0, tail_high$x))
  }
  # if the 1-level quantile is equal to the first grid value
  # or the level quantile is equal to the last grid value we
  # return the quantile. The reason is probably that there are
  # not enough point in the grid of x and y values.
  if( tail_low$idx == 1 || tail_high$idx ==  length(x)){
    return(c(tail_low$x, tail_high$x))
  }

  # if we have a well-behaved density find the HPD interval 
  # by interval halving, use zero and the 1-level quantile as start
  lower2 <- tail_low$x
  
  i<-0
  diff = 5
  while((diff^2>delta) & (i<maxiter))
  {
    low <- (lower1 + lower2)/2
    lprob<-.mycdf(low,x, y)
    idx <- which.min(abs(low - x))
    x_c <- x[idx]
    if(low > x_c){
        ldens <- y[idx] + (y[idx + 1] - y[idx])*(x_c - x[idx])/(x[idx+1] - x[idx])
    } else {
        ldens <- y[idx - 1] + (y[idx] - y[idx - 1])*(x_c - x[idx - 1])/(x[idx] - x[idx - 1])
    }
    uprob <- lprob+level
    upp <- .myquantile(uprob,x,y)
    udens <- upp$y
    diff <- (udens-ldens)/(udens+ldens)
    i <- i+1
    if (diff<0) lower2<-low else lower1<-low
  }
  result<-c(low,upp$x)
  return(result)
}



.mycdf = function(value, x, y){
  
  idx <- which(x <= value)[1] 
  x.tmp <- x[-1] - x[-length(x)]
  y.tmp <- (y[-1] + y[-length(y)])/2
  dn <- cumsum(x.tmp * y.tmp)
  dn <- dn/dn[length(dn)]

  return(dn[idx])
}


.myquantile = function(q, x, y){

    x.tmp <- x[-1] - x[-length(x)]
    y.tmp <- (y[-1] + y[-length(y)])/2
    dn <- cumsum(x.tmp * y.tmp)
    dn <- dn/dn[length(dn)]
    idx <- c()
    for(i in 1:length(q)){
        idx[i] <- which(dn >= q[i])[1] 
    }
    return(list(x=x[idx], y=y[idx], idx=idx))
}

.mydmarginal <- function(mu, y1, y2, al, bl, W, control.available, w, a, b, cons, log=F){

    if(any(mu < 0) || any(mu > 1)){
        stop("mu is not within the support of [0,1]")
    }

    if(length(w) > 1){
        prior_term <- 0
        for(i in 1:length(w)){
            prior_term <- prior_term + w[i]*dbeta(mu, shape1=a[i], shape2=b[i])
        }
    } else {
        # since mu in [0,1] we get 1 using a uniform prior
        prior_term <- dbeta(mu, a, b)
    }

    denom <- bl + cons
    if(control.available){
        denom <- denom + 1
    } 

    # res_tmp <- prior_term*mu^{y1}/W* (1- (cons*(1-mu))/(denom))^{-(al + y1 + y2)}
    # it is more stable to calculate everything on log scale
    res_tmp <- exp(log(prior_term) + y1*log(mu) - log(W) - (al + y1 + y2)*log(1-(cons*(1-mu))/denom))

    if(log){
        return(log(res_tmp))
    } else {
        return(res_tmp)
    }
}

.mydmarginalDBD <- function(mu, y1, y2, al, bl, W, control.available, w, a, b, cons, log=F){

    if(any(mu < 0) || any(mu > 1)){
        stop("mu is not within the support of [0,1]")
    }

    denom <- bl + cons
    if(control.available){
        denom <- denom + 1
    } 
    prior_term <- w[1] * as.numeric(mu==0) + w[2]*dbeta(mu, shape1=a, shape2=b) + w[3]*as.numeric(mu==1)

    # res_tmp <- prior_term*mu^{y1}/W* (1- (cons*(1-mu))/(denom))^{-(al + y1 + y2)}
    res_tmp <- exp(log(prior_term) + y1*log(mu) - log(W) - (al + y1 + y2)*log(1-(cons*(1-mu))/denom))

    if(log){
        return(log(res_tmp))
    } else {
        return(res_tmp)
    }
}


## mixture of point mass at zero, beta distribution, and point mass at one:
## w1*delta_0 + w2*Be(x,a,b) + (1-w1-w2)*delta_1
.diracBetaDirac <- function(x, w1, w2, a, b){

    y <- w2*dbeta(x, shape1=a, shape2=b)
    if(x==0){
        y <- y + w1
    }
    if(x==1){
        y <- y + (1-w1-w2)
    }
    return(y)
}


.methylEstbeta <- function(x, verbose=TRUE, controlCI=list(compute=FALSE, method="Wald", 
    level=0.95, nmarg=512, ncpu=NULL)){

    if(controlCI$compute){
        # check whether the right CI-type is chosen
        if(priorTab(x)$ncomp > 1 & controlCI$method != "quantile"){
            stop("\n\t The CI options \"Wald\" and \"HPD\" are only available when using a uniform prior for the methylation level!\n\t Please choose either the option \"quantile\" or turn the computation of CI intervals off.")
        }
    }

    f <- fOffset(x)
    paramTab <- priorTab(x)
    cpggroup <- paramTab[["CpG groups"]]
    ncomp <- paramTab[["ncomp"]]

    myMean <- myVar <- myAl <- myBl <- myW <- matrix(NA, nrow=length(x),
        ncol=nsampleInterest(x))
    ci <- list()
    control.available <- TRUE

    for(j in 1:nsampleInterest(x)){

        y1 <- as.numeric(sampleInterest(x)[,j])
        if(nrow(control(x)) == 1){
            w.na <- !is.na(y1)
            y1 <- y1[w.na]
            y2 <- 0
            control.available <- FALSE
        } else {
            if(ncontrol(x) == 1){
                y2 <- as.numeric(control(x)[,1])
            } else {
                y2 <- as.numeric(control(x)[,j])
            }
            # only compute methylation estimates for bins with 
            # sensible values for control and sample of interest
            w.na <- !is.na(y1+y2)
            y1 <- y1[w.na]
            y2 <- y2[w.na]
        }

        cpggroup <- paramTab[["CpG groups"]]
        cpggroup <- cpggroup[w.na]
        
        len <- length(y1)

        # the first two entries of paramTab save
        # cpggroup levels and number of components
        params <- paramTab[[3+j]]
        al <- params[1, cpggroup]
        bl <- params[2, cpggroup]

        if(ncomp==1){
            ## will change for different mixtures
            w <- a <- matrix(1, nrow=1, ncol=1)
            ## for b we need the full length as it directly goes to hyperg2F1_vec
            b <-  matrix(1, nrow=1, ncol=length(cpggroup))
        } else {
            if(ncomp==2){
                w1 <- params[3,cpggroup]
                w <- rbind(w1, 1-w1)
                a <- params[c(5,6),cpggroup]
                b <- params[c(7,8),cpggroup]
            } else {
                w1 <- params[3,cpggroup]
                w2 <- params[4,cpggroup]
                w <- rbind(w1, w2, 1-w1-w2)
                a <- params[c(6, 8, 7),cpggroup]
                b <- params[c(7, 8, 6),cpggroup]    
            }
        }


        cons <- f[,j]
        if(length(cons) > 1){
            cons <- cons[w.na]
        } 

        # reset for every sample!
        A <- B <- W <- 0

        if( verbose )
            message(paste0("\nSample ",j, ":\tGetting mean and variance\t"))

        if(control.available){
            z_arg <- cons/(cons + 1 + bl)
        } else {
            z_arg <- cons/(cons+bl)
        }
        a_arg <- y1 + y2 + al
        for(i in 1:ncomp){
            ab <- a[i,] + b[i,]
            y1a <- y1 + a[i,]
            y1ab <- y1a + b[i,]
            A <- A + w[i,]*exp(lgamma(ab) + lgamma(y1a + 1) - 
                (lgamma(a[i,]) + lgamma(y1ab + 1))) *
                hyperg2F1_vec(a=a_arg, b=b[i,], c=y1ab + 1, z=z_arg)
            B <- B + w[i,]*exp(lgamma(ab) + lgamma(y1a + 2) - 
                (lgamma(a[i,])+lgamma(y1ab + 2))) *
                hyperg2F1_vec(a=a_arg, b=b[i,], c=y1ab + 2, z=z_arg)
            W <- W + w[i,]*exp(lgamma(ab) + lgamma(y1a) - 
                (lgamma(a[i,])+lgamma(y1ab))) *
                hyperg2F1_vec(a=a_arg, b=b[i,], c=y1ab, z=z_arg)
        }
        tmp_mean <- A/W
        tmp_var <- B/W - (A/W)^2
        if(ncomp==1){
            rm(a_arg, ab, y1a, y1ab, z_arg, A, B, b)
            b <- matrix(1, ncol=1, nrow=1)
        } else {
            rm(a_arg, ab, y1a, y1ab, z_arg, A, B)
        }
        gc()

        if(controlCI$compute){

            if(controlCI$method=="Wald"){
                if( verbose )
                    message(paste0("Sample ",j, ":\tGetting Wald-based credible interval.\n"))
                    
                ll <- (1-controlCI$level[1])/2

                logit_mean <- .logit(tmp_mean) 
                logit_sd <- sqrt(.dlogit(tmp_mean)^2*tmp_var)

                lci <- uci <- rep(NA, length(x))
                lci[w.na] <- .invlogit(logit_mean - abs(qnorm(ll)) * logit_sd)
                uci[w.na] <- .invlogit(logit_mean + abs(qnorm(ll)) * logit_sd)

                #ci[[j]] <- cbind(lci=lci, uci=uci, width=uci-lci)
                ci[[j]] <- cbind(lci=lci, uci=uci)
            } else {
                if( verbose )
                    message(paste0("Sample ",j, ":\tGetting ", controlCI$method, " credible interval.\n"))

                if(is.null(controlCI$ncpu)){
                    controlCI$ncpu <- .getCpu(maxCPU=FALSE)
                }
                method <- controlCI$method
                level <- controlCI$level[1]
                nmarg <- controlCI$nmarg

                ci_tmp <- mclapply(1:len, .getcredible, 
                    method = method, level = level, nmarg = nmarg, 
                    y1 = y1, y2 = y2, al = al, bl = bl, W = W, 
                    control.available = control.available, w = w, 
                    a = a, b = b, cons = cons,  mc.cores=controlCI$ncpu)
                ci[[j]] <- do.call(rbind, ci_tmp)  
                colnames(ci[[j]]) <- c("lower.limit", "upper.limit")
                gc()
            }
        }
        myMean[w.na,j] <- tmp_mean
        myVar[w.na,j] <- tmp_var
        myW[w.na,j] <- W
    }
    colnames(myMean) <- colnames(myVar) <- colnames(myW) <- colnames(sampleInterest(x))
    if(controlCI$compute)
        names(ci) <- colnames(sampleInterest(x))

    methEst(x) <- list(mean=myMean, var=myVar, ci=ci, W=myW)

    return(x)
}


.methylEstDBD <- function(x, verbose=TRUE, controlCI=list(compute=FALSE, method="quantile", 
    level=0.95, nmarg=512, ncpu=NULL)){

    if(controlCI$compute){
        # check whether the right CI-type is chosen
        if(controlCI$method != "quantile"){
            stop("\n\t The CI options \"Wald\" and \"HPD\" are only available when using a uniform prior for the methylation level!\n\t Please choose either the option \"quantile\" or turn the computation of CI intervals off.")
        }
    }


    f <- fOffset(x)
    paramTab <- priorTab(x)
    cpggroup <- paramTab[["CpG groups"]]
    ncomp <- paramTab[["ncomp"]]

    myMean <- myVar <- myAl <- myBl <- myW <- matrix(NA, nrow=length(x),
        ncol=nsampleInterest(x))
    ci <- list()
    control.available <- TRUE

    for(j in 1:nsampleInterest(x)){

        y1 <- sampleInterest(x)[,j]
        if(nrow(control(x)) == 1){
            w.na <- !is.na(y1)
            y1 <- y1[w.na]
            y2 <- rep(0, length(y1))
            control.available <- FALSE
        } else {
            if(ncontrol(x) == 1){
                y2 <- control(x)[,1]
            } else {
                y2 <- control(x)[,j]
            }
            # only compute methylation estimates for bins with 
            # sensible values for control and sample of interest
            w.na <- !is.na(y1+y2)
            y1 <- y1[w.na]
            y2 <- y2[w.na]
        }

        cpggroup <- paramTab[["CpG groups"]]
        cpggroup <- cpggroup[w.na]
        
        len <- length(y1)

        # the first two entries of paramTab save
        # cpggroup levels and number of components
        params <- paramTab[[3+j]]
        al <- params[1, cpggroup]
        bl <- params[2, cpggroup]

        # beta parameters
        a <- params[3, cpggroup]
        b <- params[4, cpggroup]
    
        # weights for second and third component
        w2 <- params[5, cpggroup]
        w3 <- params[6, cpggroup]
        w <- rbind(1-w2-w3, w2, w3)

        cons <- f[,j]
        if(length(cons) > 1){
            cons <- cons[w.na]
        } else {
            cons <- rep(cons, length(y1))
        }

        # reset for every sample!
        A <- 0
        B <- 0
        W <- 0

        if( verbose )
            message(paste0("Sample ",j, ":\tGetting mean and variance\t"))

        if(control.available){
            z_arg <- cons/(cons + 1 + bl)
        } else {
            z_arg <- cons/(cons+bl)
        }
        a_arg <- y1 + y2 + al

        ab <- a + b
        y1a <- y1 + a
        y1ab <- y1a + b
        W <- w3 + w2*exp(lgamma(ab) + lgamma(y1a) - 
            (lgamma(a)+lgamma(y1ab))) * hyperg2F1_vec(a=a_arg, b=b, c=y1ab, z=z_arg)
        A <- w3 + w2*exp(lgamma(ab) + lgamma(y1a + 1) - 
            (lgamma(a) + lgamma(y1ab + 1))) * hyperg2F1_vec(a=a_arg, b=b, c=y1ab + 1, z=z_arg)
        B <- w3 + w2*exp(lgamma(ab) + lgamma(y1a + 2) - 
            (lgamma(a)+lgamma(y1ab + 2))) * hyperg2F1_vec(a=a_arg, b=b, c=y1ab + 2, z=z_arg)

        tmp_mean <- A/W
        tmp_var <- B/W - (A/W)^2

        if(ncomp==1){
            rm(a_arg, ab, y1a, y1ab, z_arg, A, B, b)
            b <- matrix(1, ncol=1, nrow=1)
        } else {
            rm(a_arg, ab, y1a, y1ab, z_arg, A, B)
        }
        gc()

        if(controlCI$compute){
               if( verbose )
                    message(paste0("Sample ",j, ":\tGetting ", controlCI$method, " credible interval.\n"))

                if(is.null(controlCI$ncpu)){
                    controlCI$ncpu <- .getCpu(maxCPU=FALSE)
                }

                method <- controlCI$method
                level <- controlCI$level[1]
                nmarg <- controlCI$nmarg

                ci_tmp <- mclapply(1:len, .getcredibleDBD, 
                    method = method, level = level, nmarg = nmarg, 
                    y1 = y1, y2 = y2, al = al, bl = bl, W = W, 
                    control.available = control.available, w = w, 
                    a = a, b = b, cons = cons,  mc.cores=controlCI$ncpu)
                ci_tmp_mclapply <- do.call(rbind, ci_tmp)  
                gc()

                ci[[j]] <- cbind(lci=ci_tmp_mclapply[1,], uci=ci_tmp_mclapply[2,])
        }

        myMean[w.na,j] <- tmp_mean
        myVar[w.na,j] <- tmp_var
        myW[w.na,j] <- W
    }
    colnames(myMean) <- colnames(myVar) <- colnames(myW) <- colnames(sampleInterest(x))
    if(controlCI$compute)
        names(ci) <- colnames(sampleInterest(x))

    methEst(x) <- list(mean=myMean, var=myVar, ci=ci, W=myW)

    return(x)
}


.belongsTo <- function(value, rangeStart, rangeEnd){
        
    res <- (value >= rangeStart) & (value <= rangeEnd)
    return(res)
}

