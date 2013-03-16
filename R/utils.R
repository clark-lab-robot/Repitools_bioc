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
        control[[i]] <- control[[i]][-w.na]
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


    ## parameters for lambda
    al <- arg[1]
    bl <- arg[2]

    ## uniform prior (alternatives in progress)
    w <- 1
    a <- 1
    b <- 1

    len <- length(y1)
    denom <- bl + 1 + cons
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


