methylEst <- function(x, controlCI=list(compute=FALSE, method="Wald", 
    level=0.95, nmarg=512, ncpu=NULL)){

    if(class(x) != "BayMethList"){
        stop("x must be a BayMethList object")
    }
    if(!is.element(controlCI$method, c("Wald"))){
        stop("The method for computing credible intervals must be \"Wald\"")
    }
    if(is.null(controlCI$compute)){
        controlCI$compute <- FALSE
    }
    if(is.null(controlCI$method)){
        controlCI$method <- "Wald"
    }
    if(is.null(controlCI$level)){
        controlCI$level <- 0.95
    }
    if(is.null(controlCI$nmarg)){
        controlCI$nmarg <- 512
    }
    if(is.null(controlCI$ncpu)){
        controlCI$ncpu <- NULL
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
            y2 <- 0
            w.na <- !is.na(y1)
            y1 <- y1[w.na]
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

        ## will change for different mixtures
        w <- matrix(1, nrow=1, ncol=length(cpggroup))
        a <- matrix(1, nrow=1, ncol=length(cpggroup))
        b <- matrix(1, nrow=1, ncol=length(cpggroup))
        
        len <- length(y1)

        # the first two entries of paramTab save
        # cpggroup levels and number of components
        params <- paramTab[[2+j]]
        al <- params[1, cpggroup]
        bl <- params[2, cpggroup]

        cons <- f[,j]
        if(length(cons) > 1){
            cons <- cons[w.na]
        }

        # reset for every sample!
        A <- 0
        B <- 0
        W <- 0

        message("\n\t Get mean and variance\n\n")
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

        ## derive Wald-based CIs on logit scale using the delta rule
        message("\n\t Get credible interval\n\n")

        ll <- (1-controlCI$level[1])/2

        logit_mean <- .logit(tmp_mean) 
        logit_sd <- sqrt(.dlogit(tmp_mean)^2*tmp_var)

        lci <- uci <- rep(NA, length(x))
        lci[w.na] <- .invlogit(logit_mean - abs(qnorm(ll)) * logit_sd)
        uci[w.na] <- .invlogit(logit_mean + abs(qnorm(ll)) * logit_sd)

        ci[[j]] <- cbind(lci=lci, uci=uci, width=uci-lci)
     
        myMean[w.na,j] <- tmp_mean
        myVar[w.na,j] <- tmp_var
        myAl[w.na,j] <- al
        myBl[w.na,j] <- bl
        myW[w.na,j] <- W
    #  }
    }
    methEst(x) <- list(mean=myMean, var=myVar, ci=ci, 
        W=myW, al=myAl, bl=myBl)

    return(x)
}
