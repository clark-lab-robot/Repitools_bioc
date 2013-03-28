methylEst <- function(x, verbose=FALSE, controlCI=list(compute=FALSE, method="Wald", 
    level=0.95, nmarg=512, ncpu=NULL)){

    if(class(x) != "BayMethList"){
        stop("x must be a BayMethList object")
    }
    if(!is.element(controlCI$method, c("Wald", "HPD", "quantile"))){
        stop("The method for computing credible intervals must be either \"Wald\", \"HPD\" or \"quantile\"")
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

                ci[[j]] <- cbind(lci=lci, uci=uci, width=uci-lci)
            } else {
                if( verbose )
                    message(paste0("Sample ",j, ":\tGetting ", controlCI$method, " credible interval.\n"))

                if(is.null(controlCI$ncpu)){
                    controlCI$ncpu <- Repitools:::.getCpu(maxCPU=FALSE)
                }

                sfInit(parallel=TRUE, cpus=controlCI$ncpu)
                sfLibrary("Repitools", character.only=TRUE )
                sfExport(".getcredible", namespace="Repitools")
                sfExport(".myquantile", namespace="Repitools")
                sfExport(".mydmarginal", namespace="Repitools")
                sfExport(".myhpd", namespace="Repitools")
                sfExport(".mysfmarginal", namespace="Repitools")
                sfExport(".mysmarginal", namespace="Repitools")
                sfExport(".myhpd", namespace="Repitools")
                sfExport(".mymarginalfix", namespace="Repitools")
                sfExport("controlCI")
                sfExport("tmp_var")
                sfExport("al")
                sfExport("bl")
                sfExport("W")
                sfExport("cons")
                sfExport("y1")
                sfExport("y2")
                        
                ci_tmp <- sfSapply(1:length(y1), .getcredible, 
                    method=controlCI$method,
                    level=controlCI$level[1], 
                    nmarg=controlCI$nmarg, 
                    y1=y1, y2=y2, al=al, bl=bl, W=W, cons=cons, my_sd= sqrt(tmp_var))

                sfStop()

                ci[[j]] <- cbind(lci=ci_tmp[1,], uci=ci_tmp[2,], width=ci_tmp[2,]-ci_tmp[1,])

            }
        }
        myMean[w.na,j] <- tmp_mean
        myVar[w.na,j] <- tmp_var
        myAl[w.na,j] <- al
        myBl[w.na,j] <- bl
        myW[w.na,j] <- W
    }
    colnames(myMean) <- colnames(myVar) <- colnames(myAl) <- colnames(myBl) <- colnames(myW) <- colnames(sampleInterest(x))
    names(ci) <- colnames(sampleInterest(x))

    methEst(x) <- list(mean=myMean, var=myVar, ci=ci, 
        W=myW, al=myAl, bl=myBl)

    return(x)
}
