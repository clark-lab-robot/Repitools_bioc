methylEst <- function(x, verbose=FALSE, controlCI=list(compute=FALSE, method="Wald", 
    level=0.95, nmarg=512, ncpu=NULL)){

    if(class(x) != "BayMethList"){
        stop("x must be a BayMethList object")
    }
    if(is.null(controlCI$method)){
        controlCI$method <- "Wald"
    }
    if(!is.element(controlCI$method, c("Wald", "HPD", "quantile"))){
        stop("The method for computing credible intervals must be either \"Wald\", \"HPD\" or \"quantile\"")
    }
    if(is.null(controlCI$compute)){
        controlCI$compute <- FALSE
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

    if(priorTab(x)$method == "beta"){
        .methylEstbeta(x=x, verbose=verbose, controlCI=controlCI)
    } else if(priorTab(x)$method == "DBD"){
        .methylEstDBD(x=x, verbose=verbose, controlCI=controlCI)
    } else {
        stop("This should not happen!")
    }
}

