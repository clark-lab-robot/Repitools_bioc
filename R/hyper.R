hyperg2F1_vec <- function(a, b, c, z)
{
    if (is.complex(c(a, b, c, z))) {
        stop("complex values of a,b,c,z not supported.")
    }
    if(any(z >= 1)){
        stop("The absolute value of z must be smaller than one!")
    }
    if(any(b <= 0) || any(c <= b)){
        stop("The condition c > b > 0 is violated.")
    }
    la <- length(a)
    lb <- length(b)
    lc <- length(c)
    lz <- length(z)
    if(length(unique(c(la,lb,lc, lz))) > 1){
        stop("All arguments must have the same length.")
    }
    
    n <- length(a)
    all <- vector(mode="numeric", length=n)
    ret <- .C("_myHyp2f1", as.double(all), as.double(a), as.double(b), 
        as.double(c), as.double(z), as.integer(n))[[1]]
    
    return(ret)
}
