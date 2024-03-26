##  VT::08.11.2023 This is an improved version of doScale() from robustbase
##
myscale <- function (x, center, scale)
{
    stopifnot(is.numeric(p <- ncol(x)))
    ## MM: follow standard R's	scale.default() as much as possible

	if(missing(center) || is.null(center) || (is.logical(center) && !center)) 
            center <- rep(0, p)
    else {
        if(is.logical(center) && center)
            center <- mean
        
        centerFn <- is.function(center)
        doIt <- if(centerFn) {
            centerName <- deparse(substitute(center)) # "median" typically
    	   center <- apply(x, 2L, center)
    	   TRUE
        } else {
        	if(length(center) == p && is.numeric(center))
        	    TRUE
        	else if(missing(center) || is.null(center)) {
                    center <- 0; FALSE
    	} else
    	    stop(gettextf("'%s' must be a function, numeric vector of length p, or NULL",
                              "center"), domain=NA)
        }
        if(doIt)
    	   x <- sweep(x, 2L, center, `-`, check.margin=FALSE)
    }
    
	if(missing(scale) || is.null(scale) || (is.logical(scale) && !scale)) 
            scale <- rep(1, p)
    else {

        if(is.logical(scale) && scale)
            scale <- sd

        scaleFn <- is.function(scale)
        doIt <- if(scaleFn) {
    	   scale <- apply(x, 2L, scale)
    	   TRUE
        } else {
        	if(length(scale) == p && is.numeric(scale))
        	    TRUE
        	else if(missing(scale) || is.null(scale)) {
        	    scale <- 1
        	    FALSE
        	} else
        	    stop(gettextf("'%s' must be a function, numeric vector of length p, or NULL",
                                  "scale"), domain=NA)
            }
            
            if(doIt) {
                if(any(is.na(scale)) || any(scale < 0))
                    stop("provide better scale; must be all positive")
                if(any(s0 <- scale == 0)) {
        ## FIXME:
        ### Better and easier alternative (and as "FAST MCD"): return "singular cov.matrix"
        ### since scale 0 ==>  more than 50% points are on hyperplane x[,j] == const.
                    ## find scale if there is any variation; otherwise use s := 1
                    S <- if(centerFn && centerName == "median")
                             abs else function(.) abs(. - median(.))
                    non0Q <- function(u) {
                        alph <- c(10:19, 19.75)/20 # not all the way to '1' {=> finite qnorm()}
                        qq <- quantile(S(u), probs=alph, names=FALSE)
                        if(any(pos <- qq != 0)) { ## the first non-0 if there is one
                            i <- which.max(pos)
                            qq[i] / qnorm((alph[i] + 1)/2)
                        } else 1
                    }
                    scale[s0] <- apply(x[,s0, drop=FALSE], 2L, non0Q)
                }
                x <- sweep(x, 2L, scale, `/`, check.margin = FALSE)
            }
        }
    
    ## return
    list(x=x, center=center, scale=scale)
}

