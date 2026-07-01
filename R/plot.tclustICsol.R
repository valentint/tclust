##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/r/tclust", load_code=roxygen2:::load_installed)

#' Plot Method for \code{tclustICsol} Objects
#'
#' The plot method for class \code{tclustICsol}.
#'
#' @name plot.tclustICsol
#' @method plot tclustICsol

#' @description Displays one of the solutions, selected by the argument 'sol'.
#'  The default display is a scatterplot matrix of the data using
#'  colors and symbols of the observations to identify the groups in the 
#'  selected solution. If the argument 'choice' is specified and its length 
#'  is two, a simple scatterplot will be shown.
#'  The function fails, if \code{store.x = FALSE} is specified in 
#'  the \code{tclustICsol()} call because the original data matrix is 
#'  required here.
#'
#' @param x The \code{tclustICsol} object to be displayed
#' @param whichIC A character value which specifies which information criteria 
#'  to use. Fpr the possible values for \code{whichIC} see the help of \code{tclustIC}.
#' @param sol Which solution to display - a number between 1 and the \code{nsol}
#'  argument of \code{tclustICsol}.
#' @param col optional colors to identify the groups. If not specified default 
#'  values will be selected.
#' @param pch optional symbols to identify the groups. If not specified default 
#'  values will be selected.
#' @param main optional title. The default title shows the information criteria, 
#'  the solution number, number of groups and restriction factor c used and 
#'  indication wheather the solution is true or spurios.
#' @param sub1 an optional subtitle. The default subtitle shows the type of the 
#'  displaued solution: 'Best in' and 'Stable in'.
#' @param choice a numeric vector of length between 2 and the number of 
#'  variables in the input data matrix. If missing, a scatterplot matrix with 
#'  all variables will be shown. If two variables are selected, a simple 
#'  scatterplot of the two selected variables is shown.
#' @param \ldots Further (optional) graphical arguments
#' 
#' @examples
#'
#'  #--- EXAMPLE 1 ------------------------------------------
#'  \donttest{
#'  data(geyser2)
#'  (out <- tclustIC(geyser2, whichIC="MIXMIX", alpha=0.1))
#'
#'  ##  Show the first two best solutions using as Information criterion MIXMIX
#'  cat("\nBest solutions using MIXMIX\n")
#'  outsol <- tclust::tclustICsol(out, whichIC="MIXMIX", nsol=2)
#'  plot(outsol)
#'  plot(outsol, choice=c(1,2))
#'  plot(outsol, choice=c(1,2), xlab="XLAB", ylab="YLAB")
#'  }
#'
plot.tclustICsol <- function(x, whichIC, sol=1, col, pch, main, sub1, choice, ...) {

    ## The requested information criterion whichIC should correspond 
    ##  to those of the 'tclustICsol' object 
    if(missing(whichIC)) {
        whichIC <- x$whichIC
        if(whichIC == "ALL")
            whichIC <- "MIXMIX"  
    } else {
        if(x$whichIC == "ALL")
            choices <- c("MIXMIX", "MIXCLA", "CLACLA")
        else if(x$whichIC == "MIXMIX")
            choices <- c("MIXMIX")
        else if(x$whichIC == "MIXCLA")
            choices <- c("MIXCLA")
        else if(x$whichIC == "CLACLA")
            choices <- c("CLACLA")
            
        whichIC <- match.arg(whichIC, choices)
    }

    df <- if(whichIC == "MIXMIX") x$MIXMIXbs else if(whichIC=="CLACLA") x$CLACLAbs else x$MIXCLAbs
    idx <- if(whichIC == "MIXMIX") x$MIXMIXbsIDX else if(whichIC=="CLACLA") x$CLACLAbsIDX else x$MIXCLAbsIDX
    
    ## The selected solution 'sol' should be between 1 and the number of solutions in 
    ##  the input 'tclustICsol' object (i.e. the number of rows of 
    ##  the chosen MIXMIXbs or CLACLAbs or MIXCLAbs matrix. Defaults to 1.
    ## The values of 'k', 'c' and IDX will be chosen according to 'sol'
    if(sol < 1 || sol > nrow(df))
        stop(paste("'sol' must be between 1 and ", nrow$df))
    k <- df[[sol,1]]
    c <- df[[sol,2]]
    id <- idx[,sol]

    ## Select colors and symbols for plotting according to the number of groups
    maxassig <- max(id)     # number of groups
    if(missing(col))
        col <- 1:(k+1)
    else
        col <- rep(col, len=maxassig + 1)  
    
    if(missing(pch))
        pch <- 1:(k+1)
    else
        pch <- rep(pch,  len=maxassig + 1)
    col <- col[id + 1]    
    pch <- pch[id + 1]

    ## Select the type of solution: 'Best in' and 'Stable in'. These will be 
    ##  used to form the subtitle.
    best <- df[[sol, 3]]
    stable <- sort(c(best, df[[sol, 4]]))
    type <- df[sol, 5]
    xbest <- if(length(best) > 1) paste("Best in c=", best[1], "-", best[length(best)], sep="") 
        else if(length(best) == 1) paste("Best in c=", best[1], sep="") else ""
    xstable <- if(length(stable) > 1) paste("Stable in c=", stable[1], "-", stable[length(stable)], sep="") 
        else if(length(stable) == 1) paste("Stable in c=", stable[1], sep="") else ""

    if(missing(main))
        main <- paste(whichIC, " (solution ", sol, ", ", type, "): k=", k, " c=", c, sep="")
    if(missing(sub1))
        sub1 <- paste(xbest, "", xstable)

    ## Chose the type of plot according to the parameter 'choice'. 
    ##  - If missing 'choice', independent of the number of variables 
    ##      a scatterplot matrix will be plotted. Same, if the length of 
    ##      'choice' is greater than 2.
    ##  - If two components are selected a simple scatterplot of the two 
    ##      selected variables will be shown
    ##  Select the corresponding columns from the data matrix.
    data <- x$x
    if(!missing(choice)) {
        if(length(choice) < 2 || length(choice) > ncol(data))
            stop(paste("The length of 'choice' must be greater than 2 and smaller than ", ncol(data)))
         if(any(choice < 0) || any(choice > ncol(data)))
            stop(paste("'chice' must  be between 1 and ", ncol(data)))
        if(length(choice) != length(unique(choice)))
            stop(paste("There are duplicated columns in 'choice'!"))
         
        data <- data[, choice]
    }

    if(missing(choice) || ncol(data) > 2) {
        ## Matrix scatterplot when p > 2
        .plot.tclustICsol.pairs(data, ngroups=k, groups=id, col=col, pch=pch, main=main, sub1=sub1, ...)
    } else
        .plot.tclustICsol.2d(data, ngroups=k, groups=id, col=col, pch=pch, main=main, sub1=sub1, ...)
    
    return(invisible(x))
}

##  Draw pairwise scatter plots for the data set 'x'
##  - upper triangle - scatter plot with classical and robust 0.975-ellipses
##  - histograms on the diagonal
##  - lower triangle - robust (MCD) and classical correlations
##
##  - x     - data
##  - main  - caption of the plot
##
.plot.tclustICsol.pairs <- function(x, ngroups, groups, main="", sub1="", ...){
    
    panel.density <- function(x,...)
    {
        usr <- par("usr"); on.exit(par(usr=usr))
        par(usr = c(usr[1:2], 0, 1.5) )

        args <- list(...)
        if("col" %in% names(args)) {
            col <- sort(unique(args[["col"]]))
            if(is.null(col))
                col <- (1:ngroups) + 1
        } else
            col <- (1:ngroups) + 1

        for(i in 1:ngroups) {
            xx <- x[groups==i]
            tryd <- try(d <- density(xx, na.rm=TRUE, bw="nrd", adjust=1.2), silent=TRUE)

            ## VT::11.08.2022: fix error "Found if() conditions comparing class() to string"
            ##  if(class(tryd) != "try-error")
            if(!is(tryd, "try-error")) {
                d$y <- d$y/max(d$y)
                lines(d, col=col[i])
            }
        }
    }

    panel.boxplot <- function(x,...)
    {
        xr <- range(x, na.rm = TRUE)
        usr <- par("usr"); on.exit(par(usr=usr))

        ff <- 0.1
        dd <- ff * (xr[2] - xr[1])
        xr[1] <- xr[1] - dd
        xr[2] <- xr[2] + dd
        par(usr=c(xr, 0, 1.5))

       
## cat("\nxr=", xr, "\n")
## cat("\nminmax=", c(min(x), max(x)), "\n")
## cat("\nusr=", par("usr"), "\n")

        args <- list(...)
        if("col" %in% names(args)) {
            col <- sort(unique(args[["col"]]))
            if(is.null(col))
                col <- (1:ngroups) + 1
        } else
            col <- (1:ngroups) + 1

##  cat("\nargs=", "col" %in% names(args), ", ngroups=", ngroups, ", col=", col, "\n")

##        if(length(groups == 0) > 0) {   # there are trimmed observations, remove them from the boxplots
##            x <- x[groups > 0]
##            groups <- groups[groups > 0]
##            if(length(col) > length(unique(groups)))
##                col <- col[-1]
##        }

##  cat("\nngroups=", ngroups, ", length(table(groups))=", length(table(groups)), "Classes: ", table(groups), "\n")
        ngroups <- length(table(groups))

        boxplot(x ~ groups, 
              add=TRUE, 
              horizontal=TRUE,
              at=seq(0.2, 0.8, length.out=ngroups), 
              axes=FALSE, 
              border=col, col=0, 
              pars=list(boxwex = 0.15)) # Adjust box thickness
    }

    panel.ellipse <- function(x, y, ...) {
        usr <- par("usr"); on.exit(par(usr=usr))

        xmin <- min(x)
        xmax <- max(x)
        ymin <- min(y)
        ymax <- max(y)

        ff <- 0.1
        xmin <- xmin - ff*(xmax-xmin); #print(ff*(xmax-xmin))
        xmax <- xmax + ff*(xmax-xmin); #print(ff*(xmax-xmin))
        ymin <- ymin - ff*(ymax-ymin); #print(ff*(ymax-ymin))
        ymax <- ymax + ff*(ymax-ymin); #print(ff*(ymax-ymin))
        par(usr = c(xmin, xmax, ymin, ymax))

        points(x,y, ...)
    }

##    cat("\nngroups=", ngroups, ", length(table(groups))=", length(table(groups)), "Classes: ", table(groups), " Not trimmed: ", length(table(groups[groups > 0])), "\n")
##    print(table(groups))

    if(ngroups > length(table(groups[groups > 0]))) {     ## there are empty clusters
        warning("There are empty clusters: ", ngroups," requested ", length(table(groups[groups > 0])), " available!")
    }
    
    pairs(x, 
        lower.panel=panel.ellipse,
        diag.panel=panel.boxplot,
        upper.panel=panel.ellipse,
        labels=names(x),
        ...)

    if(!is.null(main))
        title(main=main, line=3)
  
    if(!is.null(sub1))
        title(sub=sub1)      
}

.plot.tclustICsol.2d <- function(x, ngroups, groups, xlab, ylab, axes=3, xlim, ylim, col=col, pch=pch, main, sub1, ...){

    if(missing(xlab))
        xlab <- colnames(x)[1]
    if(missing(ylab))
        ylab <- colnames(x)[2]
        
    if(missing(xlim))
        xlim <- range(x[, 1])
    
    if(missing(ylim))
        ylim <- range(x[, 2])

    fact <- 0.04
    r <- cbind(xlim, ylim)
    rd <- apply(r, 2, diff)
    usr <- as.numeric (r + (c(-1, 1) %*% t (rd)) * fact)
    
    plot.new()
    par(usr = usr)

    points(x[,1], x[,2], pch=pch, col=col)

    if(!is.null(main))
        title(main=main)
  
    if(!is.null(sub1))
        title(sub=sub1)      

    axis.x <- axes %% 2        ## x axis 1 or 3
    axis.y <- axes >= 2        ## y axis 2 or 3
    
    cex <- par ("cex")
    if (!missing (xlab))
        mtext (side = 1, xlab, line = 1.5 + 1.5 * axis.x, cex = cex)
    if (!missing (ylab))
        mtext (side = 2, ylab, line = 1.5 + 1.5 * axis.y, cex = cex)

    box()
    if(axis.x)
        axis(1)
    if(axis.y)
        axis(2)
}
