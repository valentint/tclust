#' Plot Method for \code{tclust} and \code{tkmeans} Objects
#'
#' The plot method for classes \code{tclust} and \code{tkmeans}.
#'
#' @name plot.tclust
#' @aliases plot.tkmeans
#' @method plot tclust

#' @description One and two dimensional structures are treated separately (e.g. tolerance 
#'  intervals/ellipses are displayed). Higher dimensional structures are displayed 
#'  by plotting the two first Fisher's canonical coordinates (evaluated by 
#'  \code{tclust::discr_coords}) and derived from the final cluster assignments 
#'  (trimmed observations are not taken into account). 
#'  \code{plot.tclust.Nd} can be called with one or two-dimensional \code{tclust}- or \code{tkmeans}-objects 
#'  too. The function fails, if \code{store.x = FALSE} is specified in the \code{tclust()} or \code{tkmeans()} call, 
#'  because the original data matrix is required here.
#'
#' @param x The \code{tclust} or \code{tkmeans} object to be displayed
#' @param \ldots Further (optional) arguments which specify the details of the 
#'    resulting plot (see section "Further Arguments").
#'
#' @section Further Arguments:
#'
#' \itemize{
#'  \item \code{xlab, ylab, xlim, ylim, pch, col} Arguments passed to \code{plot()}.
#'  \item \code{main} The title of the plot. Use \code{"/p"} for displaying the chosen parameters 
#'      \code{alpha} and \code{k} or \code{"/r"} for plotting the chosen restriction.
#'  \item \code{main.pre} An optional string which is added to the plot's caption.
#'  \item \code{sub} A string specifying the subtitle of the plot. Use \code{"/p"} (default) for 
#'  displaying the chosen parameters \code{alpha} and \code{k}, \code{"/r"} for plotting 
#'  the chosen restriction and \code{"/pr"} for both. 
#'  \item \code{sub1} A secondary (optional) subtitle.
#'  \item \code{labels} A string specifying the type of labels to be drawn. Either 
#'      \code{labels="none"} (default), \code{labels="cluster"} or \code{labels="observation"} 
#'      can be specified. If specified, parameter \code{pch} is ignored.
#'  \item \code{text} A vector of length n (the number of observations) containing 
#'      strings which are used as labels for each observation. If specified, 
#'      the parameters \code{labels} and \code{pch} are ignored.
#'  \item \code{by.cluster} Logical value indicating whether parameters 
#'      \code{pch} and \code{col} refer to observations (FALSE) or clusters (TRUE).
#'  \item \code{jitter.y} Logical value, specifying whether the drawn values shall be 
#'      jittered in y-direction for better visibility of structures in 1 dimensional data.
#'  \item \code{tol} The tolerance interval. 95\% tolerance ellipsoids (assuming normality) 
#'      are plotted by default.
#'  \item \code{tol.col, tol.lty, tol.lwd} Vectors of length k or 1 containing 
#'      the \code{col}, \code{lty} and \code{lwd} arguments for the tolerance 
#'      ellipses/lines.
#' }
#'
#' @examples
#'  #--- EXAMPLE 1------------------------------
#'  sig <- diag (2)
#'  cen <- rep (1, 2)
#'  x <- rbind(MASS::mvrnorm(360, cen * 0,   sig),
#'  	       MASS::mvrnorm(540, cen * 5,   sig * 6 - 2),
#'  	       MASS::mvrnorm(100, cen * 2.5, sig * 50))
#'  # Two groups and 10\% trimming level
#'  a <- tclust(x, k = 2, alpha = 0.1, restr.fact = 12)

#'  plot (a)
#'  plot (a, labels = "observation")
#'  plot (a, labels = "cluster")
#'  plot (a, by.cluster = TRUE)

#'  #--- EXAMPLE 2------------------------------
#'  sig <- diag (2)
#'  cen <- rep (1, 2)
#'  x <- rbind(MASS::mvrnorm(360, cen * 0,   sig),
#'  	       MASS::mvrnorm(540, cen * 5,   sig),
#'  	       MASS::mvrnorm(100, cen * 2.5, sig))
#'  # Two groups and 10\% trimming level
#'  a <- tkmeans(x, k = 2, alpha = 0.1)

#'  plot (a)
#'  plot (a, labels = "observation")
#'  plot (a, labels = "cluster")
#'  plot (a, by.cluster = TRUE)
#'

plot.tclust <-
function (x,  ...) {
  if (x$int$dim[2] == 1)
    .plot.tclust.1d (x, ...)
  else if (x$int$dim[2] == 2)
    .plot.tclust.2d (x, ...)
  else
    .plot.tclust.Nd (x, ...)
}

#' @rdname plot.tclust
#' @method plot tkmeans
plot.tkmeans <- function(x,  ...) {
	plot.tclust(x, ...)
}

#######################
##  .plot.tclust.1d  ##
#######################

.plot.tclust.1d <-
function (x, xlab, ylab, xlim, ylim, tol = 0.95, tol.lwd = 1, tol.lty = 3, tol.col,
          jitter.y = FALSE, ...)
{
  if (x$int$dim[2] != 1)
    stop ("tclust object of dimension 1 expected.")
  
  if (is.null (x$par$x))
    stop ("dataset not included in tclust object - cannot plot object.")

  if (missing (xlab))
  {
    dn <- dimnames (x$par$x)
    if (is.null (dn[[2]]))
      xlab <- "x"
    else
      xlab <- dn[[2]][1]
  }

  if (missing (ylab) && jitter.y)
    ylab <- "(random jitter)"

  n <- x$int$dim[1]

  if (jitter.y)
    y <- runif (n, min = -1) / 4
  else
    y <- rep (0, n)

  if (is.numeric (tol) && length (tol) == 1 &&  0 < tol && tol < 1)
    tol.fact = sqrt (qchisq(tol, 1))
  else
    tol.fact <- NULL

  x.c = as.numeric (x$centers)
  if (!is.null (x$cov))		## tkmeans- objects don't have cov-info
  {
	  x.sd = sqrt (as.numeric (x$cov))

	  if (missing (xlim))
	  {
		xlim <- range (x$par$x)
		if (!is.null (tol.fact))
		  xlim <- range (xlim, x.c + x.sd * tol.fact, x.c - x.sd * tol.fact)
	  }
  }
  else
  {
	if (missing (xlim))
		xlim <- range (x$par$x)
  }

  if (missing (ylim))
	ylim <- c (-1, 1)

  X <- cbind (x$par$x, y)
  .plot.tclust.0 (x = x, X = X, xlab = xlab, ylab = ylab, axes = 1,
                 xlim = xlim, ylim = ylim, ...)
 
  .vline (x.c, 3, lty = 2, col = 1 + (1:x$k))    ##  cluster centers

	if (!is.null (tol.fact) &&
		!is.null (x$cov))		## tkmeans- objects don't have cov-info
  {
    #tol.fact = sqrt(qchisq(tol, 1))
    if (missing (tol.col))
      tol.col <- (1:x$k) + 1
    else
      tol.col <- rep (tol.col, x$k)
      
    tol.lty <- rep (tol.lty , x$k)
    tol.lwd <- rep (tol.lwd , x$k)   
    .vline (x.c + x.sd * tol.fact, 2, col = tol.col, lwd = tol.lwd,
            lty = tol.lty)
    .vline (x.c - x.sd * tol.fact, 2, col = tol.col, lwd = tol.lwd,
            lty = tol.lty)
  }
}

#######################
##  .plot.tclust.2d  ##
#######################

.plot.tclust.2d <-
function (x, xlab, ylab, tol = 0.95, tol.lwd = 1, tol.lty = 3, tol.col = 1, ...) {
  if (nrow (x$centers) != 2)
    stop ("tclust object of dimension 2 expected.")

  if (is.null (x$par$x))
    stop ("dataset not included in tclust object - cannot plot object.")

  dn <- dimnames (x$par$x)
  if(is.list (dn) && length (dn[[2]]) == 2) {
    if (missing (xlab))
      xlab = dn[[2]][1]
    if (missing (ylab))
      ylab = dn[[2]][2]
  }
  else {
    if (missing (xlab))
      xlab = "x1"
    if (missing (ylab))
      ylab = "x2"
  }

  X <- cbind (x$par$x[, 1:2])
  .plot.tclust.0 (x = x, X = X, xlab = xlab, ylab = ylab, axes = 3, ...)

  if (!is.null (x$cov) && is.numeric (tol) && length (tol) == 1 &&  0 < tol && tol < 1)
  {
    tol.col <- rep (tol.col, x$k)
    tol.lty <- rep (tol.lty, x$k)
    tol.lwd <- rep (tol.lwd, x$k)

    tol.fact = sqrt(qchisq(tol, 2))  
    for (k in 1:x$k)
        .doEllipses (eigen = eigen (x$cov[,,k]), center = x$centers[,k],
        lwd = tol.lwd, lty = tol.lty[k], col = tol.col[k], size = tol.fact)
  }
}

#######################
##  .plot.tclust.Nd  ##
#######################

.plot.tclust.Nd <- 
function (x, xlab, ylab, ...)
{
  if (is.null (x$par$x))
    stop ("dataset not included in tclust object - cannot plot object.")

  if (missing (xlab))
    xlab <- "First discriminant coord."
  if (missing (ylab))
    ylab <- "Second discriminant coord."

  X <- discr_coords(x, x$par$equal.weights)

  .plot.tclust.0 (x = x, X = X, xlab = xlab, ylab = ylab, axes = 0, ...)
}

######################
##  .plot.tclust.0  ##
######################

.plot.tclust.0 <-
function (x, X, labels = c ("none", "cluster", "observation"), text,
          xlab, ylab, col, pch, by.cluster = TRUE, axes = 3, xlim, ylim, ...)
{

  if(by.cluster) {
    maxassig <- max (x$cluster)
    
    if (missing (col))
      col <- 1:(x$k+1)
    else
      col <- rep (col,  len = maxassig + 1 )  

    if (missing (pch))
      pch <- 1:(x$k+1)#rep (1, x$k+1)  #
    else
      pch <- rep (pch,  len = maxassig + 1 )
    col <- col[x$cluster + 1]    
    pch <- pch[x$cluster + 1]
  }
  else {
    if (missing (col))
      col <- x$cluster + 1
    if (missing (pch))
      pch <- 1
  }

  n <- x$int$dim[1]

  if (!missing (text))
    text <- rep (text, len = n)
  else if (!missing (labels))    {
    labels <- match.arg(labels)
    if (labels == "cluster")
      text = paste (x$cluster)
    else if (labels == "observation")
      text = paste (1:nrow (X))
  }

  plot.new ()
  par (usr = .plot.tclust.calc.usr (X, xlim, ylim))

  if (missing (text))
    points (X[,1],X[,2], pch = pch, col = col)
  else
    text (X[,1], X[,2], labels = text, col = col)

  .plot.tclust.title (x, ...)

  axis.x <- axes %% 2        ## x axis 1 or 3
  axis.y <- axes >= 2        ## y axis 2 or 3

  cex <- par ("cex")
  if (!missing (xlab))
    mtext (side = 1, xlab, line = 1.5 + 1.5 * axis.x, cex = cex)
  if (!missing (ylab))
    mtext (side = 2, ylab, line = 1.5 + 1.5 * axis.y, cex = cex)

  box ()
  if (axis.x)
    axis (1)
  if (axis.y)
    axis (2)
}

##########################
##  .plot.tclust.title  ##
##########################

.plot.tclust.title <- function (x, main, main.pre, sub, sub1, ...) {
  sub.par <- TRUE
  sub.restr <- FALSE

  sub.ovr <- missing (sub)
  if (!missing (sub) && is.character (sub))   {
    sub.ovr <- TRUE
    if (sub == "/p")
      sub.par <- !(sub.restr <- FALSE)
    else if (sub ==  "/r")
      sub.par <- !(sub.restr <- TRUE)
    else if (sub ==  "/pr")
      sub.par <- sub.restr <- TRUE
    else
      sub.ovr <- FALSE
  }

	if (is.null (x$par$restr.C))
		sub.restr <- FALSE

  ralph <- round (x$par$alpha, 2)

  txt.par <- bquote(paste (k == .(x$par$k), ", ", alpha == .(ralph)))

  if (is.null(x$par$restr.C))
	txt.restr <- ""
  else if (x$par$restr.C == 2) ## this is the sigma - restriction, thus no restr.fact has to be printed.
    txt.restr <- paste ("restr = \"", x$par$restr, "\"", sep = "")
  else
    txt.restr <- paste ("restr = \"", x$par$restr, "\", restr.fact = ", x$par$restr.fact, sep = "")

  if (missing (main))
    main <- "Classification"  #"Cluster Assignment"
  else
    if (is.character (main))
      if (main == "/r" && !is.null (x$par$restr.C))
        main <- txt.restr
      if (main == "/p")
        main <- txt.par

  if (!missing (main.pre) && !is.null (main.pre))
    main <- paste (main.pre, main)

  if (sub.ovr)
    if (sub.par)
      sub  <- txt.par
    else if (!is.null (x$par$restr.C))
      sub <- txt.restr

  if (missing (sub1) && sub.restr && !(!sub.par && sub.ovr))
    sub1 <- txt.restr

  n.sub = .is.visible (sub) + .is.visible (sub1)    ##  number of subtitles to draw

  if (!is.null (main))
    title (main = main, line = ifelse (n.sub > 1, 2.3, 1.6))
  
  if (.is.visible (sub))
    mtext(sub, cex = 0.8, line = ifelse (n.sub > 1, 1.2, 0.3))

  if (.is.visible (sub1))
    mtext(sub1, cex = 0.8, line = ifelse (n.sub > 1, 0.1, 0.3))
}

.is.visible <- function (x)
{
  if (missing (x))
    return (FALSE)
  !is.null (x) && x != ""
}

.plot.tclust.calc.usr <- function (X, xlim, ylim, fact = 0.04)
{
  if (missing (xlim))
    xlim <- range (X[, 1])

  if (missing (ylim))
    ylim <- range (X[, 2])

  r <- cbind (xlim, ylim)
  rd <- apply (r, 2, diff)
#  r + rd * fact * c(-1,1)
  as.numeric (r + (c(-1, 1) %*% t (rd)) * fact)
}
