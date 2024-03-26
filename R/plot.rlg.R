#' Plot an 'rlg' object
#'
#' @description Different plots for the results of 'rlg' analysis, stored in an 
#'    \code{rlg} object, see Details.
#'
#' @param x An \code{rlg} object to plot.
#' @param which Select the required plot.
#' @param sort Whether to sort.
#' @param ask if \code{TRUE}, the user is \emph{ask}ed before each plot, see \code{par(ask=.)}. 
#'  Default is \code{ask = which=="all" && dev.interactive()}.  
#' @param \dots Other parameters to be passed to the lower level functions.
#'
#' @examples
#'  data (LG5data)
#'  x <- LG5data[, 1:10]
#'  clus <- rlg(x, d = c(2,2,2), alpha=0.1)
#'  plot(clus, which="eigenvalues") 
#'  plot(clus, which="scores") 
#'
plot.rlg <- function(x, which=c("all", "scores", "loadings", "eigenvalues"), 
    sort = TRUE,
    ask = (which=="all" && dev.interactive(TRUE)), ...){
    
    ## option "clusters" removed
    which <- match.arg(which)
      
    if(!inherits(x, "rlg"))
        stop("Use only with 'rlg' objects")
    
    op <- if (ask) par(ask = TRUE) else list()
    on.exit(par(op))
    
    if(which == "all" || which == "scores")
        .plot.scores(x)
    if(which == "all" || which == "loadings")
        .plot.loadings(x)
    if(which == "all" || which == "eigenvalues")
        .plot.eigenvalues(x)
    
    invisible(x)
}

.plot.clusters <- function(x, sort){
  index <- .sort(x, sort)
  
  p <- ncol(x$x)
  mat <- matrix(rep(x$cluster, each=p), ncol=p, byrow=T)
  
  rownames(mat) <- .get.rownames(x)
  colnames(mat) <- .get.colnames(x)
  
  heatmap(mat[index,], Rowv = NA, Colv = NA, scale = "none", col = 1:(length(x$d)+1))
}

.plot.scores <- function(x){
  
  K <- which(x$dimensions >= 2)
  
  if(length(K) > 0) {
      par(mfrow=rev(n2mfrow(length(K))))
      for(k in K){
        scores <- princomp(x$x[x$cluster == k, ],cor = TRUE)$scores
        plot(scores[,1], scores[,2],
             main=paste("Cluster", k), xlab="dim 1", ylab="dim 2", col="white")
        text(scores[,1], scores[,2],
             as.character(substr(.get.rownames(x)[x$cluster==k], 1, 3)))
      }
  }
}

.plot.loadings <- function(x){
  dim <- rev(n2mfrow(length(x$dimensions)))
  d.max <- min(max(x$dimensions),3)
  dim[1] <- dim[1]*d.max
  mat <- matrix(1:(dim[1]*dim[2]), nrow = dim[1], ncol = dim[2])
  layout(mat)
  for(k in 1:length(x$dimensions)){
    n2 <- min(x$dimensions[k], 3)
    for(i2 in 1:d.max){
      if(i2 <= n2){
        title <- ifelse(i2 == 1, paste("Cluster", k, "\nDim", i2) , paste("Dim", i2))
        barplot(x$U[[k]][,i2], main=title)
      }else{
        plot(0,type='n',axes=FALSE,ann=FALSE)
      }
    }
  }
  par(mfrow=c(1,1))
}

.plot.eigenvalues <- function(x){
  p <- ncol(x$x)
  eigenvalues <- matrix(NA,nrow=length(x$dimensions),ncol=p)
  for (k in 1:length(x$dimensions)){
    eigenvalues[k,] <- princomp(x$x[x$cluster == k, ],cor = TRUE)$sdev^2
  }
  lim.sup <- max(eigenvalues)
  plot(eigenvalues[1,],ylim=c(0,lim.sup),type="n",xlab="intrisic dimensions",ylab="eigenvalues")
  for (k in 1:length(x$dimensions)){
    lines(eigenvalues[k,],col=k+1,type="b",pch=19)
  }
  legend("topright",
         legend=paste("Cluster",as.character(1:length(x$dimensions))),
         col=1:length(x$dimensions)+1, lty=rep(1,length(x$dimensions)),
         pch=rep(19,length(x$dimensions)))
}

.sort <- function(x, sort){
  if(sort)
    index <- order(x$cluster)
  else
    index <- 1:length(x$cluster)
  
  return(index)
}

.get.rownames <- function(x){
  if(!is.null(rownames(x$x)))
    return(rownames(x$x))
  else
    return(1:nrow(x$x))
}

.get.colnames <- function(x){
  if(!is.null(colnames(x$x)))
    return(colnames(x$x))
  else
    return(1:ncol(x$x))
}
