if(FALSE)
{
    library(tclust)
    library(MASS)
    
    x <- as.matrix(iris[, -5])
    dm1 <- tclust:::dmvnrm(x, colMeans(x), cov(x))
    dm2 <- .dmnorm(x, colMeans(x), cov(x))
    dm3 <- mvtnorm::dmvnorm(x, colMeans(x), cov(x))
    cbind.data.frame(dm1=dm1, dm2=dm2, dm3=dm3, diff1=round(dm1-dm2, 12), diff2=round(dm1-dm3, 12))
}

##  Multivariate normal density
.dmnorm <- function (X, mu, sigma)
{
  ((2 * pi)^(-length(mu) / 2)) *
  (det(sigma)^(-1/ 2)) *
  exp (-0.5 * .evmaha (X, mu, sigma))
}

.evmaha <- function (X, mu, sigma)    ##  calculate mahalanobis distances 
                                      ##  using the eigenvalues and eigenvectors. 
                                      ##  thus no singularity problems arise.
{                                     ##  Mahalanobis distance == Inf is possible.
    v <- eigen (sigma)
    Zt <- t (v$vectors) %*% (t (X) - mu)
    colSums ((Zt * v$values^(-.5))^2)
}


