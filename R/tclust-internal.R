.doEllipses <- function (eigval, eigvec, eigen, center, cov, n = 100, size = 1, ...)
{
  if (!missing (cov))
    eigen <- base::eigen (cov)
  if (!missing (eigen))
  {
    eigval <- eigen$values
    eigvec <- eigen$vectors
  }

  eigval[eigval < 0] <- 0

                      ##  check dimensionality of eigenvalues
  if (!is.numeric (eigval) || !length (eigval) == 2)
    stop ("argument eigval has to be a numeric vector of length 2.")

                      ##  check dimensionality of center
  if (!is.numeric (center) || !length (center) == 2)
    stop ("argument center has to be a numeric vector of length 2.")

                      ##  check dimensionality of eigenvectors
  if (!is.matrix (eigvec) || any (dim (eigvec) != 2))
    stop ("argument eigvec has to be a numeric mamtrix of dimension 2x2.")

  r <- seq (0, 2 * pi, length.out = n)    ##  create rad rep. of circle

  uc <- rbind (sin(r), cos (r))      ##  create a unit circle
                           ##  "stretch" circle corresponding to ev & sizefact
  uc = t(uc * sqrt(eigval) * size)
                           ##  rotate resulting ellipses from PC into XY coords
  XY = uc %*% t(eigvec)          

  XY = t( t(XY) + center)  ##  move ellipses to the specified center

  lines (XY[,1], XY[,2], ...) ##  draw ellipses
}

.getsubmatrix <- function(x, idx)  
    matrix (x[drop = FALSE,,,idx], nrow = dim (x)[1])

.vline <- function (x, yfact = 2, col = 1, lty = 1, lwd = 1, ...)
{
   ylim = par ("usr")[3:4]
   ylim <- ylim + diff (ylim) * (1 - 1 / yfact) / 2 * c(1,-1)
   
   col = rep (col, length (x))
   lty = rep (lty, length (x))
   lwd = rep (lwd, length (x))
   
   for (i in 1:length (x))
     lines (rep (x[i], 2), ylim, col = col[i], lty = lty [i], lwd = lwd [i]
             , ...)
}

