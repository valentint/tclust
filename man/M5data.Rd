% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{M5data}
\alias{M5data}
\title{M5data data}
\format{
The first two columns are the two variables. The last column is the true 
 classification vector where symbol "0" stands for the contaminating data points.
}
\source{
García-Escudero, L.A.; Gordaliza, A.; Matrán, C. and Mayo-Iscar, A. (2008), 
 "A General Trimming Approach to Robust Cluster Analysis". Annals of Statistics, 
 Vol.36, pp. 1324-1345.
}
\usage{
data(M5data)
}
\description{
A bivariate data set obtained from three normal bivariate distributions with 
 different scales and proportions 1:2:2. One of the components is very overlapped 
 with another one. A 10\% background noise is added uniformly distributed in a rectangle 
 containing the three normal components and not very overlapped with the three mixture 
 components. A precise description of the M5 data set can be found in 
 García-Escudero et al. (2008).
}
\examples{
#--- EXAMPLE 1 ------------------------------------------ 
data (M5data) 
x <- M5data[, 1:2] 
clus <- tclust(x, k=3, alpha=0.1, nstart=200, niter1=3, niter2=17, 
   nkeep=10, opt="HARD", equal.weights=FALSE, restr.fact=50, trace=TRUE) 
plot (x, col=clus$cluster+1)
}
\keyword{datasets}
