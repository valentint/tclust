\name{DiscrFact}
\encoding{latin1}
\alias{DiscrFact}
\alias{print.DiscrFact}
\title{ Discriminant Factor Analysis for tclust Objects}
\description{ Analyzes a \code{tclust}-object by calculating discriminant factors and comparing the quality of the actual cluster assignments and the second best possible assignment for each observation. Discriminant factors, measuring the strength of the "trimming" decision may also be defined. % (!? what does this sentence mean?). 
Cluster assignments of observations with large discriminant factors are considered as "doubtful" decisions. Silhouette plots give a graphical overview of the discriminant factors distribution (see \code{\link{plot.DiscrFact}}). More details can be found in García-Escudero et al. (2010).}
\usage{
DiscrFact(x, threshold = 1/10)
}
\arguments{
  \item{x}{ A \code{tclust} object. }
  \item{threshold}{ %Observations with discriminant factors larger than \code{log (threshold)} are considered as to be "doubtfully" assigned to a cluster.
  A cluster assignment or a trimming decision for an observation with a discriminant factor larger than \code{log(threshold)} is considered as a "doubtful" decision. 
   }
}
\details{
This function compares the actual (best) assignment of each observation to its second best possible assignment. This comparison is based on the discriminant factors of each observation, which are calculated here. If the discriminant factor of an observation is larger than a given level (\code{log (threshold)}), the observation is considered as "doubtfully" assigned to a cluster. More information is shown when \link[=plot.DiscrFact]{plotting} the returned \code{DiscrFact} object.
}
\value{
The function returns an S3 object of type \code{DiscrFact} containing the following components:
	\item{x}{ A \code{tclust} object. }
	\item{ylimmin}{ A minimum y-limit calculated for plotting purposes. }
	\item{ind}{ The actual cluster assignment. }
	\item{ind2}{ The second most likely cluster assignment for each observation. }
	\item{disc}{ The (weighted) likelihood of the actual cluster assignment of each observation. }
	\item{disc2}{ The (weighted) likelihood of the second best cluster assignment of each observation. }
	\item{assignfact}{ The factor \code{log (disc/disc2)}. }
	\item{threshold}{ The threshold used for deciding whether \code{assignfact} indicates a "doubtful" assignment. }
	\item{mean.DiscrFact}{ A vector of length \code{k + 1} containing the mean discriminant factors for each cluster (including the outliers). }
}
\references{ 
    \enc{García}{Garcia}-Escudero, L.A.; Gordaliza, A.; \enc{Matrán}{Matran}, C. and Mayo-Iscar, A. (2011), 
    "Exploring the number of groups in robust model-based clustering." \emph{Statistics and Computing}, \bold{21} 
    pp. 585-599, <doi:10.1007/s11222-010-9194-z>
}

\author{ 
    Agustin Mayo-Iscar, Luis Angel García-Escudero, Heinrich Fritz 
}

\seealso{ \code{\link{plot.DiscrFact}}}
\examples{
sig <- diag (2)
cen <- rep (1, 2)
x <- rbind(mvtnorm::rmvnorm(360, cen * 0,   sig),
	       mvtnorm::rmvnorm(540, cen * 5,   sig * 6 - 2),
	       mvtnorm::rmvnorm(100, cen * 2.5, sig * 50)
)
clus.1 <- tclust (x, k = 2, alpha = 0.1, restr.fact = 12)

clus.2 <- tclust (x, k = 3, alpha = 0.1, restr.fact = 1)
  ##  restr.fact and k are chosen improperly for pointing out the 
  ##    difference in the plot of DiscrFact

dsc.1 <- DiscrFact (clus.1)
plot(dsc.1)

dsc.2 <- DiscrFact (clus.2)
plot (dsc.2)

}
\keyword{ hplot }
\keyword{ cluster }
\keyword{ robust }
\keyword{ multivariate }
