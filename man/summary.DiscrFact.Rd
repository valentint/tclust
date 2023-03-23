\name{summary.DiscrFact}
\encoding{latin1}
\alias{summary.DiscrFact}
\title{ summary Method for DiscrFact Objects}
\description{ The \code{\link[base]{summary}} method for class \code{DiscrFact}. }
\usage{
\method{summary}{DiscrFact} (object, hide.emtpy = TRUE, show.clust, show.alt, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class "\code{DiscrFact}" as from \code{DiscrFact ()}. }
  \item{ hide.emtpy }{ A logical value specifying whether clusters without doubtful assignment shall be hidden. }
  \item{ show.clust }{ A logical value specifying whether the number of doubtful assignments per cluster shall be displayed. }
  \item{ show.alt }{ A logical value specifying whether the alternative cluster assignment shall be displayed. }
  \item{\dots}{ Arguments passed to or from other methods. }
}
%\details{}
\value{ No return value is provided.}
\references{ García-Escudero, L.A.; Gordaliza, A.; Matrán, C. and Mayo-Iscar, A. (2009), "Exploring the number of groups in robust model-based clustering".\cr
Preprint available at www.eio.uva.es/infor/personas/langel.html. }
\author{ Agustin Mayo Iscar, Luis Angel Garcia Escudero, Heinrich Fritz }
%\note{ ~~further notes~~  ~Make other sections like Warning with \section{Warning }{....} ~}
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
summary(dsc.1)

dsc.2 <- DiscrFact (clus.2)
summary (dsc.2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ cluster }
\keyword{ robust }
\keyword{ multivariate }
