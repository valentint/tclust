\name{ctlcurves}
\encoding{latin1}
\alias{ctlcurves}
\alias{print.ctlcurves}
%\alias{summary.ctlcurves}
\title{ Classification Trimmed Likelihood Curves }
\description{
	The function applies \code{\link[tclust]{tclust}} several times on a given dataset while parameters \code{alpha} and \code{k}%and \code{restr.fact}
	 are altered. The resulting object gives an idea of the optimal trimming level and number of clusters considering a particular dataset. 
}
\usage{
ctlcurves (x, k = 1:4, alpha = seq (0, 0.2, len = 6), 
           restr.fact = 50, trace = 1, ...)
}
\arguments{
  \item{x}{ A matrix or data frame of dimension \code{n} x \code{p}, containing the observations (row-wise). }
  \item{k}{ A vector of cluster numbers to be checked. By default cluster numbers from 1 to 5 are examined. }
  \item{alpha}{ A vector containing the alpha levels to be checked. By default alpha levels from 0 to 0.2 (continuously increased by 0.01), are checked. }
%  \item{mah.alpha}{??}
%  \item{k.max}{ The maximal number of clusters to be checked (\code{k} is iterated from 1 to k.max). }
%  \item{alpha.count}{ The length of the grid of trimming levels to be examined. }
%  \item{alpha.max}{ The maximum trimming level to be checked. (trimming levels are calculated as seq(0, alpha.max, length.out = alpha.count) }
  %\item{restr.fact}{ A vector of factors constraining the ratio of clusters eigenvalues. Each single value of this vector is  passed to \code{\link{tclust}}.}
  \item{restr.fact}{ The restriction factor passed to \code{\link{tclust}}. }
  \item{\dots}{ Further arguments (as e.g. \code{restr}), passed to \code{\link[tclust]{tclust}}. }
  \item{trace}{ Defines the tracing level, which is set to \code{1} by default. Tracing level \code{2} gives additional information on the current iteration. }
%  \item{warnings}{ The warning level (0: no warnings, 1: warnings if \code{restr.fact} has been chosen too small)}
}
\details{ These curves show the values of the trimmed classification (log-)likelihoods when altering the trimming proportion \code{alpha} and the number of clusters \code{k}. The careful examination of these curves provides valuable information for choosing these parameters in a clustering problem. For instance, an appropriate \code{k} to be chosen is one that we do not observe a clear increase in the trimmed classification likelihood curve for k with respect to the k+1 curve for almost all the range of alpha values. Moreover, an appropriate choice of parameter alpha may be derived by determining where an initial fast increase of the trimmed classification likelihood curve stops for the final chosen k. A more detailed explanation can be found in García-Escudero et al. (2010). }
\value{ The function returns an S3 object of type \code{ctlcurves} with components: %Functions \code{\link[=plot.ctlcurves]{plot}}, \code{print} and \code{summary} are implemented for this type. 
	\item{par}{ A list containing all the parameters passed to this function. }
	\item{obj}{ An array containing the objective functions values of each computed cluster-solution. }
	\item{min.weights}{ An array containing the minimum cluster weight of each computed cluster-solution. }
%	\item{bestassig}{ The \code{tclust} object which provided the best cluster-solution (maximizing the target function). }
%	\item{aa}{The different levels of alpha which have been tested. }
	So far there is no output available for \code{print.ctlcurves}. Use \code{plot} on an \code{ctlcurves} object for a graphical interpretation of it.
}
\references{ %García-Escudero, L.A.; Gordaliza, A.; Matrán, C. and Mayo-Iscar, A. (2009), "Exploring the number of groups in robust model-based clustering".\cr
%Preprint available at www.eio.uva.es/infor/personas/langel.html. 
\enc{García}{Garcia}-Escudero, L.A.; Gordaliza, A.; \enc{Matrán}{Matran}, C. and Mayo-Iscar, A. (2010), "Exploring the number of groups in robust model-based clustering." Statistics and Computing, (Forthcoming). \cr
Preprint available at www.eio.uva.es/infor/personas/langel.html.
}
\author{ Agustin Mayo Iscar, Luis Angel Garcia Escudero, Heinrich Fritz }
%\note{ ~~further notes~~  ~Make other sections like Warning with \section{Warning }{....} ~ }
\seealso{ \code{\link{plot.ctlcurves}} }
\examples{
\dontshow{set.seed (0)}
\dontrun{
#--- EXAMPLE 1 ------------------------------------------

sig <- diag (2)
cen <- rep (1, 2)
x <- rbind(mvtnorm::rmvnorm(108, cen * 0,   sig),
	       mvtnorm::rmvnorm(162, cen * 5,   sig * 6 - 2),
	       mvtnorm::rmvnorm(30, cen * 2.5, sig * 50)
)

ctl <- ctlcurves (x, k = 1:4)

  ##  ctl-curves 
plot (ctl)  ##  --> selecting k = 2, alpha = 0.08

  ##  the selected model 
plot (tclust (x, k = 2, alpha = 0.08, restr.fact = 7))

#--- EXAMPLE 2 ------------------------------------------

data (geyser2)
ctl <- ctlcurves (geyser2, k = 1:5)

  ##  ctl-curves 
plot (ctl)  ##  --> selecting k = 3, alpha = 0.08

  ##  the selected model
plot (tclust (geyser2, k = 3, alpha = 0.08, restr.fact = 5))


#--- EXAMPLE 3 ------------------------------------------

data (swissbank)
ctl <- ctlcurves (swissbank, k = 1:5, alpha = seq (0, 0.3, by = 0.025))

  ##  ctl-curves 
plot (ctl)  ##  --> selecting k = 2, alpha = 0.1

  ##  the selected model
plot (tclust (swissbank, k = 2, alpha = 0.1, restr.fact = 50))
}

}
\keyword{ hplot }
\keyword{ multivariate }
\keyword{ robust }
\keyword{ cluster }
