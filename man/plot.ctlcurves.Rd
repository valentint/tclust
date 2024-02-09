\name{plot.ctlcurves}
\alias{plot.ctlcurves}
\encoding{latin1}
\title{ plot Method for ctlcurves Objects }
\description{ The \code{plot} method for class \code{ctlcurves}: This function plots a \code{ctlcurves} object, comparing the target functions values with different values of parameter \code{restr.fact}.
%These functions plot a ctlcurves object in two different styles, either classic or simple. The classical style displays three rows of diagnostic plots, analyzing different clustering results on the same data set. The first row shows the target functions value. The second row displays the increase of the target functions value gained by increasing k by one, whereas the third row shows the minimum cluster weights for each estimation. The simple plots only display the target functions values, comparing the results of different values of parameter \code{restr.fact}.
}
\usage{
\method{plot}{ctlcurves} (x, what = c ("obj", "min.weights", "doubtful"),
         main, xlab, ylab, xlim, ylim, col, lty = 1, ...)
}
\arguments{
  \item{x}{ The ctlcurves object to be printed. }
  \item{what}{A string indicating which type of plot shall be drawn. See the details section for more information.}
%  \item{simple}{ Logical value indicating whether the simple plot (\code{plot.ctlcurves.simple}) or the classical plot (\code{plot.ctlcurves.classic}) shall be displayed. If the simulation has been done for more than one factor value, by default the simple plot is displayed for each factor. If only one factor has been chosen, the classical plot is displayed by default. }
%  \item{plot.idx}{ An index vector % (\code{plot.ctlcurves.simple}) or matrix (\code{plot.ctlcurves.classic})
%  giving the indices of the plots which shall be displayed. If provided, no further plot arrangement (\code{\link[graphics]{par} (mfrow = ...)}) is done. }
%  \item{link.ylim}{ A logical value indicating whether the y-ranges of the plots shall be identical. Ignored, if parameter ylim is set. }
%  \item{ylim}{ The y-range of all plots displayed. If provided, parameter \code{link.ylim} is ignored. }
  \item{main}{ A character-string containing the title of the plot. }
  \item{xlab, ylab, xlim, ylim}{ Arguments passed to \code{plot()}. }
%  \item{set.idx}{ If the \code{ctlcurves} - object contains more than one factor, only data related to \code{factor[set.idx]} is displayed. }
%  \item{min.weights}{Locigal value indicating whether the minimum cluster weights shall be printed instead of the target functions values. By default this option is disabled}
  \item{col}{ A single value or vector of line colors passed to \code{\link[graphics]{lines}}. }
  \item{lty}{ A single value or vector of line types passed to \code{\link[graphics]{lines}}. }
   \item{\dots}{ Arguments to be passed to or from other methods. }
}
\details{ These curves show the values of the trimmed classification 
          (log-)likelihoods when altering the trimming proportion \code{alpha}
	  and the number of clusters \code{k}.
	  The careful examination of these curves provides valuable 
	  information for choosing these parameters in a clustering problem.
	  For instance, an appropriate \code{k} to be chosen is one that we do
	  not observe a clear increase in the trimmed classification likelihood
	  curve for k with respect to the k+1 curve for almost all the range of
	  alpha values. Moreover, an appropriate choice of parameter alpha may
	  be derived by determining where an initial fast increase of the
	  trimmed classification likelihood curve stops for the final chosen k.
	  A more detailed explanation can be found in García-Escudero et al. (2010).
	  
This function implements a series of plots, which display characteristic values of the each model, computed with different values for \code{k} and \code{alpha}.
The plot type is selected by setting argument \code{what} to one of the following values:
  \describe{
  \item{\code{"obj"}}{Objective function values.}
%  \item{\code{"out"}}{The number of outlying observations identified among the 
%               clusters. 
%	       Here observations are counted, which were assigned to a
%	       cluster, and within this cluster yield a squared 
%	       Mahalanobis distance larger than
%               \code{qchisq (1-mah.alpha, p)}, with
%	       \code{p} the number of dimensions of the data matrix \code{x},
%	       and \code{mah.alpha} the ?? level as chosen in 
%	       function \code{\link{ctlcurves}}.\cr
%	       A dashed line represents the chosen ?? level.
%	       }
  \item{\code{"min.weights"}}{The minimum cluster weight found for each computed model.
               This plot is intended to spot spurious clusters, which in 
	       general yield quite small weights.
	       }
  \item{\code{"doubtful"}}{The number of "doubtful" decisions identified by 
                  \code{\link{DiscrFact}}. }
  }
}
%\value{ no value is returned }
\references{ %García-Escudero, L.A.; Gordaliza, A.; Matrán, C. and Mayo-Iscar, A. (2009), "Exploring the number of groups in robust model-based clustering".\cr
%Preprint available at www.eio.uva.es/infor/personas/langel.html. 
\enc{García}{Garcia}-Escudero, L.A.; Gordaliza, A.; \enc{Matrán}{Matran}, C. and Mayo-Iscar, A. (2010), "Exploring the number of groups in robust model-based clustering." Statistics and Computing, (Forthcoming). \cr
Preprint available at www.eio.uva.es/infor/personas/langel.html.
}
\author{ Agustin Mayo Iscar, Luis Angel Garcia Escudero, Heinrich Fritz }
\examples{
\dontshow{
set.seed(0)
}
\dontrun{
sig <- diag (2)
cen <- rep (1, 2)
x <- rbind(mvtnorm::rmvnorm(108, cen * 0,   sig),
	       mvtnorm::rmvnorm(162, cen * 5,   sig * 6 - 2),
	       mvtnorm::rmvnorm(30, cen * 2.5, sig * 50)
)

ctl <- ctlcurves(x, k = 1:4)
plot(ctl)

}
}
\keyword{ hplot }
\keyword{ cluster }
\keyword{ robust }
\keyword{ multivariate }
