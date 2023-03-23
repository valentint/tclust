\name{plot.DiscrFact}
\encoding{latin1}
\alias{plot.DiscrFact}
\alias{plot_DiscrFact_p2}
\alias{plot_DiscrFact_p3}
\title{ plot Method for DiscrFact Objects }
\description{
The \code{plot} method for class \code{DiscrFact}: Next to a plot of the \code{tclust} object which has been used for creating the \code{DiscrFact} object, a silhouette plot indicates the presence of groups with a large amount of doubtfully assigned observations. A third plot similar to the standard \code{tclust} plot serves to highlight the identified doubtful observations.
}
\usage{
\method{plot}{DiscrFact} (x, enum.plots = FALSE, ...)
plot_DiscrFact_p2 (x, xlab = "Discriminant Factor", 
                   ylab = "Clusters", main, xlim, 
                   print.Discr = TRUE, main.pre, ...)
                  
                  
                  
plot_DiscrFact_p3 (x, main = "Doubtful Assignments", col, pch, 
                   col.nodoubt = grey (0.8), by.cluster = FALSE, 
                   ...)
                   

}
\arguments{
  \item{x}{ An object of class "\code{DiscrFact}" as from \code{DiscrFact ()}. }
  \item{enum.plots}{ A logical value indicating whether the plots shall be enumerated in their title ("(a)", "(b)", "(c)"). }
  \item{xlab, ylab, xlim}{ Arguments passed to funcion \code{\link[tclust]{plot.tclust}}. }
  \item{main}{ Argument passed to funcion \code{plot()}. }
%  \item{draw.legend}{ A logical value indicating whether the legend shall be displayed. }
  \item{print.Discr}{ A logical value indicating whether each clusters mean discriminant factor shall be plotted }
  \item{main.pre}{An optional string which is appended to the plot's caption.}
  \item{pch, col}{Arguments passed to function \code{plot()}.}
  \item{col.nodoubt}{Color of all observations not considered as to be assigned doubtfully.}
  \item{by.cluster}{Logical value indicating whether parameters \code{pch} and \code{col} refer to observations (\code{FALSE}) or clusters (\code{TRUE}).}
  \item{\dots}{ Arguments to be passed to or from other methods. }
}
\details{
\code{plot.DiscrFact.p2} displays a silhouette plot based on the discriminant factors of the observations. A solution with many large discriminant factors is not reliable. Such clusters can be identified with this silhouette plot. %Thus, when comparing clustering solutions, small discriminant factors are preferred. 
Thus \code{plot.DiscrFact.p3} displays the dataset, highlighting observations with discriminant factors greater than the given threshold. Function \code{plot.DiscrFact} combines the standard plot of a \code{tclust} object, and the two plots introduced here.
}
\value{ No return value is provided.}
\references{
%\enc{García}{Garcia}-Escudero, L.A.; Gordaliza, A.; \enc{Matrán}{Matran}, C. and Mayo-Iscar, A. (2009), "Exploring the number of groups in robust model-based clustering".\cr
%Preprint available at www.eio.uva.es/infor/personas/langel.html.
\enc{García}{Garcia}-Escudero, L.A.; Gordaliza, A.; \enc{Matrán}{Matran}, C. and Mayo-Iscar, A. (2010), "Exploring the number of groups in robust model-based clustering." Statistics and Computing, (Forthcoming). \cr
Preprint available at www.eio.uva.es/infor/personas/langel.html.
}
\author{ Agustin Mayo Iscar, Luis Angel \enc{García}{Garcia} Escudero, Heinrich Fritz }
%\note{ ~~further notes~~ ~Make other sections like Warning with \section{Warning }{....} ~}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
sig <- diag (2)
cen <- rep (1, 2)
x <- rbind(mvtnorm::rmvnorm(360, cen * 0,   sig),
	       mvtnorm::rmvnorm(540, cen * 5,   sig * 6 - 2),
	       mvtnorm::rmvnorm(100, cen * 2.5, sig * 50)
)

clus.1 <- tclust (x, k = 2, alpha=0.1, restr.fact=12)
clus.2 <- tclust (x, k = 3, alpha=0.1, restr.fact=1)

dsc.1 <- DiscrFact (clus.1)
plot(dsc.1)

dsc.2 <- DiscrFact (clus.2)
plot (dsc.2)

dev.off ()
plot_DiscrFact_p2 (dsc.1)
plot_DiscrFact_p3 (dsc.2)
}
\keyword{ hplot }
\keyword{ cluster }
\keyword{ robust }
\keyword{ multivariate }
