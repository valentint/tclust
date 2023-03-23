\name{plot.tclust}
\encoding{latin1}
\alias{plot.tclust}
\alias{plot.tkmeans}
\title{ plot Method for tclust Objects }
\description{
  The \code{plot} method for classes \code{tclust} and \code{tkmeans}.
}
\usage{
\method{plot}{tclust} (x, ...)
\method{plot}{tkmeans} (x, ...)
}
\arguments{
  \item{x}{
    The \code{tclust} or \code{tkmeans} object to be displayed.
  }
  \item{\dots}{
    Further (optional) arguments which specify the details of the 
    resulting plot (see section "Further Arguments").
  }
  }
\section{Further Arguments}{
 \describe{
  \item{\code{xlab, ylab, xlim, ylim, pch, col}}{
    Arguments passed to \code{plot()}.
  }
  \item{\code{main}}{
    The title of the plot.
    Use "/p" for displaying the chosen parameters \code{alpha} 
    and \code{k} or "/r" for plotting the chosen restriction (\code{tclust} only).
  }
  \item{\code{main.pre}}{
    An optional string which is added to the plot's caption.
  }
  \item{\code{sub}}{
    A string specifying the subtitle of the plot.
    Use "/p" (default) for displaying the chosen parameters \code{alpha} 
    and \code{k}, "/r" for plotting the chosen restriction (\code{tclust} only) and "/pr" for both.
  }
  \item{\code{sub1}}{
    A secondary (optional) subtitle.
  }
  \item{\code{labels}}{
    A string specifying the type of labels to be drawn. Either \code{"none"}
    (default), \code{"cluster"} or \code{"observation"} can be specified.
    If specified, parameter \code{pch} is ignored.
%    If omitted, the symbol, specified by \code{pch} is used.
  }
  \item{\code{text}}{
    A vector of length \code{n} (the number of observations) containing strings
    which are used as labels for each observation.
    If specified, the parameters \code{labels} and \code{pch} are ignored.
  }
  \item{\code{by.cluster}}{
    Logical value indicating whether parameters \code{pch} and \code{col}
    refer to observations (\code{FALSE}) or clusters (\code{TRUE}).
  }
  \item{\code{jitter.y}}{
    Logical value, specifying whether the drawn values shall be jittered in
    y-direction for better visibility of structures in 1 dimensional data.
  }
  \item{\code{tol}}{
    The tolerance interval.
    95\% tolerance ellipsoids (assuming normality) are plotted by default (\code{tclust} only).
  }
  \item{\code{tol.col, tol.lty, tol.lwd}}{
    Vectors of length \code{k} or \code{1} containing the \code{col},
    \code{lty} and \code{lwd} arguments for the tolerance ellipses/lines (\code{tclust} only).
  }
 }
}
\details{
%  The method \code{plot()} of the S3 class \code{tclust} selects the appropriate
%  function according to the dimensionality of the \code{tclust} object
%  (either\code{plot.tclust.1d},\code{plot.tclust.2d} or\code{plot.tclust.Nd}).
%  Internally functions \code{plot.tclust.0} and \code{plot.tclust.title} are 
%  called for doing the actual plotting and formatting the title.
%  
  One and two dimensional structures are treated separately
  (e.g. tolerance intervals/ellipses are displayed).
  Higher dimensional structures are displayed by plotting the two first Fisher's
  canonical coordinates (evaluated by \code{\link[tclust]{discr_coords}})
  and derived from the final cluster assignments (trimmed observations are not
  taken into account). \code{plot.tclust.Nd} can be called with one or
  two-dimensional \code{tclust}-objects too.
  The function fails, if "\code{store.x = FALSE}" is specified in
  \code{\link[tclust]{tclust}}, because the original data matrix is required
  here.
}
%\value{ No return value is provided. }

\author{ Agustin Mayo Iscar, Luis Angel Garcia Escudero, Heinrich Fritz }

\examples{
#--- EXAMPLE 1------------------------------
sig <- diag (2)
cen <- rep (1, 2)
x <- rbind(mvtnorm::rmvnorm(360, cen * 0,   sig),
	       mvtnorm::rmvnorm(540, cen * 5,   sig * 6 - 2),
	       mvtnorm::rmvnorm(100, cen * 2.5, sig * 50)
)
# Two groups and 10\% trimming level
a <- tclust (x, k = 2, alpha = 0.1, restr.fact = 12)

plot (a)
plot (a, labels = "observation")
plot (a, labels = "cluster")
plot (a, by.cluster = TRUE)

#--- EXAMPLE 2------------------------------
x <- c(	rnorm(72,	0,   1),
		rnorm(108, 10,   1),
		rnorm(20,  2.5, 10))

a <- tkmeans (x, k = 2, alpha = 0.1)
plot (a, jitter.y = TRUE)

}
\keyword{ hplot }
\keyword{ multivariate }
\keyword{ robust }
\keyword{ cluster }
