\name{discr_coords}
\encoding{latin1}
\alias{discr_coords}
\title{ Discriminant coordinates/canonical variates of tclust objects }
\description{ 
Computes the two first discriminant coordinates (canonical coordinates) directly from a \code{tclust} object to obtain a graphical representations of cluster solutions in higher dimensional (p > 2) cases. 
}

\usage{ discr_coords(x, equal.weights)}
\arguments{
  \item{x}{ A \code{tclust} object. }
  \item{equal.weights}{ A logical value, controlling whether the clusters should be considered as equal-sized (\code{TRUE}) when combining their covariance structures, or if their actual size shall be considered (\code{FALSE}). By default value \code{x$par$equal.weights} is assumed}.
}
\details{
  The functionality of \code{discr_coords} is directly derived from \code{discrcoord} as implemented in the package "fpc" by Christian Hennig. It has been adopted in order to directly use the covariance information contained in the \code{tclust}-object. The function fails, if "\code{store.x = FALSE}" is specified in \code{\link[tclust]{tclust}}, because the original data matrix is required here.
}
\value{ A two-dimensional matrix, containing the canonical coordinates of all observations given by the \code{tclust}-object. }
\references{ Hennig, C. and Christlieb, N. (2002), "Validating visual clusters in large datasets: fixed point clusters of spectral features.", Computational Statistics and Data Analysis Vol.40, 723-739. }
\author{ Agustin Mayo Iscar, Luis Angel Garcia Escudero, Heinrich Fritz }
\keyword{ hplot }
\keyword{ multivariate }
\keyword{ robust }
\keyword{ cluster }
