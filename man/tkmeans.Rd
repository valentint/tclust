\name{tkmeans}
\encoding{latin1}
\alias{ tkmeans }
\alias{ print.tkmeans }
\title{ Trimmed k-means Cluster Analysis  }
\description{ 
  \code{tkmeans} searches for \code{k} (or less) spherical clusters in a data matrix \code{x}, whereas the \code{ceiling (alpha n)} most outlying observations are trimmed.
}
\usage{
tkmeans (x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20, 
         equal.weights = FALSE, center = 0, scale = 1, store.x = TRUE,
         drop.empty.clust = TRUE, trace = 0, warnings = 2, zero.tol = 1e-16)
}
\arguments{
  \item{x}{
    A matrix or data.frame of dimension \code{n} x \code{p}, containing the
    observations (row-wise).
  }
  \item{k}{
    The number of clusters initially searched for.
  }
  \item{alpha}{
    The proportion of observations to be trimmed.
  }
  \item{nstart}{
    The number of random initializations to be performed.
  }
  \item{iter.max}{
    The maximum number of concentration steps to be performed.
    The concentration steps are stopped, whenever two consecutive steps lead
    to the same data partition.
  }
  \item{equal.weights}{
    A logical value, specifying whether equal cluster weights (\code{TRUE}) or
    not (\code{FALSE}) shall be considered in the concentration and assignment 
    steps.
  }
  \item{center, scale}{
    A center and scale vector, each of length \code{p} which can optionally 
    be specified for centering and scaling \code{x} before calculation
  % For centering and scaling \code{x} before calculation, these parameters are 
  % passed to \code{\link[pcaPP]{ScaleAdv}}. Usual choices are 
  % \code{center = mean, scale = sd} for classical, or 
  %  \code{center = median, scale = mad} for robust scaling.
  }
  \item{store.x}{
    A logical value, specifying whether the data matrix \code{x} shall be 
    included in the result structure.
    By default this value is set to \code{TRUE}, because functions
    \code{ \link{plot.tkmeans}} depends on this information.
    However, when big data matrices are handled, the result structure's size 
    can be decreased noticeably when setting this parameter to \code{FALSE}.
  }
  \item{drop.empty.clust}{
    Logical value specifying, whether empty clusters shall be omitted in the 
    resulting object.
    (The result structure does not contain center and covariance estimates of 
    empty clusters anymore.
    Cluster names are reassigned such that the first \code{l} clusters
    (\code{l <= k}) always have at least one observation.
  }
  \item{trace}{
    Defines the tracing level, which is set to \code{0} by default.
    Tracing level \code{2} gives additional information on the iteratively
    decreasing objective function's value.
  }
  \item{warnings}{
    The warning level (0: no warnings; 1: warnings on unexpected behavior.
  }
  \item{zero.tol}{ The zero tolerance used. By default set to 1e-16. }

}
%\details{ }	% 2be written..

\value{
  The function returns an S3 object of type \code{tkmeans}, containing the 
  following values:
  \item{centers}{
    A matrix of size \code{p} x \code{k} containing the centers (column-wise)
    of each cluster.
  }
  \item{cluster}{
    A numerical vector of size \code{n} containing the cluster assignment for
    each observation. Cluster names are integer numbers from \code{1} to
    \code{k}, \code{0} indicates trimmed observations.
  }
  \item{par}{
    A list, containing the parameters the algorithm has been called with
    (\code{x}, if not suppressed by \code{store.x = FALSE},\code{ k},
    \code{ alpha},\code{ restr.fact},\code{ nstart},\code{ KStep}, and
    \code{ equal.weights}).
  }
  \item{k}{
    The (final) resulting number of clusters.
    Some solutions with a smaller number of clusters might be found when using
    the option \code{equal.weights = FALSE}.
  }
  \item{obj}{
    The value of the objective function of the best (returned) solution.
  }
  \item{size}{
    An integer vector of size k, returning the number of observations contained by each cluster.
  }
  \item{weights}{
    A numerical vector of length k, containing the weights of each cluster.
  }
  \item{int}{
    A list of values internally used by function related to \code{tkmeans} objects.
  }
}

\references{ 
 Cuesta-Albertos, J. A.; Gordaliza, A. and Matrán, C. (1997), "Trimmed k-means: an attempt to robustify quantizers". Annals of Statistics, Vol. 25 (2), 553-576. 
 }
\author{ Agustin Mayo Iscar, Luis Angel Garcia Escudero, Heinrich Fritz }

\examples{
\dontshow{
set.seed (0)
}
#--- EXAMPLE 1 ------------------------------------------
sig <- diag (2)
cen <- rep (1,2)
x <- rbind(mvtnorm::rmvnorm(360, cen * 0,   sig),
            mvtnorm::rmvnorm(540, cen * 5,   sig * 6 - 2),
            mvtnorm::rmvnorm(100, cen * 2.5, sig * 50)
            )

# Two groups and 10\% trimming level
clus <- tkmeans (x, k = 2, alpha = 0.1)

plot (clus)
plot (clus, labels = "observation")
plot (clus, labels = "cluster")

#--- EXAMPLE 2 ------------------------------------------
data (geyser2)
clus <- tkmeans (geyser2, k = 3, alpha = 0.03)
plot (clus)

#--- EXAMPLE 3 ------------------------------------------
data (swissbank)
# Two clusters and 8\% trimming level
clus <- tkmeans (swissbank, k = 2, alpha = 0.08)

                            # Pairs plot of the clustering solution
pairs (swissbank, col = clus$cluster + 1)
                                  # Two coordinates
plot (swissbank[, 4], swissbank[, 6], col = clus$cluster + 1,
     xlab = "Distance of the inner frame to lower border",
     ylab = "Length of the diagonal")
plot (clus)

# Three clusters and 0\% trimming level
clus <- tkmeans (swissbank, k = 3, alpha = 0.0)

                            # Pairs plot of the clustering solution
pairs (swissbank, col = clus$cluster + 1)

                                   # Two coordinates
plot (swissbank[, 4], swissbank[, 6], col = clus$cluster + 1, 
      xlab = "Distance of the inner frame to lower border", 
      ylab = "Length of the diagonal")

plot (clus)

}

\keyword{ multivariate }
\keyword{ robust }
\keyword{ cluster }
