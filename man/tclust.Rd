\name{tclust}
\encoding{latin1}
\alias{ tclust }
%\alias{ tkmeans }
\alias{ print.tclust }
\title{ General Trimming Approach to Robust Cluster Analysis  }
\description{ 
  \code{tclust} searches for \code{k} (or less) clusters with different
  covariance structures in a data matrix \code{x}.
  Relative cluster scatter can be restricted by a constant value
  \code{restr.fact}. For robustifying the estimation, a proportion
  \code{alpha} of observations may be trimmed. 
%  \code{tkmeans} implements a robust version of the data driven trimming
%  method "k-means".
  In particular, the trimmed k-means method (\code{\link[tclust]{tkmeans}})is represented by the \code{tclust}
  method, setting parameters \code{restr = "eigen"}, \code{restr.fact = 1}
  and \code{equal.weights = TRUE}.
}
\usage{

tclust (x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20, 
        restr = c ("eigen", "deter", "sigma"), restr.fact = 12, 
        equal.weights = FALSE, center, scale, store.x = TRUE, 
        drop.empty.clust = TRUE, trace = 0, warnings = 3, 
        zero.tol = 1e-16)
%tkmeans (x, k = 3, alpha = 0.05, nstart = 50, iter.max = 20, 
%        center, scale,  store.x = TRUE, drop.empty.clust = TRUE, 
%        trace = 0, warnings = 2, zero.tol = 1e-16)
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
%  \item{fuzzy}{
%    A logical variable indicating whether the ``fuzzy´´ algorithm shall be
%    applied.
%  }
%  \item{m}{
%    A tuning parameter for the fuzzy algorithm.
%  }
  \item{nstart}{
    The number of random initializations to be performed.
  }
  \item{iter.max}{
    The maximum number of concentration steps to be performed.
    The concentration steps are stopped, whenever two consecutive steps lead
    to the same data partition.
  }
  \item{restr}{
    The type of restriction to be applied on the cluster scatter matrices.
    Valid values are \code{"eigen"} (default), \code{"deter"} and
    \code{"sigma"}.
    % , \code{"dir.eigen"}, \code{"dir.deter"} and \code{"prop"}.
    See the detail section for further explanation. }
  \item{restr.fact}{
    The constant \code{restr.fact >= 1} constrains the allowed differences
    among group scatters.
    Larger values imply larger differences of group scatters, a value of 1
    specifies the strongest restriction.
    When using \code{restr = "sigma"} this parameter is not considered, as all
    cluster variances are averaged, always implying \code{restr.fact = 1}.
    % The maximum of the quotient between the maximum and minimum eigenvalues
    % of the cluster's covariance matrices.
    % If the covariance structures differ too much, they are adjusted regarding
    % this factor.
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
    \code{ \link{plot.tclust}} and\code{ \link{DiscrFact}} depend on this 
    information.
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
%  \item{throwwarnings}{
%    Whether warnings shall be thrown due to unexpected algorithm behavior.
%    If the algorithm is used in a different context 
%    (e.g. \code{\link{graphical.tclust}}), the calling function might handle 
%    the warnings itself, not passing them to the user.
%  }
  \item{trace}{
    Defines the tracing level, which is set to \code{0} by default.
    Tracing level \code{2} gives additional information on the iteratively
    decreasing objective function's value.
  }
%  \item{iter.tune}{
%    A vector of three integer values, specifying the number of iterations for 
%    restrictions \code{"dir.eigen"}, \code{"dir.deter"} and \code{"prop"}.
%    See details for further explanation.
%  }
  \item{warnings}{
    The warning level (0: no warnings; 1: warnings on unexpected behavior;
    2: warnings if \code{restr.fact} causes artificially restricted results).
  }
  \item{zero.tol}{ The zero tolerance used. By default set to 1e-16. }
%  \item{dropclust}{
%    A logical value indicating whether too small clusters (covariance 
%    structure of rank < p) shall be dropped and the algorithm shall 
%    continue with k - 1 clusters.
%    If set to \code{FALSE} the current solution is discarded, and a new 
%    initialization is considered.
%  }
%  \item{scale}{
%    If the data shall be scaled before calculation, this parameter specifies
%    an according scale function (e.g. sd, mad,\dots) or a vector of length
%    \code{p} or 1.
%  }
}
\details{
  This iterative algorithm initializes \code{k} clusters randomly and performs
  "concentration steps" in order to improve the current cluster assignment.
  The number of maximum concentration steps to be performed is given by
  \code{iter.max}.
  For approximately obtaining the global optimum, the system is initialized
  \code{nstart} times and concentration steps are performed until convergence
  or \code{iter.max} is reached.
  When processing more complex data sets higher values of \code{nstart} and 
  \code{iter.max} have to be specified 
  (obviously implying extra computation time).
  However, if more then half of the iterations would not converge, a warning
  message is issued, indicating that \code{nstart} has to be increased.

  The parameter \code{restr} defines the cluster's shape restrictions, which 
  are applied on all clusters during each iteration.
  Options \code{"eigen"}/\code{"deter"} restrict the ratio between the maximum 
  and minimum eigenvalue/determinant of all cluster's covariance structures to 
  parameter \code{restr.fact}. Setting \code{restr.fact} to \code{1}, yields 
  the strongest restriction, forcing all eigenvalues/determinants to be equal 
  and so the method looks for similarly scattered (respectively spherical)
  clusters.
  Option \code{"sigma"} is a simpler restriction, which averages the covariance 
  structures during each iteration (weighted by cluster sizes) in order to get 
  similar (equal) cluster scatters.
  % Option \code{"dir.eigen"}/\code{"dir.deter"} is searching for clusters 
  % with equal eigenvectors whereas the ratio between maximum/minimum 
  % eigenvalue/determinant is restricted by the function argument
  % \code{restr.fact}.
  % Option \code{"prop"} searches for proportional cluster structures whereas
  % the cluster's determinants are restricted by argument \code{restr.fact}.
  % Internally the latter three restrictions implement three nested iterations.
  % The number of these iterations is controlled by the argument 
  % \code{iter.tune} which is supposed to be an integer vector of length 3.
  % By default each iterative step is carried out 20 times.
  % The parameter \code{restr.fact} restricts the relative shape of the
  % clusters. 
  % It represents the maximum possible ratio between the greatest eigenvalues
  % of the group covariance matrices of all clusters, divided by the smallest
  % of these eigenvalues. If this ratio in any concentration step is too big,
  % eigenvalues of the clusters are changed in order to get more similar group
  % scatters.
  % Setting \code{restr.fact} to \code{1}, yields the strongest restriction,
  % forcing all eigenvalues to be equal and so the method looks for similarly
  % scattered spherical clusters.
}

\value{
  The function returns an S3 object of type \code{tclust}, containing the 
  following values:
  \item{centers}{
    A matrix of size \code{p} x \code{k} containing the centers (column-wise)
    of each cluster.
  }
  \item{cov}{
    An array of size  \code{p} x \code{p} x \code{k} containing the covariance 
    matrices  of each cluster.
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
    A list of values internally used by function related to \code{tclust} objects.
  }
%  \item{iter.successful}{
%    The number of successful iterations.
%    If \code{droclust = FALSE} is specified, some iterations may fail due to 
%    too small cluster sizes.
%    When \code{k} is chosen too high, this value might decrease to zero!
%  }
%  \item{iter.converged}{
%    The number of converged iterations.
%  }
%  \item{z}{
%    If \code{"fuzzy = TRUE"} has been selected this value contains a matrix 
%    with the fuzzy cluster pertinences.
%}
}

\references{ 
Garcia-Escudero, L.A.; Gordaliza, A.; Matran, C. and Mayo-Iscar, A. (2008), "A General Trimming Approach to Robust Cluster Analysis". Annals of Statistics, Vol.36, 1324-1345. Technical Report available at www.eio.uva.es/inves/grupos/representaciones/trTCLUST.pdf\cr
Fritz, H.; Garcia-Escudero, L.A.; Mayo-Iscar, A. (2012), "tclust: An R Package for a Trimming Approach to Cluster Analysis". Journal of Statistical Software, 47(12), 1-26. URL http://www.jstatsoft.org/v47/i12/}
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
clus <- tclust (x, k = 2, alpha = 0.1, restr.fact = 8)

plot (clus)
plot (clus, labels = "observation")
plot (clus, labels = "cluster")

# Three groups (one of them very scattered) and 0\% trimming level
clus <- tclust (x, k = 3, alpha=0.0, restr.fact = 100)

plot (clus)

%#--- EXAMPLE 2 ------------------------------------------
%data (geyser2)
%clus <- tkmeans (geyser2, k = 3, alpha = 0.03)
%plot (clus)

#--- EXAMPLE 3 ------------------------------------------
data (M5data)
x <- M5data[, 1:2]

clus.a <- tclust (x, k = 3, alpha = 0.1, restr.fact =  1,
                  restr = "eigen", equal.weights = TRUE, warnings = 1)
clus.b <- tclust (x, k = 3, alpha = 0.1, restr.fact =  1,
                   equal.weights = TRUE, warnings = 1)
clus.c <- tclust (x, k = 3, alpha = 0.1, restr.fact =  1,
                  restr = "deter", equal.weights = TRUE, iter.max = 100,
		  warnings = 1)
clus.d <- tclust (x, k = 3, alpha = 0.1, restr.fact = 50,
                  restr = "eigen", equal.weights = FALSE)

pa <- par (mfrow = c (2, 2))
plot (clus.a, main = "(a) tkmeans")
plot (clus.b, main = "(b) Gallegos and Ritter")
plot (clus.c, main = "(c) Gallegos")
plot (clus.d, main = "(d) tclust")
par (pa)

#--- EXAMPLE 4 ------------------------------------------
data (swissbank)
# Two clusters and 8\% trimming level
clus <- tclust (swissbank, k = 2, alpha = 0.08, restr.fact = 50)

                            # Pairs plot of the clustering solution
pairs (swissbank, col = clus$cluster + 1)
                                  # Two coordinates
plot (swissbank[, 4], swissbank[, 6], col = clus$cluster + 1,
     xlab = "Distance of the inner frame to lower border",
     ylab = "Length of the diagonal")
plot (clus)

# Three clusters and 0\% trimming level
clus <- tclust (swissbank, k = 3, alpha = 0.0, restr.fact = 110)

                            # Pairs plot of the clustering solution
pairs (swissbank, col = clus$cluster + 1)

                                   # Two coordinates
plot (swissbank[, 4], swissbank[, 6], col = clus$cluster + 1, 
      xlab = "Distance of the inner frame to lower border", 
      ylab = "Length of the diagonal")

plot (clus)

\dontshow{
%#--- EXAMPLE 5 ------------------------------------------
%# more complex covariance structure restrictions
%data (M5data)
%x <- M5data[, 1:2]
%
%clus.a <- tclust (x, k = 3, alpha = 0.1, restr.fact = 51, restr = "dir.eigen")
%clus.b <- tclust (x, k = 3, alpha = 0.1, restr.fact = 8, restr = "dir.deter")
%clus.c <- tclust (x, k = 3, alpha = 0.1, restr.fact = 8, restr = "prop")
%
%pa <- par (mfrow = c (1, 3))
%plot (clus.a, main = "(a) dir.eigen")
%plot (clus.b, main = "(b) dir.deter")
%plot (clus.c, main = "(c) prop")
%par (pa)
}
}

\keyword{ multivariate }
\keyword{ robust }
\keyword{ cluster }
