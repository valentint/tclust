#include <RcppCommon.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// Input parameters of the procedure
struct params
{
  int n;              // Number of observations
  int p;              // Number of features
  double alpha;       // Trimming level
  int trimm;          // Number of observations to trim
  int no_trim;        // Number of observations to cluster
  int k;              // Number of clusters
  bool equal_weights; // Equal population proportion for all clusters
  double zero_tol;    // Tolerance that substitutes 0 (to avoid some numerical issues)
  int restrC;         // Restriction type (0 = eigenvalues or determinants)  
  bool deterC;        // Determinants restriction  
  double restr_fact;  // Level of constraint (of eigenvalues or determinants)
  double cshape;      // Level of constraint of the eigenvalues when deterC=true
  Rcpp::String opt;   // Model of estimate, can be "HARD" for hard assignment or "MIXTURE" for mixture assignment
};

// Cluster information
struct iteration
{
  arma::mat centers;    // Cluster centers
  arma::cube cov;       // Cluster covariance matrices
  arma::uvec cluster;   // Cluster assignation indices
  arma::vec disttom;    // Distances to cluster center
  double obj;           // Value of the objective function
  double NlogL;         // Value of the the negative of the CLASSIFICATION LOG-LIKELIHOOD  of the untrimmed units
                        // NlogL = -sum(max(ll(untrimmed units,[],2));
  arma::vec size;       // Cluster sizes
  arma::vec weights;    // Cluster weights
  int code;             // A return code signaling particular situations (e.g. data are aligned)
  arma::mat posterior;  // Cluster assignment given by 0/1 columns
};

namespace Rcpp {

// Support for wrap
template <> SEXP wrap(const iteration& iter);
template <> SEXP wrap(const params& pa);

}
