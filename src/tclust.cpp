/**
 * Tclust - Clustering with trimming: C++ Functions
 *
 */

//  Rcpp::compileAttributes("C:/users/valen/onedrive/myrepo/rrdev/robClus")
//  Rcpp::compileAttributes("C:/projects/statproj/R/tclust")

// #include <RcppCommon.h>
// #include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include "tclust_types.h"

#include <Rcpp.h>

using namespace Rcpp;

namespace Rcpp {
  template <> SEXP wrap(const iteration& iter) {
    return Rcpp::List::create(
      _["obj"] = iter.obj,
      _["cluster"] = iter.cluster,
      _["disttom"] = iter.disttom,
      _["size"] = iter.size,
      _["weights"] = iter.weights,
      _["centers"] = iter.centers,
      _["cov"] = iter.cov,
      _["code"] = iter.code,
      _["posterior"] = iter.posterior);
  }
  
  template <> SEXP wrap(const params& pa) {
    return Rcpp::List::create(
      _["n"] = pa.n,
      _["p"] = pa.p,
      _["alpha"] = pa.alpha,
      _["trimm"] = pa.trimm,
      _["no_trim"] = pa.no_trim,
      _["k"] = pa.k,
      _["equal_weights"] = pa.equal_weights,
      _["zero_tol"] = pa.zero_tol);
      _["restrC"] = pa.restrC;
      _["deterC"] = pa.deterC;
      _["restr_fact"] = pa.restr_fact;
      _["cshape"] = pa.cshape;
      _["opt"] = pa.opt;
  }
}

/**
 * Convert an iteration object to Rcpp List
 *
 * @param iter: iteration object.
 *
 * @returns Rcpp List with same fields as the iteration object.
 */
Rcpp::List iter_to_list(iteration &iter)
{
  return Rcpp::List::create(
      _["centers"] = iter.centers,
      _["cov"] = iter.cov,
      _["cluster"] = iter.cluster,
      _["disttom"] = iter.disttom,
      _["obj"] = iter.obj,
      _["size"] = iter.size,
      _["weights"] = iter.weights,
      _["code"] = iter.code,
      _["posterior"] = iter.posterior);
}

/**
 * Apply restrictions to eigenvalues.
 *
 * @param autovalues: p x k matrix containin eigenvalues.
 * @param ni_ini: current sample size of the clusters.
 * @param factor_e: the level of the constraints.
 * @param zero_tol: tolerance level.
 *
 * @returns matrix with constrained eigenvalues.
 */
arma::mat restr2Eigenv(arma::mat autovalues, arma::vec ni_ini, double factor_e, double zero_tol)
{

  // Initializations
  double c = factor_e;
  arma::mat d = autovalues.t();
  int p = autovalues.n_rows;
  int k = autovalues.n_cols;
  int n = arma::accu(ni_ini);
  arma::mat nis(k, p); 
  nis.each_col() = ni_ini;

  //VT::20.01.2024: set any negative values in d to 0
  arma::uvec idx = find(d < 0);
  d.elem(idx).fill(0);
  
  // d_ is the ordered set of values in which the restriction objective function change the definition
  // points in d_ correspond to the frontiers for the intervals in which this objective function has the same definition
  // ed is a set with the middle points of these intervals
  arma::vec d_ = arma::sort(arma::vectorise(arma::join_cols(d, d / c)));
  int dim = d_.n_elem;

  arma::vec d_1 = arma::join_cols(d_, arma::vec({d_.back() * 2}));
  arma::vec d_2 = arma::join_cols(arma::vec({0}), d_);
  arma::vec ed = (d_1 + d_2) / 2;

  dim++;

  // the only relevant eigenvalues are those belonging to a clusters with sample size greater than 0.
  // eigenvalues corresponding to a clusters with 0 individuals has no influence in the objective function
  // if all the eigenvalues are 0 during the smart initialization we assign to all the eigenvalues the value 1
  if(arma::max(d.elem(arma::find(nis > 0))) <= zero_tol)
  {
    return arma::mat(p, k, arma::fill::zeros); // solution corresponds to 0 matrix
  }

  // we check if the eigenvalues verify the restrictions
  if(std::abs(arma::max(d.elem(arma::find(nis > 0))) / arma::min(d.elem(arma::find(nis > 0)))) <= c)
  {
    d.elem(arma::find(nis == 0)).fill(arma::mean(d.elem(arma::find(nis > 0))));
    return d.t(); // the solution corresponds to the input because it verifies the constraints
  }

  arma::mat t(k, dim);
  arma::mat s(k, dim);
  arma::mat r(k, dim);

  arma::vec sol(dim);
  arma::vec sal(dim);

  for (int mp_ = 0; mp_ < dim; mp_++)
  {
    for (int i = 0; i < k; i++)
    {
      r(i, mp_) = arma::accu(d.row(i) < ed(mp_)) + arma::accu(d.row(i) > ed(mp_) * c);
      s(i, mp_) = arma::accu(d.row(i) % (d.row(i) < ed(mp_)));
      t(i, mp_) = arma::accu(d.row(i) % (d.row(i) > ed(mp_) * c));
    }

    sol(mp_) = arma::accu(ni_ini / n % (s.col(mp_) + t.col(mp_) / c)) / (arma::accu(ni_ini / n % (r.col(mp_))));

    arma::mat sol_mp_matrix(k, p, arma::fill::value(sol(mp_)));
    arma::mat c_sol_mp_matrix(k, p, arma::fill::value(c * sol(mp_)));

    arma::mat e = sol(mp_) * arma::conv_to<arma::mat>::from(d < sol(mp_)) +
                  d % arma::conv_to<arma::mat>::from((d >= sol(mp_)) % (d <= c * sol(mp_))) +
                  (c * sol(mp_)) * arma::conv_to<arma::mat>::from(d > c * sol(mp_));

    arma::mat o = -1.0 / 2.0 * nis / n % (arma::log(e) + d / e);

    sal(mp_) = arma::accu(o);
  }

  // m is the optimum value for the eigenvalues procedure
  int eo = sal.index_max();
  double m = sol(eo);

  // based on the m value we get the restricted eigenvalues
  return (
             m * arma::conv_to<arma::mat>::from(d < m) +
             d % arma::conv_to<arma::mat>::from(d >= m) % arma::conv_to<arma::mat>::from(d <= c * m) +
             (c * m) * arma::conv_to<arma::mat>::from(d > c * m))
      .t();
}

arma::mat HandleSmallEv(arma::mat autovalues, double zero_tol)
{

    // Initializations
    int p = autovalues.n_rows;
    // int k = autovalues.n_cols;       //  not used here

    // Replace all ev less than zero_tol by zero_tol    
    autovalues.elem(find(autovalues <= zero_tol)).fill(zero_tol);   // "<= zero.tol" for checking for zero
    
    arma::mat mi = arma::min(autovalues);                           // the minimum eigenvalue of each cluster
    arma::mat ma = arma::max(autovalues);                           // the maximum eigenvalue of each cluster
    
    //  (mi/ma).print("mi/ma");
    //  (mi/zero_tol).print("mi/zero_tol");
    arma::uvec idx_iter = find(mi/ma <= zero_tol);
    
    //  idx_iter.print("idx_iter");
    //  autovalues.print("autovalues");

    // for each of these clusters. set all "normal" eigenvalues to a high value
    for(arma::uword ki = 0; ki < idx_iter.n_elem; ki++) { 
        arma::uvec idx = arma::find(autovalues.cols(idx_iter(ki), idx_iter(ki)) > mi(idx_iter(ki))/zero_tol);
        //  idx.print("To change index");
        for(arma::uword ji = 0; ji < idx.n_elem; ji++)
            autovalues(idx(ji), idx_iter(ki)) = mi(idx_iter(ki))/zero_tol;
    }
    
    //  autovalues.print("Modified autovalues");
   
	arma::mat det = prod(autovalues);								// the new determinants
    //  det.print("The new determinants");

    //	autovalues_det contains the autovalues corrected by the determinant
    //	the product of the eigenvalues of each column in autovalues_det is equal to 1
    
    //  arma::diagmat(det).print("The diagonal of det");
    //  arma::diagmat(arma::pow(det, -1.0/p)).print("det^");

	arma::mat autovalues_det = autovalues * arma::diagmat(arma::pow(det, -1.0/p));

    //  autovalues_det.print("Returning autovalues corrected by the determinant");
    
    return autovalues_det;  
}

/**
 * Apply restrictions to determinants.
 *
 * @param autovalues: p x k matrix containin eigenvalues.
 * @param ni_ini: current sample size of the clusters.
 * @param factor_e: the level of the constraints.
 * @param zero_tol: tolerance level.
 *
 * @returns matrix with constrained eigenvalues.
 */
arma::mat restr2Deter_old(arma::mat autovalues, arma::vec ni_ini, double factor_e, double zero_tol)
{

    // Initializations
    double c = factor_e;
    int p = autovalues.n_rows;
    int k = autovalues.n_cols;
    //int n = arma::accu(ni_ini);
    //arma::mat nis(k, p);
    //nis.each_col() = ni_ini;
    
//    Rcout << "Entering restr2Deter" << std::endl;

    if(p == 1)
        return restr2Eigenv(autovalues, ni_ini, factor_e, zero_tol);
    
    arma::mat es = prod(autovalues, 0);
    if(arma::max(es(find(ni_ini))) <= zero_tol) {
        arma::mat aret(p, k, arma::fill::zeros);
        // aret.print("!!!! - restr2Deter() returning zeros!");
        return aret;
    }
    
    arma::mat d = es.t();
	arma::mat autovalues_det = HandleSmallEv(autovalues, zero_tol);
	
    double ma = max(d(find(ni_ini > 0)));
    double mi = min(d(find(ni_ini > 0)));
 
//    Rcout << "Max: " << ma << std::endl;
//    Rcout << "Min: " << mi << std::endl;
  
    arma::mat dfin; 		
	if(ma/mi <= c) {
//        Rcout << "ma/mi <= c " << std::endl;
//        Rcout << "Fill with " << mean(d.elem(find(ni_ini  >= 0))) << std::endl;

		d.elem(find(ni_ini <= 0)).fill(mean(d.elem(find(ni_ini  >= 0))));
		dfin = pow(d, 1.0/p);
	} else {
//        arma::pow(d, 1.0/p).print("Passing to restr2Eigenv");
//        Rcout << "pow(c, 1.0/p)" << pow(c, 1.0/p) << std::endl;
		
        dfin = restr2Eigenv(pow(d, 1.0/p).t(), ni_ini, pow(c, 1.0/p), zero_tol);
    }
    
//    dfin.print("This is dfin");
//    diagmat(dfin).print("This is the diag matrix of dfin");
    
//    Rcout << "Exiting restr2Deter" << std::endl;
	return autovalues_det * diagmat(dfin);      							
}

/**
 * Apply restrictions to determinants.
 *
 * @param autovalues: p x k matrix containin eigenvalues.
 * @param ni_ini: current sample size of the clusters.
 * @param restr_factor: constraint level for the determinant constraints
 * @param cshape: constraint level for the eigenvalues constraints, default is 1e10
 * @param zero_tol: tolerance level.
 *
 * @returns matrix with constrained eigenvalues.
 */
arma::mat restr2Deter(arma::mat autovalues, arma::vec ni_ini, double restr_factor, double cshape, double zero_tol)
{

    // Initializations
    int p = autovalues.n_rows;
    int K = autovalues.n_cols;
    arma::vec ni_ini_1 = {1};
    
    //  Rcout << "Entering restr2Deter2" << std::endl;

    arma::uvec idx = find(autovalues < 1e-16);
    autovalues.elem(idx).fill(0);
    arma::mat autovalues_ = autovalues; 
    
    if(p == 1) return restr2Eigenv(autovalues, ni_ini, restr_factor, zero_tol);
    for(int k = 0; k < K; k++)  {
        autovalues_.col(k) = restr2Eigenv(autovalues.col(k), ni_ini_1, cshape, zero_tol);
    }

    arma::mat es = prod(autovalues_, 0);                        // apply(autovalues_, 2, prod)
   
    idx = find(es == 0);                                        // es[es==0] <- 1
    es.elem(idx).fill(1);                                       // ...
    arma::mat esmat(arma::pow(es, 1.0/p));
    esmat.set_size(1, K);
    arma::mat gm = autovalues_/repmat(esmat, p, 1);
    arma::mat d = sum(autovalues/gm, 0)/p;
    
    idx = find_nan(d);                                           // d[is.nan(d)] <- 0
    d.elem(idx).fill(0);                                         // ...
    
    arma::mat dfin = restr2Eigenv(d, ni_ini, pow(restr_factor, 1.0/p), zero_tol);
    
    // gm = gm * (gm > 0) + 1 * (gm == 0)
    idx = find(gm < 0);                                         // gm[gm<0] <- 0
    gm.elem(idx).fill(0);                                       // ...
    idx = find(gm == 0);                                        // gm[gm==0] <- 1
    gm.elem(idx).fill(1);                                       // ...

    dfin.set_size(1, K);
    dfin = repmat(dfin, p, 1);
    
	return dfin % gm;
}

/**
 * Manages constraints of a clustering stage.
 *
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void fRestr(iteration &iter, params &pa) {
    arma::cube u(pa.p, pa.p, pa.k); // Eigenvectors
    arma::mat d(pa.p, pa.k);        // Eigenvalues

    for(int ki = 0; ki < pa.k; ki++) {

        arma::vec eigval;
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, iter.cov.slice(ki));
    
        u.slice(ki) = eigvec;
        d.col(ki) = eigval;
    }

    d.elem(find(d < 0)).fill(0);  // all eigenvalue < 0 are assigned to 0, this issue appears for numerical errors
    
    if(!pa.deterC)
        d = restr2Eigenv(d, iter.size, pa.restr_fact, pa.zero_tol);
    else {
        // iter.size.print("Autovalues size");
        // d.print("Autovalues to be restricted");
        d = restr2Deter(d, iter.size, pa.restr_fact, pa.cshape, pa.zero_tol);
    }
    
    // checking for singularity in all clusters.
    int code = d.max() > pa.zero_tol;
    iter.code = code;
    
    if(code)  {
        // Reconstructing the cov matrices
        // d.print();
        for (int ki = 0; ki < pa.k; ki++) { 
            // Rcout << ki << std::endl;
            iter.cov.slice(ki) = u.slice(ki) * arma::diagmat(d.col(ki)) * u.slice(ki).t();
        }
    }
 
    // Rcout << "Exiting fRestr ()" << std::endl;
}

/**
 * Calculates the initial cluster assignment and initial values for the parameters.
 *
 * @param x: matrix with observations and features.
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void initClusters(arma::mat x, iteration &iter, params &pa)
{
  arma::mat iter_center = arma::mat(pa.p, pa.k);        // Column ki stores the centers of cluster p
  arma::cube iter_sigma = arma::cube(pa.p, pa.p, pa.k); // Covariance matrix of each cluster
  arma::vec size;                                       // Cluster sizes
  arma::vec weights;                                    // Cluster weigths

  arma::ivec idx = arma::randi(pa.k * (pa.p + 1), arma::distr_param(0, pa.n - 1));

  for (int ki = 0; ki < pa.k; ki++)
  {
    // Select the p+1 next elements of idx
    arma::uvec rows_to_select = arma::conv_to<arma::uvec>::from(idx.rows(ki * (pa.p + 1), ki * (pa.p + 1) + pa.p));

    // Set cluster centers to random observation
    iter_center.col(ki) = x.row(rows_to_select(0)).t();

    arma::mat X_ini = x.rows(rows_to_select);

    // Calculate cov matrix of cluster k
    iter_sigma.slice(ki) = pa.p / (float)(pa.p + 1) * arma::cov(X_ini);
  }

  iter.centers = iter_center.t();
  iter.cov = iter_sigma;

  if (pa.equal_weights)
  {
    size = arma::vec(pa.k, arma::fill::value(pa.no_trim / pa.k));
    weights = arma::vec(pa.k, arma::fill::value(1 / (double)pa.k));
  }
  else
  {
    weights = runif(pa.k);
    weights /= arma::accu(weights);
    size = arma::round(pa.n * weights);
  }

  iter.size = size;
  iter.weights = weights;
}

/**
 * FUNCTIONS TO CALCULATE THE DENSITY OF A MULTIVARIATE NORMAL
 *
 * Taken from: Nino Hardt, Dicko Ahmadou, Benjamin Christoffersen. Faster Multivariate Normal densities with RcppArmadillo and OpenMP
 * https://gallery.rcpp.org/articles/dmvnorm_arma/
 */

static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat)
{
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;) {
    double tmp(0.);
    for (unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

arma::vec dmvnrm_arma_fast(arma::mat const &x, arma::rowvec const &mean, arma::mat const &cov, bool const logd = false)
{
  using arma::uword;
  uword const n = x.n_rows, xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(cov)));
  double const rootisum = arma::sum(log(rooti.diag())),
               constants = -(double)xdim / 2.0 * log2pi,
               other_terms = rootisum + constants;

  arma::rowvec z;
  for(uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  if(logd)
    return out;
  
  return exp(out);
}

/**
 * Calculate the objective function value.
 *
 * @param x: matrix with observations and features.
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void calcObj(arma::mat x, iteration &iter, params &pa)
{
  int n = pa.n;
  int k = pa.k;
  Rcpp::String opt = pa.opt;

  arma::vec ww(n);
  arma::vec w;

  if(opt == "HARD") {
    for (int ki = 0; ki < k; ki++) {
      w = iter.weights(ki) * dmvnrm_arma_fast(x, iter.centers.row(ki), iter.cov.slice(ki)) % arma::conv_to<arma::mat>::from(iter.cluster == ki + 1);
      ww = w % (w >= 0) + ww; //	calculates each individual contribution for the obj funct hard
    }
  } else {
    for (int ki = 0; ki < k; ki++) {
      w = iter.weights(ki) * dmvnrm_arma_fast(x, iter.centers.row(ki), iter.cov.slice(ki));
      ww = w % (w >= 0) + ww; //	calculates each individual contribution for the obj funct mixture
    }
  }
  iter.obj = arma::accu(arma::log(ww.elem(arma::find(iter.cluster > 0))));
}

/**
 * Estimate the model parameters.
 *
 * @param x: matrix with observations and features.
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void estimClustPar(arma::mat x, iteration &iter, params &pa)
{
    for (int ki = 0; ki < pa.k; ki++) {
        if (iter.size(ki) > pa.zero_tol) {
            iter.centers.row(ki) = (iter.posterior.col(ki).t() * x) / iter.size(ki);
            // x centered
            arma::mat X_c = x;
            X_c.each_row() -= iter.centers.row(ki);
            X_c.each_col() %= iter.posterior.col(ki);
            iter.cov.slice(ki) = (X_c.t() * X_c) / iter.size(ki);
        } else {
            iter.centers.row(ki) = arma::mat(1, pa.p);
            iter.cov.slice(ki) = arma::mat(pa.p, pa.p, arma::fill::eye);
        }
    }
}

/**
 * Find cluster assignment and trimming.
 *
 * @param x: matrix with observations and features.
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void findClustAssig(arma::mat x, iteration &iter, params &pa)
{
    int k = pa.k;
    int n = pa.n;
    int no_trim = pa.no_trim;
    bool equal_weights = pa.equal_weights;
    Rcpp::String opt = pa.opt;
    
    // Rcout << "Entering findClustAssig()" << std::endl;
    
    arma::mat ll(n, k);
    for(int ki = 0; ki < k; ki++)   {
        // Rcout << "ki=" << ki << "=================" << std::endl;
        // iter.cov.slice(ki).print();
        
        ll.col(ki) = iter.weights(ki) * dmvnrm_arma_fast(x, iter.centers.row(ki), iter.cov.slice(ki));
    }

    arma::uvec old_assig = iter.cluster;
    arma::vec pre_z;
    arma::uvec tc_set;
    
    if (opt == "HARD") {
        pre_z = arma::max(ll, 1);
    } else {
        pre_z = arma::sum(ll, 1);
    }
    
    // Rcout << "Before trimming ..." << std::endl;

    // Determine elements to trim
    arma::uvec sorted_index = arma::sort_index(pre_z, "descending");
    arma::uvec last_indexes = arma::linspace<arma::uvec>(no_trim, n - 1, n - no_trim);
    arma::uvec obs_to_trim = sorted_index.elem(last_indexes);
    pre_z.elem(obs_to_trim) = arma::zeros<arma::vec>(n - no_trim);
    tc_set = pre_z > 0;

    // Cluster assignment with trimming
    // Rcout << "New cluster assignment ..." << std::endl;
    iter.cluster = (arma::index_max(ll, 1) + 1) % tc_set;

    // Find assignation matrix posterior
    if(opt == "MIXT") {
        iter.posterior = ll;
        iter.posterior.each_col() /= (pre_z + (pre_z == 0));
    } else {
    
        arma::uvec one_to_n = arma::linspace<arma::uvec>(0, n - 1, n);
        arma::uvec aux_assig = iter.cluster - 1 + (iter.cluster == 0);
        
        // 2xn matrix containing (observation, cluster) pairs in each column
        arma::umat subscripts = (arma::join_rows(one_to_n, aux_assig)).t();
        
        iter.posterior = arma::mat(n, k);
        iter.posterior.elem(arma::sub2ind(arma::size(iter.posterior), subscripts)) = arma::ones<arma::vec>(n);
    }

    iter.posterior.each_col() %= arma::conv_to<arma::vec>::from(tc_set); // Set to 0 all trimmed rows
    
    if(opt == "HARD" && arma::all(old_assig == iter.cluster)) {
        iter.code = 2;
    }

    // Obtain the clusters size
    iter.size = (arma::sum(iter.posterior, 0)).t();
    
    if (!equal_weights)
        iter.weights = iter.size/no_trim;

}

/**
 * Apply niter concentration steps to initial solution given by iter and pa.
 *
 * @param niter: number of concentration steps.
 * @param x: matrix with observations and features.
 * @param iter: a reference to the cluster information. Its values are modified.
 * @param pa: a reference to the procedure parameters.
 */
void concentration_steps(int niter, arma::mat x, iteration &iter, params &pa)
{

  for(int i1 = 0; i1 < niter; i1++) {
    fRestr(iter, pa); // restricting the clusters' scatter structure (Changes the iter object)

    // Rcout << "After frestr(): iter.code=" << iter.code << std::endl; 
    
    if(iter.code == 0)  {
      if(i1 > 0) {
        calcObj(x, iter, pa);
        Rcpp::warning("Data not in general position");
      } else {
        arma::cube cov = arma::cube(pa.p, pa.p, pa.k);
        cov.each_slice() = arma::eye(pa.p, pa.p);
        iter.cov = cov;
      }
    }

    // Estimate the cluster's assigment and TRIMMING (mixture models and HARD )
    // Rcout << "Before findClustAssig():" << std::endl; 
    findClustAssig(x, iter, pa); 
    
    // Rcout << "After findClustAssig(): iter.code=" << iter.code << ", i1=" << i1 << std::endl; 
    if((int)iter.code == 2 || (i1 == niter - 1))
      break;

    // Rcout << "Estimate Cluster Par:" << std::endl; 
    estimClustPar(x, iter, pa); // estimates the cluster's parameters
  }
  
  calcObj(x, iter, pa); // calculates the objetive function value
}

// Internal function for concentration steps (refinement) in tclust2
// @name tclust_c2
// @param x Rcpp::NumericMatrix, The input data.
// @param k The number of clusters initially searched for.
// @param cluster arma::uvec A numerical vector of size \code{n} containing the 
//      cluster assignment for each observation. Cluster names are integer numbers 
//      from 1 to k, 0 indicates trimmed observations. Note that it could be empty 
//      clusters with no observations when \code{equal_weights=FALSE}.
// @param alpha double, The proportion of observations to be trimmed.
// @param restrC Restriction type (0 for restriction on eigenvalues or determinant)
// @param deterC Restriction on determinants (true or false)
// @param restr_fact The constant \code{restr_fact >= 1} constrains the allowed 
//  differences among group scatters in terms of eigenvalues ratio. Larger values 
//  imply larger differences of group scatters, a value of 1 specifies the 
//  strongest restriction.
// @param cshape Shape constraint
// @param niter2 The maximum number of concentration steps to be performed for the 
//  \code{nkeep} solutions kept for further iteration. The concentration steps are 
//  stopped, whenever two consecutive steps lead to the same data partition.
// @param opt Define the target function to be optimized. A classification likelihood 
//  target function is considered if \code{opt="HARD"} and a mixture classification 
//  likelihood if \code{opt="MIXT"}.
// @param equal_weights A logical value, specifying whether equal cluster weights 
//  shall be considered in the concentration and assignment steps.
// @param zero_tol The zero tolerance used. By default set to 1e-16.
// @export
// [[Rcpp::export]]
iteration tclust_c2(arma::mat x, int k, arma::uvec cluster, double alpha = 0.05,
                     int restrC=0, bool deterC=false, double restr_fact = 12, double cshape=1e10,
                     int niter2 = 20, Rcpp::String opt = "HARD",
                     bool equal_weights = false, double zero_tol = 1e-16)
{

    int n = x.n_rows;
    int p = x.n_cols;
    int no_trim = std::floor(n * (1 - alpha));
    
    // Build the posterior matrix
    arma::uvec tc_set;
    arma::uvec one_to_n = arma::linspace<arma::uvec>(0, n - 1, n);
    arma::uvec aux_assig = cluster - 1 + (cluster == 0);
    tc_set = cluster > 0;         // not trimmed observations
    
    // 2xn matrix containing (observation, cluster) pairs in each column
    arma::umat subscripts = (arma::join_rows(one_to_n, aux_assig)).t();
    
    arma::mat posterior(n, k);
    posterior.elem(arma::sub2ind(arma::size(posterior), subscripts)) = arma::ones<arma::vec>(n);
    
    //VT::04.05.2024 - Fix the issue reported by Domenico:
    //    all trimmed observations were assigned to the first class (0)
    posterior.each_col() %= arma::conv_to<arma::vec>::from(tc_set); // Set to 0 all trimmed rows
    
    arma::vec size = (arma::sum(posterior, 0)).t();
    
    params pa;
    pa.n = n;
    pa.p = p;
    pa.alpha = alpha;
    pa.no_trim = no_trim;
    pa.trimm = n - no_trim;
    pa.k = k;
    pa.equal_weights = equal_weights;
    pa.zero_tol = zero_tol;
    pa.restrC = restrC;
    pa.deterC = deterC;
    pa.restr_fact = restr_fact;
    pa.cshape = cshape;
    pa.opt = opt;
    
    iteration iter;
    iter.centers = arma::mat(k, p);
    iter.cov = arma::cube(p, p, k);
    iter.cluster = cluster;
    iter.obj = 0.0,
    iter.size = size;
    iter.weights = size / no_trim;
    iter.code = 0;
    iter.posterior = posterior;
    
    estimClustPar(x, iter, pa);
    
    // VT::06.05.2024 - If niter2==0 (no iterations in the second step) tclust_c2 could
    //  break with "chol(): decomposition failed" - this happens if nkeep = nstart, i.e.
    // call tclust_c2 on all tried solutions (nstart) and at leats one of the solutions
    // has a class with one member which means that its covariance is the zero matrix.
    //  To solve this, we do an exctra call to frestr(), if niter2 == 0.
    if(niter2 == 0) {
        //  iter.cov.print("Before restriction");
        //  iter.size.print("Size");
        
        fRestr(iter, pa); 
    
        //  iter.cov.print("After restriction");
    }
    
    concentration_steps(niter2, x, iter, pa);
    
    return iter;
}

// Internal function for concentration steps (initializations) in tclust2
// @name tclust_c1
// @param x Rcpp::NumericMatrix, The input data.
// @param k The number of clusters initially searched for.
// @param alpha double, The proportion of observations to be trimmed.
// @param restrC Restriction type (0 for restriction on eigenvalues or determinant)
// @param deterC Restriction on determinants (true or false)
// @param restr_fact The constant \code{restr_fact >= 1} constrains the allowed 
//  differences among group scatters in terms of eigenvalues ratio. Larger values 
//  imply larger differences of group scatters, a value of 1 specifies the 
//  strongest restriction.
// @param cshape Shape constraint
// @param niter1 int, The number of concentration steps to be performed for the 
//     nstart initializations. 
// @param opt Define the target function to be optimized. A classification likelihood 
//  target function is considered if \code{opt="HARD"} and a mixture classification 
//  likelihood if \code{opt="MIXT"}.
// @param equal_weights A logical value, specifying whether equal cluster weights 
//  shall be considered in the concentration and assignment steps.
// @param zero_tol The zero tolerance used. By default set to 1e-16.
// @export
// [[Rcpp::export]]
Rcpp::List tclust_c1(arma::mat x, int k, double alpha = 0.05,
                     int restrC=0, bool deterC=false, double restr_fact = 12, double cshape=1e10, 
                     int niter1 = 3, Rcpp::String opt = "HARD",
                     bool equal_weights = false, double zero_tol = 1e-16)
{

  int n = x.n_rows;
  int p = x.n_cols;
  int no_trim = std::floor(n * (1 - alpha));

  params pa;
  pa.n = n;
  pa.p = p;
  pa.alpha = alpha;
  pa.no_trim = no_trim;
  pa.trimm = n - no_trim;
  pa.k = k;
  pa.equal_weights = equal_weights;
  pa.zero_tol = zero_tol;
  pa.restrC = restrC;
  pa.deterC = deterC;
  pa.restr_fact = restr_fact;
  pa.cshape = cshape;
  pa.opt = opt;

  iteration iter;
  iter.obj = 0.0,
  iter.cluster = arma::uvec(n);
  iter.size = arma::vec(k);
  iter.weights = arma::vec(k);
  iter.centers = arma::mat(k, p);
  iter.cov = arma::cube(p, p, k);
  iter.code = 0;
  iter.posterior = arma::mat(n, p);

  initClusters(x, iter, pa);                // Cluster random initialization
  concentration_steps(niter1, x, iter, pa); // Apply niter1 concentration steps

/*
  iteration iter1;
  Rcout << "Before tclust_c2():" << std::endl; 
  iter1 = tclust_c2(x, k, iter.cluster, alpha,
                     restrC, deterC, restr_fact, cshape,
                     0, opt, equal_weights, zero_tol);
  Rcout << "After tclust_c2():" << std::endl; 
  iter.centers.print("C1");
  iter1.centers.print("C2");
*/
  
  return Rcpp::List::create(
      _["obj"] = iter.obj,
      _["cluster"] = iter.cluster);
}


// [[Rcpp::export]]
arma::mat tclust_restr2_eigenv(arma::mat autovalues, arma::vec ni_ini, double factor_e=12, double zero_tol=1e-16) {
    return restr2Eigenv(autovalues, ni_ini, factor_e, zero_tol);
}

// [[Rcpp::export]]
arma::mat tclust_restr2_deter_old(arma::mat autovalues, arma::vec ni_ini, double factor_e=12, double zero_tol=1e-16) {
    return restr2Deter_old(autovalues, ni_ini, factor_e, zero_tol);
}

// [[Rcpp::export]]
arma::mat tclust_restr2_deter(arma::mat autovalues, arma::vec ni_ini, double restr_factor=12, double cshape=1e10, double zero_tol=1e-16) {
    return restr2Deter(autovalues, ni_ini, restr_factor, cshape, zero_tol);
}

// [[Rcpp::export]]
arma::mat tclust_HandleSmallEv(arma::mat autovalues, double zero_tol=1e-16) {
    return HandleSmallEv(autovalues, zero_tol);
}

// [[Rcpp::export]]
arma::vec dmvnrm(arma::mat x, arma::rowvec mean, arma::mat cov) {
    return dmvnrm_arma_fast(x, mean, cov, false);
}
