// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// This allows row-wise lambda functions in arma
// [[Rcpp::plugins("cpp11")]]

using namespace arma;

//' Solve regression problem to fit b-splines
//'
// [[Rcpp::export(fit_bsplines_cpp)]]
mat fit_bsplines_cpp(mat x_basis, mat y_basis, mat x_penalty, mat y_penalty,
                     vec x_obs, vec y_obs, double lambda = 1e-8){

  int num_basis_fns = x_basis.n_cols;
  rowvec kroenecker_ones;
  kroenecker_ones.ones(num_basis_fns);
  
  mat x_basis_kroenecker = kron(x_basis,kroenecker_ones);
  mat y_basis_kroenecker = kron(kroenecker_ones, y_basis);
  mat cross_basis =  x_basis_kroenecker % y_basis_kroenecker;
  
  mat penalty_id = eye( num_basis_fns, num_basis_fns );
  mat x_penalty_kroenecker = kron(x_penalty,penalty_id);
  mat y_penalty_kroenecker = kron(penalty_id, y_penalty);
  mat cross_penalty =  x_penalty_kroenecker + y_penalty_kroenecker;
  
  mat coefs_x = inv(cross_basis.t() * cross_basis + lambda*cross_penalty) *
    cross_basis.t() * x_obs;
  mat coefs_y = inv(cross_basis.t() * cross_basis + lambda*cross_penalty) *
    cross_basis.t() * y_obs;
  
  mat coeff_mat = join_rows(coefs_x, coefs_y);

  return coeff_mat;
}






