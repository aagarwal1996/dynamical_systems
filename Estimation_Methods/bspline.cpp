// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// This allows row-wise lambda functions in arma
// [[Rcpp::plugins("cpp11")]]

using namespace arma;

//' Solve regression problem to fit b-splines
//'
// [[Rcpp::export(fit_bsplines_cpp)]]
mat fit_bsplines_cpp(mat x_basis, mat y_basis, mat x_penalty, mat y_penalty,
                     vec x_obs, vec y_obs, double lambda){

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

// [[Rcpp::export(get_roughness_pen)]]
mat get_roughness_pen(mat x_basis, mat y_basis, mat x_penalty, mat y_penalty,
                     vec x_obs, vec y_obs, double lambda){
  
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
  
  mat roughness_x = coefs_x.t() * cross_penalty * coefs_x;
  mat roughness_y = coefs_y.t() * cross_penalty * coefs_y;
  mat roughness_total = roughness_x + roughness_y;
  return roughness_total;
}

// [[Rcpp::export(get_roughness_pen_redux)]]
mat get_roughness_pen_redux(mat coefs_x, mat coefs_y, mat cross_penalty){
  
  mat roughness_x = coefs_x.t() * cross_penalty * coefs_x;
  mat roughness_y = coefs_y.t() * cross_penalty * coefs_y;
  mat roughness_total = roughness_x + roughness_y;
  return roughness_total;
  
}

// [[Rcpp::export(get_cross_penalty)]]
mat get_cross_penalty(mat x_basis, mat y_basis, mat x_penalty, mat y_penalty){
  
  int num_basis_fns = x_basis.n_cols;
  
  mat penalty_id = eye( num_basis_fns, num_basis_fns );
  mat x_penalty_kroenecker = kron(x_penalty,penalty_id);
  mat y_penalty_kroenecker = kron(penalty_id, y_penalty);
  mat cross_penalty =  x_penalty_kroenecker + y_penalty_kroenecker;
  return cross_penalty;
  
}


// [[Rcpp::export(get_cross_basis)]]
mat get_cross_basis(mat x_basis, mat y_basis){
  
  int num_basis_fns = x_basis.n_cols;
  rowvec kroenecker_ones;
  kroenecker_ones.ones(num_basis_fns);
  
  mat x_basis_kroenecker = kron(x_basis,kroenecker_ones);
  mat y_basis_kroenecker = kron(kroenecker_ones, y_basis);
  mat cross_basis =  x_basis_kroenecker % y_basis_kroenecker;
  mat penalty_id = eye( num_basis_fns, num_basis_fns );
  return cross_basis;
  
}

// [[Rcpp::export(get_coefs)]]
mat get_coefs(mat cross_penalty, mat cross_basis, vec obs, double lambda){
  
  mat coefs = inv(cross_basis.t() * cross_basis + lambda*cross_penalty) * cross_basis.t() * obs;
  
  return coefs;
  
}




