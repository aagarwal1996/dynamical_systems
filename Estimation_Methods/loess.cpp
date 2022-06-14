// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// This allows row-wise lambda functions in arma
// [[Rcpp::plugins("cpp11")]]

using namespace arma;

// From: https://gist.github.com/boennecd/db36747096e0d8a996a06e571f181ede
static double const log2pi = std::log(2.0 * M_PI);
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat const &x,
                      arma::rowvec const &mean,
                      arma::mat const &sigma,
                      bool logd = false) {
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())),
    constants = -xdim/2.0 * log2pi,
    other_terms = rootisum + constants;

  arma::rowvec z;
  for (int i=0; i < n; i++) {
    z      = (x.row(i) - mean) * rooti;
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  if (logd == false)
    out = exp(out);
  return(out);
}


//' Evaluates LOESS at a single point
//' 
// [[Rcpp::export(get_loess_pred)]]

rowvec get_loess_pred(rowvec query_point, mat data, double bw){
  
  vec weights(data.n_rows);
  uvec data_coord_index = {0,1};
  uvec data_grad_index = {2,3};
  mat data_coord = data.cols(data_coord_index);
  mat data_grad = data.cols(data_grad_index);
  
  mat bw_matrix = mat(2,2,fill::eye) * bw;
  vec gauss_kernel_weights = dmvnrm_arma(data_coord, query_point, bw_matrix);
  vec normalized_weights = gauss_kernel_weights/sum(gauss_kernel_weights);
  mat weight_matrix = mat(data.n_rows,data.n_rows,fill::eye).each_col() % normalized_weights;
  
  mat inv_kernel_mat = inv(data_coord.t() * weight_matrix * data_coord);
  rowvec loess_prediction = query_point * inv_kernel_mat * data_coord.t() * weight_matrix * data_grad;
  
  return(loess_prediction);
}

//' Estimates LOESS over an input grid
//' 
// [[Rcpp::export(eval_loess_fit)]]

mat eval_loess_fit(mat eval_coords, mat data, double bw){
  mat estimated_f = mat(eval_coords.n_rows, eval_coords.n_cols);
  for(int i = 0; i < eval_coords.n_rows; i++){
    estimated_f.row(i) = get_loess_pred(eval_coords.row(i), data, bw);
    
  }
  return(estimated_f);
}

// [[Rcpp::export(loess_bw_cv)]]
double loess_bw_cv(mat grid, mat data, vec bw_grid){
  double min_bw;
  min_bw = 0.5;
  return(min_bw);
}