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


//' Evaluate NW kernel at coordinates `pred_coords`
//' 
// [[Rcpp::export(eval_gaussian_NW)]]
rowvec eval_gaussian_NW(rowvec pred_coords, mat data, mat bw_matrix){
  
  vec weights(data.n_rows);
  uvec data_coord_index = {0,1};
  uvec data_grad_index = {2,3};
  mat data_coord = data.cols(data_coord_index);
  mat data_grad = data.cols(data_grad_index); 
    
  vec gauss_kernel_weights = dmvnrm_arma(data_coord, pred_coords, bw_matrix);
  double total_weight = sum(gauss_kernel_weights);
  mat gauss_kernel_preds = data_grad.each_col() % gauss_kernel_weights;
  rowvec nw_estimate = sum(gauss_kernel_preds,0) / total_weight;
  
  return nw_estimate;
  
  // WT: not sure why `each_row()` did not work... eg:
  //data.each_row( [](rowvec& row){ getconstant(row); } ).print();
  
}

//' Estimate NW kernel for given coordinates conditional on data
//' 
// [[Rcpp::export(NW_regression_cpp)]]
mat NW_regression_cpp(mat eval_coords, mat data, mat bw_matrix){
  mat NW_estimate = mat(eval_coords.n_rows, eval_coords.n_cols);
  
  for (int i = 0; i < eval_coords.n_rows; i++){
    NW_estimate.row(i) = eval_gaussian_NW(eval_coords.row(i), data, bw_matrix);
  }
  
  return NW_estimate;
}
