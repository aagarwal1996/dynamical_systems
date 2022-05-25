library(fda)
library(deSolve)
library(Rcpp)
sourceCpp(here::here('Estimation_Methods/bspline.cpp'))

generate_bspline_basis <- function(data, x_grid, y_grid, norder = 4, nbasis = 12,
                                   penalty_order = 2){
  # number of breaks equals `nbasis` - `norder` + 2; 10 for default parameters (8 interior knots)
  xbasis = create.bspline.basis(rangeval=c(min(x_grid),max(x_grid)),norder=norder,nbasis=nbasis)
  ybasis = create.bspline.basis(rangeval=c(min(y_grid),max(y_grid)),norder=norder,nbasis=nbasis)
  xbasis.vals = eval.basis(data$x,xbasis)
  ybasis.vals = eval.basis(data$y,ybasis)
  
  # get roughness penalty for each axis
  xPen = eval.penalty(xbasis, penalty_order)
  yPen = eval.penalty(ybasis, penalty_order)
  
  spline_fit_list <- list(xbasis = xbasis, ybasis = ybasis, 
                          xbasis.eval = xbasis.vals, ybasis.eval = ybasis.vals,
                          xpenalty = xPen, ypenalty = yPen)
  return(spline_fit_list)
}
  
calculate_spline_gradient_field <- function(data, x_grid, y_grid, norder = 4, nbasis = 12){
  # evaluate b-spline basis functions at coordinates in data (N x 2 matrix)
  bspline_basis_fns <- generate_bspline_basis(data, x_grid, y_grid, norder = 4, nbasis = 12)
  
  bspline_fit_coeffs <- fit_bsplines_cpp(bspline_basis_fns$xbasis.eval, bspline_basis_fns$ybasis.eval,
                                  bspline_basis_fns$xpenalty, bspline_basis_fns$ypenalty,
                                  data$f_x, data$f_y)
  
  # create bivariate functional data objects for our fit splines 
  spline.fd_x <- bifd(t(matrix(bspline_fit_coeffs[,1],12,12)), 
                      bspline_basis_fns$xbasis,  bspline_basis_fns$ybasis)
  spline.fd_y <- bifd(t(matrix(bspline_fit_coeffs[,2],12,12)),
                      bspline_basis_fns$ xbasis,  bspline_basis_fns$ybasis)
  
  # evaluate over user-specified grid
  spline_grid_x <- eval.bifd(x_grid,y_grid,spline.fd_x)
  spline_grid_y <- eval.bifd(x_grid,y_grid,spline.fd_y)
  
  # return as |Grid| x 2 matrix
  spline_field <- matrix(c(c(spline_grid_x),c(spline_grid_y)),ncol=2)

  return(spline_field)
}
  

bifd_spline_gradient <- function(t,x,p,sfd_x,sfd_y){
  # Helper function for `lsoda` which evaluates derivatives of a spline fit
  dx = eval.bifd(x[1],x[2],sfd_x)
  dy = eval.bifd(x[1],x[2],sfd_y)
  return(list(as.vector(c(dx,dy))))
}