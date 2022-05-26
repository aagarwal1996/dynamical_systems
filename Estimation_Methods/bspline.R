library(fda)
library(deSolve)
library(Rcpp)
library(ggplot2)
library(superheat)
sourceCpp(here::here('Estimation_Methods/bspline.cpp'))

calculate_spline_gradient_field <- function(data, x_grid, y_grid, norder = 4, nbasis = 12,
                                            penalty_order = 2, plot_penalty = FALSE){
  # evaluate b-spline basis functions at coordinates in data (N x 2 matrix)
  bspline_basis_fns <- generate_bspline_basis(data, x_grid, y_grid, norder = norder, 
                                              nbasis = nbasis, penalty_order = penalty_order)
  
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
  
  if (plot_penalty){
    plot_spline_penalty(x_grid,y_grid, spline.fd_x, spline.fd_y, penalty_order)
  }

  return()
  #return(spline_field)
}
  
plot_spline_penalty <- function(x_grid, y_grid, spline_bifd_x, spline_bifd_y, penalty_order = 2){
  # plots nth-derivative penalty at each point in eval grid
  # get penalty at each point in grid
  
  #bspline_penalty <- get_bspline_penalty_cpp(xbasis, ybasis, xpen, ypen, xdata, ydata)
  #print(bspline_penalty)
  # build plots
  
  # eval.basis(x_grid,spline_bifd_x$sbasis)
  sbasis_penalty_x <- getbasismatrix(x_grid,spline_bifd_x$sbasis, penalty_order)
  tbasis_penalty_x <- getbasismatrix(y_grid,spline_bifd_x$tbasis, penalty_order)
  grid_penalty_x <- abs(sbasis_penalty_x %*% spline_bifd_x$coef %*% t(tbasis_penalty_x))
  # prep for plotting
  grid_penalty_x <- t(grid_penalty_x)
  rownames(grid_penalty_x) <- round(y_grid,2)
  colnames(grid_penalty_x) <- round(x_grid,2)
  superheat(grid_penalty_x,
            bottom.label.text.angle = 270,
            title = "Second Derivative Penalty for x Over Grid")
  
  sbasis_penalty_y <- getbasismatrix(x_grid,spline_bifd_y$sbasis, penalty_order)
  tbasis_penalty_y <- getbasismatrix(y_grid,spline_bifd_y$tbasis, penalty_order)
  grid_penalty_y <- abs(sbasis_penalty_y %*% spline_bifd_y$coef %*% t(tbasis_penalty_y))
  grid_penalty_y <- t(grid_penalty_y)
  rownames(grid_penalty_y) <- round(y_grid,2)
  colnames(grid_penalty_y) <- round(x_grid,2)
  superheat(grid_penalty_y,
            bottom.label.text.angle = 270,
            title = "Second Derivative Penalty for y Over Grid")
  
  eval_grid <- unname(as.matrix(expand.grid(x_grid,y_grid)))
  long_penalty_x <- cbind(eval_grid, c(t(grid_penalty_x)))
  colnames(long_penalty_x) <- c("x", "y", "penalty")
  x_penalty_plot <- ggplot(data.frame(long_penalty_x), aes(x, y, fill=penalty)) + 
    geom_raster() +
    labs(title="Second Derivative Penalty for x Over Grid")
  print(x_penalty_plot)
  long_penalty_y <- cbind(eval_grid, c(t(grid_penalty_y)))
  colnames(long_penalty_y) <- c("x", "y", "penalty")
  y_penalty_plot <- ggplot(data.frame(long_penalty_y), aes(x, y, fill=penalty)) + 
    geom_raster() +
    labs(title = "Second Derivative Penalty for y Over Grid")
  print(y_penalty_plot)
  
  return()
}

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

bifd_spline_gradient <- function(t,x,p,sfd_x,sfd_y){
  # Helper function for `lsoda` which evaluates derivatives of a spline fit
  dx = eval.bifd(x[1],x[2],sfd_x)
  dy = eval.bifd(x[1],x[2],sfd_y)
  return(list(as.vector(c(dx,dy))))
}

### Testing

abhi_data <- generate_limit_cycle_data("abhi", c())
calculate_spline_gradient_field(abhi_data,
                                seq(-3,3,length.out=20), seq(-9,9,length.out=30),
                                plot_penalty = TRUE)
