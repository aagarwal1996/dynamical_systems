##
## Imports
##
library(fda)
library(ggplot2)
library(deSolve)
library(mvtnorm) # for Multivariate Normal density estimates in KDE
source(here::here('data_generation.R')) # functions to generate data from specific DS are in another file

##
## Data Generation
##

generate_data <- function(model, params, num_samples = 1000, 
                          add_obs_noise = F, add_grad_nose = F,
                          use_seed = F, save_csv = F){
  # This is a wrapper for data generation under a number of DS models
  # See `data_generation.R` for individual models
  #
  ## Inputs:
  # model (string): name of the model to sample from the limit cycle of
  #   implemented: van_der_pol, abhi
  #   to-implement: unstable_spiral, morris_lecar
  # params (vector): vector of parameters for the specified model
  # add_obs_noise (logical)[optional]: whether to add random noise to each observation (x,y); TODO: implement
  # add_grad_nose (logical)[optional]: whether to add random noise to each gradient (f_x, f_y); TODO: implement
  # use_seed (logical)[optional]: whether to set a seed before drawing data
  # save_csv (logical)[optional]: whether to save the CSV of the data; TODO: implement
  #
  ## Outputs:
  # sampled_data (data.frame): `num_samples` sampled points from the limit cycle; columns in [x,y,f_x,f_y]
  
  if (use_seed){set.seed(2022)}
  
  if (model == "van_der_pol"){
    sampled_data <- generate_van_der_pol(params)
  }
  else if(model == "unstable_spiral"){
    stop('Error: Unstable sprial not yet implemented.')
  }
  else if(model == "abhi"){
    sampled_data <- get_abhi_data()
  }
  else{
    stop('Error: ', model, ' model not implemented.')
  }
  
  if (save_csv){
    # TODO: Implement
  }
  
  return(sampled_data)
}

##
## Gradient Field Approximation
##

evaluate_gradient_methods <- function(data, tail_n = 500, 
                            x_grid_size = 24, y_grid_size = 24, extrapolation_size = 2){
  # Function which runs multiple methods to estimate the gradient field
  #   implemented: splines and kernel regression
  #   to-implement: local linear regression
  #
  ## Inputs:
  # data (data.frame): contains data with columns in [x, y, f_x, f_y]
  # tail_n (integer)[optional]: we only use the last `tail_n` samples of the data, ideally from the limit cycle
  # x_grid_size (integer)[optional]: size of x-axis grid to extrapolate over
  # y_grid_size (integer)[optional]: size of y-axis grid to extrapolate over
  # extrapolation_size (numeric)[optional]: grid will be from max and min of observations +/- this variable
  #
  # Outputs:
  # None (Displays various plots)

  limit_cycle.data <- tail(data, n = tail_n)
  
  # get extreme points of the data
  xmin <- floor(min(limit_cycle.data$x))
  xmax <- ceiling(max(limit_cycle.data$x))
  ymin <- floor(min(limit_cycle.data$y))
  ymax <- ceiling(max(limit_cycle.data$y))
  
  # grid to evaluate gradient extrapolation over
  x_grid <- seq(xmin - extrapolation_size , xmax + extrapolation_size, len = x_grid_size)
  y_grid <- seq(ymin - extrapolation_size, ymax + extrapolation_size, len = y_grid_size)
  
  
  spline_result <- spline_gradient(limit_cycle.data, x_grid, y_grid)
  NW_regression_result <- NW_regression(data,x_grid,y_grid) # takes some time to run
   
  plot_field(data, NW_regression_result, x_grid, y_grid, title="NW Kernel Regression")
  ## Uncomment line below for original plots from Abhi's .Rmd
  #generate_spline_plots(data, spline_result, x_grid, y_grid)
  
  # TODO: change output format of splines to use `plot_field function`
  plot(data$x,data$y,type='l',xlim=c(min(x_grid),max(x_grid)),ylim=c(min(y_grid),max(y_grid)),
       main = "Spline Estimate", xlab = "x", ylab = "y")
  eta <- 0.1
  for(i in seq(1,length(x_grid), by = 2)){
    for(j in seq(1,length(y_grid), by = 2)){
      arrows(x_grid[i],y_grid[j],x_grid[i]+eta*spline_result$x_grad_eval[i,j],y_grid[j]+ eta*spline_result$y_grad_eval[i,j])
    }
  }
  return(NW_regression_result)
}

##
## Spline fitting
##

spline_gradient <- function(data, x_grid, y_grid){
  # Estimates the gradient field after fitting b-splines to the limit cycle data
  #
  ## Inputs:
  # data (data.frame): contains data with columns in [x, y, f_x, f_y]
  # x_grid (numeric; vector): x-axis grid points to evaluate extrapolation over
  # y_grid (numeric; vector): y-axis grid points to evaluate extrapolation over
  #
  ## Outputs:
  # output_list (named list): see construction for contents; TODO: update documentaiton
  
  # create basis functions with 10 knots
  # TODO: adapt number of basis functions to the size of the limit cycle?
  xbasis = create.bspline.basis(range=c(min(x_grid),max(x_grid)),norder=4,nbasis=12)
  ybasis = create.bspline.basis(range=c(min(y_grid),max(y_grid)),norder=4,nbasis=12)
  xbasis.vals = eval.basis(data$x,xbasis)
  ybasis.vals = eval.basis(data$y,ybasis)
  
  # What we need is the evaluation of phi_j(x_i)*psi_k(x_j) for
  # each x and y. This creates 144 rows. I'll produce this using
  # Kronecker products
  Xmat = (xbasis.vals%x%matrix(1,1,12)) * (matrix(1,1,12)%x%ybasis.vals)
  #print(dim(Xmat))
  # Here the columns of xbasis.vals are repeated 12 times in sequence
  # while each column of ybasis.vals is repeated 12 times together.
  
  # Now we need a penalty matrix. We can get the penalty for one
  # basis from
  xPen = eval.penalty(xbasis, 2)
  yPen = eval.penalty(ybasis, 2)
  
  # to create the combined penalty we take the same Kronecker
  # product form
  
  allPen = xPen%x%diag(12) + diag(12)%x%yPen
  
  # (note that this penalizes the sum of squared second derivative,
  # without the cross term that would go into a thin plate spline
  # penalty)
  
  # And we can put it all together as
  
  lambda = 1e-8
  coefs_x =  solve( t(Xmat)%*%Xmat  + lambda*allPen, t(Xmat)%*%data$f_x)  
  coefs_y =  solve( t(Xmat)%*%Xmat  + lambda*allPen, t(Xmat)%*%data$f_y)
  
  # We'll reshape the coefficients into a matrix and put it in
  # a bivariate functional data object
  sfd_x = bifd(t(matrix(coefs_x,12,12)), xbasis, ybasis)
  sfd_y = bifd(t(matrix(coefs_y,12,12)), xbasis, ybasis)
  
  smat_x = eval.bifd(x_grid,y_grid,sfd_x)
  smat_y = eval.bifd(x_grid,y_grid,sfd_y)
  
  traj = lsoda(c(data$x[1],data$y[1]),seq(0,100,by=0.1),sVdP,0,sfd_x=sfd_x,sfd_y=sfd_y)
  
  output_list <- list(x_grad_bifd = sfd_x, x_grad_eval = smat_x, y_grad_bifd = sfd_y, y_grad_eval = smat_y, 
                      trajectory = traj, second.deriv_penalty = allPen)
  return(output_list)
}

sVdP = function(t,x,p,sfd_x,sfd_y){
  # Helper function for `lsoda` which evaluates derivatives of a spline fit
  dx = eval.bifd(x[1],x[2],sfd_x)
  dy = eval.bifd(x[1],x[2],sfd_y)
  return(list(as.vector(c(dx,dy))))
}

##
## Kernel regression
##

NW_regression <- function(data, x_grid, y_grid){
  nw_results_array <- array(-5, c(length(x_grid),length(y_grid),2))

  bw_matrix <- matrix(c(.1,0,0,0.1),2)
  for (i in 1:length(x_grid)){
    for (j in 1:length(y_grid)){
      nw_results_array[i,j,] <- eval_gaussian_NW(x_grid[i], y_grid[j], data, bw_matrix)
    }
  }
  return(nw_results_array)
}

eval_gaussian_NW <- function(x_val, y_val, data, bw_matrix){
  eval_point <- c(x_val, y_val)
  evaluation_results <- data
  evaluation_results$kernel_distance <- apply(evaluation_results[,1:2],1,  function(x) dmvnorm(x = x, mean = eval_point, sigma = bw_matrix))

  evaluation_results$numerator_x_terms <- evaluation_results$f_x * evaluation_results$kernel_distance
  evaluation_results$numerator_y_terms <- evaluation_results$f_y * evaluation_results$kernel_distance
  
  gradient_estimate_x <- sum(evaluation_results$numerator_x_terms)/sum(evaluation_results$kernel_distance)
  gradient_estimate_y <- sum(evaluation_results$numerator_y_terms)/sum(evaluation_results$kernel_distance)
  
  return(c(gradient_estimate_x, gradient_estimate_y))
  
}

##
## Local linear regression
##

# TODO: Implement!

##
## Plotting functions
##

# TODO: Add better documentation

plot_observations <- function(data){
  # helper function which generates a scatter plot of observations 
  ggplot(data,aes(x=x,y=y)) +
    geom_point() + 
    labs(x="x",y="y",title="Observations from a Dynamical System")
}

plot_field <- function(ls_samples, gradient_data, x_grid,y_grid,title_str=""){
  # plot a gradient field
  eta = 0.1
  plot(x_grid,y_grid,type='l',
       xlim=c(min(x_grid)-1,max(x_grid)+1),
       ylim=c(min(y_grid)-1,max(y_grid)+1), xlab = "", ylab = "")
  plot(ls_samples$x, ls_samples$y,
       xlim=c(min(x_grid)-1,max(x_grid)+1),
       ylim=c(min(y_grid)-1,max(y_grid)+1), 
       main = title_str, xlab = "x", ylab = "y")
  for(i in seq(1,length(x_grid), by = 2)){
    for(j in seq(1,length(y_grid), by = 2)){
      arrows(x_grid[i],y_grid[j],x_grid[i]+eta*gradient_data[i,j,1],y_grid[j]+ eta*gradient_data[i,j,2])
    }
  }
}

generate_spline_plots <- function(data, spline_output, x_grid, y_grid){
  # Generate the plots from Abhi's original code

  contour(x_grid,y_grid,spline_output$x_grad_eval)
  contour(x_grid,y_grid,spline_output$y_grad_eval)
  eta = 0.1
  plot(data$x,data$y,type='l',xlim=c(min(x_grid),max(x_grid)),ylim=c(min(y_grid),max(y_grid)))
  for(i in seq(1,length(x_grid),by=5)){
    for(j in seq(1,length(y_grid),by=5)){
      arrows(x_grid[i],y_grid[j],x_grid[i]+eta*spline_output$x_grad_eval[i,j],y_grid[j]+ eta*spline_output$y_grad_eval[i,j])
    }
  }
  
  #  Check how well we interpolate
  try1 = eval.bifd(data$x,data$y,spline_output$x_grad_bifd)
  try2 = eval.bifd(data$x,data$y,spline_output$y_grad_bifd)
  
  plot(data$f_x,diag(try1))
  plot(data$f_y,diag(try2))
  abline(c(0,1))
}

##
## Evaluation
##

abhi_data <- generate_data("abhi", c())
saved_abhi_NW <- evaluate_gradient_methods(abhi_data, extrapolation_size = 2)

