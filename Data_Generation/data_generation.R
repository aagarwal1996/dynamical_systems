##
## Data Generation
##

# TODO: Add observation noise

get_abhi_data <- function(){
  # Returns Abhi's simulated data
  # 
  ## Inputs:
  # None
  #
  ## Outputs:
  # sampled_data (data.frame): 1000 samples [x,y,f_x,f_y]
  
  # load data from csv into df
  x_df = read.csv("Saved_Data/Abhi/x.csv", header = FALSE, col.names = c("x"))
  y_df = read.csv("Saved_Data/Abhi/y.csv", header = FALSE, col.names = c("y"))
  fx_df = read.csv("Saved_Data/Abhi/f_x.csv", header = FALSE, col.names = c("f_x"))
  fy_df = read.csv("Saved_Data/Abhi/f_y.csv", header = FALSE, col.names = c("f_y"))
  
  # convert to consistent numeric encoding
  x =  sapply(x_df[, c(1)], as.numeric)
  y =  sapply(y_df[, c(1)], as.numeric)
  
  f_x =  sapply(fx_df[, c(1)], as.numeric)
  f_y =  sapply(fy_df[, c(1)], as.numeric)
  
  # return as data.frame
  sampled_data <- data.frame(x = x,y = y,f_x, f_y)
  return(sampled_data)
  
}

generate_van_der_pol <- function(params, num_samples = 1500, sample_density = 0.1){
  # function which returns samples from a Van der Pol system specified by \mu
  # initial condition is randomly selected
  # 
  ## Inputs:
  # params (numeric): mu parameter for Van der pol system
  # num_samples (integer)[optional]: number of samples to generate
  # sample_density (numeric)[optional]: `lsoda` Delta t between samples
  #
  ## Outputs:
  # sampled_data (data.frame): num_samples samples; columns in [x,y]
  
  # the only parameter for this model is mu
  if (length(params)>1){stop('Error: Too many Van Der Pol parameters')}
  mu <- params$mu
  
  # randomly pick the IC in [-10, 10] x [-10, 10]
  initial_condition <- runif(2, min = -10, max = 10)
  names(initial_condition) <- c('x','y')
  
  # compute the prediction and its gradient at each sample
  traj <- data.frame(lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density), eval_van_der_pol_gradient, mu))
  derivative_evals <- t(apply(traj[,2:3], 1, eval_van_der_pol_gradient, t = 0, mu = mu, list = T))
  traj['f_x'] <- derivative_evals[,1] # TODO: is there a way to get the gradient with `lsoda` as opposed to manually?
  traj['f_y'] <- derivative_evals[,2]
  
  return(traj[,-1])
  
}

eval_van_der_pol_gradient <- function(t, v, mu, list = F){
  # helper function to evaluate gradient of van der pol at point v = (x,y)
  x <- v[1]
  y <- v[2]
  
  dx <- mu*(x - x^3/3 - y)
  dy <- (1/mu)*x
  
  ## this form allows sampling when mu = 0 <-> uniform circular motion
  ## however, it doesn't work with `lsoda` for unknown reasons
  #dx <- y
  #dy <- mu*(1-x^2)*y - x
  
  # different return format if used in `lsoda` or not
  if(list){return(c(dx,dy))}
  return(list(as.vector(c(dx,dy))))
}

van_der_pol_gradient_helper <- function(v, mu){

  # redcued form not for lsoa
  x <- v[1]
  y <- v[2]
  
  dx <- mu*(x - x^3/3 - y)
  dy <- (1/mu)*x
  
  ## this form allows sampling when mu = 0 <-> uniform circular motion
  ## however, it doesn't work with `lsoda` for unknown reasons
  #dx <- y
  #dy <- mu*(1-x^2)*y - x
  
  # different return format if used in `lsoda` or not
  return(matrix(c(dx,dy), nrow = 1))
}

##
## Generate Solution Paths
##

## Truth

generate_true_path_vdp <- function(data, mu, initial_condition, num_samples = 500, sample_density = 0.1){
  # function which generates a random trajectory along the true gradient field
  # NOTE currently for VdP ONLY
  #
  ## Inputs:
  # data (matrix): matrix of limit cycle data to regress on; columns (x,y,f_x,f_y)
  # num_samples (integer)[optional]: number of samples to generate
  #
  ## Outputs:
  # sampled_path (data.frame): num_samples samples; columns in [x,y]
  
  initial_condition <- unlist(initial_condition)
  # compute the prediction and its gradient at each sample
  traj <- data.frame(lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density), eval_van_der_pol_gradient, mu))
  return(traj[,-1])
  
}

## Nadaraya Watson

nw_path_helper <- function(t,y,parms){
  # helper function to get nw gradient in lsoda
  data <- parms[[1]]
  bw_matrix <- parms[[2]]
  grad_pred <- NW_regression_cpp(unname(matrix(y,nrow = 1)), data, bw_matrix)
  grad_pred.lsoda <- list(c(x = grad_pred[1,1], y = grad_pred[1,2])) # reformat as list with named vector
  return(grad_pred.lsoda)
}

generate_nw_path <- function(data, estimator, initial_condition, num_samples = 500, sample_density = 0.1){
  # function which generates a random trajectory along the gradient field
  # of a NW kernel regression fit. The IC is randomized
  # 
  ## Inputs:
  # data (matrix): matrix of limit cycle data to regress on; columns (x,y,f_x,f_y)
  # bw (double): parameter which controls Kernel bandwidth (isotropic Normal variance)
  # num_samples (integer)[optional]: number of samples to generate
  #
  ## Outputs:
  # sampled_path (data.frame): num_samples samples; columns in [x,y]
  
  # the only parameter for this model is mu
  # randomly pick the IC
  h <- estimator$params$h
  initial_condition <- unlist(initial_condition)
  bw_matrix <- h*diag(2)
  # compute the prediction and its gradient at each sample
  traj <- data.frame(lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density),
                           nw_path_helper,list(as.matrix(data),bw_matrix)))
  return(traj[,-1])
  
}

## LOESS

loess_path_helper <- function(t,y,parms){
  # helper function to get nw gradient in lsoda
  data <- parms[[1]]
  h <- parms[[2]]
  grad_pred <- get_loess_pred(unname(matrix(y,nrow = 1)), data, h)
  grad_pred.lsoda <- list(c(x = grad_pred[1,1], y = grad_pred[1,2])) # reformat as list with named vector
  return(grad_pred.lsoda)
}

generate_loess_path <- function(data, estimator, initial_condition, num_samples = 500, sample_density = 0.1){
  # function which generates a random trajectory along the gradient field
  # of the LOESS Solution. The IC is randomized

  h <- estimator$params$h
  initial_condition <- unlist(initial_condition)
  # compute the prediction and its gradient at each sample
  traj <- data.frame(lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density),
                           loess_path_helper,list(as.matrix(data),h)))
  return(traj[,-1])
  
}

## b-splines

generate_spline_path <- function(estimator, initial_condition, num_samples = 500, sample_density = 0.1){
  initial_condition <- unlist(initial_condition)
  
  # compute the prediction and its gradient at each sample
  traj <- data.frame(lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density),
                           bifd_spline_gradient,estimator$sfd_list))
  return(traj[,-1])
  
}

## k-NN

knn_path_helper <- function(t,y,parms){
  # helper function to get nw gradient in lsoda
  data <- parms[[1]]
  k <- parms[[2]]
  grad_pred <- eval_knn_fit(unname(matrix(y,nrow = 1)), data, k)
  grad_pred.lsoda <- list(c(x = grad_pred[1,1], y = grad_pred[1,2])) # reformat as list with named vector
  return(grad_pred.lsoda)
}

generate_knn_path <- function(data, estimator, initial_condition, num_samples = 500, sample_density = 0.1){
  # function which generates a random trajectory along the gradient field
  # of the k-NN Solution. The IC is randomized
  # 
  ## Inputs:
  # data (matrix): matrix of limit cycle data to regress on; columns (x,y,f_x,f_y)
  # k (integer): parameter which controls number of neighbors
  # num_samples (integer)[optional]: number of samples to generate
  #
  ## Outputs:
  # sampled_path (data.frame): num_samples samples; columns in [x,y]
  
  # the only parameter for this model is mu
  # randomly pick the 
  # compute the prediction and its gradient at each sample
  k <- estimator$params$k
  initial_condition <- unlist(initial_condition)
  
  traj <- data.frame(lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density),
                           knn_path_helper,list(as.matrix(data),k)))
  return(traj[,-1])
  
}


#' This function evaluates the gradient field over a user-specified grid
#'
#' @param model (string): name of the model to generate the gradient field of; see `generate_limit_cycle_data()`
#' @param params (vector): vector of parameters for the specified model; see `generate_limit_cycle_data()`
#' @param eval_grid (matrix): n x 2 matrix of grid points to calculate gradient at
#'
#' @return grid_gradient (matrix): n x 2 matrix of gradient evaluated at `eval_grid`
#'
#' @examples
#' x_grid <- seq(-1, 1, len = x_grid_size)
#' y_grid <- seq(-1, 1, len = y_grid_size)
#' eval_grid <- unname(as.matrix(expand.grid(x_grid,y_grid))) # convert to N x 2 matrix
#' true_field <- generate_grid_data("van_der_pol", c(2), eval_grid)
generate_true_grid_data <- function(model, params, eval_grid){
  if (model == "van_der_pol"){
    grid_gradient <- t(apply(eval_grid, 1, van_der_pol_gradient_helper, mu = params$mu))
  }
  else if(model == "asymmetric_circle"){
    stop('Error: Asymmetric circle gradient grid not yet implemented.')
  }
  else if(model == "abhi"){
    # true field is unknown
    grid_gradient <- matrix(0, nrow = nrow(eval_grid), ncol = ncol(eval_grid))
  }
  else{
    stop('Error: ', model, ' model not implemented.')
  }
  
  return(grid_gradient)
}

###################
### Generalized ###
###################

get_gradient_field <- function(data, estimator){
  
  grid <- data$grid
  
  if (estimator$method == "spline"){
      spline_fit <- calculate_spline_gradient_field(data$limit_cycle_tail,
                 grid$x_grid, grid$y_grid, lambda = estimator$params$lambda)
      sfd_list <- list(fdx = spline_fit$fdx, fdy = spline_fit$fdy)
      return_list <- list(field = spline_fit$field, sfd_list = sfd_list)
      return(return_list)
  }
  else if (estimator$method == "truth"){
    evaluated_field <- generate_true_grid_data(data$system, data$params, grid$eval_grid)
  }
  else if (estimator$method == "knn"){
    evaluated_field <- eval_knn_fit(grid$eval_grid, as.matrix(data$limit_cycle_tail), estimator$params$k)
  }
  else if (estimator$method == "nw"){
    nw_bw_matrix <- estimator$params$h*diag(2)
    evaluated_field <- NW_regression_cpp(grid$eval_grid, as.matrix(data$limit_cycle_tail), nw_bw_matrix)
  }
  else if (estimator$method == "loess"){
    evaluated_field <- eval_loess_fit(grid$eval_grid, as.matrix(data$limit_cycle_tail), estimator$params$h)
  }
  else {
    evaluated_field <- NA
  }

  return(evaluated_field)
}

generate_solution_path <- function(data, estimator, initial_condition){
  evaluated_field <- NA
  
  if (estimator$method == "truth"){
    if (estimator$params$name == "van_der_pol"){
      evaluated_field <- generate_true_path_vdp(data, estimator$params$mu, initial_condition) 
      }
    }
  else if (estimator$method == "knn"){
    evaluated_field <- generate_knn_path(data, estimator, initial_condition)
  }
  else if (estimator$method == "nw"){
    evaluated_field <- generate_nw_path(data, estimator, initial_condition)
  }
  else if (estimator$method == "loess"){
    evaluated_field <- generate_loess_path(data, estimator, initial_condition)
  }
  else if (estimator$method == "spline"){
    evaluated_field <- generate_spline_path(estimator, initial_condition)
  }

  return(evaluated_field)
}
