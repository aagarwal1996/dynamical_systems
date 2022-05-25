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
  mu <- params[1]
  
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

## Nadaraya Watson

nw_path_helper <- function(t,y,parms){
  # helper function to get nw gradient in lsoda
  data <- parms[[1]]
  bw_matrix <- parms[[2]]
  grad_pred <- NW_regression_cpp(unname(matrix(y,nrow = 1)), data, bw_matrix)
  grad_pred.lsoda <- list(c(x = grad_pred[1,1], y = grad_pred[1,2])) # reformat as list with named vector
  return(grad_pred.lsoda)
}

generate_nw_path <- function(data, bw, num_samples = 2000, sample_density = 0.1){
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
  # randomly pick the IC in [-10, 10] x [-10, 10]
  initial_condition <- runif(2, min = -6, max = 6)
  names(initial_condition) <- c('x','y')
  bw_matrix <- bw*diag(2)
  # compute the prediction and its gradient at each sample
  traj <- data.frame(lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density),nw_path_helper,list(data,bw_matrix)))
  return(traj[,-1])
  
}

## LOESS

loess_path_helper <- function(t,y,parms){
  # helper function to get nw gradient in lsoda
  data <- parms[[1]]
  bw_matrix <- parms[[2]]
  grad_pred <- get_loess_pred(unname(matrix(y,nrow = 1)), data, bw_matrix)
  grad_pred.lsoda <- list(c(x = grad_pred[1,1], y = grad_pred[1,2])) # reformat as list with named vector
  return(grad_pred.lsoda)
}

generate_loess_path <- function(data, bw, num_samples = 2000, sample_density = 0.1){
  # function which generates a random trajectory along the gradient field
  # of the LOESS Solution. The IC is randomized
  # 
  ## Inputs:
  # data (matrix): matrix of limit cycle data to regress on; columns (x,y,f_x,f_y)
  # bw (double): parameter which controls Kernel bandwidth (isotropic Normal variance)
  # num_samples (integer)[optional]: number of samples to generate
  #
  ## Outputs:
  # sampled_path (data.frame): num_samples samples; columns in [x,y]
  
  # the only parameter for this model is mu
  # randomly pick the IC in [-10, 10] x [-10, 10]
  initial_condition <- runif(2, min = -6, max = 6)
  names(initial_condition) <- c('x','y')
  bw_matrix <- bw*diag(2)
  # compute the prediction and its gradient at each sample
  traj <- data.frame(lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density),nw_path_helper,list(data,bw_matrix)))
  return(traj[,-1])
  
}

## b-splines

generate_spline_path <- function(data, bw, num_samples = 2000, sample_density = 0.1){
  # function which generates a random trajectory along the gradient field
  # of the LOESS Solution. The IC is randomized
  # 
  ## Inputs:
  # data (matrix): matrix of limit cycle data to regress on; columns (x,y,f_x,f_y)
  # bw (double): parameter which controls Kernel bandwidth (isotropic Normal variance)
  # num_samples (integer)[optional]: number of samples to generate
  #
  ## Outputs:
  # sampled_path (data.frame): num_samples samples; columns in [x,y]
  
  # the only parameter for this model is mu
  # randomly pick the IC in [-10, 10] x [-10, 10]
  initial_condition <- runif(2, min = -6, max = 6)
  names(initial_condition) <- c('x','y')
  bw_matrix <- bw*diag(2)
  # compute the prediction and its gradient at each sample
  traj <- data.frame(lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density),nw_path_helper,list(data,bw_matrix)))
  return(traj[,-1])
  
}