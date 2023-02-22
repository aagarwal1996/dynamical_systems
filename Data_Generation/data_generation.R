##
## Data Generation
##

library(deSolve)

#####################
# Smooth Noisy Data #
#####################

spline_smooth_samples <- function(data, impute_grad = F, nbasis = 48, max_t = 1,
								  lambda = 1e-10, norder = 4, penalty_order = 2,
								  smooth_type = "bspline"){
	if (smooth_type == "none"){
		return(data)
	}
	
	# initialize FDA objects
	penalty.Lfd = int2Lfd(penalty_order)
	xbasis = create.bspline.basis(rangeval=c(0,max_t),norder=norder,nbasis=nbasis)
	xbasis.fdPar = fdPar(xbasis,penalty.Lfd,lambda)
	ybasis = create.bspline.basis(rangeval=c(0,max_t),norder=norder,nbasis=nbasis)
	ybasis.fdPar = fdPar(ybasis,penalty.Lfd,lambda)

	lambdas = 10^seq(-6,4,by=0.25)    # lambdas to look over
	x.mean.gcv = rep(0,length(lambdas)) # store mean gcv
	y.mean.gcv = rep(0,length(lambdas))

	for(ilam in 1:length(lambdas)){
		# Set lambda
		xbasis.fdPar_i = xbasis.fdPar
		xbasis.fdPar_i$lambda = lambdas[ilam]
		ybasis.fdPar_i = ybasis.fdPar
		ybasis.fdPar_i$lambda = lambdas[ilam]

		# Smooth
		x_smooth_i <- smooth.basis(argvals = seq(0,max_t,len=length(data$x)), y=data$x, fdParobj=xbasis.fdPar_i)
		y_smooth_i <- smooth.basis(argvals = seq(0,max_t,len=length(data$y)), y=data$y, fdParobj=ybasis.fdPar_i)

		# Record average gcv
		x.mean.gcv[ilam] = mean(x_smooth_i$gcv)
		y.mean.gcv[ilam] = mean(y_smooth_i$gcv)
	}

	# pick the best smoothing parameters from CV
	best_x = which.min(x.mean.gcv)
	lambdabest_x = lambdas[best_x]
	best_y = which.min(y.mean.gcv)
	lambdabest_y = lambdas[best_y]
	xbasis.fdPar$lambda = lambdabest_x
	ybasis.fdPar$lambda = lambdabest_y
	
	x_smooth <- smooth.basis(argvals = seq(0,max_t,len=length(data$x)), y=data$x, fdParobj=xbasis.fdPar)
	y_smooth <- smooth.basis(argvals = seq(0,max_t,len=length(data$y)), y=data$y, fdParobj=ybasis.fdPar)
	smoothed_samples <- cbind(eval.fd(seq(0,max_t,len=length(data$x)),x_smooth$fd),
							  eval.fd(seq(0,max_t,len=length(data$y)),y_smooth$fd))
	colnames(smoothed_samples) <- c("x","y")
	smoothed_samples <- as.data.frame(smoothed_samples)

	if(impute_grad){
		x_basis_deriv <- deriv.fd(x_smooth$fd, 1)
		x_basis_deriv_eval <- eval.fd(evalarg = seq(0,max_t,len=length(data$x)), fdobj=x_basis_deriv)
		y_basis_deriv <- deriv.fd(y_smooth$fd, 1)
		y_basis_deriv_eval <- eval.fd(evalarg = seq(0,max_t,len=length(data$y)), fdobj=y_basis_deriv)
		smoothed_samples <- smoothed_samples %>% mutate("f_x" = x_basis_deriv_eval, "f_y" = y_basis_deriv_eval)
	} else{
		smoothed_samples$f_x <- data$f_x
		smoothed_samples$f_y <- data$f_y
	}

	return(smoothed_samples)
}

########################
# Abhi's Original Data #
########################

#' Get Abhi's VdP Data
#' 
#' This function reads the saved `.csv` files for the first example in the project
#'
#' @return sampled_data (data.frame) 1000 samples [x,y,f_x,f_y]
#' @export
get_abhi_data <- function(){
  
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

##############
# Van der Pol #
##############

generate_van_der_pol_samples <- function(system_params, num_samples = 1500, sample_density = 0.1){
  # function which returns samples from a Van der Pol system specified by \mu
  # initial condition is randomly selected
  # 
  ## Inputs:
  # system_params (numeric): mu parameter for Van der pol system
  # num_samples (integer)[optional]: number of samples to generate
  # sample_density (numeric)[optional]: `lsoda` Delta t between samples
  #
  ## Outputs:
  # sampled_data (data.frame): num_samples samples; columns in [x,y]
  
  # the only parameter for this model is mu
  if (length(system_params)>1){stop('Error: Too many Van Der Pol parameters')}

  # randomly pick the IC in [-10, 10] x [-10, 10]
  initial_condition <- runif(2, min = -10, max = 10)
  names(initial_condition) <- c('x','y')

  # compute the prediction and its gradient at each sample
  # unlist needed because of lsoda-friendly return format in eval_van_der_pol_gradient
  traj <- data.frame(deSolve::lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density), eval_van_der_pol_gradient, system_params))
  derivative_evals <- matrix(
    unlist(apply(traj[,2:3], 1, eval_van_der_pol_gradient, t = 0, system_params = system_params),use.names=F),
    ncol=(ncol(traj)-1),byrow=T)

  traj['f_x'] <- derivative_evals[,1] # TODO: is there a way to get the gradient with `lsoda` as opposed to manually?
  traj['f_y'] <- derivative_evals[,2]
  
  return(traj[,-1])
  
}

eval_van_der_pol_gradient <- function(t, v, system_params){
  # helper function to evaluate gradient of van der pol at point v = (x,y)
  
  mu <- system_params$mu
  x <- v[1]
  y <- v[2]
  
  dx <- mu*(x - x^3/3 - y)
  dy <- (1/mu)*x
  
  return(list(as.vector(c(dx,dy))))
}

##################
# Lotka-Volterra #
##################

generate_lv_samples <- function(system_params, num_samples = 1500, sample_density = 0.1){
  # randomly pick the IC
  initial_condition <- runif(2, min = 0, max = 10)
  names(initial_condition) <- c('x','y')
  
  # compute the prediction and its gradient at each sample
  traj <- data.frame(deSolve::lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density), 
                           eval_lv_gradient, system_params))
  derivative_evals <- matrix(
    unlist(apply(traj[,2:3], 1, eval_lv_gradient, t = 0, system_params = system_params),use.names=F),
    ncol=(ncol(traj)-1),byrow=T)
  
  traj['f_x'] <- derivative_evals[,1]
  traj['f_y'] <- derivative_evals[,2]
  
  return(traj[,-1])
  
}

eval_lv_gradient <- function(t, v, system_params){
  # helper function to evaluate gradient of van der pol at point v = (x,y)
  x <- v[1]
  y <- v[2]
  
  dx <- system_params$alpha*x - system_params$beta*x*y
  dy <- system_params$delta*x*y - system_params$gamma*y
  
  return(list(as.vector(c(dx,dy))))
}

##################################
# Log-Transformed Lotka-Volterra #  
##################################

generate_log_lv_samples <- function(system_params, num_samples = 1500, sample_density = 0.1){
  # randomly pick the IC
  initial_condition <- runif(2, min = 5, max = 10)
  names(initial_condition) <- c('x','y')
  
  # compute the prediction and its gradient at each sample
  traj <- data.frame(deSolve::lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density), 
                           eval_log_lv_gradient, system_params))
  derivative_evals <- matrix(
    unlist(apply(traj[,2:3], 1, eval_log_lv_gradient, t = 0, system_params = system_params),use.names=F),
    ncol=(ncol(traj)-1),byrow=T)
  
  traj['f_x'] <- derivative_evals[,1]
  traj['f_y'] <- derivative_evals[,2]
  
  return(traj[,-1])
  
}

eval_log_lv_gradient <- function(t, v, system_params){
  # helper function to evaluate gradient of van der pol at point v = (x,y)
  x <- v[1]
  y <- v[2]
  
  d_x <- system_params$alpha - system_params$beta*exp(y)
  d_y <- system_params$delta*exp(x) - system_params$gamma
  
  return(list(as.vector(c(d_x,d_y))))
}

########################
# Rosenzweig-MacArthur #  
########################

generate_rzma_samples <- function(system_params, num_samples = 1500, sample_density = 0.1){
  a <- system_params$a # saturation point
  b <- system_params$b #  half-saturation constant
  c <- system_params$c # predator growth rate
  r <- system_params$r # prey growth rate
  K <- system_params$K # carrying capacity
  
  # system only defined in Q1
  # randomly pick the IC in [0, 10] x [0, 10]
  initial_condition <- runif(2, min = 3, max = 7)
  names(initial_condition) <- c('x','y')
  
  # compute the prediction and its gradient at each sample
  traj <- data.frame(deSolve::lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density), eval_rzma_gradient, system_params))
  
  derivative_evals <- matrix(
    unlist(apply(traj[,2:3], 1, eval_rzma_gradient, t = 0, params = system_params),use.names=F),
    ncol=(ncol(traj)-1),byrow=T)
  #t(apply(traj[,2:3], 1, eval_rzma_gradient, t = 0, params = system_params))
  traj['f_x'] <- derivative_evals[,1] 
  traj['f_y'] <- derivative_evals[,2]
  
  return(traj[,-1])
  
}

eval_rzma_gradient <- function(t, v, params){
  # helper function to evaluate gradient of van der pol at point v = (x,y)
  x <- v[1]
  y <- v[2]
  
  f_x <- params$r * (1 - x/params$K)
  phi <- params$a / (params$b + x)
  
  dx <- x*f_x - y*phi*x
  dy <- y*(-params$c + params$K*x*phi)

  return(list(as.vector(c(dx,dy))))
}


###########################
# Generate Solution Paths #
###########################

# TODO: Avoid dulicated code

generate_true_path_vdp <- function(data, system_params, initial_condition, num_samples = 500, sample_density = 0.1){
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
  traj <- data.frame(deSolve::lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density), eval_van_der_pol_gradient, system_params))
  return(traj[,-1])
  
}

generate_true_path_rzma <- function(data, system_params, initial_condition, num_samples = 500, sample_density = 0.1){
  # function which generates a random trajectory along the true gradient field
  #
  ## Inputs:
  # data (matrix): matrix of limit cycle data to regress on; columns (x,y,f_x,f_y)
  # num_samples (integer)[optional]: number of samples to generate
  #
  ## Outputs:
  # sampled_path (data.frame): num_samples samples; columns in [x,y]
  initial_condition <- unlist(initial_condition)
  # compute the prediction and its gradient at each sample
  traj <- data.frame(deSolve::lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density), eval_rzma_gradient, system_params))
  return(traj[,-1])
  
}

generate_true_path_lv <- function(data, system_params, initial_condition, num_samples = 500, sample_density = 0.1){
  initial_condition <- unlist(initial_condition)
  traj <- data.frame(deSolve::lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density), eval_lv_gradient, system_params))
  return(traj[,-1])
  
}

generate_true_path_log_lv <- function(data, system_params, initial_condition, num_samples = 500, sample_density = 0.1){
  initial_condition <- unlist(initial_condition)
  traj <- data.frame(deSolve::lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density), eval_log_lv_gradient, system_params))
  return(traj[,-1])
  
}

########################
# Classical Estimators #
########################

## Nadaraya Watson

nw_path_helper <- function(t,y,params){
  # helper function to get nw gradient in deSolve::lsoda
  data <- params[[1]]
  bw_matrix <- params[[2]]
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
  traj <- data.frame(deSolve::lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density),
                           nw_path_helper,list(as.matrix(data),bw_matrix)))
  return(traj[,-1])
  
}

## LOESS

loess_path_helper <- function(t,y,params){
  # helper function to get nw gradient in lsoda
  data <- params[[1]]
  h <- params[[2]]
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
  traj <- data.frame(deSolve::lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density),
                           loess_path_helper,list(as.matrix(data),h)))
  return(traj[,-1])
  
}

## b-splines

generate_spline_path <- function(estimator, initial_condition, num_samples = 3000, sample_density = 0.1){
  initial_condition <- unlist(initial_condition)
  
  # compute the prediction and its gradient at each sample
  traj <- data.frame(deSolve::lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density),
                           bifd_spline_gradient,estimator$sfd_list))
  return(traj[,-1])
  
}

## k-NN

knn_path_helper <- function(t,y,params){
  # helper function to get nw gradient in lsoda
  data <- params[[1]]
  k <- params[[2]]
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
  
  traj <- data.frame(deSolve::lsoda(initial_condition,seq(0,num_samples*sample_density,by=sample_density),
                           knn_path_helper,list(as.matrix(data),k)))
  return(traj[,-1])
  
}

########################
# Simulation Functions #
########################

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
generate_true_grid_data <- function(model, system_params, eval_grid){
  if (model == "van_der_pol"){
    grid_gradient <- matrix(unlist(apply(eval_grid, 1, eval_van_der_pol_gradient, t = 0, system_params = system_params),
                                   use.names=F),ncol=ncol(eval_grid),byrow=T)
  } else if(model == "rzma"){
    grid_gradient <- matrix(unlist(apply(eval_grid, 1, eval_rzma_gradient, t = 0, system_params = system_params),
                                   use.names=F),ncol=ncol(eval_grid),byrow=T)
  }
  else if(model == "lotka_volterra"){
    grid_gradient <- matrix(unlist(apply(eval_grid, 1, eval_lv_gradient, t = 0, system_params = system_params),
                                   use.names=F),ncol=ncol(eval_grid),byrow=T)
  }  
  else if(model == "log_lotka_volterra"){
    grid_gradient <- matrix(unlist(apply(eval_grid, 1, eval_log_lv_gradient, t = 0, system_params = system_params),
                                   use.names=F),ncol=ncol(eval_grid),byrow=T)
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

get_gradient_field <- function(data, estimator){
  
  grid <- data$grid
  
  if (estimator$method == "spline"){
      spline_fit <- calculate_spline_gradient_field(data$limit_cycle_tail, norder = estimator$params$norder, nbasis = estimator$params$nbasis, 
                 grid$x_grid, grid$y_grid, lambda = estimator$params$lambda, side_info = estimator$params$side_info)
      sfd_list <- list(fdx = spline_fit$fdx, fdy = spline_fit$fdy)
      return_list <- list(field = spline_fit$field, sfd_list = sfd_list, params = estimator$params)
      return(return_list)
  }
  else if (estimator$method == "gd_spline"){
    spline_fit <- calculate_gd_spline_gradient_field(data$limit_cycle_tail, grid$x_grid, grid$y_grid, gd_params=estimator$params$gd_params,
                           side_info = estimator$params$side_info, norder = estimator$params$norder, nbasis = estimator$params$nbasis, lambda = estimator$params$lambda)
    sfd_list <- list(fdx = spline_fit$fdx, fdy = spline_fit$fdy)
    return_list <- list(field = spline_fit$field, sfd_list = sfd_list, params = estimator$params)
    return(return_list)
  }
  else if (estimator$method == "truth"){
    evaluated_field <- generate_true_grid_data(data$system_name, data$system_params, grid$eval_grid)
  }
  else if (estimator$method == "knn"){
    evaluated_field <- eval_knn_fit(grid$eval_grid, as.matrix(data$limit_cycle_tail), estimator$params$k)
  }
  else if (estimator$method == "nw"){
    nw_bw_matrix <- estimator$params$h*diag(2)
    evaluated_field <- NW_regression_cpp(grid$eval_grid, as.matrix(data$limit_cycle_tail), nw_bw_matrix)
  }
  else if (estimator$method == "loess"){
    if (estimator$params$h == "gcv"){
      bw_grid <- seq(0.05, 2, length.out = 10) # TODO: Allow user to specify
      estimator$params$h <- loess_bw_cv(grid$eval_grid, as.matrix(data$limit_cycle_tail), bw_grid)
    }
    evaluated_field <- eval_loess_fit(grid$eval_grid, as.matrix(data$limit_cycle_tail), estimator$params$h)
  }
  else {
    evaluated_field <- NA
  }

  return_list = list(evaluated_field = evaluated_field, params = estimator$params)
  return(return_list)
}

generate_solution_path <- function(data, estimator, initial_condition){
  evaluated_field <- NA
  if (estimator$method == "truth"){
    if (estimator$params$name == "van_der_pol"){
      evaluated_field <- generate_true_path_vdp(data, estimator$params, initial_condition) 
    } else if (estimator$params$name == "lotka_volterra") {
      evaluated_field <- generate_true_path_lv(data, estimator$params, initial_condition) 
    } else if (estimator$params$name == "log_lotka_volterra") {
      evaluated_field <- generate_true_path_log_lv(data, estimator$params, initial_condition) 
    } else if (estimator$params$name == "rzma") {
      evaluated_field <- generate_true_path_rzma(data, estimator$params, initial_condition)
    }
  } else if (estimator$method == "knn"){
    evaluated_field <- generate_knn_path(data, estimator, initial_condition)
  } else if (estimator$method == "nw"){
    evaluated_field <- generate_nw_path(data, estimator, initial_condition)
  } else if (estimator$method == "loess"){
    evaluated_field <- generate_loess_path(data, estimator, initial_condition)
  } else if ((estimator$method == "spline") | (estimator$method == "gd_spline")){
    evaluated_field <- generate_spline_path(estimator, initial_condition)
  }
  
  return(evaluated_field)
}

return_true_gradient_fn <- function(system_name){

	if (system_name == "van_der_pol"){
		gradient_fn <- eval_van_der_pol_gradient
	} else if (system_name == "lotka_volterra") {
		gradient_fn <- eval_lv_gradient
	} else if (system_name == "log_lotka_volterra") {
		gradient_fn <- eval_log_lv_gradient
	} else if (system_name == "rzma") {
		gradient_fn <- eval_rzma_gradient
	} else {
		stop('Error: True gradient not available for the ', estimator_object$method, ' system')
	}

	return(gradient_fn)
}