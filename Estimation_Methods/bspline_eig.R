# IMPORTANT: You must source `run_simulation_redux.R` first
source(here::here('Estimation_Methods','bspline.R'))
sourceCpp(here::here('Estimation_Methods/bspline.cpp'))

get_eigenconstrained_matrix <- function(dim, eigs = c(-10,-12,-15,-16, -18, -21)){
  stopifnot(dim == length(eigs))
  random_mat <- matrix(runif(dim^2,-10,10),nrow=dim)
  random_mat_sym <- random_mat + t(random_mat)
  eigenvectors <- eigen(random_mat_sym)$vectors
  eigenvalues <- diag(eigs)
  eig_sample <- eigenvectors %*% eigenvalues %*% t(eigenvectors)
  return(eig_sample)
}

calculate_eigpen_spline_gradient_field <- function(data, x_grid, y_grid, init_coefs_x = NA, init_coefs_y = NA, norder = 6, nbasis = 12, penalty_order = 2, plot_penalty = FALSE){
  
  # evaluate b-spline basis functions at coordinates in data (N x 2 matrix)
  bspline_basis_fns <- generate_bspline_basis(data, x_grid, y_grid, norder = norder, 
                                              nbasis = nbasis, penalty_order = penalty_order)
  # optimizes over log-values to enforce positivity
  cross_penalty_matrix <- get_cross_penalty(bspline_basis_fns$xbasis.eval, bspline_basis_fns$ybasis.eval,
                                           bspline_basis_fns$xpenalty, bspline_basis_fns$ypenalty)
  cross_basis <- get_cross_basis(bspline_basis_fns$xbasis.eval, bspline_basis_fns$ybasis.eval)
  
  #init_coefs <- fit_bsplines_cpp(bspline_basis_fns$xbasis.eval, bspline_basis_fns$ybasis.eval,
  #                               bspline_basis_fns$xpenalty, bspline_basis_fns$ypenalty,
  #                               data$f_x, data$f_y, 1e-2)
  
  # init_coefs_x.eig <-eigen(matrix(init_coefs[,1], nrow = nbasis))
  # init_coefs_x <- Re(init_coefs_x.eig$vectors %*% (diag(init_coefs_x.eig$values) * -1 * sign(Re(diag(init_coefs_x.eig$values)))) %*% solve(init_coefs_x.eig$vectors))
  # init_coefs_y.eig <-eigen(matrix(init_coefs[,2], nrow = nbasis))
  # init_coefs_y <- Re(init_coefs_y.eig$vectors %*% (diag(init_coefs_y.eig$values) * -1 * sign(Re(diag(init_coefs_y.eig$values)))) %*% solve(init_coefs_y.eig$vectors))
  
  #init_coefs_x <- get_eigenconstrained_matrix(nbasis) #init_coefs[,1]
  #init_coefs_x <- init_coefs[,1]
  #init_coefs_y <- get_eigenconstrained_matrix(nbasis) #init_coefs[,2]#get_eigenconstrained_matrix(nbasis) # # runif(nbasis^2,-5,5) #init_coefs[,2]
  #init_coefs_y <- init_coefs[,2]
  
  #init_coefs_x = runif(144,-5,5)
  #init_coefs_y = runif(144,-5,5)

  optim_fit <- optim(c(0.1,0.1, init_coefs_x,init_coefs_y), optim_penalized_fit, data = data, bspline_basis_fns = bspline_basis_fns, 
                     cross_penalty_matrix = cross_penalty_matrix, method = "BFGS", control = list(trace = 3, maxit = 1000))
  #View(optim_fit$par)
  fit_result <- get_penalized_fit(optim_fit$par, data, bspline_basis_fns, x_grid, y_grid, nbasis)
  #fit_result <- get_penalized_fit(c(0.1,0.1, init_coefs_x,init_coefs_y), data, bspline_basis_fns, x_grid, y_grid, nbasis)
  return(fit_result)
}

optim_penalized_fit <- function(params, data, bspline_basis_fns, cross_penalty_matrix){
  lambda <- 1e0 #exp(params[1])
  eig_pen <- 1e4 #exp(params[2])
  
  nbasis = sqrt((length(params)-2)/2)
  x_coeff <- params[3:(2+nbasis^2)]
  y_coeff <- params[(3+nbasis^2):(2+2*nbasis^2)]

  spline.fd_x <- bifd(t(matrix(x_coeff,nbasis,nbasis)), 
                      bspline_basis_fns$xbasis,  bspline_basis_fns$ybasis)
  spline.fd_y <- bifd(t(matrix(y_coeff,nbasis,nbasis)),
                      bspline_basis_fns$xbasis,  bspline_basis_fns$ybasis)

  num_artifical <- 200
  set.seed(2)
  x_samples <- runif(num_artifical, bspline_basis_fns$xbasis$rangeval[1], bspline_basis_fns$xbasis$rangeval[2])
  y_samples <- runif(num_artifical, bspline_basis_fns$ybasis$rangeval[1], bspline_basis_fns$ybasis$rangeval[2])
  artificial_samples <- matrix(c(x_samples,y_samples,rep(NA,num_artifical),rep(NA,num_artifical)),ncol=4)
  colnames(artificial_samples) <- c("x","y","f_x","f_y")

  sse <- sum((diag(eval.bifd(data[,1],data[,2],spline.fd_x))-data[,3])^2) + sum(diag((eval.bifd(data[,1],data[,2],spline.fd_y))-data[,4])^2)
  
  extended_data <- rbind(data,artificial_samples)
  sum_max_eig <- max(apply(extended_data, 1, get_max_eig, fd_x = spline.fd_x, fd_y = spline.fd_y))
  
  roughness_pen <- get_roughness_pen_redux(matrix(spline.fd_x$coefs,ncol=1),matrix(spline.fd_y$coefs,ncol=1),cross_penalty_matrix)[1,1]
  loss <- (1+(sum_max_eig > 0))*(sse + lambda*roughness_pen) + (eig_pen*sum_max_eig)*(sum_max_eig > 0) # indicator penalty
  cat("Loss: ",loss," Max Eig Pen: ", sum_max_eig," Roughness:",round(roughness_pen,3)," SSE:",round(sse,3),"\n")
  return(loss)
  
}

get_penalized_fit <- function(params, data, bspline_basis_fns, x_grid, y_grid, nbasis){
  
  nbasis = sqrt((length(params)-2)/2)
  x_coeff <- params[3:(2+nbasis^2)]
  y_coeff <- params[(3+nbasis^2):(2+2*nbasis^2)]

  spline.fd_x <- bifd(t(matrix(x_coeff,nbasis,nbasis)), 
                      bspline_basis_fns$xbasis,  bspline_basis_fns$ybasis)
  spline.fd_y <- bifd(t(matrix(y_coeff,nbasis,nbasis)),
                      bspline_basis_fns$xbasis,  bspline_basis_fns$ybasis)

  spline_grid_x <- eval.bifd(x_grid,y_grid,spline.fd_x)
  spline_grid_y <- eval.bifd(x_grid,y_grid,spline.fd_y)
  
  # return as |Grid| x 2 matrix
  spline_field <- matrix(c(c(spline_grid_x),c(spline_grid_y)),ncol=2)
  
  #if (plot_penalty){
  #  plot_spline_penalty(x_grid,y_grid, spline.fd_x, spline.fd_y, penalty_order)
  #}
  
  return(list(field = spline_field, fdx = spline.fd_x, fdy = spline.fd_y))
}



# https://stackoverflow.com/questions/64749481/how-to-solve-a-quadratic-equation-in-r
quadratic_real <- function(a, b, c)
{
  a <- as.complex(a)
  answer <- c((-b + sqrt(b^2 - 4 * a * c)) / (2 * a),
              (-b - sqrt(b^2 - 4 * a * c)) / (2 * a))
  answer <- Re(answer) # get only real component
  return(answer)
}

get_max_eig <- function(eval_vec, fd_x, fd_y){
  # why do I need to add artificial data here?
  eval_result_x <- eval.bifd(sevalarg = c(unname(eval_vec[1]),1), tevalarg = c(unname(eval_vec[2]),1), bifd = fd_x,
                 sLfdobj = 1, tLfdobj = 1)[1,]
  eval_result_y <- eval.bifd(sevalarg = c(unname(eval_vec[1]),1), tevalarg = c(unname(eval_vec[2]),1), bifd = fd_y,
                             sLfdobj = 1, tLfdobj = 1)[1,]
  local_jacobian <- unname(rbind(eval_result_x, eval_result_y))
  
  eig_penalty <- max(quadratic_real(1, -sum(diag(local_jacobian)), det(local_jacobian)))
  if(is.na(eig_penalty)){ eig_penalty <- 0}
  
  return(eig_penalty)
  
}


## Eval

some_noise <- list(name = "low stiff; sigma^2 = 0.02", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = 0.02, var_y = 0.02,
                   lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 8, noise_seed = 8)
experiment_list = list(some_noise)
experiment_data <- generate_data_object(experiment_list, noisy_smooth_basis = 96)

spline_control <- list(method = "spline",  params = list(lambda = 1e0, norder = 4, nbasis = 12, side_info = list()))
experiment_estimators <- list(spline_control) # optim uses first spline fit as a hot start
experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)
plot1 <- list(type = "field", experiments = list(c(data = 1, estimator = 1)))
visualize_results(experiment_results, list(plot1))

optim_fit <- calculate_eigpen_spline_gradient_field(experiment_data[[1]]$limit_cycle_tail,
                                                    experiment_data[[1]]$grid$x_grid, experiment_data[[1]]$grid$y_grid,
                                                    experiment_results[[1]]$estimates[[1]]$sfd_list$fdx$coefs,
                                                    experiment_results[[1]]$estimates[[1]]$sfd_list$fdy$coefs,
                                                    norder = 4, nbasis = 12)

generate_optim_plots <- function(experiment_data, optim_fit){
  optim_field <- ggplot_field(experiment_data[[1]]$limit_cycle_tail, experiment_data[[1]]$grid$eval_grid,
                              optim_fit$field,
                              title="")
  print(optim_field)
  
  outcome_list <- experiment_data
  
  outcome_list[[1]]$estimates <- list(list(method = "spline", sfd_list = optim_fit, estimated_field = optim_fit$field))
  
  plot_solution_paths(list(c(1,1)), outcome_list)
}

generate_optim_plots(experiment_data, optim_fit)





