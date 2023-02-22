NW_regression <- function(eval_grid, data){
  x_grid <- eval_grid[,1]
  y_grid <- eval_grid[,2]
  nw_results_array <- array(-5, c(nrow(eval_grid), 2))
  
  bw_matrix <- matrix(c(.1,0,0,0.1),2)
  for (i in 1:nrow(eval_grid)){
    nw_results_array[i,] <-  eval_gaussian_NW(eval_grid[i,1], eval_grid[i,2], data, bw_matrix)
  }
  # for (i in 1:length(x_grid)){
  #   for (j in 1:length(y_grid)){
  #     nw_results_array[i,j,] <- eval_gaussian_NW(x_grid[i], y_grid[j], data, bw_matrix)
  #   }
  # }
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

plot_field <- function(ls_samples, eval_grid, gradient_data,title_str=""){
  # plot a gradient field
  
  eta = 0.1
  x_grid = eval_grid[,1]
  y_grid = eval_grid[,2]
  
  plot(ls_samples$x, ls_samples$y,
       xlim=c(min(x_grid)-1,max(x_grid)+1),
       ylim=c(min(y_grid)-1,max(y_grid)+1), 
       main = title_str, xlab = "x", ylab = "y")
  
  for (i in 1:nrow(eval_grid)){
    # prevent over plotting
    if (i %% 2 == 0){
      Arrows(x_grid[i], y_grid[i], x_grid[i] + eta*gradient_data[i,1], y_grid[i] + eta*gradient_data[i,2], lwd=0.75) 
    }
  }
}

plot_gradient_path_single <- function(ls_samples, eval_grid, gradient_data, bw, num_reps = 9, title =""){
  # plots NW estimated trajectories
  nw_trajectory <- generate_nw_path(ls_samples, bw, num_samples = 500)
  trajectory_tibble <- tibble(x = nw_trajectory[,1], y = nw_trajectory[,2], u = NA, v = NA, run = 1)
  
  if (num_reps > 1){
    for (i in 2:num_reps){
      nw_trajectory <- generate_nw_path(ls_samples, bw, num_samples = 500)
      new_trajectory_tibble <- tibble(x = nw_trajectory[,1], y = nw_trajectory[,2], u = NA, v = NA, run = i)
      trajectory_tibble <- bind_rows(trajectory_tibble, new_trajectory_tibble)
    }
  }
  
  sample_tibble <- tibble(x = ls_samples[,1], y = ls_samples[,2], u = NA, v = NA)
  gradient_tibble <- tibble(x = eval_grid[,1], y = eval_grid[,2], u = gradient_data[,1], v = gradient_data[,2])
  
  path_plot <- ggplot(gradient_tibble, aes(x = x, y = y, u = u, v = v)) +
    geom_quiver(color = "#003262") +
    geom_path(data = trajectory_tibble, aes(x = x, y = y)) +
    geom_point(data = sample_tibble, aes(x = x, y = y), color = "#FDB515") +
    geom_point(x = as.numeric(trajectory_tibble[1,1]), y =as.numeric(trajectory_tibble[1,2]), color = "red", size = 3) +
    labs(title = title) +
    facet_wrap(~run)
  
  return(path_plot)
  
}

####### Old Code

spline_gradient <- function(data, x_grid, y_grid){
  # Estimates the gradient field after fitting b-splines to the limit cycle data
  #
  ## Inputs:
  # data (data.frame): contains data with columns in [x, y, f_x, f_y]
  # x_grid (numeric; vector): x-axis grid points to evaluate extrapolation over
  # y_grid (numeric; vector): y-axis grid points to evaluate extrapolation over
  #
  ## Outputs:
  # output_list (named list): see construction for contents; TODO: update documentation
  
  # create basis functions with (deault = 10) knots
  # TODO: adapt number of basis functions to the size of the limit cycle?
  xbasis = create.bspline.basis(range=c(min(x_grid),max(x_grid)),norder=4,nbasis=12)
  ybasis = create.bspline.basis(range=c(min(y_grid),max(y_grid)),norder=4,nbasis=12)
  xbasis.vals = eval.basis(data$x,xbasis)
  ybasis.vals = eval.basis(data$y,ybasis)
  
  # What we need is the evaluation of phi_j(x_i)*psi_k(x_j) for
  # each x and y. This creates 144 rows. I'll produce this using
  # Kronecker products
  Xmat = (xbasis.vals%x%matrix(1,1,12)) * (matrix(1,1,12)%x%ybasis.vals)
  
  # Here the columns of xbasis.vals are repeated 12 times in sequence
  # while each column of ybasis.vals is repeated 12 times together.
  
  # Now we need a penalty matrix. We can get the penalty for one
  # basis from
  xPen = eval.penalty(xbasis, 2)
  yPen = eval.penalty(ybasis, 2)
  
  # to create the combined penalty we take the same Kronecker
  # product form
  
  allPen = xPen%x%diag(12) + diag(12)%x%yPen
  #View(allPen)
  
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
  
  # build a penalty matrix to plot...
  #sfd_x.penalty = bifd(t(matrix(coefs_x,12,12)), xbasis, ybasis)
  #sfd_y.penalty = bifd(t(matrix(coefs_y,12,12)), xbasis, ybasis)
  #smat_x.penalty = eval.bifd(x_grid,y_grid,sfd_x)
  #smat_y.penalty = eval.bifd(x_grid,y_grid,sfd_y)
  #View(allPen)
  #print(dim(xPen))
  
  
  traj = lsoda(c(data$x[1],data$y[1]),seq(0,100,by=0.1),bifd_spline_gradient,0,sfd_x=sfd_x,sfd_y=sfd_y)
  
  output_list <- list(x_grad_bifd = sfd_x, x_grad_eval = smat_x, y_grad_bifd = sfd_y, y_grad_eval = smat_y, 
                      trajectory = traj, second.deriv_penalty = allPen)
  return(coefs_y)
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


plot_observations <- function(data){
  # helper function which generates a scatter plot of observations 
  ggplot(data,aes(x=x,y=y)) +
    geom_point() + 
    labs(x="x",y="y",title="Observations from a Dynamical System")
}

#' This function estimates the gradient field over a grid given limit cycle data
#' 4 different estimators are compared: NW KDE, LOESS, kNN, and bSplines
#' The output is a series of plots allowing direct comparison between methods
#'
#' @param data (data.frame): sequential solution trajectory data with columns in [x, y, f_x, f_y]
#' @param model_str (string): name of the model to generate the gradient field of; see `generate_limit_cycle_data()`
#'     This is used to generate the true gradient field
#' @param model_params (vector): vector of parameters for the specified model
#' @param tail_n (integer)[optional]: only the last `tail_n` samples are considered as the limit cycle
#' @param x_grid_size (integer)[optional]: number of x-axis samples in extrapolation grid
#' @param y_grid_size integer)[optional]: number of y-axis samples in extrapolation grid
#' @param extrapolation_size (numeric)[optional]: axis limits will be the extreme observations +/- this parameter
#' @param nw_bandwidth (numeric)[optional]: Scalar bandwidth to use in NW kernel regression
#' @param loess_bandwidth (numeric)[optional]: Scalar bandwidth to use in LOESS
#' @param model_title (string)[optional]: Name of ODE model to put on plots e.g. "Stiff VdP"
#'
#' @return
#'
#' @examples
#' evaluate_gradient_methods(vp_data, model_str = "van_der_pol", model_params = c(0.5),
#'  extrapolation_size = 1, model_title = "Van der Pol; mu = 0.5")
evaluate_gradient_methods <- function(data, model_str, model_params, tail_n = 700, 
                                      x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5,
                                      nw_bandwidth = 0.1, loess_bandwidth = 0.1,
                                      model_title = ""){
  limit_cycle.data <- tail(data, n = tail_n)
  # convert data to matrix (no labels) for C++ implementations
  limit_cycle.data.matrix <- unname(as.matrix(limit_cycle.data))
  
  # get extreme points of the data
  xmin <- floor(min(limit_cycle.data$x))
  xmax <- ceiling(max(limit_cycle.data$x))
  ymin <- floor(min(limit_cycle.data$y))
  ymax <- ceiling(max(limit_cycle.data$y))
  
  # grid to evaluate gradient extrapolation over
  x_grid <- seq(xmin - extrapolation_size , xmax + extrapolation_size, len = x_grid_size)
  y_grid <- seq(ymin - extrapolation_size, ymax + extrapolation_size, len = y_grid_size)
  eval_grid <- unname(as.matrix(expand.grid(x_grid,y_grid)))
  
  # Generate plots
  
  # get truth over evaluation grid
  true_field <- generate_grid_data(model_str, model_params, eval_grid)
  #true_field_plot <- ggplot_field(limit_cycle.data, eval_grid, true_field, title="Truth")
  
  # kNN
  nn_fit <- eval_knn_fit(eval_grid, limit_cycle.data.matrix,1)
  #nn_field_plot <- ggplot_field(limit_cycle.data, eval_grid, nn_fit, title="1-NN")
  
  multi_nn_fit <- eval_knn_fit(eval_grid, limit_cycle.data.matrix,400)
  #multi_nn_field_plot <- ggplot_field(limit_cycle.data, eval_grid, multi_nn_fit, title="400-NN")
  
  #nn_gradient_plot <- plot_gradient_path_multi(limit_cycle.data.matrix, eval_grid,
  #                                               model_str, "knn", nn_fit, 1,
  #                                               title = paste(model_title, "lsoda Solutions Along 1-NN Gradient Field", sep = ": "))
  
  # Truth
  #true_gradient_plot <- plot_gradient_path_multi(limit_cycle.data.matrix, eval_grid,
  #                                                 model_str, "truth", true_field, model_params,
  #                                                 title = paste(model_title, "lsoda Solutions Along True Gradient Field", sep = ": "))
  
  # b-splines (mix of R and C++)
  # to avoid simulated solution path out of support issues in `lsoa`, use a larger support
  #x_grid_spline <- seq(xmin - spline_support_extension , xmax + spline_support_extension, len = x_grid_size)
  #y_grid_spline <- seq(ymin - spline_support_extension, ymax + spline_support_extension, len = y_grid_size)
  spline_result <- calculate_spline_gradient_field(limit_cycle.data, x_grid, y_grid, plot_penalty = FALSE)
  #spline_field_plot <- ggplot_field(limit_cycle.data, eval_grid, spline_result$field, title="Spline")
  #spline_gradient_plot <- plot_gradient_path_multi(limit_cycle.data.matrix, eval_grid,
  #                                                 model_str, "spline", spline_result$field, list(fdx = spline_result$fdx, fdy = spline_result$fdy),
  #                                             title = paste(model_title, "lsoda Solutions Along Spline Gradient Field", sep = ": "))
  
  # NW
  nw_bw_matrix <- nw_bandwidth*diag(2) # TODO: Vary bandwidth (okay as is in no noise setting)
  NW_regression_result <- NW_regression_cpp(eval_grid, limit_cycle.data.matrix, nw_bw_matrix)
  #nw_field_plot <- ggplot_field(limit_cycle.data, eval_grid, NW_regression_result, title="NW Kernel Regression")
  #nw_gradient_plot <- plot_gradient_path_multi(limit_cycle.data.matrix, eval_grid,
  #                                             model_str, "nw", NW_regression_result, nw_bandwidth,
  #                                             title = paste(model_title, "lsoda Solutions Along NW Gradient Field", sep = ": "))
  
  # LOESS
  loess_fit <- eval_loess_fit(eval_grid, limit_cycle.data.matrix, loess_bandwidth)
  #loess_field_plot <- ggplot_field(limit_cycle.data, eval_grid, loess_fit, title="LOESS")
  #loess_gradient_plot <- plot_gradient_path_multi(limit_cycle.data.matrix, eval_grid, 
  #                                                model_str, "loess", loess_fit, loess_bandwidth,
  #                                                title = paste(model_title, "lsoda Solutions Along LOESS Gradient Field", sep = ": "))
  
  # Output plots
  
  #field_plots <- ggarrange(true_field_plot, nn_field_plot, multi_nn_field_plot,
  #                         spline_field_plot, nw_field_plot, loess_field_plot,
  #                         ncol=3, nrow=2)
  #field_plots <- annotate_figure(field_plots, top = text_grob(label = model_title))
  #print(field_plots)
  plot_solution_paths(limit_cycle.data.matrix, eval_grid, NA, NA, NA, 
                      samples_per_path = 500,
                      methods_toplot = c("truth"),
                      title ="")
  #print(true_gradient_plot)
  #print(nn_gradient_plot)
  #print(nw_gradient_plot)
  #print(loess_gradient_plot)
  #print(spline_gradient_plot)
}


plot_gradient_path_multi <- function(ls_samples, eval_grid, model_str, method_str, gradient_data, params, num_reps = 9, samples_per_path = 500, title =""){
  # generate random solution paths using the NW or LOESS estimated gradient field
  plot_list <- list()
  sample_tibble <- tibble(x = ls_samples[,1], y = ls_samples[,2], u = NA, v = NA)
  gradient_tibble <- tibble(x = eval_grid[,1], y = eval_grid[,2], u = gradient_data[,1], v = gradient_data[,2])
  
  
  for (i in 1:num_reps){
    if ( (method_str == "truth") & (model_str == "van_der_pol") ){
      trajectory <- generate_true_path_vdp(ls_samples, params, num_samples = samples_per_path)
    }
    else if (method_str == "loess"){
      trajectory <- generate_loess_path(ls_samples, params, num_samples = samples_per_path)
    }
    else if (method_str == "nw"){
      trajectory <- generate_nw_path(ls_samples, params, num_samples = samples_per_path)
    }
    else if (method_str == "spline"){
      trajectory <- generate_spline_path(ls_samples, params, num_samples = samples_per_path)
    }
    else if (method_str == "knn"){
      trajectory <- generate_knn_path(ls_samples, params, num_samples = samples_per_path)
    }
    else{
      stop("Method not implemented")
    }
    
    trajectory_tibble <- tibble(x = trajectory[,1], y = trajectory[,2], u = NA, v = NA, run = i)
    path_plot <- ggplot(gradient_tibble, aes(x = x, y = y, u = u, v = v)) +
      geom_quiver(color = "#003262") +
      geom_point(data = sample_tibble, aes(x = x, y = y), color = "#FDB515") +
      geom_path(data = trajectory_tibble, aes(x = x, y = y)) +
      geom_point(x = as.numeric(trajectory_tibble[1,1]), y =as.numeric(trajectory_tibble[1,2]), color = "red", size = 3) +
      facet_wrap(~run)
    plot_list <- append(plot_list, list(path_plot))
  }
  
  path_plots <- ggarrange(plotlist = plot_list, nrow = ceiling(length(plot_list)/3), ncol = 3)
  path_plots <- annotate_figure(path_plots, top = text_grob(label = title))
  return(path_plots)
  
}


# # DEPRICATED BELOW:
# #abhi_data <- generate_limit_cycle_data("abhi", c())
# #evaluate_gradient_methods(abhi_data, "abhi", c(), extrapolation_size = 0.5, model_title = "Abhi's Van der Pol")
# 
# stiff_mu <- 20
# vp_stiff_samples <- generate_limit_cycle_data("van_der_pol", c(stiff_mu))
# vp_stiff_list <- list(system = "van_der_pol", params = c(stiff_mu), limit_cycle_samples = vp_stiff_samples)
# evaluate_gradient_methods(vp_data.stiff, model_str = "van_der_pol", model_params = c(stiff_mu), extrapolation_size = 1, model_title = "Van der Pol; mu = 20")
#  
# #nonstiff_mu <- 1.5
# #vp_data.nonstiff <- generate_limit_cycle_data("van_der_pol", c(nonstiff_mu))
# #evaluate_gradient_methods(vp_data.nonstiff, model_str = "van_der_pol", model_params = c(nonstiff_mu), extrapolation_size = 1, model_title = "Van der Pol; mu = 1.5")


#' Smooth Noisy Samples
#' 
#' This function takes noisy data as an input and returns a smoothed version using splines
#' If gradient data is not observed, this function can also impute that after smoothing
#' Additional control is given over the splines through optional parameters
#'
#' TODO: Complete
#' @param orig_data 
#' @param data 
#' @param impute_grad 
#' @param lambda 
#' @param norder 
#' @param nbasis 
#' @param penalty_order 
#' @param max_t 
#' @param smooth_type 
#' @param title 
#'
#' @return
#' @export
#'
#' @examples
spline_smooth_noisy_samples <- function(orig_data, data, impute_grad = F, 
										lambda = 1e-10, norder = 4, nbasis = 48, penalty_order = 2, max_t = 1,
										smooth_type = "bspline", title = ""){
	
	if (smooth_type == "orig"){
		return(orig_data)
	}
	
	time_grid <- seq(0,max_t,length.out=nrow(data))
	penalty.Lfd = int2Lfd(penalty_order)
	
	xbasis = create.bspline.basis(rangeval=c(0,max_t),norder=norder,nbasis=nbasis)
	xbasis.fdPar = fdPar(xbasis,penalty.Lfd,lambda)
	ybasis = create.bspline.basis(rangeval=c(0,max_t),norder=norder,nbasis=nbasis)
	ybasis.fdPar = fdPar(ybasis,penalty.Lfd,lambda)
	
	# code taken from Giles' tutorial
	
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
	
	# We can plot what we have
	
	#plot(lambdas,x.mean.gcv,type='b',log='x')
	#plot(lambdas,y.mean.gcv,type='b',log='x')
	# Lets select the lowest of these and smooth
	
	best_x = which.min(x.mean.gcv)
	lambdabest_x = lambdas[best_x]
	best_y = which.min(y.mean.gcv)
	lambdabest_y = lambdas[best_y]
	
	xbasis.fdPar$lambda = lambdabest_x
	ybasis.fdPar$lambda = lambdabest_y
	cat("After GCV, Lambda x:", lambdabest_x, " and Lambda y:", lambdabest_y)
	x_smooth <- smooth.basis(argvals = seq(0,max_t,len=length(data$x)), y=data$x, fdParobj=xbasis.fdPar)
	y_smooth <- smooth.basis(argvals = seq(0,max_t,len=length(data$y)), y=data$y, fdParobj=ybasis.fdPar)
	
	smoothed_samples <- cbind(eval.fd(seq(0,max_t,len=length(data$x)),x_smooth$fd),
							  eval.fd(seq(0,max_t,len=length(data$y)),y_smooth$fd))
	colnames(smoothed_samples) <- c("x","y")
	smoothed_samples <- as.data.frame(smoothed_samples)
	
	verbose = FALSE
	if (verbose){
		tps_x <- Tps(time_grid, data$x)$fitted.values
		tps_y <- Tps(time_grid, data$y)$fitted.values
		
		plotting_df <- rbind(tibble(x = orig_data$x, y = orig_data$y, pair = 1:nrow(orig_data), label = "Truth"),
							 tibble(x = data$x, y = data$y, pair = 1:nrow(orig_data), label = "Noisy Samples"),
							 tibble(x = smoothed_samples$x, y = smoothed_samples$y, pair = 1:nrow(orig_data), label = "b-Spline"),
							 tibble(x = tps_x, y = tps_y, pair = 1:nrow(orig_data), label = "TPS"))
		
		# plot all smoothers separately
		smooth_side_by_side <- plotting_df %>%
			ggplot(aes(x=x, y=y)) +
			geom_point() +
			labs(title=paste0(title,"; : ",norder,"; # of b-spline basis: ", nbasis,"; x GCV lambda: ", round(lambdabest_x,2),"; y GCV lambda: ", round(lambdabest_y,2))) +
			facet_wrap(~label)
		#print(smooth_side_by_side)
		file_name = paste0("comparison_",norder,"_",nbasis,"_mu20.png")
		file_path = paste0("Result_Images/2022-07-18/",file_name)
		ggsave(file_path,smooth_side_by_side,width=14,height=7)
		
		# plot all smoothers together
		connected_smooth <- plotting_df %>%
			filter(label != "TPS") %>%
			ggplot(aes(x=x,y=y, color=label)) +
			geom_point(aes(fill=label),size=3) +
			geom_line(aes(group = pair),color="grey")
		#print(connected_smooth)
		
		# visualize derivative of b-spline smooth
		
		time_axis_tibble <- rbind(tibble(t=seq(0,max_t,len=length(data$x)),pos=smoothed_samples$x,estimate="b-spline",axis="x"),
								  tibble(t=seq(0,max_t,len=length(data$x)),pos=orig_data$x,estimate="truth",axis="x"),
								  tibble(t=seq(0,max_t,len=length(data$x)),pos=data$x,estimate="noisy",axis="x"),
								  tibble(t=seq(0,max_t,len=length(data$y)),pos=smoothed_samples$y,estimate="b-spline",axis="y"),
								  tibble(t=seq(0,max_t,len=length(data$y)),pos=orig_data$y,estimate="truth",axis="y"),
								  tibble(t=seq(0,max_t,len=length(data$y)),pos=data$y,estimate="noisy",axis="y"))
		b_spline_position_plot <- time_axis_tibble %>%
			ggplot(aes(x=t, y=pos, color = estimate))+
			geom_point(data = . %>% filter(estimate %in% c("truth")))+
			geom_point(data = . %>% filter(estimate %in% c("noisy")))+
			geom_line(data = . %>% filter(estimate %in% c("b-spline")))+
			labs(title=paste0("Position; ",title,"; Order: ",norder,"; # of b-spline basis: ", nbasis,"; x GCV lambda: ", round(lambdabest_x,2),"; y GCV lambda: ", round(lambdabest_y,2))) +
			facet_wrap(~axis, scales = "free")
		#print(b_spline_position_plot)
		file_name = paste0("position_",norder,"_",nbasis,"_mu20.png")
		file_path = paste0("Result_Images/2022-07-18/",file_name)
		ggsave(file_path,b_spline_position_plot,width=14,height=7)
		
		# visualize derivative of b-spline smooth
		x_basis_deriv <- deriv.fd(x_smooth$fd, 1)
		x_basis_deriv_eval <- eval.fd(evalarg = seq(0,max_t,len=length(data$x)), fdobj=x_basis_deriv)
		y_basis_deriv <- deriv.fd(y_smooth$fd, 1)
		y_basis_deriv_eval <- eval.fd(evalarg = seq(0,max_t,len=length(data$y)), fdobj=y_basis_deriv)
		
		derivative_tibble <- rbind(tibble(t=seq(0,max_t,len=length(data$x)),grad=x_basis_deriv_eval,estimate="b-spline",axis="x"),
								   tibble(t=seq(0,max_t,len=length(data$x)),grad=data$f_x,estimate="truth",axis="x"),
								   tibble(t=seq(0,max_t,len=length(data$y)),grad=y_basis_deriv_eval,estimate="b-spline",axis="y"),
								   tibble(t=seq(0,max_t,len=length(data$y)),grad=data$f_y,estimate="truth",axis="y"))
		b_spline_deriv_plot <- derivative_tibble %>%
			ggplot(aes(x=t, y=grad, color = estimate))+
			geom_point(data = . %>% filter(estimate %in% c("truth")))+
			geom_line(data = . %>% filter(estimate %in% c("b-spline")))+
			labs(title=paste0("Gradient; ",title,"; Order: ",norder,"; # of b-spline basis: ", nbasis,"; x GCV lambda: ", round(lambdabest_x,2),"; y GCV lambda: ", round(lambdabest_y,2))) +
			facet_wrap(~axis, scales = "free")
		#print(b_spline_deriv_plot)
		file_name = paste0("gradient_",norder,"_",nbasis,"_mu20.png")
		file_path = paste0("Result_Images/2022-07-18/",file_name)
		ggsave(file_path,b_spline_deriv_plot,width=14,height=7)
	}
	
	if(impute_grad){
		x_basis_deriv <- deriv.fd(x_smooth$fd, 1)
		x_basis_deriv_eval <- eval.fd(evalarg = seq(0,max_t,len=length(data$x)), fdobj=x_basis_deriv)
		y_basis_deriv <- deriv.fd(y_smooth$fd, 1)
		y_basis_deriv_eval <- eval.fd(evalarg = seq(0,max_t,len=length(data$y)), fdobj=y_basis_deriv)
		
		smoothed_samples <- smoothed_samples %>% mutate("f_x" = x_basis_deriv_eval, "f_y" = y_basis_deriv_eval)
	}
	
	if (smooth_type == "bspline"){
		return(smoothed_samples)
	} else if (smooth_type == "tps"){
		smoothed_samples <- matrix(c(tps_x,tps_y), ncol = 2)
		return(smoothed_samples)
	} else {
		stop("Spline smooth not implemented")
	}
}