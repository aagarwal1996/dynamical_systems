#############
## Imports ##
#############

# Packages
library(fda)
library(tidyverse)
library(paletteer)

# Local Files
# functions to generate data from specific system
source(here::here('Data_Generation','data_generation.R'))
source(here::here('fit_helpers.R'))
source(here::here('Estimation_Methods','bspline_gd_update.R'))

###################
# Data Generation #
###################

# This is simple enough to be built by hand for real world data
generate_experiment_trial <- function(experiment_name, system_name = NULL, system_parameters = list(),
									  random_seed = NULL, num_samples = 1000, delta_t = 0.1,
									  replicates = 1, var_x = 0, var_y = 0){
	# build a list of lists to hold system information and samples
	trial_object <- list()
	trial_object$experiment_name <- experiment_name
	trial_object$system_name <- system_name
	trial_object$system_parameters <- system_parameters

	# generate noisy samples
	if (is.numeric(random_seed)){ set.seed(random_seed) }
	trial_object$replicates <- list()
	if (replicates > 0){
		for (i in 1:replicates){
			trial_object$replicates[[i]] <- list(samples = generate_limit_cycle_data(system_name, system_parameters,
																					 var_x, var_y, num_samples, delta_t))
		}
	}

	return(trial_object)
}

generate_limit_cycle_data <- function(system_name, system_parameters, var_x, var_y,
									  num_samples = 1000, delta_t = 0.1){
	if (system_name == "van_der_pol"){
		sampled_data <- generate_van_der_pol_samples(system_parameters, num_samples = num_samples,
			sample_density = delta_t)
	} else if (system_name == "rzma"){
		sampled_data <- generate_rzma_samples(system_parameters, num_samples = num_samples,
			sample_density = delta_t)
	} else if(system_name == "lotka_volterra"){
		sampled_data <- generate_lv_samples(system_parameters, num_samples = num_samples,
			sample_density = delta_t)
	} else if(system_name == "log_lotka_volterra"){
		sampled_data <- generate_log_lv_samples(system_parameters, num_samples = num_samples,
			sample_density = delta_t)
	} else if(system_name == "abhi"){
		sampled_data <- get_abhi_data()
	} else{
		stop('Error: ', system_name, ' data generation not implemented.')
	}

	if ((var_x != 0) | (var_y != 0)){
		# add Gaussian noise to position of observations
		noise_matrix <- mvrnorm(nrow(sampled_data), c(0,0,0,0), diag(c(var_x,var_y,0,0)))
		noisy_samples <- sampled_data + noise_matrix
	}

	colnames(sampled_data) <- c("x","y","f_x","f_y")
	return(sampled_data)
}

smooth_trial <- function(experiment_trial, smooth = "bspline",
						 smooth_lambda = 1e-12, num_basis_fn = 48, impute_grad = FALSE){
	for (i in 1:length(experiment_trial$replicates)){
		# specify smooth = "none" to fit estimators to original data downstream
		experiment_trial$replicates[[i]]$smooth <- spline_smooth_samples(experiment_trial$replicates[[i]]$samples,
													impute_grad, num_basis_fn, max_t = 1,
													smooth_lambda, norder = 6, penalty_order = 2,
													smooth_type = smooth)
		# TODO: Detect tail automatically
		n_samples <- dim(experiment_trial$replicates[[i]]$samples)[1]
		experiment_trial$replicates[[i]]$smooth_tail <- tail(experiment_trial$replicates[[i]]$smooth, 0.7*n_samples)
		experiment_trial$replicates[[i]]$smooth_params <- list(method = smooth, n_basis_fn = num_basis_fn)
	}
	return(experiment_trial)
}

#######################
# Dynamics Estimation #
#######################

apply_estimators <- function(experiment_trial, estimators_list){
	for (i in 1:length(experiment_trial$replicates)){
		experiment_trial$replicates[[i]]$estimators <- list()
		for (j in 1:length(estimators_list)){
			estimator <- estimators_list[[j]]
			field_fn <- get_field_fn(experiment_trial$replicates[[i]]$smooth_tail, estimator,
										experiment_trial$system_name, experiment_trial$system_parameters)
			experiment_trial$replicates[[i]]$estimators[[j]] <- list(estimator = estimator, fn = field_fn)
		}
	}
	return(experiment_trial)
}

get_field_fn <- function(data, estimator, sys_name, sys_params){
	field_fn <- NA
	if (estimator$method == "truth"){
		if(sys_name == "van_der_pol"){
			field_fn <- function(v){ eval_van_der_pol_gradient(0,v,sys_params) }
		} else if(sys_name == "rzma"){
			field_fn <- function(v){ eval_rzma_gradient(0,v,sys_params) }
		} else if(sys_name == "lotka_volterra"){
			field_fn <- function(v){ eval_lv_gradient(0,v,sys_params) }
		} else if(sys_name == "log_lotka_volterra"){
			field_fn <- function(v){ eval_log_lv_gradient(0,v,sys_params) }
		} else {
			warning("Valid model not specified for true field. Zero field returned.")
			field_fn <- function(v){ return(list(rep(0,length(v)))) }
		}
	} else if (estimator$method == "vanilla_spline"){
		fda_fit <- fit_vanilla_spline_field(data, estimator$side_info, estimator$norder,
			estimator$nbasis,estimator$penalty_order,estimator$lambda)
		field_fn <- function(v){ return(bifd_spline_gradient(0,v,fda_fit)) }
	} else if (estimator$method == "gd_spline"){
		gd_fit <- calculate_gd_spline_gradient_field(data, estimator$gd_params, estimator$side_info,
													 estimator$norder,estimator$nbasis,estimator$penalty_order,estimator$lambda)
		field_fn <- function(v){ return(bifd_spline_gradient(0,v,gd_fit)) }
	} else if (estimator$method == "nw"){
		nw_params <- list(data, estimator$bandwidth_matrix)
		field_fn <- function(v){ return(nw_get_grad(0,v,nw_params)) }
	} else if (estimator$method == "loess"){
		loess_params <- list(data, estimator$bandwidth_scalar)
		field_fn <- function(v){ return(loess_get_grad(0,v,loess_params)) }
	}else if (estimator$method == "knn"){
		knn_params <- list(data, estimator$k)
		field_fn <- function(v){ return(knn_get_grad(0,v,knn_params)) }
	} else{ warning("Invalid estimator specified.") }
	return(field_fn)
}

##############
# Evaluation #
##############

# Grid MSE

get_grid_mse <- function(fit_samples, x_grid_count = 24, y_grid_count = 24, frame_frac = 0.2){
	# TODO: Replace MSE with a lambda cost function
	x_range = max(fit_samples$smooth$x) - min(fit_samples$smooth$x)
	x_grid <- seq(min(fit_samples$smooth$x) - frame_frac*x_range,
				  max(fit_samples$smooth$x) + frame_frac*x_range, length.out =  x_grid_count)
	y_range = max(fit_samples$smooth$y) - min(fit_samples$smooth$y)
	y_grid <- seq(min(fit_samples$smooth$y) - frame_frac*y_range,
				  max(fit_samples$smooth$y) + frame_frac*y_range, length.out =  y_grid_count)
	eval_grid <- unname(as.matrix(expand.grid(x_grid,y_grid))) # convert to N x 2 matrix
	
	predictions <- list()
	mse_results <- c()
	for (i in 1:length(fit_samples$estimators)){
		preds <- apply(eval_grid,1,fit_samples$estimators[[i]]$fn)
		preds <- do.call(rbind,lapply(preds,unlist))
		predictions[[i]] <- preds
		delta <- predictions[[i]] - predictions[[1]]
		mse <- mean(apply(delta,1,function(x){ sqrt(sum(x^2) + 1e-10)}))
		mse_results <- c(mse_results, mse)
	}
	
	
	return(mse_results)
}

# Delta Tube MSE

get_dt_mse <- function(fit_samples, radius_val = 0.2, sample_frac = 0.5, plot = FALSE){
	data <- fit_samples$smooth_tail
	delta_ring_samples <- data[sort(sample(nrow(data),size=floor(nrow(data)*sample_frac),replace=FALSE)),]
	delta_ring_sideinfo <- t(apply(delta_ring_samples,1,get_normal_vectors, radius = radius_val))
	side_information_pos <- rbind(delta_ring_sideinfo[,c(1,2)],delta_ring_sideinfo[,c(3,4)])
	colnames(side_information_pos) <- c("x","y")
	
	predictions <- list()
	mse_results <- c()
	for (i in 1:length(fit_samples$estimators)){
		preds <- apply(side_information_pos,1,fit_samples$estimators[[i]]$fn)
		preds <- do.call(rbind,lapply(preds,unlist))
		predictions[[i]] <- preds
		delta <- predictions[[i]] - predictions[[1]]
		mse <- mean(apply(delta,1,function(x){ sqrt(sum(x^2) + 1e-10)}))
		mse_results <- c(mse_results, mse)
	}
	
	# format for plotting
	delta_list <- lapply(predictions, function(x){x - predictions[[1]]})[-1] # remove truth
	delta_list <- lapply(delta_list, as_tibble)
	delta_list <- lapply(delta_list, function(x){names(x) = c("f_x","f_y");return(x)})
	delta_list <- lapply(seq_along(delta_list), function(i){delta_list[[i]] %>%
			mutate(estimator=as.factor(paste0(fit_samples$estimators[[i+1]]$estimator$method,":",round(mse_results[i+1],3)))) %>%
			mutate(x = as.numeric(side_information_pos[,1]), y = as.numeric(side_information_pos[,2]))})
	plot_tibble <- do.call(rbind,delta_list)
	plot_tibble <- plot_tibble %>%
		mutate(mse = sqrt(f_x^2 + f_y^2 + 1e-12))
	
	if (plot){
		field_plot <- ggplot(plot_tibble,aes(x=x,y=y,color=mse)) +
			geom_point(alpha=0.3) +
			scale_colour_gradient2() +
			geom_point(data = fit_samples$smooth_tail, aes(x = x, y = y), color = "#FDB515", alpha = 0.3) +
			facet_wrap(~estimator) +
			labs(title = "Delta Ring MSE Across Methods")
		print(field_plot)
	}
	
	return(mse_results)
}

get_normal_vectors <- function(sample, radius){
	rotation_matrix <- matrix(c(0,-1,1,0),ncol=2)
	magnitude = sqrt(sample["f_x"]^2 + sample["f_y"]^2)
	point <- c(sample["x"], sample["y"])
	normal_vec_out <- (rotation_matrix %*% c(sample["f_x"], sample["f_y"]))/magnitude
	
	ring_sample <- c(point + radius*normal_vec_out, point - radius*normal_vec_out) 
	# normal out x, y; normal in x, y
	return(ring_sample)
}

# Solution Path Integral

	
apply_transient_error <- function(experiment_trial,
								  ic_counts, ic_radii, ic_sigma, ic_filter_fn,
								  integral_N){
	# for each estimator and replicate, get the integrated error
	# IMPORTANT: This implicitly assumes for all relative errors that the truth is the first estimator
	for (i in 1:length(experiment_trial$replicates)){
		experiment_trial$replicates[[i]]$ic <- sample_ic(experiment_trial$replicates[[i]]$smooth_tail,
			ic_counts, ic_radii, ic_sigma, filter_fn=ic_filter_fn)
		estimators_list <- experiment_trial$replicates[[i]]$estimators
		for (j in 1:length(estimators_list)){
			t_star <- get_lc_cutoff(experiment_trial$replicates[[i]]$smooth)
			experiment_trial$replicates[[i]]$estimators[[j]]$transient_mse <- apply(
				experiment_trial$replicates[[i]]$ic, 1,  calc_transient_error,
				experiment_trial$replicates[[i]]$estimators[[1]]$fn,
				experiment_trial$replicates[[i]]$estimators[[j]]$fn,
				t_star, integral_N)
		}
	}
}

calc_transient_error <- function(initial_condition, true_grad, est_grad, t_star, N, losses, plot = FALSE, return_transients = FALSE){
	# Step 1: Fit estimated transient from initial condition
	time_vec.true <- seq(0,t_star,length.out=N)
	true_grad_augment <- function(t,v,parms){return(true_grad(v))} # play nice with deSolve
	true_transient <- deSolve::lsoda(initial_condition, time_vec.true, true_grad_augment)
	if(is.null(nrow(true_transient))){
		true_transient <- matrix(c(true_transient),nrow=1,dimnames=list(NULL,names(true_transient)))
		true_transient <- cbind(true_transient,source = 1)
	} else {true_transient <- cbind(true_transient,source = 1)}
	true_transient_true_grad <- matrix(unlist(apply(true_transient,1,function(row){true_grad(c(x=row[2],y=row[3]))})),ncol=2,byrow=TRUE)
	
	time_vec.est <- seq(0,10*t_star,length.out=10*N) # solve for too long
	est_grad_augment <- function(t,v,parms){return(est_grad(v))}
	estimated_path <- deSolve::lsoda(initial_condition, time_vec.est, est_grad_augment)
	est_transient.cutoff <- get_lc_cutoff(estimated_path)
	est_transient <- estimated_path[1:est_transient.cutoff,]
	if(is.null(nrow(est_transient))){
		est_transient <- matrix(c(est_transient),nrow=1,dimnames=list(NULL,names(est_transient)))
		est_transient <- cbind(est_transient,source = 2)
	} else {est_transient <- cbind(est_transient,source = 2)}
	est_transient_est_grad <- matrix(unlist(apply(est_transient,1,function(row){est_grad(c(x=row[2],y=row[3]))})),ncol=2,byrow=TRUE)
	
	# Step 2: Calculate mean nearest-neighbor distance
	computed_losses <- c()
	if ("squared_model" %in% losses){
		distances <- (RANN::nn2(matrix(est_transient[,c("x","y")],ncol=2),
										matrix(true_transient[,c("x","y")],ncol=2),
										k=1)$nn.dists)^2
		transient_mse <- apply(distances,2,mean)
		computed_losses <- c(computed_losses, c("squared_model" = transient_mse))
	}
	if ("squared_estimator" %in% losses){
		distances <- (RANN::nn2(matrix(true_transient[,c("x","y")],ncol=2),
										matrix(est_transient[,c("x","y")],ncol=2),
										k=1)$nn.dists)^2
		transient_mse <- apply(distances,2,mean)
		computed_losses <- c(computed_losses, c("squared_estimator" = transient_mse))
	}
	if ("trig_model" %in% losses){
		true_transient_est_grad <- matrix(unlist(apply(true_transient,1,function(row){est_grad(c(x=row[2],y=row[3]))})),ncol=2,byrow=TRUE)
		true_transient_grad_mat <- cbind(true_transient_true_grad,true_transient_est_grad)
		distances <- apply(true_transient_grad_mat,1,function(row){
			sin(acos((row[1]*row[3] + row[2]*row[4]) / (sqrt(row[1]^2 + row[2]^2)*sqrt(row[3]^2 + row[4]^2)) )/2 )
		})
		transient_mse <- mean(distances, na.rm = TRUE)
		computed_losses <- c(computed_losses, c("trig_model" = transient_mse))
	}
	if ("trig_estimator" %in% losses){
		est_transient_true_grad <- matrix(unlist(apply(est_transient,1,function(row){true_grad(c(x=row[2],y=row[3]))})),ncol=2,byrow=TRUE)
		est_transient_grad_mat <- cbind(est_transient_true_grad,est_transient_est_grad)
		distances <- apply(est_transient_grad_mat,1,function(row){
			sin(acos((row[1]*row[3] + row[2]*row[4]) / (sqrt(row[1]^2 + row[2]^2)*sqrt(row[3]^2 + row[4]^2)) )/2 )
		})
		transient_mse <- mean(distances, na.rm = TRUE)
		computed_losses <- c(computed_losses, c("trig_estimator" = transient_mse))
	}

	# Step 3: Diagnostics
	if (plot){
		full_data <- rbind(true_transient,est_transient)
		transient_plot <- ggplot(as_tibble(full_data), aes(x=x, y=y)) +
			geom_point(aes(color=factor(source))) +
			geom_path(aes(color=factor(source))) +
			geom_point(x = initial_condition[1], y = initial_condition[2] , color = "red", size = 3) + 
			labs(title="Transient")
		print(transient_plot)
	}
	
	if (return_transients){
		true_transient <- cbind(true_transient, true_transient_true_grad)
		est_transient <- cbind(est_transient, est_transient_est_grad)
		output_list <- list(losses = computed_losses,
							true_transient = true_transient,
							est_transient = est_transient
							)
		return(output_list)
	}
	return(computed_losses)
}

sample_ic <- function(data, counts, radial_vec, sigma, filter_fn = NULL, plot = FALSE){
	# verify length of radii matches counts
	stopifnot(length(counts)==length(radial_vec))

	position_data <- data[,c("x","y")]
	center <- apply(position_data,2,mean)
	position_data.centered <- sweep(position_data,2,center,FUN="-")
	
	sampled_points <- position_data.centered[sample(1:nrow(position_data.centered),sum(counts),replace=TRUE),]
	radii <- rep(radial_vec,counts)
	# add some Gaussian noise to the radius if not equal to 1
	radii[-which(radii==1)] <- radii[-which(radii==1)] + rnorm(length(radii[-which(radii==1)]),sd=sigma)
	initial_conditions <- sweep(sampled_points*radii,2,center,FUN="+")
	
	# apply condition filtering if applicable
	if(!is.null(filter_fn)){
		filter_mask <-apply(initial_conditions,1,filter_fn)
		orig_n <- nrow(initial_conditions)
		initial_conditions <- initial_conditions[filter_mask,]
		if(nrow(initial_conditions) != orig_n){
			warning(paste0(orig_n - nrow(initial_conditions), " initial conditions filtered. ",
						   nrow(initial_conditions), " remaining."))
		}
	}

	if (plot){
		ic_plot <- ggplot(position_data, aes(x=x, y=y)) +
			geom_point(color = "#FDB515",alpha=0.15) +
			geom_point(data = initial_conditions, aes(x=x,y=y),color = "#ADC900") +
			labs(title="Initial Conditions")
		print(ic_plot)
	}

	initial_conditions <- apply(initial_conditions,2,as.numeric) # play nice with deSolve
	return(initial_conditions)
}

get_lc_cutoff <- function(full_samples,epsilon=NULL,gradient=FALSE){
	# this naive algorithm delineates the limit cycle from transient
	# it places a nearest neighbor ball around the final sample and starts
	#   the signal at the earliest observation within that ball 
	x_var <- ifelse(gradient,"f_x","x")
	y_var <- ifelse(gradient,"f_y","y")
	
	final_observation <- full_samples[nrow(full_samples),c(x_var,y_var)]
	if(is.null(nrow(final_observation))){final_observation <- matrix(final_observation,ncol=length(final_observation))}
	other_observations <- full_samples[-nrow(full_samples),c(x_var,y_var)]
	if(is.null(nrow(other_observations))){other_observations <- matrix(other_observations,ncol=length(other_observations))}
	
	cutoff_index <- min(RANN::nn2(other_observations,final_observation,
								  k=nrow(other_observations))$nn.idx[1:(nrow(other_observations)/50)])
	
	if (cutoff_index == nrow(full_samples)){warning("No transient cutoff found. This system may be unstable.")}
	return(cutoff_index)
}

################
# Visualizaton #
################

plot_field_vectors <- function(samples, grad_fn,
					   x_grid_count = 24, y_grid_count = 24, frame_frac = 0.2,
					   title="", rescale = 1, sample_alpha=1){
	# plot a gradient field with limit cycle samples overlayed

	x_range = max(samples$x) - min(samples$x)
	x_grid <- seq(min(samples$x) - frame_frac*x_range,
				  max(samples$x) + frame_frac*x_range, length.out =  x_grid_count)
	y_range = max(samples$y) - min(samples$y)
	y_grid <- seq(min(samples$y) - frame_frac*y_range,
				  max(samples$y) + frame_frac*y_range, length.out =  y_grid_count)
	eval_grid <- unname(as.matrix(expand.grid(x_grid,y_grid))) # convert to N x 2 matrix'
	
	preds <- apply(eval_grid,1,grad_fn)
	preds <- do.call(rbind,lapply(preds,unlist))
	preds_tibble <- as_tibble(preds)
	colnames(preds_tibble) <- c("u","v")
	
	eval_tibble <- as_tibble(eval_grid)
	colnames(eval_tibble) <- c("x","y")
	plot_tibble <- cbind(eval_tibble, preds_tibble)
	plot_tibble <- plot_tibble %>% mutate(u = rescale*u, v = rescale*v)
	colnames(samples) <- c("x","y","u","v")
	
	field_plot <- ggplot(plot_tibble, aes(x = x, y = y, u = u, v = v)) +
		geom_quiver(color = "#003262",vecsize=rescale) +
		geom_point(data = samples, aes(x = x, y = y), color = "#FDB515", alpha = sample_alpha) +
		labs(title = title) + 
		xlim(min(plot_tibble$x), max(plot_tibble$x)) + 
		ylim(min(plot_tibble$y), max(plot_tibble$y))
	
	return(field_plot)
}

plot_field_heatmap <- function(samples, grad_fn,
					   x_grid_count = 24, y_grid_count = 24, frame_frac = 0.2,
					   title="", rescale = 1, sample_alpha=1){
	# plot a gradient field with limit cycle samples overlayed
	
	x_range = max(samples$x) - min(samples$x)
	x_grid <- seq(min(samples$x) - frame_frac*x_range,
				  max(samples$x) + frame_frac*x_range, length.out =  x_grid_count)
	y_range = max(samples$y) - min(samples$y)
	y_grid <- seq(min(samples$y) - frame_frac*y_range,
				  max(samples$y) + frame_frac*y_range, length.out =  y_grid_count)
	eval_grid <- unname(as.matrix(expand.grid(x_grid,y_grid))) # convert to N x 2 matrix'
	
	preds <- apply(eval_grid,1,grad_fn)
	preds <- do.call(rbind,lapply(preds,unlist))
	preds_tibble <- as_tibble(preds)
	colnames(preds_tibble) <- c("f_x","f_y")
	
	eval_tibble <- as_tibble(eval_grid)
	colnames(eval_tibble) <- c("x","y")
	plot_tibble <- cbind(eval_tibble, preds_tibble)
	plot_tibble <- plot_tibble %>% mutate(f_x = rescale*f_x, f_y = rescale*f_y)
	plot_tibble <- plot_tibble %>% 
		pivot_longer(cols=c("f_x","f_y"),names_to="axis",values_to="magnitude") %>%
		mutate(sign = ifelse(magnitude >= 0,"positive","negative")) %>%
		mutate(magnitude = abs(magnitude))

	colnames(samples) <- c("x","y","u","v")

	field_plot <- ggplot(plot_tibble,aes(x=x,y=y)) +
		geom_tile(aes(fill=sign,alpha=magnitude),color="white") +
		scale_fill_manual(values=c(positive="salmon", negative="steelblue")) +
		geom_point(data = samples, aes(x = x, y = y), color = "#FDB515", alpha = sample_alpha) +
		facet_wrap(~axis) +
		labs(title = title)
	return(field_plot)
}