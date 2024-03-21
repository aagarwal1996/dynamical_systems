source(here::here('pipeline.R'))

# Define Functions

generate_trajectory_data <- function(ic_matrix, fit_estimators, t_star = 10, N = 300){
	return_list <- list()
	for (ic_index in seq_len(nrow(ic_matrix))){
		path_list <- list()
		for (estimator_index in seq_along(fit_estimators)){
			time_vec.est <- seq(0,5*t_star,length.out=5*N)
			est_grad_augment <- function(t,v,parms){return(fit_estimators[[estimator_index]]$fn(v))}
			estimated_path <- deSolve::lsoda(ic_matrix[ic_index,], time_vec.est, est_grad_augment)
			est_transient.cutoff <- get_lc_cutoff(estimated_path)
			est_lc <- estimated_path[est_transient.cutoff:nrow(estimated_path),]
			path_list[[estimator_index]] <- est_lc
		}
		return_list[[ic_index]] <- list(ic = ic_matrix[ic_index,], paths = path_list)
	}
	
	return(return_list)
}

eval_lc_fit <- function(ic_trajectories, fit_estimators){
	labeled_paths <- lapply(seq_along(ic_trajectories$paths), function(i){
		path_mat <- ic_trajectories$paths[[i]]
		estimator <- rep(i, nrow(path_mat))
		path_mat <- cbind(path_mat, estimator)
		return(path_mat) })
	plot_object <- as.data.frame(do.call(rbind,labeled_paths))
	plot_object$estimator <- as.factor(plot_object$estimator)
	
	x_plot <- ggplot2::ggplot(plot_object, aes(x=time,y=x,color=estimator)) +
		ggplot2::geom_point()
	print(x_plot)
	y_plot <- ggplot2::ggplot(plot_object, aes(x=time,y=y,color=estimator)) +
		ggplot2::geom_point()
	print(y_plot)
	summary_plot <- ggplot2::ggplot(plot_object, aes(x=x,y=y,color=estimator)) +
		ggplot2::geom_point()
	print(summary_plot)
	return()
}


# Run Example

vdp_trial <- generate_experiment_trial("example_vdp", system_name = "van_der_pol", system_parameters = list(mu=1.5),
									   random_seed = 2, num_samples = 1000, delta_t = 0.1,
									   replicates = 1, var_x = 0, var_y = 0)
vdp_trial_smooth <- smooth_trial(vdp_trial)

estimator_list <- list(
	list(method="truth"),
	list(method="vanilla_spline",side_info=list(),norder=6,nbasis=36,penalty_order=2,lambda=1e-2),
	list(method="nw",bandwidth_matrix = 0.25*diag(2))
)

vdp_trial_fit <- apply_estimators(vdp_trial_smooth,estimator_list)

ic_counts <- c(0,0,0,1,1)
ic_rads <- c(0.25,0.75,1,1.25,1.5)
sampled_ics <- sample_ic(vdp_trial_smooth$replicates[[1]]$smooth_tail, ic_counts, ic_rads, .1, plot = TRUE)

traj_data <- generate_trajectory_data(sampled_ics, vdp_trial_fit$replicates[[1]]$estimators)
eval_lc_fit(traj_data[[1]], vdp_trial_fit$replicates[[1]]$estimators)