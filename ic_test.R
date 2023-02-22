source("pipeline.R")

vdp_trial <- generate_experiment_trial("example_vdp", system_name = "van_der_pol", system_parameters = list(mu=1.5),
									   random_seed = 2, num_samples = 1000, delta_t = 0.1,
									   replicates = 1, var_x = 0, var_y = 0)
vdp_trial_smooth <- smooth_trial(vdp_trial)

some_data <- vdp_trial_smooth$replicates[[1]]$smooth_tail
first_counts <- c(50,50,50,50,50)
first_rads <- c(0.25,0.75,1,1.25,1.5)

fn <- function(v){v[1] > 0 && v[2] > 0}
sample_ic(some_data, first_counts, first_rads, .1, fn, plot = TRUE)

estimator_list <- list(
 	list(method="truth"),
 	list(method="gd_spline",side_info=list(),norder=6,nbasis=36,penalty_order=2,lambda=1e-2,
 		 gd_params = list(algorithm="Vanilla",eig_rank=1,eta=0.02,dt_radius=0.2,
 		 				 batching=list(num_iterations=1,batches_per=25,batch_size=1,skip_negative=T))),
 	list(method="vanilla_spline",side_info=list(),norder=6,nbasis=36,penalty_order=2,lambda=1e-2),
 	list(method="nw",bandwidth_matrix = 0.25*diag(2))
)


vdp_trial <- generate_experiment_trial("example_vdp", system_name = "van_der_pol", system_parameters = list(mu=5),
									   random_seed = 31, num_samples = 5000, delta_t = 0.01,
									   replicates = 1, var_x = 0.05, var_y = 0.05)
vdp_trial_smooth <- smooth_trial(vdp_trial)
vdp_trial_fit <- apply_estimators(vdp_trial_smooth,estimator_list)
get_dt_mse(vdp_trial_fit$replicates[[1]],plot=TRUE)

plot_field_heatmap(vdp_trial_fit$replicates[[1]]$smooth_tail, vdp_trial_fit$replicates[[1]]$estimators[[2]]$fn,
		   rescale = 1, sample_alpha = 0.3)
plot_field_heatmap(vdp_trial_fit$replicates[[1]]$smooth_tail, vdp_trial_fit$replicates[[1]]$estimators[[3]]$fn,
				   rescale = 1, sample_alpha = 0.3)
plot_field_heatmap(vdp_trial_fit$replicates[[1]]$smooth_tail, vdp_trial_fit$replicates[[1]]$estimators[[4]]$fn,
				   rescale = 1, sample_alpha = 0.3)

t_star <- get_lc_cutoff(vdp_trial_fit$replicates[[1]]$smooth)

ic <- sample_ic(vdp_trial_fit$replicates[[1]]$smooth_tail, c(1), c(1.3), 0)
calc_transient_error(ic, vdp_trial_fit$replicates[[1]]$estimators[[1]]$fn, 
					vdp_trial_fit$replicates[[1]]$estimators[[2]]$fn, t_star, 1000, plot = TRUE,int_pred=TRUE)
calc_transient_error(ic, vdp_trial_fit$replicates[[1]]$estimators[[1]]$fn, 
					 vdp_trial_fit$replicates[[1]]$estimators[[3]]$fn, t_star, 1000, plot = TRUE)
calc_transient_error(ic, vdp_trial_fit$replicates[[1]]$estimators[[1]]$fn, 
					 vdp_trial_fit$replicates[[1]]$estimators[[4]]$fn, t_star, 1000, plot = TRUE,int_pred=TRUE)

# start.time <- proc.time()[3]
# j <- apply_transient_error(vdp_trial_fit,c(1,1),c(.9,1.1),0,NULL,100)
# end.time <- proc.time()[3]
# delta.time <- end.time - start.time
# paste0("Runtime: ", delta.time)