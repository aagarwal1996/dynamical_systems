# Summarize Fit
source("pipeline.R")

get_estimated_transient <- function(initial_condition,est_grad,t_star,N){
	time_vec.est <- seq(0,3*t_star,length.out=3*N) # solve for too long
	est_grad_augment <- function(t,v,parms){return(est_grad(v))}
	estimated_path <- deSolve::lsoda(initial_condition, time_vec.est, est_grad_augment)
	est_transient <- estimated_path
	est_transient_est_grad <- matrix(unlist(apply(est_transient,1,function(row){est_grad(c(x=row[2],y=row[3]))})),ncol=2,byrow=TRUE)
	
	est_transient <- cbind(est_transient, est_transient_est_grad)
	return(est_transient)
}

plot_single_knee <- function(data_object, initial_conditions, replicate_idx = 1){
	path_tibble_list <- list()
	selected_replicate <- data_object$replicates[[replicate_idx]]
	smoothed_samples <- selected_replicate$smooth
	t_star <- get_lc_cutoff(smoothed_samples)
	
	global_index <- 1
	for (ic_index in 1:nrow(initial_conditions)){
		for (estimator_index in 1:length(selected_replicate$estimators)){
			est_transient <- get_estimated_transient(initial_conditions[ic_index,],
													 selected_replicate$estimators[[estimator_index]]$fn,
												 t_star, 300)

			path_tibble <- as_tibble(est_transient) %>% 
				mutate(Estimator = factor(selected_replicate$estimators[[estimator_index]]$estimator$method),IC=factor(ic_index))
			colnames(path_tibble) <- c("time","x","y","f_x","f_y","Estimator","IC")

			path_tibble_list[[global_index]] <- path_tibble
			global_index <- global_index + 1
			
		}
	}
	path_tibble_long <- dplyr::bind_rows(path_tibble_list)
	
	transient_plot <- ggplot(path_tibble_long, aes(x=x, y=y)) +
		geom_point(aes(shape=Estimator,color=IC)) +
		geom_path(aes(linetype=Estimator,color=IC)) +
		geom_point(data = as_tibble(initial_conditions), aes(x=x, y=y), color = "red", size = 3) + 
		geom_point(data = as_tibble(smoothed_samples), aes(x = x, y = y), color = "#FDB515", alpha = 0.3) +
		facet_wrap(~Estimator)
	
	print(transient_plot)
}

# Run 

# Z = angle with respect to forward and back extension/flexion (pitch)
# Y = angle with respect to rotation (rudder)
# X = angle with respect to vertical pivot (roll)

younger_tibble <- read_table("gait/gait_data/WBDSascii/WBDS20walkT08ang.txt") %>% mutate(unit="younger")
older_tibble <- read_table("gait/gait_data/WBDSascii/WBDS41walkT08ang.txt") %>% mutate(unit="older")
full_tibble <- rbind(younger_tibble, older_tibble) %>% mutate(unit = as.factor(unit))

x_covar = 'LKneeAngleZ'
y_covar = 'RKneeAngleZ'

younger_samples <- younger_tibble %>% select(x = !!x_covar, y = !!y_covar) %>% mutate(f_x = 0, f_y = 0)
younger_object <- generate_experiment_trial(experiment_name="younger_cohort",replicates=0)
younger_object$replicates <- list(list(samples=younger_samples))
younger_object <- smooth_trial(younger_object, impute_grad=TRUE)

initial_conditions <-  sample_ic(younger_object$replicates[[1]]$smooth, c(1,1,1), c(.95,1,1.05), 1)
estimator_list <- list(
	list(method="loess",bandwidth_scalar=1),
	list(method="vanilla_spline",side_info=list(),norder=6,nbasis=36,penalty_order=2,lambda=1e-2),
	list(method="gd_spline",side_info=list(),norder=6,nbasis=36,penalty_order=2,lambda=1e-2,
		 gd_params = list(algorithm="Projection_Motion",eig_rank=1,eta=1,dt_radius=5,
		 				 batching=list(num_iterations=2,batches_per=100,batch_size=1,skip_negative=T)))
)
younger_object <- apply_estimators(younger_object,estimator_list)
plot_single_knee(younger_object,initial_conditions)
plot_field_heatmap(younger_object$replicates[[1]]$smooth,
				   younger_object$replicates[[1]]$estimators[[2]]$fn)
plot_field_heatmap(younger_object$replicates[[1]]$smooth,
				   younger_object$replicates[[1]]$estimators[[3]]$fn)


