library(patchwork)
library(gridExtra)
source("pipeline.R")
source("plot_eigen_decomp.R")

build_lc_scorecard <- function(data_object,initial_conditions,t_star,true_grad_index,estimators,losses,N=1000,title_str=""){
	loss_tibble_list <- list()
	path_tibble_list <- list()
	global_index <- 1
	for (ic_index in 1:nrow(initial_conditions)){
		new_ic = TRUE # to avoid duplicate plotting of the ground truth
		for (estimator_index in estimators){
			outcome_list <- calc_lc_error(initial_conditions[ic_index,],data_object$estimators[[true_grad_index]]$fn,
												 data_object$estimators[[estimator_index]]$fn,
												 t_star, N, losses, return_transients = TRUE)
			# CREATE LOSS LONG 
			l <- length(losses)
			loss_row <- tibble(IC_x=rep(initial_conditions[ic_index,1],l),IC_y=rep(initial_conditions[ic_index,2],l),
							   Estimator=factor(data_object$estimators[[estimator_index]]$estimator$method),
							   Metric=names(outcome_list$losses),Loss=unname(outcome_list$losses))
			loss_tibble_list[[global_index]] <- loss_row
			
			# CREATE PATHS LONG
			path_tibble <- as_tibble(outcome_list$est_lc) %>% 
				select(-source) %>%
				mutate(Estimator = factor(data_object$estimators[[estimator_index]]$estimator$method),IC=factor(ic_index))
			print(paste0("Est LC n:",nrow(path_tibble)))
			print(head(path_tibble))
			colnames(path_tibble) <- c("time","x","y","f_x","f_y","Estimator","IC")
			if (new_ic){
				true_tibble <- as_tibble(outcome_list$true_lc) %>% 
					select(-source) %>%
					mutate(Estimator = factor("Truth"),IC=factor(ic_index))
				colnames(true_tibble) <- c("time","x","y","f_x","f_y","Estimator","IC")
				path_tibble <- rbind(path_tibble, true_tibble)
			}
			new_ic = FALSE 
			
			global_index <- global_index + 1
			path_tibble_list[[global_index]] <- path_tibble
		}
	}
	
	loss_tibble_long <- dplyr::bind_rows(loss_tibble_list)
	loss_tibble_wide <- pivot_wider(loss_tibble_long,
									names_from=c(Estimator,Metric),names_sep="*",values_from=Loss) %>% 
		mutate(across(where(is.numeric), round, 3))
	path_tibble_long <- dplyr::bind_rows(path_tibble_list)
	
	transient_plot <- ggplot(path_tibble_long, aes(x=x, y=y)) +
		geom_point(aes(shape=Estimator,color=IC)) +
		geom_path(aes(linetype=Estimator,color=IC)) +
		geom_point(data = as_tibble(initial_conditions), aes(x=x, y=y), color = "red", size = 3) + 
		geom_point(data = as_tibble(data_object$smooth_tail), aes(x = x, y = y), color = "#FDB515", alpha = 0.3) +
		labs(title=title_str) +
		facet_wrap(~Estimator)
	
	if (TRUE){
		mytheme <- gridExtra::ttheme_default(
			colhead = list(fg_params=list(cex = 0.5)))
		output_table <- gridExtra::tableGrob(loss_tibble_wide, theme = mytheme)
		# Set widths/heights to 'fill whatever space I have'
		output_table$widths <- unit(rep(1, ncol(output_table)), "null")
		output_table$heights <- unit(rep(1, nrow(output_table)), "null")
		
		# Format table as plot
		p3 <- ggplot() +
			annotation_custom(output_table)
	}
	scorecard_plot <- transient_plot + p3 + plot_layout(ncol = 1)
	print(scorecard_plot)
	return(list(tibble = loss_tibble_long, plot = scorecard_plot))
}

set.seed(22)
estimator_list <- list(
	list(method="truth"),
	list(method="gd_spline",side_info=list(),norder=6,nbasis=36,penalty_order=2,lambda=1e-2,
		 gd_params = list(algorithm="Vanilla",eig_rank=1,eta=0.1,dt_radius=0.2,
		 				 batching=list(num_iterations=4,batches_per=25,batch_size=1,skip_negative=T))),
	list(method="vanilla_spline",side_info=list(),norder=6,nbasis=36,penalty_order=2,lambda=1e-2)
)

vdp_trial <- generate_experiment_trial("example_vdp", system_name = "van_der_pol", system_parameters = list(mu=0.5),
									   random_seed = 31, num_samples = 1500, delta_t = 0.05,
									   replicates = 1, var_x = 0.05, var_y = 0.05)
vdp_trial_smooth <- smooth_trial(vdp_trial)
vdp_trial_fit <- apply_estimators(vdp_trial_smooth,estimator_list)

initial_conditions <-  sample_ic(vdp_trial_fit$replicates[[1]]$smooth_tail, c(1,1), c(0.9,1.1), 1)
losses <- c("squared_model","squared_estimator")
t_star <- get_lc_cutoff(vdp_trial_fit$replicates[[1]]$smooth)

# Compute
result <- build_lc_scorecard(vdp_trial_fit$replicates[[1]],initial_conditions,
	t_star,1,c(2,3),losses, N = 6000, title_str = "VdP with mu = 3")

# Compare Eigenvectors
vdp_limits <- list(x_lower = -0.1, x_upper = 0.1, y_lower = -1, y_upper = -0.95)
possible_data <- vdp_trial_fit$replicates[[1]]$smooth_tail
plot_data <- possible_data %>% filter(x > vdp_limits$x_lower & 
	x < vdp_limits$x_upper & 
	y > vdp_limits$y_lower & 
	y < vdp_limits$y_upper)
eigen_plot_vanilla <- plot_eigen_decomp(vdp_trial_fit$replicates[[1]]$estimators[[2]],
	vdp_limits, vector_density = 10, sample_df = plot_data, title = "Vanilla")
eigen_plot_gd <- plot_eigen_decomp(vdp_trial_fit$replicates[[1]]$estimators[[3]],
	vdp_limits, vector_density = 10, sample_df = plot_data, title = "Penalized")
eigen_combined <- eigen_plot_vanilla + eigen_plot_gd + plot_layout(ncol = 2)
eigen_combined
eigen_delta <- plot_eigen_decomp_delta(vdp_trial_fit$replicates[[1]]$estimators[[2]],
	vdp_trial_fit$replicates[[1]]$estimators[[3]],
	vdp_limits, vector_density = 10, sample_df = plot_data, title = "Penalized - Truth")
eigen_delta