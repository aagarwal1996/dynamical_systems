# this file contains functions to help compare different loss functions to use in transient fitting
library(patchwork)
library(gridExtra)
source("pipeline.R")

build_scorecard <- function(data_object,initial_conditions,t_star,true_grad,estimators,losses,N=1000,title_str=""){
	loss_tibble_list <- list()
	path_tibble_list <- list()
	
	global_index <- 1
	for (ic_index in 1:nrow(initial_conditions)){
		new_ic = TRUE # to avoid duplicate plotting of the ground truth
		for (estimator_index in estimators){
			outcome_list <- calc_transient_error(initial_conditions[ic_index,],data_object$estimators[[true_grad]]$fn,
												 data_object$estimators[[estimator_index]]$fn,
												 t_star, N, losses, return_transients = TRUE)
			# CREATE LOSS LONG 
			l <- length(losses)
			loss_row <- tibble(IC_x=rep(initial_conditions[ic_index,1],l),IC_y=rep(initial_conditions[ic_index,2],l),
				   Estimator=factor(data_object$estimators[[estimator_index]]$estimator$method),
				   Metric=names(outcome_list$losses),Loss=unname(outcome_list$losses))
			loss_tibble_list[[global_index]] <- loss_row
			
			# CREATE PATHS LONG
			
			path_tibble <- as_tibble(outcome_list$est_transient) %>% 
				select(-source) %>%
				mutate(Estimator = factor(data_object$estimators[[estimator_index]]$estimator$method),IC=factor(ic_index))
			colnames(path_tibble) <- c("time","x","y","f_x","f_y","Estimator","IC")
			if (new_ic){
				true_tibble <- as_tibble(outcome_list$true_transient) %>% 
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
	 	geom_point(data = as_tibble(data_object$smooth), aes(x = x, y = y), color = "#FDB515", alpha = 0.3) +
	 	labs(title=title_str)
	
	mytheme <- gridExtra::ttheme_default(
		colhead = list(fg_params=list(cex = 0.5)))
	output_table <- gridExtra::tableGrob(loss_tibble_wide, theme = mytheme)
	# Set widths/heights to 'fill whatever space I have'
	output_table$widths <- unit(rep(1, ncol(output_table)), "null")
	output_table$heights <- unit(rep(1, nrow(output_table)), "null")
	
	# Format table as plot
	p3 <- ggplot() +
		annotation_custom(output_table)
	scorecard_plot <- transient_plot + p3 + plot_layout(ncol = 1)
	print(scorecard_plot)

	return(list(tibble = loss_tibble_long, plot = scorecard_plot))
}

run_replicated_score_card <- function(data_object,true_grad,estimators,losses,ic_counts,ic_radii,ic_noise,N=1000,title_str=""){
	run_outcome_list <- list()
	for (i in 1:length(data_object$replicates)){
		initial_conditions <-  sample_ic(data_object$replicates[[i]]$smooth_tail, ic_counts, ic_radii, ic_noise)
		t_star <- get_lc_cutoff(data_object$replicates[[i]]$smooth)
		run_result <- build_scorecard(data_object$replicates[[i]],initial_conditions,t_star,true_grad,estimators,losses, N = N, title_str = title_str)
		run_outcome_list <- append(run_outcome_list,run_result$tibble)
	}
	loss_tibble_long <- dplyr::bind_rows(loss_tibble_list)
	return(loss_tibble_long)
}

##### TESTER ##### 

estimator_list <- list(
	list(method="truth"),
	list(method="gd_spline",side_info=list(),norder=6,nbasis=36,penalty_order=2,lambda=1e-2,
		 gd_params = list(algorithm="Vanilla",eig_rank=1,eta=0.02,dt_radius=0.2,
		 				 batching=list(num_iterations=1,batches_per=25,batch_size=1,skip_negative=T))),
	list(method="vanilla_spline",side_info=list(),norder=6,nbasis=36,penalty_order=2,lambda=1e-2)
)

vdp_trial <- generate_experiment_trial("example_vdp", system_name = "van_der_pol", system_parameters = list(mu=10),
									   random_seed = 31, num_samples = 1500, delta_t = 0.05,
									   replicates = 1, var_x = 0.05, var_y = 0.05)
vdp_trial_smooth <- smooth_trial(vdp_trial)
vdp_trial_fit <- apply_estimators(vdp_trial_smooth,estimator_list)

initial_conditions <-  sample_ic(vdp_trial_fit$replicates[[1]]$smooth_tail, c(1,1,1), c(0.5,1,1.5), 1)
losses <- c("squared_model","squared_estimator", "trig_model", "trig_estimator")
t_star <- get_lc_cutoff(vdp_trial_fit$replicates[[1]]$smooth)
build_scorecard(vdp_trial_fit$replicates[[1]],initial_conditions,t_star,1,c(2,3),losses, N = 1000, title_str = "Scorecard for VdP mu = 10")

run_replicated_score_card(vdp_trial_fit,1,c(2,3),losses,c(1,1,1),c(0.5,1,1.5),1)