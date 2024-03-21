#source(here::here('..','dynamical_systems_simulations','run_simulation_redux.R'))
source(here::here('pipeline.R'))
source(here::here('transient_metrics.R'))

estimator_list <- list(
	list(method="loess",bandwidth_scalar=1),
	list(method="gd_spline",side_info=list(),norder=6,nbasis=36,penalty_order=2,lambda=1e-2,
		 gd_params = list(algorithm="Vanilla",eig_rank=1,eta=0.02,dt_radius=0.2,
		 				 batching=list(num_iterations=1,batches_per=25,batch_size=1,skip_negative=T))),
	list(method="vanilla_spline",side_info=list(),norder=6,nbasis=36,penalty_order=2,lambda=1e-2)
)


# Z = angle with respect to forward and back extension/flexion (pitch)
# Y = angle with respect to rotation (rudder)
# X = angle with respect to vertical pivot (roll)

younger_tibble <- read_table("gait/gait_data/WBDSascii/WBDS21walkT08ang.txt") %>% mutate(unit="younger")
older_tibble <- read_table("gait/gait_data/WBDSascii/WBDS42walkT08ang.txt") %>% mutate(unit="older")
full_tibble <- rbind(younger_tibble, older_tibble) %>% mutate(unit = as.factor(unit))

x_covar = 'RHipAngleZ'
y_covar = 'LKneeAngleZ'

younger_samples <- younger_tibble %>% select(x = !!x_covar, y = !!y_covar) %>% mutate(f_x = 0, f_y = 0)
younger_object <- generate_experiment_trial(experiment_name="younger_cohort",replicates=0)
younger_object$replicates <- list(list(samples=younger_samples))
younger_object <- smooth_trial(younger_object, impute_grad=TRUE)
younger_object <- apply_estimators(younger_object,estimator_list)

initial_conditions <-  sample_ic(younger_object$replicates[[1]]$smooth, c(2,0,0), c(0.9,1,1.1), 1)
losses <- c("squared_model","squared_estimator", "trig_model", "trig_estimator")
t_star <- get_lc_cutoff(younger_object$replicates[[1]]$smooth)
build_scorecard(younger_object$replicates[[1]],initial_conditions,t_star,1,c(1,2,3),losses, N = 1000, title_str = "Scorecard for Younger Knee Cohort")

older_samples <- older_tibble %>% select(x = !!x_covar, y = !!y_covar) %>% mutate(f_x = 0, f_y = 0)
older_object <- generate_experiment_trial(experiment_name="older_cohort",replicates=0)
older_object$replicates <- list(list(samples=older_samples))
older_object <- smooth_trial(older_object, impute_grad=TRUE)
older_object <- apply_estimators(older_object,estimator_list)

initial_conditions <-  sample_ic(older_object$replicates[[1]]$smooth, c(1,1,1), c(0.5,1,1.5), 1)
losses <- c("squared_model","squared_estimator", "trig_model", "trig_estimator")
t_star <- get_lc_cutoff(older_object$replicates[[1]]$smooth)
build_scorecard(older_object$replicates[[1]],initial_conditions,t_star,1,c(1,2),losses, N = 1000, title_str = "Scorecard for VdP mu = 10")


# 
# ## EDA
# 
# ggplot(full_tibble, aes(x = LHipAngleZ, y= LKneeAngleZ)) +
#   geom_point() +
#   facet_wrap(~unit)
# 
# # https://stackoverflow.com/questions/9282726/r-how-can-i-use-a-variable-to-define-scale-colour-manual-and-key-label-colours
# #legend_palette <-  scale_colour_manual("Legend", values = c(x_covar = "red", y_covar = "blue"))
# #names(legend_palette$values) <- c(x_covar,y_covar)
# 
# full_tibble %>%
#    ggplot() +
#    geom_point(aes(x=Time, y = LHipAngleZ), color = "blue") +
#    geom_point(aes(x=Time, y = LKneeAngleZ), color = "red")  +
#    labs(title="Angle over Time",x="Time",y="Angle") +
#    facet_wrap(~unit) +
#    labs(caption=paste0("Blue = ",x_covar,"; Red = ",y_covar))
# 
# # Spline lambda
# 
# ## Gradient Field Fitting
# spline <- list(method = "spline",  params = list(lambda = 1e-2, norder = 4, nbasis = 12, side_info = list()))
# spline_box <- list(method = "spline",  params = list(lambda = 1e-2, norder = 4, nbasis = 12, side_info = list(list(name="boundary_box"))))
# spline_dt <- list(method = "spline",  params = list(lambda = 1e-2, norder = 4, nbasis = 12, side_info = list(list(name="delta_ring",
#                                                                       sample_frac = .75, radius=0.01, shift = 0, si_magnitude = .05, lc_params = NA))))
# spline_dt_copy <- spline_dt
# younger_experiment <-  list(name = "younger", system = "knee", limit_cycle_samples = younger_smooth,
#                             lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 4)
# older_experiment <- list(name = "older", system = "knee", limit_cycle_samples = older_smooth,
#                          lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 4)
# 
# experiment_data <- generate_data_object_obs(list(younger_experiment, older_experiment))
# experiment_estimators <- list(spline_box, spline_dt)
# experiment_results_knee <- evaluate_gradient_methods(experiment_data, experiment_estimators)
# 
# younger_fields <- list(type = "field", experiments = list(c(data = 1, estimator = 1)))#, c(data = 1, estimator = 2), c(data = 1, estimator = 3)))
# younger_paths <- list(type = "solution_path", experiments = list(c(data = 1, estimator = 1),c(data = 1, estimator = 2)))#, c(data = 1, estimator = 2), c(data = 1, estimator = 3)))
# 
# older_fields <- list(type = "field", experiments = list(c(data = 2, estimator = 1)))#, c(data = 2, estimator = 2), c(data = 2, estimator = 3)))
# older_paths <- list(type = "solution_path", experiments = list(c(data = 2, estimator = 1),c(data = 2, estimator = 2)))#, c(data = 2, estimator = 2), c(data = 2, estimator = 3)))
# 
# #visualize_results(experiment_results_knee, list(younger_fields, older_fields))
# visualize_results(experiment_results_knee, list(younger_paths))
# visualize_results(experiment_results_knee, list(older_paths))



# spline_16_box <- list(method = "spline",  params = list(lambda = 1e-16, side_info = list(list(name="boundary_box"))))
# spline_16 <- list(method = "spline",  params = list(lambda = 1e-16, side_info = list()))
# spline_8_box <- list(method = "spline",  params = list(lambda = 1e-8, side_info = list(list(name="boundary_box"))))
# spline_8 <- list(method = "spline",  params = list(lambda = 1e-8, side_info = list()))
# spline_4_box <- list(method = "spline",  params = list(lambda = 1e-4, side_info = list(list(name="boundary_box"))))
# spline_4 <- list(method = "spline",  params = list(lambda = 1e-4, side_info = list()))
# spline_2_box <- list(method = "spline",  params = list(lambda = 1e-2, side_info = list(list(name="boundary_box"))))
# spline_2 <- list(method = "spline",  params = list(lambda = 1e-2, side_info = list()))
# spline_1_box <- list(method = "spline",  params = list(lambda = 1e-1, side_info = list(list(name="boundary_box"))))
# spline_1 <- list(method = "spline",  params = list(lambda = 1e-1, side_info = list()))
# spline_0_box <- list(method = "spline",  params = list(lambda = 1, side_info = list(list(name="boundary_box"))))
# spline_0 <- list(method = "spline",  params = list(lambda = 1, side_info = list()))


# younger_fields <- list(type = "field", experiments = list(c(data = 1, estimator = 1), c(data = 1, estimator = 2),
#                                                           c(data = 1, estimator = 3), c(data = 1, estimator = 4),
#                                                           c(data = 1, estimator = 5), c(data = 1, estimator = 6),
#                                                           c(data = 1, estimator = 7), c(data = 1, estimator = 8),
#                                                           c(data = 1, estimator = 9), c(data = 1, estimator = 10),
#                                                           c(data = 1, estimator = 11), c(data = 1, estimator = 12)))#,
# older_fields <- list(type = "field", experiments = list(c(data = 2, estimator = 1), c(data = 2, estimator = 2),
#                                                         c(data = 2, estimator = 3), c(data = 2, estimator = 4),
#                                                         c(data = 2, estimator = 5), c(data = 2, estimator = 6),
#                                                         c(data = 2, estimator = 7), c(data = 2, estimator = 8),
#                                                         c(data = 2, estimator = 9), c(data = 2, estimator = 10),
#                                                         c(data = 2, estimator = 11), c(data = 2, estimator = 12)))
# 
# younger_paths <- list(type = "solution_path", experiments = list(c(data = 1, estimator = 1), c(data = 1, estimator = 2),
#                                                           c(data = 1, estimator = 3), c(data = 1, estimator = 4),
#                                                           c(data = 1, estimator = 9), c(data = 1, estimator = 10)))
# older_paths <- list(type = "solution_path", experiments = list(c(data = 2, estimator = 1), c(data = 2, estimator = 2),
#                                                         c(data = 2, estimator = 3), c(data = 2, estimator = 4),
#                                                         c(data = 2, estimator = 9), c(data = 2, estimator = 10)))

## Gradient Field Fitting
# spline_box <- list(method = "spline",  params = list(lambda = 1e-2, side_info = list(list(name="boundary_box"))))
# spline_nobox <- list(method = "spline",  params = list(lambda = 1e-2, side_info = list()))
# one_nn <- list(method = "knn",  params = list(k = 1))
# nw_base <- list(method = "nw",  params = list(h = 0.01))
# 
# younger_experiment <-  list(name = "younger", system = "knee", limit_cycle_samples = younger_smooth,
#                             lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 4)
# older_experiment <- list(name = "older", system = "knee", limit_cycle_samples = older_smooth,
#                          lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 4)
# 
# experiment_data <- generate_data_object_obs(list(younger_experiment, older_experiment))
# experiment_estimators <- list(one_nn, nw_base, spline_nobox, spline_box)
# experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)
# 
# younger_fields <- list(type = "field", experiments = list(c(data = 1, estimator = 1), c(data = 1, estimator = 2),
#                                                           c(data = 1, estimator = 3), c(data = 1, estimator = 4)))#,
# older_fields <- list(type = "field", experiments = list(c(data = 2, estimator = 1), c(data = 2, estimator = 2),
#                                                         c(data = 2, estimator = 3), c(data = 2, estimator = 4)))#,
# plot2 <- list(type = "solution_path", experiments = list(c(data = 1, estimator = 1), c(data = 1, estimator = 2)))#,
# visualize_results(experiment_results, list(younger_fields, older_fields))
