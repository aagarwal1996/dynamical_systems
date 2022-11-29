#############
## Testing ##
#############

## Van der Pol
truth <- list(method = "truth",  params = list())
spline <- list(method = "spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20, side_info = list()))
spline_ring <- list(method = "spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20,
                                                      side_info = list(list(name="delta_ring",sample_frac = .2, radius=0.2, shift = 0, si_magnitude = .25,
                                                                            lc_params = list(system = "van_der_pol", params = list(mu = 3))))
))

gd_spline_vanilla2_baseline <- list(method = "gd_spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20,
                                                                         side_info = list(),
                                                                         gd_params = list(algorithm = "Vanilla", eig_rank = 2,eta=0.02, dt_radius=0, 
                                                                                          batching = list(num_iterations=1,batches_per=25,batch_size=1,skip_negative=T))
))

gd_spline_vanilla <- list(method = "gd_spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20,
                                                               side_info = list(),
                                                               gd_params = list(algorithm = "Vanilla", eig_rank = 1,eta=0.02, dt_radius=0, 
                                                                                batching = list(num_iterations=1,batches_per=25,batch_size=1,skip_negative=T))
))

gd_spline_vanilla2 <- list(method = "gd_spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20,
                                                                side_info = list(),
                                                                gd_params = list(algorithm = "Vanilla", eig_rank = 2,eta=0.02, dt_radius=0, 
                                                                                 batching = list(num_iterations=1,batches_per=25,batch_size=1,skip_negative=T))
))

gd_spline_vanilla3 <- list(method = "gd_spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20,
                                                                side_info = list(),
                                                                gd_params = list(algorithm = "Vanilla", eig_rank = 1,eta=0.02, dt_radius=0, 
                                                                                 batching = list(num_iterations=1,batches_per=50,batch_size=1,skip_negative=T))
))

gd_spline_vanilla4 <- list(method = "gd_spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20,
                                                                side_info = list(),
                                                                gd_params = list(algorithm = "Vanilla", eig_rank = 1,eta=0.02, dt_radius=0, 
                                                                                 batching = list(num_iterations=1,batches_per=25,batch_size=1,skip_negative=T))
))

gd_spline_vanilla_proj <- list(method = "gd_spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20,
                                                                    side_info = list(),
                                                                    gd_params = list(algorithm = "Projection", eig_rank = 1,eta=0.02, dt_radius=0.2, 
                                                                                     batching = list(num_iterations=5,batches_per=300,batch_size=1,skip_negative=T))
))

gd_spline_vanilla_proj2 <- list(method = "gd_spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20,
                                                                     side_info = list(),
                                                                     gd_params = list(algorithm = "Projection", eig_rank = 2,eta=0.02, dt_radius=0.2, 
                                                                                      batching = list(num_iterations=5,batches_per=300,batch_size=1,skip_negative=T))
))

gd_spline_vanilla_proj2_bigR <- list(method = "gd_spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20,
                                                                          side_info = list(),
                                                                          gd_params = list(algorithm = "Projection", eig_rank = 2,eta=0.02, dt_radius=0.5, 
                                                                                           batching = list(num_iterations=5,batches_per=300,batch_size=1,skip_negative=T))
))

gd_spline_vanilla_proj2_NS <- list(method = "gd_spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20,
                                                                        side_info = list(),
                                                                        gd_params = list(algorithm = "Projection", eig_rank = 2,eta=0.02, dt_radius=0.2, 
                                                                                         batching = list(num_iterations=5,batches_per=300,batch_size=1,skip_negative=F))
))

gd_spline_vanilla_proj2_bigR_NS <- list(method = "gd_spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20,
                                                                             side_info = list(),
                                                                             gd_params = list(algorithm = "Projection", eig_rank = 2,eta=0.02, dt_radius=0.5, 
                                                                                              batching = list(num_iterations=5,batches_per=300,batch_size=1,skip_negative=F))
))

gd_spline_vanilla_motion <- list(method = "gd_spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20,
                                                                      side_info = list(),
                                                                      gd_params = list(algorithm = "Projection_Motion", eig_rank = 1,eta=0.02, dt_radius=0, 
                                                                                       batching = list(num_iterations=5,batches_per=300,batch_size=1,skip_negative=T))
))

gd_spline_vanilla_motion_skip <- list(method = "gd_spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20,
                                                                           side_info = list(),
                                                                           gd_params = list(algorithm = "Projection_Motion", eig_rank = 1,eta=0.02, dt_radius=0, 
                                                                                            batching = list(num_iterations=5,batches_per=300,batch_size=1,skip_negative=F))
))


some_noise <- list(name = "VDP", system = "van_der_pol", params = list(mu=3), n = 1000, sample_density = 0.1, var_x = 0.05, var_y = 0.05,
                   lc_tail_n = 500, x_grid_size = 36, y_grid_size = 36, extrapolation_size = 0.5, smoother = "bspline", data_seed = 2, noise_seed = 2)
experiment_list <- list(some_noise)
experiment_data <- generate_data_object(experiment_list)
experiment_estimators <- list(truth, spline)
experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)

# plot_alpha <- list(type = "solution_path", experiments = list(c(data = 1, estimator = 1), c(data = 1, estimator = 2),
#                                                               c(data = 1, estimator = 3),  c(data = 1, estimator = 4),
#                                                               c(data = 1, estimator = 5),  c(data = 1, estimator = 6)))
# visualize_results(experiment_results, list(plot_alpha))

plot_delta <- list(type = "field_delta_paths", experiments = list(c(data = 1, estimator = 1, ref = 1),c(data = 1, estimator = 2, ref = 1)))
visualize_results(experiment_results, list(plot_delta))

## Lotka-Volterra
first_lv <- list(name = "LV", system = "lotka_volterra", params = list(alpha=2/3,beta=4/3,delta=1,gamma=1), n = 10000, sample_density = 0.01, var_x = 0.0, var_y = 0.0,
                 lc_tail_n = 10000, x_grid_size = 36, y_grid_size = 36, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 2)
experiment_list <- list(first_lv)
experiment_data <- generate_data_object(experiment_list)
experiment_estimators <- list(truth, spline,gd_spline_vanilla,gd_spline_vanilla_proj)
experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)

plot_lv <- list(type = "field_delta_paths", experiments = list(c(data = 1, estimator = 1, ref = 1),c(data = 1, estimator = 2, ref = 1),
                                                               c(data = 1, estimator = 3, ref = 1),c(data = 1, estimator = 4, ref = 1)))
visualize_results(experiment_results, list(plot_lv))

## Log Lotka-Volterra
lg_lv <- list(name = "LV", system = "log_lotka_volterra", params = list(alpha=2/3,beta=2/3,delta=1,gamma=1), n = 10000, lc_length = 3000, sample_density = 0.25, var_x = 0.0, var_y = 0.0,
              lc_tail_n = 10000, x_grid_size = 36, y_grid_size = 36, extrapolation_size = 0.5, smoother = "bspline", data_seed = 2, noise_seed = 2)
experiment_list <- list(lg_lv)
experiment_data <- generate_data_object(experiment_list)
experiment_estimators <- list(truth)#, spline,gd_spline_vanilla,gd_spline_vanilla_proj)
experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)

plot_lg_lv <- list(type = "field_delta_paths", experiments = list(c(data = 1, estimator = 1, ref = 1)))#,c(data = 1, estimator = 2, ref = 1),
#                                                               c(data = 1, estimator = 3, ref = 1),c(data = 1, estimator = 4, ref = 1)))
visualize_results(experiment_results, list(plot_lg_lv))

## Gause
first_gause <- list(name = "RZ Test", system = "gause", params = list(a = 2.8, b = 0.7, c = 1.35, r = 3.5, k = 1.5, K = 1.4), n = 10000, sample_density = 0.1, var_x = 0.03, var_y = 0.03,
                    lc_tail_n = 10000, x_grid_size = 36, y_grid_size = 36, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 2)
experiment_list <- list(first_gause)
experiment_data <- generate_data_object(experiment_list)
experiment_estimators <- list(truth, spline,gd_spline_vanilla,gd_spline_vanilla2)
experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)


plot_gause <- list(type = "field_delta_paths", experiments = list(c(data = 1, estimator = 1, ref = 1),c(data = 1, estimator = 2, ref = 1),
                                                                  c(data = 1, estimator = 3, ref = 1),c(data = 1, estimator = 4, ref = 1)))#,
#c(data = 1, estimator = 5, ref = 1),c(data = 1, estimator = 6, ref = 1)))
visualize_results(experiment_results, list(plot_gause))


# plot1 <- list(type = "field", experiments = list(c(data = 1, estimator = 1), c(data = 1, estimator = 2),
#                                                  c(data = 1, estimator = 3), c(data = 1, estimator = 4),
#                                                  c(data = 1, estimator = 5), c(data = 1, estimator = 6)))
# visualize_results(experiment_results, list(plot1))
# plot2 <- list(type = "solution_path", experiments = list(c(data = 1, estimator = 1), c(data = 1, estimator = 2),
#                                                          c(data = 1, estimator = 3),c(data = 1, estimator = 4),
#                                                          c(data = 1, estimator = 5), c(data = 1, estimator = 6)))
# 
# visualize_results(experiment_results, list(plot2))