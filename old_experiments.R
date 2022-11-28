
#############
# #NOT Testing ##
#############
# #Variable Speed Circle:
# #some_noise <- list(name = "test", system = "var_circle", params = NA, n = 1000, sample_density = 0.1, var_x = 0, var_y = 0,
# #                   lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 2)
# 
# 
# # Loop 
# mu_vec <- c(1.5)
# var_vec <- c(0.01, 0.05, 0.25)
# n_basis_vec <- c(48, 96, 192, 384)
# 
# for(i in 1:length(mu_vec)){
#   for(j in 1:length(var_vec)){
#     for(k in 1:length(n_basis_vec)){
#       mu_val <- mu_vec[i]
#       var_val <- var_vec[j]
#       n_basis_val <- n_basis_vec[k]
#       name_str <- paste0("mu = ",mu_val,"; sigma^2 = ",var_val)
#       #data_list <- list(name = name_str, system = "van_der_pol", params = list(mu = mu_val), n = 1000, sample_density = 0.1, var_x = var_val, var_y = var_val,
#       #                  lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 10, noise_seed = 20)
#       data_list <- list(name = "test", system = "var_circle", params = NA, n = 1000, sample_density = 0.1, var_x = var_val, var_y = var_val,
#                       lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 2)
#       experiment_list = list(data_list)
#       experiment_estimators <- list(truth, spline_nobox)
#       experiment_data <- generate_data_object(experiment_list, noisy_smooth_basis = n_basis_val)
#       experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)
#     }
#   }
# }
# 

# set.seed(202)
# some_noise <- list(name = "low stiff; sigma^2 = 0.02", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = 0.02, var_y = 0.02,
#                    lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 10, noise_seed = 20)
# n_basis_vec <- c(24, 48, 96, 192, 384, 500)
# 
# experiment_list = list(some_noise)
# #experiment_data <- generate_data_object(experiment_list, noisy_smooth_basis = 96)
# for(i in 1:length(n_basis_vec)){
#   n_basis <- n_basis_vec[i]
#   experiment_data <- generate_data_object(experiment_list, noisy_smooth_basis = n_basis)
# }
# nw_larger <- list(method = "nw",  params = list(h = 1))
# 
# experiment_estimators <- list(truth, spline_box, spline_nobox)
# # 
# experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)
# plot1 <- list(type = "field", experiments = list(c(data = 1, estimator = 1), c(data = 1, estimator = 2),  c(data = 1, estimator = 3)))#,
# # #                                                 c(data = 1, estimator = 4), c(data = 1, estimator = 5), c(data = 1, estimator = 6)))
# plot2 <- list(type = "solution_path", experiments = list(c(data = 1, estimator = 1), c(data = 1, estimator = 2),  c(data = 1, estimator = 3)))#,
# # #                                                         c(data = 1, estimator = 4), c(data = 1, estimator = 5), c(data = 1, estimator = 6)))
# visualize_results(experiment_results, list(plot1))

# Noisy sample smoothing:
# some_noise <- list(name = "low stiff; sigma^2 = 0.01", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = 0.01, var_y = 0.01,
#                     lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 2)
# experiment_list <- list(some_noise)
# experiment_data <- generate_data_object(experiment_list)
# truth <- list(method = "truth",  params = list())
# spline <- list(method = "spline",  params = list(lambda = 1e-2, side_info = list()))
# spline_bounded <- list(method = "spline",  params = list(lambda = 1e-2, side_info = list("boundary_box")))
# experiment_estimators <- list(truth, spline, spline_bounded)
# experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)
# plot1 <- list(type = "field", experiments = list(c(data = 1, estimator = 1), c(data = 1, estimator = 2),  c(data = 1, estimator = 3)))
# visualize_results(experiment_results, list(plot1))

# Other

# # specify data
# no_noise <- list(name = "low stiff", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = 0, var_y = 0,
#                  lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 2)
#some_noise <- list(name = "low stiff; sigma^2 = 0.01", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = 0.01, var_y = 0.01,
#                    lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 2)
# some_noise_tps <- list(name = "low stiff; sigma^2 = 0.01; TPS", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = 0.01, var_y = 0.01,
#                    lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "tps", data_seed = 1, noise_seed = 2)
# more_noise <- list(name = "low stiff; sigma^2 = 0.05", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = 0.05, var_y = 0.05,
#                    lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 2)
# some_noise_hs <- list(name = "high stiff; sigma^2 = 0.01", system = "van_der_pol", params = list(mu = 20), n = 1000, sample_density = 0.1, var_x = 0.01, var_y = 0.01,
#                    lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 2)
# more_noise_hs <- list(name = "high stiff; sigma^2 = 0.05", system = "van_der_pol", params = list(mu = 20), n = 1000, sample_density = 0.1, var_x = 0.05, var_y = 0.05,
#                       lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 2)
# 
# experiment_list <- list(no_noise, some_noise, more_noise, some_noise_hs, more_noise_hs, some_noise_tps)
# experiment_data <- generate_data_object(experiment_list)
# 
# # specify estimators
# truth <- list(method = "truth",  params = list())
#one_nn <- list(method = "knn",  params = list(k = 1))
# many_nn <- list(method = "knn",  params = list(k = 400))
#nw_base <- list(method = "nw",  params = list(h = 0.01))
# loess <- list(method = "loess",  params = list(h = 0.1))
# loess_gcv <- list(method = "loess",  params = list(h = "gcv"))
# spline0 <- list(method = "spline",  params = list(lambda = 1e-16, side_info = list(list(name="boundary_box"))))
# spline1 <- list(method = "spline",  params = list(lambda = 1e-8, side_info = list(list(name="boundary_box"))))
# spline2 <- list(method = "spline",  params = list(lambda = 1e-6, side_info = list(list(name="boundary_box"))))
# spline3 <- list(method = "spline",  params = list(lambda = 1e-4, side_info = list(list(name="boundary_box"))))
# spline4 <- list(method = "spline",  params = list(lambda = 1e-2, side_info = list(list(name="boundary_box"))))
# spline5 <- list(method = "spline",  params = list(lambda = 1e-1, side_info = list(list(name="boundary_box"))))
# spline6 <- list(method = "spline",  params = list(lambda = 1, side_info = list(list(name="boundary_box"))))
# spline7 <- list(method = "spline",  params = list(lambda = 2, side_info = list(list(name="boundary_box"))))
# 
# experiment_estimators <- list(truth, spline0, spline1, spline2, spline3, spline4, spline5, spline6, spline7)
# experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)
# 
# plot1 <- list(type = "field", experiments = list(c(data = 1, estimator = 1), c(data = 1, estimator = 2),  c(data = 1, estimator = 3),
#                                                  c(data = 1, estimator = 4), c(data = 1, estimator = 5), c(data = 1, estimator = 6),
#                                                  c(data = 1, estimator = 7), c(data = 1, estimator = 8), c(data = 1, estimator = 9)))
# plot2 <- list(type = "solution_path", experiments = list(c(data = 1, estimator = 1), c(data = 1, estimator = 2),  c(data = 1, estimator = 3),
#                                                         c(data = 1, estimator = 4), c(data = 1, estimator = 5), c(data = 1, estimator = 6),
#                                                         c(data = 1, estimator = 7), c(data = 1, estimator = 8), c(data = 1, estimator = 9)))
# visualize_results(experiment_results, list(plot2))

# spline_eps <- list(method = "spline",  params = list(lambda = 1e-32))
# 
# experiment_estimators <- list(truth, one_nn, spline_eps)
# experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)
# plot1 <- list(type = "field", experiments = list(c(data = 2, estimator = 1), c(data = 2, estimator = 2),  c(data = 2, estimator = 3)))
# visualize_results(experiment_results, list(plot1))
# # pick estimators to run on each data set
# #experiment_estimators <- list(truth, loess, loess_gcv)
# #experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)
# 
# #plot1 <- list(type = "field", experiments = list(c(data = 2, estimator = 1), c(data = 2, estimator = 2),  c(data = 2, estimator = 3)))
# #plot2 <- list(type = "solution_path", experiments = list(c(data = 2, estimator = 1), c(data = 2, estimator = 2),  c(data = 2, estimator = 3)))
# #to_plot <- list(plot1, plot2)
# #visualize_results(experiment_results, to_plot)
# 
# experiment_estimators <- list(truth, one_nn, spline1, spline2, spline3, spline4, spline5, spline6, spline7)
# experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)
# plot3 <- list(type = "field", experiments = list(c(data = 2, estimator = 1), 
#                                                  c(data = 2, estimator = 2),
#                                                  c(data = 2, estimator = 3),
#                                                  c(data = 2, estimator = 4),
#                                                  c(data = 2, estimator = 5),
#                                                  c(data = 2, estimator = 6)
#                                                  ))
# plot4 <- list(type = "field", experiments = list(c(data = 3, estimator = 1), 
#                                                  c(data = 3, estimator = 2),
#                                                  c(data = 3, estimator = 3),
#                                                  c(data = 3, estimator = 4),
#                                                  c(data = 3, estimator = 5),
#                                                  c(data = 3, estimator = 6)
# ))
# 
# plot5 <- list(type = "solution_path", experiments = list(
#                                                  c(data = 6, estimator = 1), 
#                                                  c(data = 6, estimator = 2),
#                                                  c(data = 6, estimator = 3),
#                                                  c(data = 6, estimator = 4),
#                                                  c(data = 6, estimator = 5),
#                                                  c(data = 6, estimator = 6),
#                                                  c(data = 6, estimator = 7),
#                                                  c(data = 6, estimator = 8),
#                                                  c(data = 6, estimator = 9)
#                                                  ))
# 
# #to_plot <- list(plot3, plot4)
# #visualize_results(experiment_results, to_plot)
# 
# to_plot <- list(plot5)
# visualize_results(experiment_results, to_plot)
# 
# plot6 <- list(type = "solution_path", experiments = list(
#  c(data = 2, estimator = 4),
#  c(data = 2, estimator = 5),
#  c(data = 2, estimator = 6)
# ))
# 
# to_plot <- list(plot6)
# visualize_results(experiment_results, to_plot)
# 
# 
# # Experiment which varies noise
# varx <- 0.01
# vary <- 0.01
# noise_experiment_1 <- list(name = "low stiff; sigma^2 = 0.01", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = varx, var_y = vary,
#                            lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 2)
# noise_experiment_2 <- list(name = "low stiff; sigma^2 = 0.01", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = varx, var_y = vary,
#                            lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 3)
# noise_experiment_3 <- list(name = "low stiff; sigma^2 = 0.01", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = varx, var_y = vary,
#                            lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 4)
# noise_experiment_4 <- list(name = "low stiff; sigma^2 = 0.01", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = varx, var_y = vary,
#                            lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 5)
# noise_experiment_5 <- list(name = "low stiff; sigma^2 = 0.01", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = varx, var_y = vary,
#                            lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 6)
# noise_experiment_6 <- list(name = "low stiff; sigma^2 = 0.01", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = varx, var_y = vary,
#                            lc_tail_n = 500, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5, smoother = "bspline", data_seed = 1, noise_seed = 7)
# 
# experiment_list <- list(noise_experiment_1, noise_experiment_2, noise_experiment_3, noise_experiment_4, noise_experiment_5, noise_experiment_6)
# experiment_data <- generate_data_object(experiment_list)
# experiment_estimators <- list(truth, one_nn, spline1, spline2, spline3, spline4, spline5, spline6)
# experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)
# 
# 
# noise_comp_paths <- list(type = "solution_path", experiments = list(
#   c(data = 1, estimator = 6), 
#   c(data = 2, estimator = 6),
#   c(data = 3, estimator = 6),
#   c(data = 4, estimator = 6),
#   c(data = 5, estimator = 6),
#   c(data = 6, estimator = 6)
# ))
# 
# to_plot <- list(noise_comp_paths)
# visualize_results(experiment_results, to_plot)
# 
# 
