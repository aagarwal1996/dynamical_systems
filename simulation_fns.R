#############
## Imports ##
#############

library(fda)
library(tidyverse)
library(mvtnorm) # for Multivariate Normal density estimates in KDE
library(shape) # for pretty arrows
library(ggquiver) # for vector field plots
library(ggpubr) # to generate side-by-side plots
library(fields) # for thin-plate splines

source(here::here('Data_Generation','data_generation.R')) # functions to generate data from specific DS are in another file
source(here::here('Estimation_Methods','bspline.R'))
source(here::here('Estimation_Methods','bspline_gd_redux.R'))
source(here::here('Estimation_Methods','knn.R'))
library(Rcpp)
sourceCpp(here::here('Estimation_Methods','nw_regression.cpp'))
sourceCpp(here::here('Estimation_Methods','loess.cpp'))

#####################
## Data Generation ##
#####################

generate_limit_cycle_data <- function(system, params, var_x, var_y, data_seed = 2022, noise_seed = 302,
                                      num_samples = 1500, sample_density = 0.1,
                                      save_csv = F, smooth = "bspline", title = "",
                                      noisy_smooth_basis = 48){  
  set.seed(data_seed)
  
  if (system == "van_der_pol"){
    sampled_data <- generate_van_der_pol_samples(params, num_samples=num_samples,sample_density=sample_density)
  } else if (system == "rzma"){
    sampled_data <- generate_rzma_samples(params, num_samples=num_samples,sample_density=sample_density)
  } else if(system == "lotka_volterra"){
    sampled_data <- generate_lv_samples(params, num_samples=num_samples,sample_density=sample_density)
  } else if(system == "log_lotka_volterra"){
    sampled_data <- generate_log_lv_samples(params, num_samples=num_samples,sample_density=sample_density)
  } else if(system == "abhi"){
    sampled_data <- get_abhi_data()
  } else{
    stop('Error: ', system, ' system not implemented.')
  }
  
  if (save_csv){
    # TODO: improve file name
    file_name <- paste0("Saved_Data/",system,"-",format(Sys.time(), "%m_%d_%Y-%H_%M_%S"),".csv")
    write_csv(sampled_data, file_name)
  }
  
  if ((var_x != 0) | (var_y != 0)){
    # add Gaussian noise to position of observations
    set.seed(noise_seed)
    noise_matrix <- mvrnorm(nrow(sampled_data), c(0,0,0,0), diag(c(var_x,var_y,0,0)))
    noisy_samples <- sampled_data + noise_matrix 
    smoothed_positions <- spline_smooth_noisy_samples(sampled_data, noisy_samples, nbasis=noisy_smooth_basis, max_t = num_samples * sample_density,
                                                      lambda = 1e-12, return_type = smooth, title = title)
    #smoothed_positions <- loess_smooth_noisy_samples(noisy_samples,h=25)
    sampled_data <- cbind(smoothed_positions, noisy_samples[,c(3,4)])
  }
  colnames(sampled_data) <- c("x","y","f_x","f_y")
  return(sampled_data)
}

generate_data_object_model <- function(experiment_list, noisy_smooth_basis = 48, save_csv = F){
  # copy to modify
  data_list <- experiment_list
  
  for (i in 1:length(data_list)){
    # extract params
    experiment_name <- data_list[[i]]$name
    system_name <- data_list[[i]]$system
    system_params <- data_list[[i]]$params
    num_samples <- data_list[[i]]$n
    sample_density <- data_list[[i]]$sample_density
    var_x <- data_list[[i]]$var_x
    var_y <- data_list[[i]]$var_y
    smoother <- data_list[[i]]$smoother
    data_seed <- data_list[[i]]$data_seed
    noise_seed <- data_list[[i]]$noise_seed
    lc_tail_n <- data_list[[i]]$lc_tail_n
    
    # generate data
    data_list[[i]]$limit_cycle_samples <- generate_limit_cycle_data(system_name, system_params, 
                                   var_x = var_x, var_y = var_y, data_seed = data_seed, noise_seed = noise_seed,
                                   num_samples = num_samples, sample_density = sample_density,
                                   smooth = smoother, noisy_smooth_basis = noisy_smooth_basis, title = experiment_name,
                                   save_csv = save_csv)
    
    # truncate to tail
    data_list[[i]]$limit_cycle_tail <- tail(data_list[[i]]$limit_cycle_samples, n = data_list[[i]]$lc_tail_n)
    
    # build grid to evaluate over
    grid_object <- list(xmin = min(data_list[[i]]$limit_cycle_tail$x),
                        xmax = max(data_list[[i]]$limit_cycle_tail$x),
                        ymin = min(data_list[[i]]$limit_cycle_tail$y),
                        ymax = max(data_list[[i]]$limit_cycle_tail$y)) 
    
    grid_object$x_grid = seq(floor(grid_object$xmin) - data_list[[i]]$extrapolation_size , ceiling(grid_object$xmax) + data_list[[i]]$extrapolation_size, 
                             len = data_list[[i]]$x_grid_size)
    grid_object$y_grid = seq(floor(grid_object$ymin) - data_list[[i]]$extrapolation_size, ceiling(grid_object$ymax) + data_list[[i]]$extrapolation_size, 
                             len = data_list[[i]]$y_grid_size)
    grid_object$eval_grid = unname(as.matrix(expand.grid(grid_object$x_grid,grid_object$y_grid)))
    
    data_list[[i]]$grid <-grid_object
  }
  
  return(data_list)
}

generate_data_object_obs <- function(experiment_list, noisy_smooth_basis = 48, save_csv = F){
  # No explicit model
  data_list <- experiment_list
  
  for (i in 1:length(data_list)){
    # extract params
    experiment_name <- data_list[[i]]$name
    system_name <- data_list[[i]]$system
    system_params <- NA
    num_samples <- NA
    sample_density <- NA
    var_x <- NA
    var_y <- NA
    smoother <- NA
    data_seed <- NA
    noise_seed <- NA
    
    # truncate to tail
    data_list[[i]]$limit_cycle_tail <- tail(data_list[[i]]$limit_cycle_samples, n = data_list[[i]]$lc_tail_n)
    
    # build grid to evaluate over
    grid_object <- list(xmin = min(data_list[[i]]$limit_cycle_tail$x),
                        xmax = max(data_list[[i]]$limit_cycle_tail$x),
                        ymin = min(data_list[[i]]$limit_cycle_tail$y),
                        ymax = max(data_list[[i]]$limit_cycle_tail$y)) 
    
    grid_object$x_grid = seq(floor(grid_object$xmin) - data_list[[i]]$extrapolation_size , ceiling(grid_object$xmax) + data_list[[i]]$extrapolation_size, 
                             len = data_list[[i]]$x_grid_size)
    grid_object$y_grid = seq(floor(grid_object$ymin) - data_list[[i]]$extrapolation_size, ceiling(grid_object$ymax) + data_list[[i]]$extrapolation_size, 
                             len = data_list[[i]]$y_grid_size)
    grid_object$eval_grid = unname(as.matrix(expand.grid(grid_object$x_grid,grid_object$y_grid)))
    
    data_list[[i]]$grid <-grid_object
  }
  
  return(data_list)
}

################
## Evaluation ##
################


evaluate_gradient_single <- function(data_object, estimators_list, 
                                      lc_tail_n = 700, 
                                      x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5){
  for (i in 1:length(estimators_list)){
    if (estimators_list[[i]]$method == "truth"){
        estimators_list[[i]]$params <- data_object$params
        estimators_list[[i]]$params$name <- data_object$system
        estimators_list[[i]]$estimated_field <- get_gradient_field(data_object,estimators_list[[i]])$evaluated_field
    }
  
    else if (estimators_list[[i]]$method == "spline"){
      evaluation_result <-  get_gradient_field(data_object, estimators_list[[i]])
      estimators_list[[i]]$params <- evaluation_result$params
      estimators_list[[i]]$estimated_field = evaluation_result$field
      estimators_list[[i]]$sfd_list = evaluation_result$sfd_list
    }
    else if (estimators_list[[i]]$method == "gd_spline"){
      evaluation_result <-  get_gradient_field(data_object, estimators_list[[i]])
      estimators_list[[i]]$params <- evaluation_result$params
      estimators_list[[i]]$estimated_field = evaluation_result$field
      estimators_list[[i]]$sfd_list = evaluation_result$sfd_list
    }
    else{
      evaluation_result <-  get_gradient_field(data_object, estimators_list[[i]])
      estimators_list[[i]]$estimated_field = evaluation_result$evaluated_field
      estimators_list[[i]]$params <- evaluation_result$params
    }
    
    if (TRUE){
      title_str = title_str = paste("Data Set: ", data_object$name,
                                    "\nEstimator: ", estimators_list[[i]]$method,
                                    "\nParams: ", paste(names(estimators_list[[i]]$params),
                                                        estimators_list[[i]]$params, sep = ":", collapse = "\n"))
      estimators_list[[i]]$field_plot = ggplot_field(data_object$limit_cycle_tail, data_object$grid$eval_grid, 
                                                     estimators_list[[i]]$estimated_field, rescale = .025,
                                                     title=title_str) + theme(text = element_text(size = 4))
      }
  }
  
  return(estimators_list)
}

evaluate_gradient_methods <- function(data_list, estimator_list){
  for (i in 1:length(data_list)){
    estimation_results <- evaluate_gradient_single(data_list[[i]], estimator_list)  
    data_list[[i]]$estimates <- estimation_results
  }
  
  return(data_list)

}

###################
## Visualization ##
###################

ggplot_field <- function(ls_samples, eval_grid, gradient_data, title="", rescale = 1, sample_alpha=1){
  # plot a gradient field with limit cycle samples overlayed
  
  sample_tibble <- tibble(x = ls_samples[,1], y = ls_samples[,2], u = NA, v = NA)
  gradient_tibble <- tibble(x = eval_grid[,1], y = eval_grid[,2], 
                            u = rescale*gradient_data[,1], v =rescale*gradient_data[,2])
  field_plot <- ggplot(gradient_tibble, aes(x = x, y = y, u = u, v = v)) +
    geom_quiver(color = "#003262",vecsize=rescale) +
    geom_point(data = sample_tibble, aes(x = x, y = y), color = "#FDB515", alpha = sample_alpha) +
    labs(title = title) + 
    xlim(min(gradient_tibble$x)-.5, max(gradient_tibble$x)+.5) + #TODO: Make proportional
    ylim(min(gradient_tibble$y)-.5, max(gradient_tibble$y)+.5)
  
  return(field_plot)
}

plot_estimated_fields <- function(plot_specifications, experiment_outcome){
  n_col = min(2, length(plot_specifications))
  plot_list <- list()
  for (i in 1:length(plot_specifications)){
    data_num <- plot_specifications[[i]][1]
    estimator_num <- plot_specifications[[i]][2]
    plot_list[[i]] <- experiment_outcome[[data_num]]$estimates[[estimator_num]]$field_plot
  }
  
  field_plots <- cowplot::plot_grid(plotlist = plot_list, ncol = n_col)
  print(field_plots)
}

plot_field_delta <- function(plot_specifications, experiment_outcome, sol_paths = F){
  plot_list <- list()
  if (sol_paths){
    solution_graphs <- list()
    # get shared ICs based on limit cycle data from first experiment
    shared_ic <- get_shared_ic(experiment_outcome[[1]]$limit_cycle_tail)
    
    for (i in 1:length(plot_specifications)){
      data_num <- plot_specifications[[i]][1]
      estimator_num <- plot_specifications[[i]][2]
      
      ic_results_list <- list()
      
      for (j in 1:nrow(shared_ic)){
        ic_results_list[[j]] <- generate_solution_path(experiment_outcome[[data_num]]$limit_cycle_tail,
                                                       experiment_outcome[[data_num]]$estimates[[estimator_num]], 
                                                       shared_ic[j,])
      }
      names(ic_results_list) <- rownames(shared_ic)
      title_str = paste("Data Set: ", experiment_outcome[[data_num]]$name,
                        "\nEstimator: ", experiment_outcome[[data_num]]$estimates[[estimator_num]]$method,
                        "\nParams: ", paste(names(experiment_outcome[[data_num]]$estimates[[estimator_num]]$params),
                                            experiment_outcome[[data_num]]$estimates[[estimator_num]]$params, sep = ":", collapse = "\n"))
      gradient_list <- list(estimated_field = experiment_outcome[[data_num]]$estimates[[estimator_num]]$estimated_field,
                            eval_grid = experiment_outcome[[data_num]]$grid$eval_grid)
      
      solution_path_plot <- plot_single_solution_path(ic_results_list, gradient_list, experiment_outcome[[data_num]]$limit_cycle_tail, title_str) + 
        theme(text = element_text(size = 4))
      
      solution_graphs[[i]] <- solution_path_plot
    }
  }
  
  for (i in 1:length(plot_specifications)){
    data_num <- plot_specifications[[i]][1]
    estimator_num <- plot_specifications[[i]][2]
    ref_num <- plot_specifications[[i]][3]
    ls_samples <- experiment_outcome[[data_num]]$limit_cycle_samples
    eval_grid <- experiment_outcome[[data_num]]$grid$eval_grid
    reference_data <- experiment_outcome[[data_num]]$estimates[[ref_num]]$estimated_field
    gradient_data <- experiment_outcome[[data_num]]$estimates[[estimator_num]]$estimated_field
    delta_grad <- gradient_data - reference_data
    if (sol_paths){
      plot_list[[(2*i - 1)]] <- solution_graphs[[i]]
    } else{ plot_list[[(2*i - 1)]] <- ggplot_field(ls_samples, eval_grid, reference_data, title="Reference", rescale = 1)} # truth
    
    title_str = paste("Alg: ",experiment_outcome[[data_num]]$estimates[[estimator_num]]$params$gd_params$algorithm," Skip:",
                      experiment_outcome[[data_num]]$estimates[[estimator_num]]$params$gd_params$batching$skip_negative)
    plot_list[[(2*i)]] <- ggplot_field(ls_samples, eval_grid, delta_grad, title=title_str, rescale = 50, sample_alpha = 0.05) # delta
  }
  
  n_col = length(plot_specifications)
  delta_plots <- cowplot::plot_grid(plotlist = plot_list, ncol = n_col,byrow=F)
  print(delta_plots)
}

get_shared_ic <- function(ls_samples, delta_ring = FALSE){
  # samples random ICs from within, on, outside close, and outside far from the LC
  if (!delta_ring){
    # sample points in a radial fashion, for symmetric LCs
    ic_init_samples <- ls_samples[runif(4, 1, nrow(ls_samples)),c(1,2)]
    ic_radial <- c(0.5, 1, 1.25, 1.4)
    ic_samples <- ic_init_samples * ic_radial
    
  } else if (delta_ring){
    set.seed(20)
    # for the knee data, use a perturbation approach
    delta_ring_samples <- ls_samples[sample(nrow(ls_samples),size=2,replace=FALSE),]
    colnames(delta_ring_samples) <- c("x","y","f_x","f_y")
    delta_ring_sideinfo <- t(apply(delta_ring_samples,1,get_normal_vectors, radius = .75)) # TODO: functionalize
    ic_samples <- rbind(delta_ring_sideinfo[,c(1,2)],delta_ring_sideinfo[,c(3,4)])
  }
  # format
  rownames(ic_samples) <- c("Within", "On", "OutsideClose", "OutsideFar")
  colnames(ic_samples) <- c("x", "y")
  
  return(ic_samples)
}

plot_single_solution_path <- function(solution_path_list, gradient_list, lc_data, title_str = ""){
  within_tibble <- tibble(x = solution_path_list$Within$x, y = solution_path_list$Within$y)
  on_tibble <- tibble(x = solution_path_list$On$x, y = solution_path_list$On$y)
  outside_close_tibble <- tibble(x = solution_path_list$OutsideClose$x, y = solution_path_list$OutsideClose$y)
  outside_far_tibble <- tibble(x = solution_path_list$OutsideFar$x, y = solution_path_list$OutsideFar$y)

  gradient_matrix <- matrix(cbind(gradient_list$eval_grid, gradient_list$estimated_field), ncol=4)
  colnames(gradient_matrix) <- c("x","y","u","v")
  gradient_tibble <- as_tibble(gradient_matrix)
  
  path_plot <- ggplot(gradient_tibble, aes(x = x, y = y, u = u, v = v)) + 
    geom_quiver(color = "#003262") +
    geom_point(data = lc_data, aes(x = x, y = y, u = NA, v = NA), color = "#FDB515") +
    geom_path(data = within_tibble, aes(x = x, y = y, u = NA, v = NA), color = "red") +
    geom_point(x = as.numeric(within_tibble[1,1]), y =as.numeric(within_tibble[1,2], u = NA, v = NA), color = "red", size = 3) + 
    geom_path(data = on_tibble, aes(x = x, y = y, u = NA, v = NA), color = "blue") +
    geom_point(x = as.numeric(on_tibble[1,1]), y =as.numeric(on_tibble[1,2], u = NA, v = NA), color = "blue", size = 3) +
    geom_path(data = outside_close_tibble, aes(x = x, y = y, u = NA, v = NA), color = "green") +
    geom_point(x = as.numeric(outside_close_tibble[1,1]), y =as.numeric(outside_close_tibble[1,2], u = NA, v = NA), color = "green", size = 3) + 
    geom_path(data = outside_far_tibble, aes(x = x, y = y, u = NA, v = NA), color = "purple") +
    geom_point(x = as.numeric(outside_far_tibble[1,1]), y =as.numeric(outside_far_tibble[1,2], u = NA, v = NA), color = "purple", size = 3) +
    labs(title=title_str) 
  
  return(path_plot)
}

plot_solution_paths <- function(plot_specifications, experiment_outcome, 
                                samples_per_path = 500){
  
  solution_graphs <- list()
  # get shared ICs based on limit cycle data from first experiment
  shared_ic <- get_shared_ic(experiment_outcome[[1]]$limit_cycle_tail)
  
  for (i in 1:length(plot_specifications)){
    data_num <- plot_specifications[[i]][1]
    estimator_num <- plot_specifications[[i]][2]
    
    ic_results_list <- list()
    
    for (j in 1:nrow(shared_ic)){
        ic_results_list[[j]] <- generate_solution_path(experiment_outcome[[data_num]]$limit_cycle_tail,
                                                       experiment_outcome[[data_num]]$estimates[[estimator_num]], 
                                                       shared_ic[j,])
    }
    names(ic_results_list) <- rownames(shared_ic)
    title_str = paste("Data Set: ", experiment_outcome[[data_num]]$name,
                      "\nEstimator: ", experiment_outcome[[data_num]]$estimates[[estimator_num]]$method,
                      "\nParams: ", paste(names(experiment_outcome[[data_num]]$estimates[[estimator_num]]$params),
                                          experiment_outcome[[data_num]]$estimates[[estimator_num]]$params, sep = ":", collapse = "\n"))
    gradient_list <- list(estimated_field = experiment_outcome[[data_num]]$estimates[[estimator_num]]$estimated_field,
                 eval_grid = experiment_outcome[[data_num]]$grid$eval_grid)
    
    solution_path_plot <- plot_single_solution_path(ic_results_list, gradient_list, experiment_outcome[[data_num]]$limit_cycle_tail, title_str) + 
      theme(text = element_text(size = 4))
    
    solution_graphs[[i]] <- solution_path_plot
  }

  n_col = min(2, length(plot_specifications))
  sol_plot <- cowplot::plot_grid(plotlist = solution_graphs, ncol = n_col)
  print(sol_plot)
}



visualize_results <- function(experiment_outcome, plot_list){
  for (plot in plot_list){
    if (plot$type == "field"){
      plot_estimated_fields(plot$experiments, experiment_outcome)
    }  else if (plot$type == "field_delta"){
      plot_field_delta(plot$experiments, experiment_outcome)
    }  else if (plot$type == "field_delta_paths"){
      plot_field_delta(plot$experiments, experiment_outcome, sol_paths = T)
    }
    else if (plot$type == "solution_path"){
      plot_solution_paths(plot$experiments, experiment_outcome)
    }
    else {
      print("Plot not implemented")
    }
  }
} 