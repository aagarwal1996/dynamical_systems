#############
## Imports ##
#############

library(fda)
library(deSolve)
library(tidyverse)
library(mvtnorm) # for Multivariate Normal density estimates in KDE
library(shape) # for pretty arrows
library(ggquiver) # for vector field plots
library(ggpubr) # to generate side-by-side plots

source(here::here('Data_Generation/data_generation.R')) # functions to generate data from specific DS are in another file
source(here::here('Estimation_Methods/bspline.R'))
source(here::here('Estimation_Methods/knn.R'))
library(Rcpp)
sourceCpp(here::here('Estimation_Methods/nw_regression.cpp'))
sourceCpp(here::here('Estimation_Methods/loess.cpp'))

#####################
## Data Generation ##
#####################

generate_limit_cycle_data <- function(system, params, var_x, var_y, sample_seed = 2022,
                                      num_samples = 1500, sample_density = 0.1,
                                      use_seed = F, save_csv = F){  
  set.seed(sample_seed)
  
  if (system == "van_der_pol"){
    sampled_data <- generate_van_der_pol(params, num_samples=num_samples,sample_density=sample_density)
  }
  else if(system == "abhi"){
    sampled_data <- get_abhi_data()
  }
  else{
    stop('Error: ', system, ' system not implemented.')
  }
  
  if (save_csv){
    # TODO: improve file name
    file_name <- paste0("Saved_Data/",system,"-",format(Sys.time(), "%m_%d_%Y-%H_%M_%S"),".csv")
    write_csv(sampled_data, file_name)
  }
  
  # add Gaussian noise to position of observations
  noise_matrix <- mvrnorm(nrow(sampled_data), c(0,0,0,0), diag(c(var_x,var_y,0,0)))
  noisy_samples <- sampled_data + noise_matrix
  
  return(noisy_samples)
}

generate_data_object <- function(experiment_list,
                               fixed_seed = F, save_csv = F){
  # copy to modify
  data_list <- experiment_list

  sample_seed <- 1
  # all limit cycle samples are the same before adding varying levels of noise
  ifelse(fixed_seed, sample_seed <- 2022, sample_seed <- round(runif(1,1,10000)))
  
  for (i in 1:length(data_list)){
    # extract params
    experiment_name <- data_list[[i]]$name
    system_name <- data_list[[i]]$system
    system_params <- data_list[[i]]$params
    num_samples <- data_list[[i]]$n
    sample_density <- data_list[[i]]$sample_density
    var_x <- data_list[[i]]$var_x
    var_y <- data_list[[i]]$var_y
    
    # generate data
    data_list[[i]]$limit_cycle_samples <- generate_limit_cycle_data(system_name, system_params, 
                                   var_x = var_x, var_y = var_y, sample_seed = sample_seed,
                                   num_samples = num_samples, sample_density = sample_density,
                                   save_csv = save_csv)
    
    # truncate to tail
    data_list[[i]]$limit_cycle_tail <- tail(data_list[[i]]$limit_cycle_samples, n = data_list[[i]]$lc_tail_n)
    
    # build grid to evaluate over
    grid_object <- list(xmin = min(data_list[[i]]$limit_cycle_tail$x),
                        xmax = max(data_list[[i]]$limit_cycle_tail$x),
                        ymin = min(data_list[[i]]$limit_cycle_tail$y),
                        ymax = max(data_list[[i]]$limit_cycle_tail$y)) 
    #View(data_list[[i]])
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
    }
    
    if (estimators_list[[i]]$method == "spline"){
      spline_result <- get_gradient_field(data_object, estimators_list[[i]])
      estimators_list[[i]]$estimated_field = spline_result$field
      estimators_list[[i]]$sfd_list = spline_result$sfd_list
    }
    
    else{
      estimators_list[[i]]$estimated_field = get_gradient_field(data_object, estimators_list[[i]])
    }
    
    estimators_list[[i]]$field_plot = ggplot_field(data_object$limit_cycle_tail, data_object$grid$eval_grid, 
                                                   estimators_list[[i]]$estimated_field, rescale = .025,
                                                   title=paste(estimators_list[[i]]$method, 
                                                               paste(names(estimators_list[[i]]$params), estimators_list[[i]]$params,
                                                                     sep = ":", collapse = ",")))
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

ggplot_field <- function(ls_samples, eval_grid, gradient_data, title="", rescale = 0.005){
  # plot a gradient field with limit cycle samples overlayed
  sample_tibble <- tibble(x = ls_samples[,1], y = ls_samples[,2], u = NA, v = NA)
  gradient_tibble <- tibble(x = eval_grid[,1], y = eval_grid[,2], 
                            u = rescale*gradient_data[,1], v =rescale*gradient_data[,2])
  
  field_plot <- ggplot(gradient_tibble, aes(x = x, y = y, u = u, v = v)) +
    geom_quiver(color = "#003262", vecsize = 0) +
    geom_point(data = sample_tibble, aes(x = x, y = y), color = "#FDB515") +
    labs(title = title)
  
  return(field_plot)
}

plot_estimated_fields <- function(plot_specifications, experiment_outcome){
  n_col = min(3, length(plot_specifications))
  plot_list <- list()
  for (i in 1:length(plot_specifications)){
    data_num <- plot_specifications[[i]][1]
    estimator_num <- plot_specifications[[i]][2]
    plot_list[[i]] <- experiment_outcome[[data_num]]$estimates[[estimator_num]]$field_plot
  }
  
  field_plots <- cowplot::plot_grid(plotlist = plot_list, ncol = n_col)
  print(field_plots)
}

get_shared_ic <- function(ls_samples){
  # samples random ICs from within, on, outside close, and outside far from the LC
  
  # sample random LC points
  ic_init_samples <- ls_samples[runif(4, 1, nrow(ls_samples)),c(1,2)]
  ic_radial <- c(0.5, 1, 1.25, 1.4)
  ic_samples <- ic_init_samples * ic_radial
  
  # format
  rownames(ic_samples) <- c("Within", "On", "OutsideClose", "OutsideFar")
  colnames(ic_samples) <- c("x", "y")
  
  return(ic_samples)
}

plot_single_solution_path <- function(solution_path_list, lc_data){

  within_tibble <- tibble(x = solution_path_list$Within$x, y = solution_path_list$Within$y)
  on_tibble <- tibble(x = solution_path_list$On$x, y = solution_path_list$On$y)
  outside_close_tibble <- tibble(x = solution_path_list$OutsideClose$x, y = solution_path_list$OutsideClose$y)
  outside_far_tibble <- tibble(x = solution_path_list$OutsideFar$x, y = solution_path_list$OutsideFar$y)

  gradient_tibble <- tibble(x = c(), y = c())
  # TODO: Plot field
  #gradient_tibble <- cbind(grid$eval_grid, estimators[[i]]$estimated_field)
  #colnames(gradient_tibble) <- c("x","y","u","v")
  
  path_plot <- ggplot(as_tibble(lc_data),aes(x = x, y = y)) +
    geom_point(data = lc_data, aes(x = x, y = y), color = "#FDB515") +
    geom_path(data = within_tibble, aes(x = x, y = y), color = "red") +
    geom_point(x = as.numeric(within_tibble[1,1]), y =as.numeric(within_tibble[1,2]), color = "red", size = 3) + 
    geom_path(data = on_tibble, aes(x = x, y = y), color = "blue") +
    geom_point(x = as.numeric(on_tibble[1,1]), y =as.numeric(on_tibble[1,2]), color = "blue", size = 3) +
    geom_path(data = outside_close_tibble, aes(x = x, y = y), color = "green") +
    geom_point(x = as.numeric(outside_close_tibble[1,1]), y =as.numeric(outside_close_tibble[1,2]), color = "green", size = 3) + 
    geom_path(data = outside_far_tibble, aes(x = x, y = y), color = "purple") +
    geom_point(x = as.numeric(outside_far_tibble[1,1]), y =as.numeric(outside_far_tibble[1,2]), color = "purple", size = 3) #+
    #labs(title=paste(estimators[[i]]$method, 
    #                 paste(names(estimators[[i]]$params), estimators[[i]]$params,
    #                       sep = ":", collapse = ",")))
  
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
    
    solution_path_plot <- plot_single_solution_path(ic_results_list, experiment_outcome[[data_num]]$limit_cycle_tail)
    solution_graphs[[i]] <- solution_path_plot
  }

  n_col = min(3, length(plot_specifications))
  sol_plot <- cowplot::plot_grid(plotlist = solution_graphs, ncol = n_col)
  print(sol_plot)
}



visualize_results <- function(experiment_outcome, plot_list){
  
  for (plot in plot_list){
    if (plot$type == "field"){
      plot_estimated_fields(plot$experiments, experiment_outcome)
    }  
    else if (plot$type == "solution_path"){
      plot_solution_paths(plot$experiments, experiment_outcome)
    }
    else {
      print("Plot not implemented")
    }
  }
} 

#############
## Testing ##
#############

# specify data
no_noise <- list(name = "test", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = 0, var_y = 0,
                 lc_tail_n = 700, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5)
some_noise <- list(name = "noisy_y", system = "van_der_pol", params = list(mu = 1.5), n = 1000, sample_density = 0.1, var_x = 1, var_y = 1,
                   lc_tail_n = 700, x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5)
experiment_list <- list(no_noise, some_noise)
experiment_data <- generate_data_object(experiment_list)

# specify estimators
truth <- list(method = "truth",  params = list())
one_nn <- list(method = "knn",  params = list(k = 1))
many_nn <- list(method = "knn",  params = list(k = 400))
nw_base <- list(method = "nw",  params = list(h = 0.01))
loess <- list(method = "loess",  params = list(h = 0.1))
spline1 <- list(method = "spline",  params = list(lambda = 1e-16))
spline2 <- list(method = "spline",  params = list(lambda = 1e-8))
spline3 <- list(method = "spline",  params = list(lambda = 1))
# pick estimators to run on each data set
experiment_estimators <- list(truth, spline1, loess)

experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)

plot1 <- list(type = "field", experiments = list(c(data = 1, estimator = 3), c(data = 2, estimator = 2)))
plot2 <- list(type = "solution_path", experiments = list(c(data = 1, estimator = 1), c(data = 1, estimator = 2), c(data = 1, estimator = 3)))
to_plot <- list(plot1, plot2)
visualize_results(experiment_results, to_plot)

