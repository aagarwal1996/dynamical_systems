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

#' This is a wrapper for data generation under a number of DS models
#' See `data_generation.R` for individual models
#'
#' @param model (string): name of the model to sample from the limit cycle of
#'     implemented: "van_der_pol", "abhi"   
#'     to-implement: asymmetrically sampled unit circle, others... 
#' @param params (vector): vector of parameters for the specified model
#'     Will vary my model, see `data_generation.R` for specific cases
#' @param num_samples (integer)[optional] Number of samples to generate, 
#'     starting from some initial condition
#' @param sample_density (numeric)[optional] Delta t time-step used in `lsoda`
#' @param add_obs_noise (logical)[optional]: whether to add random noise to each observation
#'     TODO: implement, will involve additional parameters
#' @param use_seed (logical)[optional]: whether to set a seed before drawing data
#' @param save_csv (logical)[optional]: whether to save the CSV of the data
#'
#' @return sampled_data (data.frame): `num_samples` sampled points from the limit cycle
#'     columns in c(x,y,f_x,f_y)
#'
#' @examples
#' generate_limit_cycle_data("van_der_pol", c(.5)) # non-stiff VdP
#' generate_limit_cycle_data("abhi", c()) # Abhi's .csv data with no params
generate_limit_cycle_data <- function(model, params, 
                                      num_samples = 1500, sample_density = 0.1,
                                      add_obs_noise = F,
                                      use_seed = F, save_csv = F){  
  if (use_seed){set.seed(2022)}
  
  if (model == "van_der_pol"){
    sampled_data <- generate_van_der_pol(params,num_samples=num_samples,sample_density=sample_density)
  }
  else if (model == "asymmetric_circle"){
    sampled_data <- generate_asymmetric_circle(params,num_samples=num_samples,sample_density=sample_density)
  }
  else if(model == "abhi"){
    sampled_data <- get_abhi_data()
  }
  else{
    stop('Error: ', model, ' model not implemented.')
  }
  
  if (save_csv){
    file_name <- paste0("Saved_Data/",model,"-",format(Sys.time(), "%m_%d_%Y-%H_%M_%S"),".csv")
    write_csv(sampled_data, file_name)
  }
  
  return(sampled_data)
}

##################################
## Gradient Field Approximation ##
##################################

evaluate_gradient_methods <- function(data_list, estimators_list, 
                                      lc_tail_n = 700, 
                                      x_grid_size = 24, y_grid_size = 24, extrapolation_size = 0.5){
  
  # truncate LC samples
  data_list$limit_cycle_tail <- tail(data_list$limit_cycle_samples, n = lc_tail_n)
  
  # build list which holds information about eval grid
  grid_object <- list(xmin = min(data_list$limit_cycle_tail$x),
                      xmax = max(data_list$limit_cycle_tail$x),
                      ymin = min(data_list$limit_cycle_tail$y),
                      ymax = max(data_list$limit_cycle_tail$y)) 
  
  grid_object$x_grid = seq(floor(grid_object$xmin) - extrapolation_size , ceiling(grid_object$xmax) + extrapolation_size, len = x_grid_size)
  grid_object$y_grid = seq(floor(grid_object$ymin) - extrapolation_size, ceiling(grid_object$ymax) + extrapolation_size, len = y_grid_size)
  grid_object$eval_grid = unname(as.matrix(expand.grid(grid_object$x_grid,grid_object$y_grid)))
    
  for (i in 1:length(estimators_list)){

    if (estimators_list[[i]]$method == "spline"){
      spline_result <- get_gradient_field(data_list, grid_object, estimators_list[[i]])
      estimators_list[[i]]$estimated_field = spline_result$field
      estimators_list[[i]]$sfd_list = spline_result$sfd_list
    }
    
    else{
      estimators_list[[i]]$estimated_field = get_gradient_field(data_list, grid_object, estimators_list[[i]])
    }
    estimators_list[[i]]$field_plot = ggplot_field(data_list$limit_cycle_tail, grid_object$eval_grid, 
                 estimators_list[[i]]$estimated_field, rescale = .025,
                 title=paste(estimators_list[[i]]$method, 
                             paste(names(estimators_list[[i]]$params), estimators_list[[i]]$params,
                                   sep = ":", collapse = ",")))
  }
  all_fields <- plot_estimated_fields(data_list, estimators_list)
  print(all_fields)
  
  all_paths <- plot_solution_paths(data_list, grid_object, estimators_list)
}

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

plot_estimated_fields <- function(data, estimator_results){
  n_col = min(3, length(estimator_results))
  plot_list <- list()
  for (i in 1:length(estimator_results)){
    plot_list[[i]] <- estimator_results[[i]]$field_plot
  }

  field_plots <- cowplot::plot_grid(plotlist = plot_list, ncol = n_col)
  return(field_plots)
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

plot_single_solution_path <- function(data, grid, estimators, i){
  solution_path_list <- estimators[[i]]$solution_paths
  
  within_tibble <- tibble(x = solution_path_list$Within$x, y = solution_path_list$Within$y)
  on_tibble <- tibble(x = solution_path_list$On$x, y = solution_path_list$On$y)
  outside_close_tibble <- tibble(x = solution_path_list$OutsideClose$x, y = solution_path_list$OutsideClose$y)
  outside_far_tibble <- tibble(x = solution_path_list$OutsideFar$x, y = solution_path_list$OutsideFar$y)
  gradient_tibble <- cbind(grid$eval_grid, estimators[[i]]$estimated_field)
  colnames(gradient_tibble) <- c("x","y","u","v")
  
  # TODO: Plot field
  path_plot <- ggplot(as_tibble(gradient_tibble),aes(x = x, y = y)) +
    geom_point(data = data$limit_cycle_tail, aes(x = x, y = y), color = "#FDB515") +
    geom_path(data = within_tibble, aes(x = x, y = y), color = "red") +
    geom_point(x = as.numeric(within_tibble[1,1]), y =as.numeric(within_tibble[1,2]), color = "red", size = 3) + 
    geom_path(data = on_tibble, aes(x = x, y = y), color = "blue") +
    geom_point(x = as.numeric(on_tibble[1,1]), y =as.numeric(on_tibble[1,2]), color = "blue", size = 3) +
    geom_path(data = outside_close_tibble, aes(x = x, y = y), color = "green") +
    geom_point(x = as.numeric(outside_close_tibble[1,1]), y =as.numeric(outside_close_tibble[1,2]), color = "green", size = 3) + 
    geom_path(data = outside_far_tibble, aes(x = x, y = y), color = "purple") +
    geom_point(x = as.numeric(outside_far_tibble[1,1]), y =as.numeric(outside_far_tibble[1,2]), color = "purple", size = 3) +
    labs(title=paste(estimators[[i]]$method, 
                paste(names(estimators[[i]]$params), estimators[[i]]$params,
                      sep = ":", collapse = ",")))
 
  return(path_plot)
}

plot_solution_paths <- function(data, grid, estimators, 
                                samples_per_path = 500){
  
  shared_ic <- get_shared_ic(data$limit_cycle_tail)
  
  sample_tibble <- tibble(x = data$limit_cycle_tail[,1], y = data$limit_cycle_tail[,2], u = NA, v = NA)
  ic_tibble <- tibble(x = shared_ic[,1], y = shared_ic[,2], u = NA, v = NA)
  
  for (i in 1:length(estimators)){
    
    # get the parameters for the appropriate method
    
    estimator_name <- estimators[[i]]$method
    estimator_params <- estimators[[i]]$params
    
    results_list <- NA
    for (ic_type in rownames(shared_ic)){
      
      if (typeof(results_list) == "logical"){
        results_list <- list(generate_solution_path(data, grid, estimators[[i]], shared_ic[ic_type,]))
      }
      else{
        results_list <- append(results_list,  list(generate_solution_path(data, grid, estimators[[i]], shared_ic[ic_type,])))
      }
    }
    names(results_list) <- rownames(shared_ic)
    estimators[[i]]$solution_paths <- results_list
    
    estimators[[i]]$solution_graph <- plot_single_solution_path(data, grid, estimators, i)
  }
  
  n_col = min(3, length(estimators))
  plot_list <- list()
  for (i in 1:length(estimators)){
    plot_list[[i]] <- estimators[[i]]$solution_graph
  }
  
  sol_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = n_col)
  print(sol_plot)
  
  return()
}

nonstiff_mu <- 20
truth_nonstiff <- list(method = "truth",  params = list(mu = nonstiff_mu))
one_nn <- list(method = "knn",  params = list(k = 1))
nw_base <- list(method = "nw",  params = list(h = 0.01))
many_nn <- list(method = "knn",  params = list(k = 400))
loess <- list(method = "loess",  params = list(h = 0.1))
zero_spline <- list(method = "spline",  params = list(lambda = 1e-16))
spline <- list(method = "spline",  params = list(lambda = 1e-8))
spline2 <- list(method = "spline",  params = list(lambda = 1))
spline3 <- list(method = "spline",  params = list(lambda = 5))
spline_bigpen <- list(method = "spline",  params = list(lambda = 30))
#nonstiff_estimator_list <- list(truth_nonstiff, one_nn, nw_base, many_nn, loess, spline)
nonstiff_estimator_list <- list(truth_nonstiff, zero_spline, spline)

 
#################
## Evaluation ##
################

nonstiff_mu <- 20
vp_nonstiff_samples <- generate_limit_cycle_data("van_der_pol", c(nonstiff_mu))
vp_nonstiff_list <- list(system = "van_der_pol", params = list(mu = nonstiff_mu), limit_cycle_samples = vp_nonstiff_samples)
evaluate_gradient_methods(vp_nonstiff_list, nonstiff_estimator_list)
