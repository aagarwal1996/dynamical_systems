##
## Imports
##
library(fda)
library(deSolve)
library(tidyverse)
library(mvtnorm) # for Multivariate Normal density estimates in KDE
library(shape) # for pretty arrows
library(ggquiver) # for vector field plots
library(ggpubr) # to generate side-by-side plots

source(here::here('Data_Generation/data_generation.R')) # functions to generate data from specific DS are in another file
source(here::here('Estimation_Methods/bspline.R'))
library(Rcpp)
sourceCpp(here::here('Estimation_Methods/nw_regression.cpp'))
sourceCpp(here::here('Estimation_Methods/loess.cpp'))
sourceCpp(here::here('Estimation_Methods/bspline.cpp'))
##
## Data Generation
##

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

generate_grid_data <- function(model, params, eval_grid){
  # This is a wrapper for data generation under a number of DS models
  # See `data_generation.R` for individual models
  #
  ## Inputs:
  # model (string): name of the model to sample from the limit cycle of
  #   implemented: van_der_pol, abhi
  #   to-implement: unstable_spiral, morris_lecar
  # params (vector): vector of parameters for the specified model
  # eval_grid (matrix): n x 2 matrix of grid points to calcualte gradient at
  #
  ## Outputs:
  # grid_gradient (data.frame): sampled grid from the gradient field; columns in [x,y,f_x,f_y]
  
  if (model == "van_der_pol"){
    grid_gradient <- t(apply(eval_grid, 1, van_der_pol_gradient_helper, mu = params[1]))
  }
  else if (model == "sinusoidal_van_der_pol"){
    grid_gradient <- t(apply(eval_grid, 1, van_der_pol_gradient_helper, mu = params[1], sinusoid = T))
  }
  else if(model == "unstable_spiral"){
    stop('Error: Unstable sprial not yet implemented.')
  }
  else if(model == "abhi"){
    # true field is unknown
    grid_gradient <- matrix(0, nrow = nrow(eval_grid), ncol = ncol(eval_grid))
  }
  else{
    stop('Error: ', model, ' model not implemented.')
  }
  
  return(grid_gradient)
}

##
## Gradient Field Approximation
##

evaluate_gradient_methods <- function(data, tail_n = 700, 
                            x_grid_size = 24, y_grid_size = 24, extrapolation_size = 2,
                            nw_bandwidth = 0.1, loess_bandwidth = 0.1,
                            model_title = "",
                            model_str = "", method_params = NA){
  # Function which runs multiple methods to estimate the gradient field
  #   implemented: splines and kernel regression
  #   to-implement: local linear regression
  #
  ## Inputs:
  # data (data.frame): contains data with columns in [x, y, f_x, f_y]
  # tail_n (integer)[optional]: we only use the last `tail_n` samples of the data, ideally from the limit cycle
  # x_grid_size (integer)[optional]: size of x-axis grid to extrapolate over
  # y_grid_size (integer)[optional]: size of y-axis grid to extrapolate over
  # extrapolation_size (numeric)[optional]: grid will be from max and min of observations +/- this variable
  #
  # Outputs:
  # None (Displays various plots)

  limit_cycle.data <- tail(data, n = tail_n)
  
  # get extreme points of the data
  xmin <- floor(min(limit_cycle.data$x))
  xmax <- ceiling(max(limit_cycle.data$x))
  ymin <- floor(min(limit_cycle.data$y))
  ymax <- ceiling(max(limit_cycle.data$y))
  
  # grid to evaluate gradient extrapolation over
  x_grid <- seq(xmin - extrapolation_size , xmax + extrapolation_size, len = x_grid_size)
  y_grid <- seq(ymin - extrapolation_size, ymax + extrapolation_size, len = y_grid_size)
  eval_grid <- unname(as.matrix(expand.grid(x_grid,y_grid)))
  
  # get truth over evaluation grid
  true_field <- generate_grid_data(model_str, c(2),eval_grid)
  true_field_plot <- ggplot_field(limit_cycle.data, eval_grid, true_field, title="Truth")
  
  # Splines
  spline_result <- spline_gradient(limit_cycle.data, x_grid, y_grid)
  spline_result_toplot <- matrix(c(c(spline_result$x_grad_eval),c(spline_result$y_grad_eval)),ncol=2)
  spline_field_plot <- ggplot_field(limit_cycle.data, eval_grid, spline_result_toplot, title="Spline")
  ## Uncomment line below for original plots from Abhi's .Rmd
  #generate_spline_plots(data, spline_result, x_grid, y_grid)

  # convert data to matrix (no labels) for CPP
  limit_cycle.data.matrix <- unname(as.matrix(limit_cycle.data))
  # NW
  nw_bw = .1
  nw_bw_matrix <- nw_bw*diag(2) # TODO: Vary bandwidth (okay as is in no noise setting)
  NW_regression_result <- NW_regression_cpp(eval_grid, limit_cycle.data.matrix, nw_bw_matrix)
  nw_field_plot <- ggplot_field(limit_cycle.data, eval_grid, NW_regression_result, title="NW Kernel Regression")
  nw_gradient_plot <- plot_gradient_path_multi(limit_cycle.data.matrix, eval_grid, NW_regression_result, nw_bw, title = paste(model_title, "lsoda Solutions Along NW Gradient Field", sep = ": "))
  
  # LOESS
  loess_bw = .1
  loess_fit <- eval_loess_fit(eval_grid, limit_cycle.data.matrix, loess_bw)
  loess_field_plot <- ggplot_field(limit_cycle.data, eval_grid, loess_fit, title="LOESS")
  loess_gradient_plot <- plot_gradient_path_multi(limit_cycle.data.matrix, eval_grid, 
                                                  loess_fit, loess_bw, loess = T, title = paste(model_title, "lsoda Solutions Along LOESS Gradient Field", sep = ": "))
  
  # Output plots
  field_plots <- ggarrange(true_field_plot, nw_field_plot, loess_field_plot, ncol=3, nrow=1) # spline_field_plot
  field_plots <- annotate_figure(field_plots, top = text_grob(label = model_title))
  print(field_plots)
  print(nw_gradient_plot)
  print(loess_gradient_plot)
}

##
## Spline fitting
##

# See: `bspline.R` and `bspline.cpp`

##
## Kernel regression
##

# See: `nw_regression.cpp`

##
## Local linear regression
##

# See: `loess.cpp`

##
## Plotting functions
##

# TODO: Add better documentation and move to separate file

plot_observations <- function(data){
  # helper function which generates a scatter plot of observations 
  ggplot(data,aes(x=x,y=y)) +
    geom_point() + 
    labs(x="x",y="y",title="Observations from a Dynamical System")
}

ggplot_field <- function(ls_samples, eval_grid, gradient_data,title=""){
  # plot a gradient field

  sample_tibble <- tibble(x = ls_samples[,1], y = ls_samples[,2], u = NA, v = NA)
  gradient_tibble <- tibble(x = eval_grid[,1], y = eval_grid[,2], u = gradient_data[,1], v = gradient_data[,2])
  
  field_plot <- ggplot(gradient_tibble, aes(x = x, y = y, u = u, v = v)) +
    geom_quiver(color = "#003262") +
    geom_point(data = sample_tibble, aes(x = x, y = y), color = "#FDB515") +
    labs(title = title)
  
  return(field_plot)
}

plot_gradient_path_multi <- function(ls_samples, eval_grid, gradient_data, bw, loess = F, num_reps = 9, title =""){
  
  plot_list <- list()
  sample_tibble <- tibble(x = ls_samples[,1], y = ls_samples[,2], u = NA, v = NA)
  gradient_tibble <- tibble(x = eval_grid[,1], y = eval_grid[,2], u = gradient_data[,1], v = gradient_data[,2])
  

  for (i in 1:num_reps){
      if (loess){
        trajectory <- generate_loess_path(ls_samples, bw, num_samples = 200)
      }
      else{
        trajectory <- generate_nw_path(ls_samples, bw, num_samples = 200)
      }
      trajectory_tibble <- tibble(x = trajectory[,1], y = trajectory[,2], u = NA, v = NA, run = i)
      path_plot <- ggplot(gradient_tibble, aes(x = x, y = y, u = u, v = v)) +
        geom_quiver(color = "#003262") +
        geom_point(data = sample_tibble, aes(x = x, y = y), color = "#FDB515") +
        geom_path(data = trajectory_tibble, aes(x = x, y = y)) +
        geom_point(x = as.numeric(trajectory_tibble[1,1]), y =as.numeric(trajectory_tibble[1,2]), color = "red", size = 3) +
        facet_wrap(~run)
      plot_list <- append(plot_list, list(path_plot))
  }
  
  path_plots <- ggarrange(plotlist = plot_list, nrow = ceiling(length(plot_list)/3), ncol = 3)
  path_plots <- annotate_figure(path_plots, top = text_grob(label = title))
  return(path_plots)
  
}


##
## Evaluation
##

abhi_data <- generate_limit_cycle_data("abhi", c())
evaluate_gradient_methods(abhi_data, extrapolation_size = 1, model_title = "Abhi's Van der Pol", model_str = "abhi", method_params = c())

#vp_data <- generate_limit_cycle_data("van_der_pol", c(20))
#evaluate_gradient_methods(vp_data, extrapolation_size = 1, model_title = "Van der Pol; mu = 20", model_str = "van_der_pol", method_params = c(20))
 
#vp_data <- generate_limit_cycle_data("van_der_pol", c(.5))
#evaluate_gradient_methods(vp_data, extrapolation_size = 1, model_title = "Van der Pol; mu = 0.5",model_str = "van_der_pol", method_params = c(0.5))

# svp_data <- generate_limit_cycle_data("sinusoidal_van_der_pol", c(2.05))
# evaluate_gradient_methods(svp_data, extrapolation_size = 1)