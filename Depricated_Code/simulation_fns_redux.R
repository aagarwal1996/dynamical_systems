#############
## Imports ##
#############

# Packages
library(fda)
library(tidyverse)
library(mvtnorm) # for Multivariate Normal density estimates in KDE
library(ggquiver) # for vector field plots
library(shape) # for pretty vector arrows in plots
library(ggpubr) # to generate side-by-side plots
library(fields) # for thin-plate splines # TODO: Remove?

# Local Files
# functions to generate data from specific system
source(here::here('Data_Generation','data_generation.R'))
# wrappers for fast Rcpp fit of gradient fields
source(here::here('Estimation_Methods','bspline.R'))
# eigenvalue-penalized gradient descent fit
source(here::here('Estimation_Methods','bspline_gd.R'))
source(here::here('Estimation_Methods','knn.R'))
library(Rcpp)
sourceCpp(here::here('Estimation_Methods','nw_regression.cpp'))
sourceCpp(here::here('Estimation_Methods','loess.cpp'))

#####################
## Data Generation ##
#####################

#' Generate data object for parameterized systems
#' 
#' This function generates samples from specified models and returns a structured list 
#'
#' @param data_parameters (list of lists) Each sublist contains the parameterization details for a single data draw
#'
#' @return data_list (list of lists) Each sublist contains samples
#'
#' @examples
#' experiment_details <- list(
#'     list(experiment_name = "Example", system_name = "van_der_pol",
#'         system_params = list(mu = 1.5), var_x = 0.1, var_y = 0.1,
#'         data_seed = 2, noise_seed = 2023, num_samples = 1000, sample_density = 0.1,
#'         smooth_type = "bspline", smooth_lambda = 1e-12, num_basis_fn = 48)
#' )
#' example_data <- generate_data_object_model(experiment_details)
generate_data_object_model <- function(data_parameters){

	data_list <- data_parameters # copy to modify
	# iterate over each requested dataset
	for (i in 1:length(data_list)){
		experiment_name <- data_list[[i]]$experiment_name
		# extract sampling parameters. See `generate_limit_cycle_data()` for details
		system_name <- data_list[[i]]$system_name
		system_params <- data_list[[i]]$system_params
		var_x <- data_list[[i]]$var_x
		var_y <- data_list[[i]]$var_y
		data_seed <- data_list[[i]]$data_seed
		noise_seed <- data_list[[i]]$noise_seed
		num_samples <- data_list[[i]]$num_samples
		sample_density <- data_list[[i]]$sample_density
		smooth_type <- data_list[[i]]$smooth_type
		smooth_lambda <- data_list[[i]]$smooth_lambda
		num_basis_fn <- data_list[[i]]$num_basis_fn
		
		# generate data.frame of (noisy) samples
		data_list[[i]]$limit_cycle_samples <- generate_limit_cycle_data(system_name, system_params,
			var_x, var_y, data_seed, noise_seed,
			num_samples, sample_density,
			smooth_type, smooth_lambda, num_basis_fn,
			plot_title = experiment_name)
	}
	return(data_list)
}

#' Generate limit cycle data for a given system
#'
#' This function returns (noisy) samples along one trajectory of a system parameterized by the inputs.
#'
#' @param system_name (string) The name of the system to generate data from
#' @param system_params (list of numeric) A list setting all parameters for the system
#' @param var_x (numeric) The variance of white noise added to the x-coordinate of samples
#' @param var_y (numeric) The variance of white noise added to the y-coordinate of samples
#' @param data_seed (integer) A seed for generating samples from the system
#' @param noise_seed (integer) A seed for positional noise added to the samples
#' @param num_samples (integer) The number of samples along a solution path to generate
#' @param sample_density (integer) The time elapsed between each sample
#' @param smooth_type (string; "bspline","tps", or "orig") Type of smoothing to apply to noisy samples
#' @param smooth_lambda (numeric) Smoothing penalty weight, if applicable
#' @param num_basis_fn (integer) The number of spline basis functions for each coordinate
#' @param plot_title (string) A title for a plot of original and smoothed samples
#'
#' @return sampled_data (data.frame) Samples of gradient and (noisy) position
#'
#' @examples
#' vdp_df <- generate_limit_cycle_data("van_der_pol", list(mu=1.5), 0, 0)
generate_limit_cycle_data <- function(system_name, system_params, var_x, var_y,
	data_seed = 2022, noise_seed = 302,
	num_samples = 1500, sample_density = 0.1,
	smooth_type = "bspline", smooth_lambda = 1e-12, num_basis_fn = 48,
	plot_title = ""){

	# TODO: Add gradient noise; include time as a covariate in returned data.frame
	# TODO: Is it useful to return the original and noisy samples?
	# get samples without noise
	set.seed(data_seed)
	if (system_name == "van_der_pol"){
		sampled_data <- generate_van_der_pol_samples(system_params, num_samples = num_samples,
			sample_density = sample_density)
	} else if (system_name == "rzma"){
		sampled_data <- generate_rzma_samples(system_params, num_samples = num_samples,
			sample_density = sample_density)
	} else if(system_name == "lotka_volterra"){
		sampled_data <- generate_lv_samples(system_params, num_samples = num_samples,
			sample_density = sample_density)
	} else if(system_name == "log_lotka_volterra"){
		sampled_data <- generate_log_lv_samples(system_params, num_samples = num_samples,
			sample_density = sample_density)
	} else if(system_name == "abhi"){
		sampled_data <- get_abhi_data()
	} else{
		stop('Error: ', system_name, ' system not implemented.')
	}

	if ((var_x != 0) | (var_y != 0)){
		# add Gaussian noise to position of observations
		set.seed(noise_seed)
		noise_matrix <- mvrnorm(nrow(sampled_data), c(0,0,0,0), diag(c(var_x,var_y,0,0)))
		noisy_samples <- sampled_data + noise_matrix

		# apply smoothing to the noisy positions
		smoothed_positions <- spline_smooth_noisy_samples(sampled_data, noisy_samples,
			nbasis = num_basis_fn, max_t = num_samples*sample_density,
			lambda = smooth_lambda, smooth_type = smooth_type, title = plot_title)
		sampled_data <- cbind(smoothed_positions, noisy_samples[,c(3,4)])
	}

	colnames(sampled_data) <- c("x","y","f_x","f_y")
	return(sampled_data)
}

#' Fill skeleton data object with user-provided data
#' 
#' This allows for interface with other simulation functions
#' For now this does not allow additional noise or smoothing
#' 
#' If a data.frame of samples with columns in ("x","y","f_x","f_y") is available,
#' this function generates a list with that data and the slots from the output of 
#' `generate_data_object_model()` 
#'
#' @param experiment_list (list of lists) Each expected to have a experiment name, system name, and `limit_cycle_samples` data.frame
#'
#' @return data_list (list of lists) Fit into a consistent format for analysis
#'
#' @examples
#' user_data <- read.csv(example.csv) # loads data.frame with appropriate columns
#' experiment_details <- list(
#'     list(experiment_name = "Example", system_name = "My System", limit_cycle_samples = user_data)
#' )
generate_data_object_obs <- function(experiment_list){

	# TODO: Allow for user to add their own noise and smoothing; check format of user data
	data_list <- experiment_list

	for (i in 1:length(data_list)){
		experiment_name <- data_list[[i]]$experiment_name
		# extract sampling parameters. See `generate_limit_cycle_data()` for details
		system_name <- data_list[[i]]$system_name
		system_params <- NA
		var_x <- NA
		var_y <- NA
		data_seed <- NA
		noise_seed <- NA
		num_samples <- NA
		sample_density <- NA
		smooth_type <- NA
		smooth_lambda <- NA
		num_basis_fn <- NA
	}

	return(data_list)
}

#######################
## Estimator Fitting ##
#######################

get_gradient_predictor(data_object,estimator_object){

	if (estimator_object$method == "truth"){
		gradient_predictor <- return_true_gradient_fn(estimator_object$system_name)
	} else if (estimator_object$method == "spline"){
		gradient_predictor <- return_spline_gradient(data_object, estimator_object)
	} else {
		stop('Error: ', estimator_object$method, ' estimator not implemented.')
	}
	return(gradient_predictor)
}
