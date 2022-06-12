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
    data_list[[i]]$data <- generate_limit_cycle_data(system_name, system_params, 
                                   var_x = var_x, var_y = var_y, sample_seed = sample_seed,
                                   num_samples = num_samples, sample_density = sample_density,
                                   save_csv = save_csv)
  }
  
  return(data_list)
}

# sample_experiment <- list(name = "test", system = "van_der_pol", params = c(1.5), n = 1000, sample_density = 0.1, var_x = 0, var_y = 1)
# y_noisy <- list(name = "noisy_y", system = "van_der_pol", params = c(1.5), n = 1000, sample_density = 0.1, var_x = 0, var_y = 100)
# experiment_list <- list(sample_experiment, y_noisy)
# test_result <- generate_data_object(experiment_list)
# View(test_result)