##
## Data Generation
##

# TODO: Add observation noise? 
# Should this been done uniquely for each system or as a general function?

get_abhi_data <- function(){
  # Returns Abhi's simulated data
  # 
  ## Inputs:
  # None
  #
  ## Outputs:
  # sampled_data (data.frame): 1000 samples [x,y,f_x,f_y]
  
  # load data from csv into df
  x_df = read.csv("Saved_Data/Abhi/x.csv", header = FALSE, col.names = c("x"))
  y_df = read.csv("Saved_Data/Abhi/y.csv", header = FALSE, col.names = c("y"))
  fx_df = read.csv("Saved_Data/Abhi/f_x.csv", header = FALSE, col.names = c("f_x"))
  fy_df = read.csv("Saved_Data/Abhi/f_y.csv", header = FALSE, col.names = c("f_y"))
  
  # convert to consistent numeric encoding
  x =  sapply(x_df[, c(1)], as.numeric)
  y =  sapply(y_df[, c(1)], as.numeric)
  
  f_x =  sapply(fx_df[, c(1)], as.numeric)
  f_y =  sapply(fy_df[, c(1)], as.numeric)
  
  # return as data.frame
  sampled_data <- data.frame(x = x,y = y,f_x, f_y)
  return(sampled_data)
  
}

generate_van_der_pol <- function(params, num_samples = 1000){
  # function which returns samples from a Van der Pol system specified by \mu
  # initial condition is randomly selected
  # 
  ## Inputs:
  # params (numeric): mu parameter for Van der pol system
  # num_samples (integer)[optional]: number of samples to generate
  #
  ## Outputs:
  # sampled_data (data.frame): num_samples samples; columns in [x,y]
  
  # the only parameter for this model is mu
  if (length(params)>1){stop('Error: Too many Van Der Pol parameters')}
  mu <- params[1]
  
  # randomly pick the IC in [-10, 10] x [-10, 10]
  initial_condition <- runif(2, min = -10, max = 10)
  names(initial_condition) <- c('x','y')
  
  # compute the prediction and its gradient at each sample
  traj <- data.frame(lsoda(initial_condition,seq(0,num_samples,by=0.1), eval_van_der_pol_gradient, mu))
  derivative_evals <- t(apply(traj[,2:3], 1, eval_van_der_pol_gradient, t = 0, mu = mu, list = T))
  traj['f_x'] <- derivative_evals[,1] # TODO: is there a way to get the gradient with `lsoda` as opposed to manually?
  traj['f_y'] <- derivative_evals[,2]
  
  return(traj)
  
}

eval_van_der_pol_gradient <- function(t, v, mu, list = F){
  # helper function to evaluate gradient of van der pol at point v = (x,y)
  x <- v[1]
  y <- v[2]
  
  dx <- mu*(x - x^3/3 - y)
  dy <- (1/mu)*x
  
  ## this form allows sampling when mu = 0 <-> uniform circular motion
  ## however, it doesn't work with `lsoda` for unknown reasons
  #dx <- y
  #dy <- mu*(1-x^2)*y - x
  
  # different return format if used in `lsoda` or not
  if(list){return(c(dx,dy))}
  return(list(as.vector(c(dx,dy))))
}