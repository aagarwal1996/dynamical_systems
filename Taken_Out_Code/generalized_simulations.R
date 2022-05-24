library(fda)
library(ggplot2)
source(here::here('data_generation.R')) # functions to generate data are in another file
## Data Generation

get_abhi_data <- function(){
  # function which returns Abhi's simulated data
  # 
  ## Inputs:
  # None
  #
  ## Outputs:
  # sampled_data (data.frame): 1000 samples [x,y,f_x,f_y]
  
  # load data from csv into df
  x_df = read.csv("x.csv", header = FALSE, col.names = c("x"))
  y_df = read.csv("y.csv", header = FALSE, col.names = c("y"))
  fx_df = read.csv("f_x.csv", header = FALSE, col.names = c("f_x"))
  fy_df = read.csv("f_y.csv", header = FALSE, col.names = c("f_y"))
  
  # convert to consistent numeric encoding
  x =  sapply(x_df[, c(1)], as.numeric)
  y =  sapply(y_df[, c(1)], as.numeric)
  
  f_x =  sapply(fx_df[, c(1)], as.numeric)
  f_y =  sapply(fy_df[, c(1)], as.numeric)
  
  # return as data.frame
  sampled_data <- data.frame(x = x,y = y,f_x, f_y)
  return(sampled_data)
  
}

get_van_der_pol_grad <- function(mu, x, y, subsitution_form = F){
  # evaluate the gradient of the Van der Pol oscillator at (x,y) for given mu
  #
  ## Inputs:
  # mu (numeric): mu parameter of DS
  # x (numeric): x coord
  # y (numeric): y coord
  # substitution_form (logical)[optional]: use a form with x = y'; avoids degeneracy issues near mu = 9
  #
  ## Outputs:
  # grad (numeric): list of (f_x, f_y) evaluated at input
  
  # equations from Wikipedia
  if(subsitution_form){
    f_x <- y
    f_y <- x #mu*(1 - x^2)*y - x
  }
  else{
    f_x <- mu*(x - (1/3)*x^3 - y)
    f_y <- (1/mu)*x
  }
  grad <- c(f_x,f_y)
  return(grad)
}

generate_van_der_pol <- function(params, num_samples = 5000, 
                                 x_step_size = 0.002, y_step_size = 0.01){
  # function which returns Euler steps from a specified Van der Pol system
  # initial condition is randomly selected
  # 
  ## Inputs:
  # params (numeric): mu parameter for Van der pol system
  # num_samples (integer)[optional]: number of smaples to generate
  # x_step_size (numeric)[optional]: x_new = x_old + x_step_size + f_x_old
  # y_step_size (numeric)[optional]: ibid.
  #
  ## Outputs:
  # sampled_data (data.frame): num_samples samples; columns in [x,y,f_x,f_y]
  
  # TODO: Add noise?
  # the only parameter for this model is mu
  if (length(params)>1){stop('Error: Too many Van Der Pol parameters')}
  mu <- params[1]
  
  # randomly pick the IC in [-10, 10] x [-10, 10]
  initial_condition <- runif(2, min = -10, max = 10)
  
  # init output
  x <- c(initial_condition[1])
  y <- c(initial_condition[2])
  f_x <- c()
  f_y <- c()
  
  # iterate to generate more data
  for (i in 1:num_samples){
    grad <- get_van_der_pol_grad(mu, x[i],y[i], subsitution_form = T)
    x <- c(x, x[i] + x_step_size * grad[1])
    y <- c(y, y[i] + y_step_size * grad[2])
    f_x <- c(f_x, grad[1])
    f_y <- c(f_y, grad[2])
  }

  sampled_data <- data.frame(x = head(x,num_samples),y = head(y,num_samples),f_x,f_y)
  return(sampled_data)
}

generate_spiral <- function(num_samples = 5000, 
                                 x_step_size = 0.002, y_step_size = 0.01){
  initial_condition <- runif(2, min = -1, max = 1)
  # function which generates data from a DS Giles described. This spirals
  # inwards if IC within unit circle; otherwise outward. Due to degeneracy 
  # issues, IC always within unit cicrle for now
  # 
  ## Inputs:
  # num_samples (integer)[optional]: number of smaples to generate
  # x_step_size (numeric)[optional]: x_new = x_old + x_step_size + f_x_old
  # y_step_size (numeric)[optional]: ibid.
  #
  ## Outputs:
  # sampled_data (data.frame): num_samples samples; columns in [x,y,f_x,f_y]
  
  # init output
  x <- c(initial_condition[1])
  y <- c(initial_condition[2])
  f_x <- c()
  f_y <- c()
  
  for (i in 1:num_samples){
    # update rule
    grad <- c(y[i] + (x[i]^2 + y[i]^2-1)*x[i],x[i] + (x[i]^2 + y[i]^2-1)*y[i])
    x <- c(x, x[i] + x_step_size * grad[1])
    y <- c(y, y[i] + y_step_size * grad[2])
    f_x <- c(f_x, grad[1])
    f_y <- c(f_y, grad[2])
  }
  
  sampled_data <- data.frame(x = head(x,num_samples),y = head(y,num_samples),f_x,f_y)
  return(sampled_data)
}

generate_morris_lecar <- function(){
  stop("Function not implemented")
}

generate_data <- function(model, params, num_samples = 1000, 
                          add_obs_noise = F, add_grad_nose = F,
                          use_seed = F, save_csv = F){
  # This is a wrapper function for the data generation under a number of models
  #
  ## Inputs:
  # model (string): name of the model to sample from the limit cycle of
  #   implemented: van der pol oscillator, Abhi
  #   to implement: morris-lecar, spiral
  # params (vector): vector of parameters for the specified model
  # add_obs_noise (logical)[optional]: whether to add random noise to each observation (x,y); TODO: implement
  # add_grad_nose (logical)[optional]: whether to add random noise to each gradient (f_x, f_y); TODO: implement
  # use_seed (logical)[optional]: whether to set a seed before drawing data
  # save_csv (logical)[optional]: whether to save the CSV of the data
  #
  ## Outputs:
  # sampled_data (data.frame): sampled data; columns in [x,y,f_x,f_y]
  
  if (use_seed){set.seed(2022)}
  
  if (model == "van_der_pol"){
    sampled_data <- generate_van_der_pol(params)
  }
  else if(model == "spiral"){
    sampled_data <- generate_spiral()
  }
  else if(model == "abhi"){
    sampled_data <- get_abhi_data()
  }
  else{
    stop('Error: ', model, ' model not implemented.')
  }
  
  # apply Abhi's spline method to estimate the gradient field
  spline_fit(sampled_data)
  
  if (save_csv){
    # TODO: Implement
  }
  return(sampled_data)
}

# estimate parameters

spline_fit <- function(data){
  x <- data$x
  y <- data$y
  f_x <- data$f_x
  f_y <- data$f_y

  xmin <- floor(min(x))
  xmax <- ceiling(max(x))
  ymin <- floor(min(y))
  ymax <- ceiling(max(y))
  
  # create basis functions with 10 knots
  # TODO: adapt number of basis functions?
  xbasis = create.bspline.basis(range=c(xmin,xmax),norder=4,nbasis=12)
  ybasis = create.bspline.basis(range=c(ymin,ymax),norder=4,nbasis=12)

  xbvals = eval.basis(x,xbasis)
  ybvals = eval.basis(y,ybasis)
  
  # What we need is the evaluation of phi_j(x_i)*psi_k(x_j) for
  # each x and y. This creates 144 rows. I'll produce this using
  # kronecker products
  
  Xmat = (xbvals%x%matrix(1,1,12)) * (matrix(1,1,12)%x%ybvals)
  print(dim(Xmat))
  # Here the columns of xbvals are repeated 12 times in sequence
  # while each column of ybvals is repeated 12 times togehter.
  
  # Now we need a penalty matrix. We can get the penalty for one
  # basis from
  
  xPen = eval.penalty(xbasis,2)
  yPen = eval.penalty(ybasis,2)
  
  # to create the combined penalty we take the same kronecker
  # product form
  
  allPen = xPen%x%diag(12) + diag(12)%x%yPen
  
  # (note that this penalizes the sum of squared second derivative,
  # without the cross term that would go into a thin plate spline
  # penalty)
  
  # And we can put it all together as
  
  lambda = 1e-8
  coefs_1 =  solve( t(Xmat)%*%Xmat  + lambda*allPen, t(Xmat)%*%f_x)  
  coefs_2 =  solve( t(Xmat)%*%Xmat  + lambda*allPen, t(Xmat)%*%f_y)
  
  # We'll reshape the coefficients into a matrix and put it in
  # a bivariate functional data object
  
  sfd_1 = bifd(t(matrix(coefs_1,12,12)), xbasis,ybasis)
  sfd_2 = bifd(t(matrix(coefs_2,12,12)), xbasis,ybasis)
  # and evaluate
  
  xpts = seq(xmin,xmax,len=41)
  ypts = seq(ymin,ymax,len=51)
  
  smat1 = eval.bifd(xpts,ypts,sfd_1)
  smat2 = eval.bifd(xpts,ypts,sfd_2)
  
  contour(xpts,ypts,smat1)
  contour(xpts,ypts,smat2)
  
  
  eta = 0.1
  plot(x,y,type='l',xlim=c(xmin,xmax),ylim=c(ymin,ymax))
  
  for(i in seq(1,41,by=5)){
    for(j in seq(1,51,by=5)){
      arrows(xpts[i],ypts[j],xpts[i]+eta*smat1[i,j],ypts[j]+ eta*smat2[i,j])
    }
  }
  
  traj = lsoda(c(x[1],y[1]),seq(0,100,by=0.1),sVdP,0,sfd_1=sfd_1,sfd_2=sfd_2)
  
  #  Check how well we interpolate
  
  try1 = eval.bifd(x,y,sfd_1)
  try2 = eval.bifd(x,y,sfd_2)
  
  plot(f_x,diag(try1))
  plot(f_y,diag(try2))
  abline(c(0,1))
}

sVdP = function(t,x,p,sfd_1,sfd_2){
  dx = eval.bifd(x[1],x[2],sfd_1)
  dy = eval.bifd(x[1],x[2],sfd_2)
  
  return(list(as.vector(c(dx,dy))))
}

local_linear_regression <- function(){
  stop("Function not implemented")
}

## analyze results

plot_observations <- function(data){
  # helper function which generates a scatter plot of observations 
  ggplot(data,aes(x=x,y=y)) +
    geom_point() + 
    labs(x="x",y="y",title="Observations from a Dynamical System")
}

