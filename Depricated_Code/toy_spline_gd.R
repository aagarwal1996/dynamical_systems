library(fda)
library(tidyverse)
library(here)

here::here()
source(here::here("Data_Generation","data_generation.R"))
source(here::here("run_simulation_redux.R"))

## Generate data
mu_param = 5
vdp_data <- generate_limit_cycle_data("van_der_pol",list(mu=mu_param),var_x=0.25,var_y=0.25)

plot_samples <- function(sampled_gradients, title_str = ""){
  colnames(sampled_gradients) <- c("x","y","u","v")
  gradient_tibble <- as_tibble(sampled_gradients)
  path_plot <- ggplot(gradient_tibble, aes(x = x, y = y, u = u, v = v)) + 
    geom_quiver(color = "#003262") +
    geom_point(aes(x = x, y = y, u = NA, v = NA), color = "#FDB515") +
    labs(title=title_str) 
  return(path_plot)
}
plot_samples(vdp_data,title_str=paste0("Sampled Van der Pol; mu = ",mu_param))

## Fit splines
lambda = 1
spline_order = 6
penalty_order = 2
nbasis = 20

x_grid <- seq(min(vdp_data$x) - 0.5, max(vdp_data$x) + 0.5, length.out = 24)
y_grid <- seq(min(vdp_data$y) - 0.5, max(vdp_data$y) + 0.5, length.out = 24)

generate_bspline_basis <- function(data, x_grid, y_grid, norder = 4, nbasis = 12,
                                   penalty_order = 2){
  # number of breaks equals `nbasis` - `norder` + 2; 10 for default parameters (8 interior knots)
  xbasis = create.bspline.basis(rangeval=c(min(x_grid),max(x_grid)),norder=norder,nbasis=nbasis)
  ybasis = create.bspline.basis(rangeval=c(min(y_grid),max(y_grid)),norder=norder,nbasis=nbasis)
  xbasis.eval = eval.basis(data$x,xbasis)
  ybasis.eval = eval.basis(data$y,ybasis)
  xbasis.deriv = eval.basis(data$x,xbasis, Lfdobj = 1) # get derivative of bassis functions
  ybasis.deriv = eval.basis(data$y,ybasis, Lfdobj = 1) 
  
  #View(solve(t(xbasis.vals) %*% xbasis.vals) %*% t(xbasis.vals) %*% data$x)
  # get roughness penalty for each axis
  xPen = eval.penalty(xbasis, penalty_order)
  yPen = eval.penalty(ybasis, penalty_order)
  
  spline_fit_list <- list(xbasis = xbasis, ybasis = ybasis, 
                          xbasis.eval = xbasis.eval, ybasis.eval = ybasis.eval,
                          xbasis.deriv = xbasis.deriv, ybasis.deriv = ybasis.deriv,
                          xpenalty = xPen, ypenalty = yPen)
  return(spline_fit_list)
}

bspline_basis_fns <- generate_bspline_basis(vdp_data, x_grid, y_grid, norder = spline_order, 
                                            nbasis = nbasis, penalty_order = penalty_order)
# lambda != 0
init_bspline_fit_coeffs <- fit_bsplines_cpp(bspline_basis_fns$xbasis.eval, bspline_basis_fns$ybasis.eval,
                                         bspline_basis_fns$xpenalty, bspline_basis_fns$ypenalty,
                                       vdp_data$f_x, vdp_data$f_y, lambda)

# create bivariate functional data objects for our fit splines 
init_spline.fd_x <- bifd(t(matrix(init_bspline_fit_coeffs[,1],nbasis,nbasis)), 
                    bspline_basis_fns$xbasis,  bspline_basis_fns$ybasis)
init_spline.fd_y <- bifd(t(matrix(init_bspline_fit_coeffs[,2],nbasis,nbasis)),
                    bspline_basis_fns$xbasis,  bspline_basis_fns$ybasis)

# evaluate over user-specified grid
init_spline_grid_x <- eval.bifd(x_grid,y_grid,init_spline.fd_x)
init_spline_grid_y <- eval.bifd(x_grid,y_grid,init_spline.fd_y)

init_spline_field <- matrix(c(c(init_spline_grid_x),c(init_spline_grid_y)),ncol=2)

## Visualize

# Plot raw field

ggplot_field <- function(ls_samples, eval_grid, gradient_data, title="", rescale = 0.005){
  # plot a gradient field with limit cycle samples overlayed
  sample_tibble <- tibble(x = ls_samples[,1], y = ls_samples[,2], u = NA, v = NA)
  gradient_tibble <- tibble(x = eval_grid[,1], y = eval_grid[,2], 
                            u = rescale*gradient_data[,1], v =rescale*gradient_data[,2])
field_plot <- ggplot(gradient_tibble, aes(x = x, y = y, u = u, v = v)) +
    geom_quiver(color = "#003262") +
    geom_point(data = sample_tibble, aes(x = x, y = y), color = "#FDB515") +
    labs(title = title)
  return(field_plot)
}

joint_grid = expand.grid(x = x_grid, y = y_grid)
ggplot_field(vdp_data,joint_grid,init_spline_field)

# Plot some solution paths

sol_paths_ic <- get_shared_ic(vdp_data,delta_ring=F)
init_sol_paths <- list()
init_spline_experiment = list(method="spline",sfd_list=list(fdx = init_spline.fd_x, fdy = init_spline.fd_y))
for (j in 1:nrow(sol_paths_ic)){
  init_sol_paths[[j]] <- generate_solution_path(vdp_data,init_spline_experiment,sol_paths_ic[j,])
}

plot_solution_paths <- function(solution_path_list, eval_grid, estimated_field, lc_data, title_str = ""){
  
  within_tibble <- tibble(x = solution_path_list[[1]]$x, y = solution_path_list[[1]]$y)
  on_tibble <- tibble(x = solution_path_list[[2]]$x, y = solution_path_list[[2]]$y)
  outside_close_tibble <- tibble(x = solution_path_list[[3]]$x, y = solution_path_list[[3]]$y)
  outside_far_tibble <- tibble(x = solution_path_list[[4]]$x, y = solution_path_list[[4]]$y)
  
  gradient_matrix <- matrix(cbind(eval_grid, estimated_field), ncol=4)
  colnames(gradient_matrix) <- c("x","y","u","v")
  gradient_tibble <- as_tibble(gradient_matrix)
  
  colnames(lc_data) <- c("x","y","u","v")
  lc_data <- as_tibble(lc_data)

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

plot_solution_paths(init_sol_paths, as.matrix(joint_grid), init_spline_field, vdp_data, title_str = "Initialization")

## Run Gradient Descent

get_normal_vectors <- function(sample, radius){
  rotation_matrix <- matrix(c(0,-1,1,0),ncol=2)
  magnitude = sqrt(sample["f_x"]^2 + sample["f_y"]^2)
  point <- c(sample["x"], sample["y"])
  normal_vec_out <- (rotation_matrix %*% c(sample["f_x"], sample["f_y"]))/magnitude
  
  ring_sample <- c(point + radius*normal_vec_out, point - radius*normal_vec_out) 
  # normal out x, y; normal in x, y
  return(ring_sample)
}

get_gd_delta_tube <- function(data, radius = 0.2, sample_frac = 1){
  delta_ring_samples <- data[sort(sample(nrow(data),size=floor(nrow(data)*sample_frac),replace=FALSE)),]
  delta_ring_sideinfo <- t(apply(delta_ring_samples,1,get_normal_vectors, radius = radius))
  side_information_pos <- rbind(delta_ring_sideinfo[,c(1,2)],delta_ring_sideinfo[,c(3,4)])
  colnames(side_information_pos) <- c("x","y")
  return(side_information_pos)
}

gd_points <- get_gd_delta_tube(vdp_data)

random_point_jacobian_gd <- function(init_spline.fd_x,init_spline.fd_y,gd_points,max_steps=5000,eta=0.002){
  # assumes same basis for f_x and f_y at the moment
  xbasis.eval = eval.basis(gd_points[,1],init_spline.fd_x$sbasis)
  ybasis.eval = eval.basis(gd_points[,2],init_spline.fd_x$tbasis)
  xbasis.deriv = eval.basis(gd_points[,1],init_spline.fd_x$sbasis, Lfdobj = 1) # get derivative of basis functions
  ybasis.deriv = eval.basis(gd_points[,2],init_spline.fd_x$tbasis, Lfdobj = 1) 
  
  spline.fd_x = init_spline.fd_x
  spline.fd_y = init_spline.fd_y
  is_neg <- c()
  for(i in 1:max_steps){
    sample_index <- sample(nrow(gd_points),1)
    linearization_point <- gd_points[sample_index,]
    dx_y = c(outer(xbasis.deriv[sample_index,],ybasis.eval[sample_index,]))
    x_dy = c(outer(xbasis.eval[sample_index,], ybasis.deriv[sample_index,]))
    x_coefs <- c(spline.fd_x$coefs)
    y_coefs <- c(spline.fd_y$coefs)
    jacobian <- matrix(c(x_coefs%*%dx_y,y_coefs%*%dx_y,x_coefs%*%x_dy,y_coefs%*%x_dy),nrow=2)
    jacobian_eigen <- eigen(jacobian)
    eig_order <- order(Re(jacobian_eigen$values),decreasing=T)
    ordered_values = jacobian_eigen$values[eig_order]
    ifelse(Re(ordered_values[1])<0,is_neg <- append(is_neg,1),is_neg <- append(is_neg,0))
    ordered_vectors = jacobian_eigen$vectors[,eig_order]
    max_eig_outer <- outer(ordered_vectors[,2],ordered_vectors[,2]) # Edited for min!
    lambda_step.x <- Re(max_eig_outer[1,1]) * dx_y + Re(max_eig_outer[1,2]) * x_dy
    lambda_step.y <- Re(max_eig_outer[2,1]) * dx_y + Re(max_eig_outer[2,2]) * x_dy
    spline.fd_x$coefs <- matrix(spline.fd_x$coefs - eta*lambda_step.x, nrow = nrow(spline.fd_x$coefs))
    spline.fd_y$coefs <- matrix(spline.fd_y$coefs - eta*lambda_step.y, nrow = nrow(spline.fd_y$coefs))
  }
  
  updated_fd <- list(fdx = spline.fd_x, fdy = spline.fd_y, is_neg = is_neg)
  return(updated_fd)
}

num_steps = 20000
eta=0.002
gd_result <- random_point_jacobian_gd(init_spline.fd_x,init_spline.fd_y,gd_points,max_steps=num_steps,eta=eta)

rolling_frac <- zoo::rollsum(gd_result$is_neg,1000)
ggplot(data=tibble(x=1000:(999+length(rolling_frac)),y=rolling_frac),aes(x=x,y=rolling_frac/1000)) +
  geom_line() + 
  labs(title="Precision of Gradient Descent",y="% of Last 1000 Gradient Descent Samples with Negative Maximum Real Eigenvalue",x="Iteration")

## Explore Results
gd_sol_paths <- list()
gd_spline_experiment = list(method="spline",sfd_list=list(fdx = gd_result$fdx, fdy = gd_result$fdy))
for (j in 1:nrow(sol_paths_ic)){
  gd_sol_paths[[j]] <- generate_solution_path(vdp_data,gd_spline_experiment,sol_paths_ic[j,])
}
gd_spline_grid_x <- eval.bifd(x_grid,y_grid,gd_result$fdx)
gd_spline_grid_y <- eval.bifd(x_grid,y_grid,gd_result$fdy)
gd_spline_field <- matrix(c(c(gd_spline_grid_x),c(gd_spline_grid_y)),ncol=2)

plot_solution_paths(gd_sol_paths, as.matrix(joint_grid), gd_spline_field, vdp_data,title_str = paste0("Gradient Descent, Iterations:",num_steps," Eta:",eta))
