NW_regression <- function(eval_grid, data){
  x_grid <- eval_grid[,1]
  y_grid <- eval_grid[,2]
  nw_results_array <- array(-5, c(nrow(eval_grid), 2))
  
  bw_matrix <- matrix(c(.1,0,0,0.1),2)
  for (i in 1:nrow(eval_grid)){
    nw_results_array[i,] <-  eval_gaussian_NW(eval_grid[i,1], eval_grid[i,2], data, bw_matrix)
  }
  # for (i in 1:length(x_grid)){
  #   for (j in 1:length(y_grid)){
  #     nw_results_array[i,j,] <- eval_gaussian_NW(x_grid[i], y_grid[j], data, bw_matrix)
  #   }
  # }
  return(nw_results_array)
}

eval_gaussian_NW <- function(x_val, y_val, data, bw_matrix){
  eval_point <- c(x_val, y_val)
  evaluation_results <- data
  evaluation_results$kernel_distance <- apply(evaluation_results[,1:2],1,  function(x) dmvnorm(x = x, mean = eval_point, sigma = bw_matrix))
  
  evaluation_results$numerator_x_terms <- evaluation_results$f_x * evaluation_results$kernel_distance
  evaluation_results$numerator_y_terms <- evaluation_results$f_y * evaluation_results$kernel_distance
  
  gradient_estimate_x <- sum(evaluation_results$numerator_x_terms)/sum(evaluation_results$kernel_distance)
  gradient_estimate_y <- sum(evaluation_results$numerator_y_terms)/sum(evaluation_results$kernel_distance)
  
  return(c(gradient_estimate_x, gradient_estimate_y))
  
}

plot_field <- function(ls_samples, eval_grid, gradient_data,title_str=""){
  # plot a gradient field
  
  eta = 0.1
  x_grid = eval_grid[,1]
  y_grid = eval_grid[,2]
  
  plot(ls_samples$x, ls_samples$y,
       xlim=c(min(x_grid)-1,max(x_grid)+1),
       ylim=c(min(y_grid)-1,max(y_grid)+1), 
       main = title_str, xlab = "x", ylab = "y")
  
  for (i in 1:nrow(eval_grid)){
    # prevent over plotting
    if (i %% 2 == 0){
      Arrows(x_grid[i], y_grid[i], x_grid[i] + eta*gradient_data[i,1], y_grid[i] + eta*gradient_data[i,2], lwd=0.75) 
    }
  }
}

plot_gradient_path_single <- function(ls_samples, eval_grid, gradient_data, bw, num_reps = 9, title =""){
  # plots NW estimated trajectories
  nw_trajectory <- generate_nw_path(ls_samples, bw, num_samples = 500)
  trajectory_tibble <- tibble(x = nw_trajectory[,1], y = nw_trajectory[,2], u = NA, v = NA, run = 1)
  
  if (num_reps > 1){
    for (i in 2:num_reps){
      nw_trajectory <- generate_nw_path(ls_samples, bw, num_samples = 500)
      new_trajectory_tibble <- tibble(x = nw_trajectory[,1], y = nw_trajectory[,2], u = NA, v = NA, run = i)
      trajectory_tibble <- bind_rows(trajectory_tibble, new_trajectory_tibble)
    }
  }
  
  sample_tibble <- tibble(x = ls_samples[,1], y = ls_samples[,2], u = NA, v = NA)
  gradient_tibble <- tibble(x = eval_grid[,1], y = eval_grid[,2], u = gradient_data[,1], v = gradient_data[,2])
  
  path_plot <- ggplot(gradient_tibble, aes(x = x, y = y, u = u, v = v)) +
    geom_quiver(color = "#003262") +
    geom_path(data = trajectory_tibble, aes(x = x, y = y)) +
    geom_point(data = sample_tibble, aes(x = x, y = y), color = "#FDB515") +
    geom_point(x = as.numeric(trajectory_tibble[1,1]), y =as.numeric(trajectory_tibble[1,2]), color = "red", size = 3) +
    labs(title = title) +
    facet_wrap(~run)
  
  return(path_plot)
  
}

####### Old Code

spline_gradient <- function(data, x_grid, y_grid){
  # Estimates the gradient field after fitting b-splines to the limit cycle data
  #
  ## Inputs:
  # data (data.frame): contains data with columns in [x, y, f_x, f_y]
  # x_grid (numeric; vector): x-axis grid points to evaluate extrapolation over
  # y_grid (numeric; vector): y-axis grid points to evaluate extrapolation over
  #
  ## Outputs:
  # output_list (named list): see construction for contents; TODO: update documentation
  
  # create basis functions with (deault = 10) knots
  # TODO: adapt number of basis functions to the size of the limit cycle?
  xbasis = create.bspline.basis(range=c(min(x_grid),max(x_grid)),norder=4,nbasis=12)
  ybasis = create.bspline.basis(range=c(min(y_grid),max(y_grid)),norder=4,nbasis=12)
  xbasis.vals = eval.basis(data$x,xbasis)
  ybasis.vals = eval.basis(data$y,ybasis)
  
  # What we need is the evaluation of phi_j(x_i)*psi_k(x_j) for
  # each x and y. This creates 144 rows. I'll produce this using
  # Kronecker products
  Xmat = (xbasis.vals%x%matrix(1,1,12)) * (matrix(1,1,12)%x%ybasis.vals)
  
  # Here the columns of xbasis.vals are repeated 12 times in sequence
  # while each column of ybasis.vals is repeated 12 times together.
  
  # Now we need a penalty matrix. We can get the penalty for one
  # basis from
  xPen = eval.penalty(xbasis, 2)
  yPen = eval.penalty(ybasis, 2)
  
  # to create the combined penalty we take the same Kronecker
  # product form
  
  allPen = xPen%x%diag(12) + diag(12)%x%yPen
  #View(allPen)
  
  # (note that this penalizes the sum of squared second derivative,
  # without the cross term that would go into a thin plate spline
  # penalty)
  
  # And we can put it all together as
  
  lambda = 1e-8
  coefs_x =  solve( t(Xmat)%*%Xmat  + lambda*allPen, t(Xmat)%*%data$f_x)  
  coefs_y =  solve( t(Xmat)%*%Xmat  + lambda*allPen, t(Xmat)%*%data$f_y)
  
  
  # We'll reshape the coefficients into a matrix and put it in
  # a bivariate functional data object
  sfd_x = bifd(t(matrix(coefs_x,12,12)), xbasis, ybasis)
  sfd_y = bifd(t(matrix(coefs_y,12,12)), xbasis, ybasis)
  
  smat_x = eval.bifd(x_grid,y_grid,sfd_x)
  smat_y = eval.bifd(x_grid,y_grid,sfd_y)
  
  # build a penalty matrix to plot...
  #sfd_x.penalty = bifd(t(matrix(coefs_x,12,12)), xbasis, ybasis)
  #sfd_y.penalty = bifd(t(matrix(coefs_y,12,12)), xbasis, ybasis)
  #smat_x.penalty = eval.bifd(x_grid,y_grid,sfd_x)
  #smat_y.penalty = eval.bifd(x_grid,y_grid,sfd_y)
  #View(allPen)
  #print(dim(xPen))
  
  
  traj = lsoda(c(data$x[1],data$y[1]),seq(0,100,by=0.1),bifd_spline_gradient,0,sfd_x=sfd_x,sfd_y=sfd_y)
  
  output_list <- list(x_grad_bifd = sfd_x, x_grad_eval = smat_x, y_grad_bifd = sfd_y, y_grad_eval = smat_y, 
                      trajectory = traj, second.deriv_penalty = allPen)
  return(coefs_y)
}

generate_spline_plots <- function(data, spline_output, x_grid, y_grid){
  # Generate the plots from Abhi's original code
  
  contour(x_grid,y_grid,spline_output$x_grad_eval)
  contour(x_grid,y_grid,spline_output$y_grad_eval)
  eta = 0.1
  plot(data$x,data$y,type='l',xlim=c(min(x_grid),max(x_grid)),ylim=c(min(y_grid),max(y_grid)))
  for(i in seq(1,length(x_grid),by=5)){
    for(j in seq(1,length(y_grid),by=5)){
      arrows(x_grid[i],y_grid[j],x_grid[i]+eta*spline_output$x_grad_eval[i,j],y_grid[j]+ eta*spline_output$y_grad_eval[i,j])
    }
  }
  
  #  Check how well we interpolate
  try1 = eval.bifd(data$x,data$y,spline_output$x_grad_bifd)
  try2 = eval.bifd(data$x,data$y,spline_output$y_grad_bifd)
  
  plot(data$f_x,diag(try1))
  plot(data$f_y,diag(try2))
  abline(c(0,1))
}
