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