source(here::here("estimation_methods","bspline.R"))

library(MASS)

# new implementation. only supports stochastic batching, with and without a check of eigenvectors
calculate_gd_spline_gradient_field <- function(data, x_grid, y_grid, gd_params, side_info = list(), 
                                               norder = 6, nbasis = 24, penalty_order = 2, lambda = 1e-8){
  init_spline_fit <- calculate_spline_gradient_field(data,x_grid, y_grid,side_info=side_info,
                                                     norder=norder,nbasis=nbasis,penalty_order=penalty_order,lambda=lambda)
  # radius for delta-tube evaluation points. radius = 0 uses samples instead
  if (!(gd_params$dt_radius == 0)){           
    gd_points <- get_gd_delta_tube(data, radius = gd_params$dt_radius)
  } else{
    gd_points <- data
  }

  gd_update_results <- batch_jacobian_gd(gd_points, gd_params$eig_rank, init_spline_fit$fdx, init_spline_fit$fdy,
                                            gd_params$batching, gd_params$eta, gd_params$algorithm)

  gd_spline_grid_x <- eval.bifd(x_grid,y_grid,gd_update_results$fdx)
  gd_spline_grid_y <- eval.bifd(x_grid,y_grid,gd_update_results$fdy)
  
  # return as |Grid| x 2 matrix
  spline_field <- matrix(c(c(gd_spline_grid_x),c(gd_spline_grid_y)),ncol=2)
  
  gd_params$batching$algorithm <- gd_params$algorithm
  analyze_gd_convergence(gd_points, gd_params$batching, gd_update_results, 
                         list(x=x_grid, y=y_grid),
                         list(x=init_spline_fit$fdx,y=init_spline_fit$fdy),
                         gd_update_results$basis_info) # fda objects

  return(list(field = spline_field, fdx = gd_update_results$fdx, fdy = gd_update_results$fdy, 
              gd_log = gd_update_results$batch_results))
}

batch_jacobian_gd <- function(gd_data, eig_rank, init_fdx, init_fdy, batch_params, eta, algorithm){
  # evaluate basis function and derivatives at knots
  # IMPORTANT: assumes identical basis functions for both directions, hence repeated use of fdx
  xbasis.eval = eval.basis(gd_data[,1],init_fdx$sbasis)
  ybasis.eval = eval.basis(gd_data[,2],init_fdx$tbasis)
  xbasis.deriv = eval.basis(gd_data[,1],init_fdx$sbasis, Lfdobj = 1)
  ybasis.deriv = eval.basis(gd_data[,2],init_fdx$tbasis, Lfdobj = 1) 
  basis_list = list(x0 = xbasis.eval, y0 = ybasis.eval, x1 = xbasis.deriv, y1 = ybasis.deriv)
  
  # for analysis, each set of batches is aggregated into an iteration
  total_batches <- batch_params$num_iterations*batch_params$batches_per
  eigen_results <- list() # to be filled with eigenvalues
  coefficient_results_x <- matrix(-5, nrow = total_batches, ncol = nrow(init_fdx$coefs)*nrow(init_fdy$coefs))
  coefficient_results_y <- matrix(-5, nrow = total_batches, ncol = nrow(init_fdx$coefs)*nrow(init_fdy$coefs))
  gd_results <- list(fdx = init_fdx, fdy = init_fdy, results = list(eigs = eigen_results, coeff_x = coefficient_results_x, coeff_y = coefficient_results_y))
  
  for (batch_num in 1:total_batches){
    gd_results <- run_gd(gd_data, eig_rank, basis_list, gd_results, batch_params, eta, batch_num, algorithm)
  }
  gd_results$basis_info <- basis_list
  return(gd_results)
}

run_gd <- function(gd_data, eig_rank, basis_list, prior_results, batch_params, eta, batch_num, algorithm){ 
  init_fdx <- prior_results$fdx
  init_fdy <- prior_results$fdy
  
  # get estimated field at samples using initial spline fit
  init_fit_fx <- eval.bifd(gd_data[,1],gd_data[,2],init_fdx)
  init_fit_fy <- eval.bifd(gd_data[,1],gd_data[,2],init_fdy)
  
  # this algorithm calculates the eigenvalue penalty of the Jacobian of residuals from the largest eigenvector 
  fdx_updates <- matrix(0, nrow = nrow(init_fdx$coefs)*nrow(init_fdy$coefs), ncol = batch_params$batch_size)
  fdy_updates <- matrix(0, nrow = nrow(init_fdx$coefs)*nrow(init_fdy$coefs), ncol = batch_params$batch_size)

  sampled_indices <- c()
  eigen_results <- c() # keep track of eigenvalues at sample indices
  batch_index <- 1
  current_run <- 0
  while (length(sampled_indices) < batch_params$batch_size){
    
    sample_index <- sample(setdiff(1:nrow(gd_data),sampled_indices),1)
    
    # TODO: replace with Jacobian helper function
    x_y_matrix <- outer(basis_list$x0[sample_index,],basis_list$y0[sample_index,])
    dx_y_matrix <- outer(basis_list$x1[sample_index,],basis_list$y0[sample_index,])
    x_dy_matrix <- outer(basis_list$x0[sample_index,], basis_list$y1[sample_index,])
    x_y_vec <- c(x_y_matrix)
    dx_y_vec <- c(dx_y_matrix)
    x_dy_vec <- c(x_dy_matrix)
    x_coefs <- c(init_fdx$coefs)
    y_coefs <- c(init_fdy$coefs)
    
    # get eigen-stuff from the original Jacobian
    jacobian <- matrix(c(x_coefs%*%c(dx_y_matrix),y_coefs%*%c(dx_y_matrix),x_coefs%*%c(x_dy_matrix),y_coefs%*%c(x_dy_matrix)),nrow=2)
    jacobian_eigen <- eigen(jacobian)
    jacobian_eig_order <- order(Re(jacobian_eigen$values),decreasing=T)
    jacobian_eig_odered_values = jacobian_eigen$values[jacobian_eig_order]
    jacobian_eig_ordered_vectors = jacobian_eigen$vectors[,jacobian_eig_order]
    
    if(algorithm == "Vanilla"){
      # if desired, skip points with a currently negative max real eigenvalue
      if (batch_params$skip_negative){
        if (Re(jacobian_eig_odered_values[1]) < 0){
          current_run <- current_run + 1
          if (current_run > nrow(gd_data)){ break }
          next  
        }
      }
      
      eig_outer <- outer(jacobian_eig_ordered_vectors[,eig_rank],jacobian_eig_ordered_vectors[,eig_rank])
      lambda_step.x <- Re(eig_outer[1,1]) * dx_y_vec + Re(eig_outer[1,2]) * x_dy_vec
      lambda_step.y <- Re(eig_outer[2,1]) * dx_y_vec + Re(eig_outer[2,2]) * x_dy_vec
      fdx_updates[,batch_index] <- lambda_step.x
      fdy_updates[,batch_index] <- lambda_step.y
    } else if (algorithm == "Projection"){
      # IMPORTANT: for now we project off estimated gradient, not observed
      true_sample_gradient <- gd_data[sample_index, c(3,4)]
      estimated_sample_gradient <- c(init_fit_fx[sample_index,sample_index],init_fit_fy[sample_index,sample_index])
      # TODO: check discrepancy is tolerable at initialization
  
      gradient_norm <- norm(estimated_sample_gradient, type="2")
      projection_matrix <- outer(estimated_sample_gradient,estimated_sample_gradient) / gradient_norm
      jacobian_residuals <- (diag(2) -  projection_matrix) %*% jacobian
      jacobian_residuals_eigen <- eigen(jacobian_residuals)
      jacobian_residuals_eig_order <- order(Re(jacobian_residuals_eigen$values),decreasing=T)
      jacobian_residuals_eig_ordered_values = jacobian_residuals_eigen$values[jacobian_residuals_eig_order]
      jacobian_residuals_eig_ordered_vectors = jacobian_residuals_eigen$vectors[,jacobian_residuals_eig_order]
      
      # makes more sense to skip negative eigenvalues in the projeciton sense
      if (batch_params$skip_negative){
        if (Re(jacobian_residuals_eig_ordered_values[1]) < 0){
          current_run <- current_run + 1
          if (current_run > nrow(gd_data)){ break }
          next  
        }
      }
      
      # the hermitian transpose is required for Thm. 2 of Magnus (1985)
      jacobian_residuals_conjugate <- Conj(t(jacobian_residuals))
      jacobian_residuals_conjugate_eigen <- eigen(jacobian_residuals_conjugate)
      jacobian_residuals_conjugate_eig_order <- order(Re(jacobian_residuals_conjugate_eigen$values),decreasing=T)
      jacobian_residuals_conjugate_eig_ordered_values = jacobian_residuals_conjugate_eigen$values[jacobian_residuals_conjugate_eig_order]
      jacobian_residuals_conjugate_eig_ordered_vectors = jacobian_residuals_conjugate_eigen$vectors[,jacobian_residuals_conjugate_eig_order]
      
    
      # these lists will contain the matrix of partial derivatives wrt each spline coefficient per direction
      jacobian_partial_list_x <- list()
      projection_partial_list_x <- list()
      jacobian_partial_list_y <- list()
      projection_partial_list_y <- list()
      
      # get partials of eigen vectors wrt Jacobian
      for (i in 1:length(x_coefs)){
        dJ_dc_x = matrix(c(dx_y_vec[i],0,x_dy_vec[i],0),2,2)
        jacobian_partial_list_x[[i]] <- dJ_dc_x
        # again assuming same bivariate basis for fdx and fdy
        dJ_dc_y = matrix(c(0,dx_y_vec[i],0,x_dy_vec[i]),2,2)
        jacobian_partial_list_y[[i]] <- dJ_dc_y
      }
      
      for (i in 1:length(x_coefs)){
        estimated_sample_gradient_x_dc <- c(x_y_matrix[i], 0)
        dP_dc_x = (1/gradient_norm) * ( outer(estimated_sample_gradient_x_dc, estimated_sample_gradient) +
                                        outer(estimated_sample_gradient, estimated_sample_gradient_x_dc) -
                                        (as.vector(estimated_sample_gradient %*% estimated_sample_gradient_x_dc)/gradient_norm^2 ) * 
                                        outer(estimated_sample_gradient,estimated_sample_gradient))
        projection_partial_list_x[[i]] <- dP_dc_x
        
        # again assuming same bivariate basis for fdx and fdy
        estimated_sample_gradient_y_dc <- c(0, x_y_matrix[i])
        dP_dc_y = (1/gradient_norm) * ( outer(estimated_sample_gradient_y_dc, estimated_sample_gradient) +
                                          outer(estimated_sample_gradient, estimated_sample_gradient_y_dc) -
                                          (as.vector(estimated_sample_gradient %*% estimated_sample_gradient_y_dc)/gradient_norm^2) *
                                          outer(estimated_sample_gradient, estimated_sample_gradient))
        projection_partial_list_y[[i]] <- dP_dc_y
      }
      
      dPJ_dC_x <- Map(function(x,y) x %*% jacobian + projection_matrix %*% y, projection_partial_list_x, jacobian_partial_list_x)
      dPJ_dC_y <- Map(function(x,y) x %*% jacobian + projection_matrix %*% y, projection_partial_list_y, jacobian_partial_list_y)
      
      dMJ_dC_x <- Map(function(x,y) x - y, jacobian_partial_list_x, dPJ_dC_x)
      dMJ_dC_y <- Map(function(x,y) x - y, jacobian_partial_list_y, dPJ_dC_y)
      
      # apply formula of partial eigenvalue wrt coefficients
      w_0 <- jacobian_residuals_conjugate_eig_ordered_vectors[,eig_rank]
      v_0 <- jacobian_residuals_eig_ordered_vectors[,eig_rank]
      dL_dC_x <-unlist(Map(function(x) (1/sqrt(w_0%*%v_0) * (t(w_0) %*% x %*% v_0) ),dMJ_dC_x))
      dL_dC_y <- unlist(Map(function(x) (1/sqrt(w_0%*%v_0) * (t(w_0) %*% x %*% v_0) ),dMJ_dC_y))
  
      # coefficient update step may be complex
      fdx_updates[,batch_index] <- Re(dL_dC_x)
      fdy_updates[,batch_index] <- Re(dL_dC_y)
    } else {stop("Invalid SGD algorithm")}
    
    eigen_results <- c(eigen_results, jacobian_eig_odered_values)
    sampled_indices <- c(sampled_indices, sample_index)
    batch_index <- batch_index + 1
  }
  
  init_fdx$coefs <- matrix(init_fdx$coefs - eta*apply(fdx_updates,1,mean), nrow = nrow(init_fdx$coefs))
  init_fdy$coefs <- matrix(init_fdy$coefs - eta*apply(fdy_updates,1,mean), nrow = nrow(init_fdy$coefs))
  
  names(eigen_results) <- rep(sampled_indices, each=2)
  prior_results$results$eigs[[batch_num]] <- eigen_results
  prior_results$results$coeff_x[batch_num,] <- c(init_fdx$coefs)
  prior_results$results$coeff_y[batch_num,] <- c(init_fdy$coefs)

  return(list(fdx = init_fdx, fdy = init_fdy, results = prior_results$results))
}

#################
# Visualization #
#################
analyze_gd_convergence <- function(gd_points, batch_params, results, grid_list, fda_list, basis_list){
  #browser()
  results_list <- results$results
  # Code to get eigen values at points we optimized WRT to. Not used
  batch_indices <- unlist(lapply(results_list$eigs,function(x){as.numeric(unique(names(x)))}))
  sampled_points <- gd_points[batch_indices,]
  # num_iterations=10,batches_per=50,batch_size
  sampled_points$iteration <- rep(1:batch_params$num_iterations, each = batch_params$batches_per*batch_params$batch_size)
  sampled_points$batch <- rep(rep(1:batch_params$batches_per,each=batch_params$batch_size),batch_params$num_iterations)
  sampled_points$lambda1 <- unlist(lapply(results_list$eigs,function(x){as.numeric(x[c(TRUE,FALSE)])}))
  sampled_points$lambda2 <- unlist(lapply(results_list$eigs,function(x){as.numeric(x[c(FALSE,TRUE)])}))
  sampled_points <- mutate(sampled_points, neg_eig = (Re(lambda1) < 0))
  
  title_str <- paste0("Alg:",batch_params$algorithm,
                      " Iterations:",batch_params$num_iterations," Batches Per:",batch_params$batches_per,
                      " Batch Size:",batch_params$batch_size)
  
  # plot the first batch from each iteration and final iteration
  plot_indices <- c(batch_params$batches_per*c(0:(batch_params$num_iterations-1)) + 1,
                    batch_params$batches_per*batch_params$num_iterations)

  grid_matrix <- unname(as.matrix(expand.grid(grid_list$x,grid_list$y)))
  spline_result_list <- list()
  
  for (iteration_index in 1:length(plot_indices)){
    spline_coeff_x <- results_list$coeff_x[plot_indices[iteration_index],]
    fx_spline <- fda_list$x
    fx_spline$coefs <- matrix(spline_coeff_x,ncol=ncol(fx_spline$coefs))
    
    spline_coeff_y <- results_list$coeff_y[plot_indices[iteration_index],]
    fy_spline <- fda_list$y
    fy_spline$coefs <- matrix(spline_coeff_y,ncol=ncol(fy_spline$coefs))
    
    # get field
    spline_grid_x <- eval.bifd(grid_list$x,grid_list$y,fx_spline)
    spline_grid_y <- eval.bifd(grid_list$x,grid_list$y,fy_spline)
    spline_field <- matrix(c(c(spline_grid_x),c(spline_grid_y)),ncol=2)
    spline_field <- cbind(spline_field, plot_indices[iteration_index]) # add indicator of iteration for later ggplot
    #field_plot <- ggplot_field(gd_points, grid_matrix, spline_field)
    
    # get eigenvalues at LC samples

    lc_eval_x <- diag(eval.bifd(gd_points$x,gd_points$y,fx_spline))
    lc_eval_y <- diag(eval.bifd(gd_points$x,gd_points$y,fy_spline))
    eig_mat <- do.call(rbind,lapply(c(1:nrow(gd_points)),get_jacobian_eig, basis_list, fx_spline, fy_spline))
    projected_eig_mat <- do.call(rbind,lapply(c(1:nrow(gd_points)),get_projected_jacobian_eig, basis_list, 
                                              fx_spline, fy_spline, lc_eval_x, lc_eval_y))
    
    prediction_tibble <- tibble(x = gd_points$x, y = gd_points$y,
                                f_x = lc_eval_x, f_y = lc_eval_y,
                                lambda1 = eig_mat[,1], lambda2 = eig_mat[,2],
                                lambda1P = projected_eig_mat[,1], lambda2P = projected_eig_mat[,2],
                                iteration = as.factor(plot_indices[iteration_index]))
    
    spline_result_list[[iteration_index]] <- list(x_spline = fx_spline,
                                                                 y_spline = fy_spline,
                                                                 field = spline_field,
                                                                 predictions = prediction_tibble)
  }
  
  full_eig_tibble <- do.call(rbind, lapply(spline_result_list, function(x){x$predictions}))
  full_eig_tibble <- do.call(rbind, lapply(spline_result_list, function(x){x$predictions}))
  
  eig_scatter1 <- ggplot(full_eig_tibble, aes(x=x,y=y,color=Re(lambda1))) +
     geom_point() +
     scale_color_viridis()  +
     labs(title="1st Eigenvalue", x = "x", y = "y", subtitle = title_str) +
     facet_wrap(~iteration)
  print(eig_scatter1)
  
  eig_scatter1_bool <- ggplot(full_eig_tibble, aes(x=x,y=y,color=(Re(lambda1)<0))) +
    geom_point() +
    labs(title="1st Eigenvalue Sign", x = "x", y = "y", subtitle = title_str) +
    facet_wrap(~iteration)
  print(eig_scatter1_bool)
  
  eig_violin1 <- ggplot(full_eig_tibble, aes(x=as.factor(iteration),y=Re(lambda1))) +
    geom_violin(draw_quantiles = c(0.5)) +
    labs(title="Distribution of Re(Lambda1) per Iteration",
         x = "Iteration", y = "Re(Lambda1)", subtitle = title_str)
  print(eig_violin1)
  
  eig_scatter2 <- ggplot(full_eig_tibble, aes(x=x,y=y,color=Re(lambda2))) +
    geom_point() +
    scale_color_viridis()  +
    labs(title="2nd Eigenvalue", x = "x", y = "y", subtitle = title_str) +
    facet_wrap(~iteration)
  print(eig_scatter2)
  
  eig_scatter2_bool <- ggplot(full_eig_tibble, aes(x=x,y=y,color=(Re(lambda2)<0))) +
    geom_point() +
    labs(title="2nd Eigenvalue Sign", x = "x", y = "y", subtitle = title_str) +
    facet_wrap(~iteration)
  print(eig_scatter2_bool)
  
  eig_violin2 <- ggplot(full_eig_tibble, aes(x=as.factor(iteration),y=Re(lambda2))) +
    geom_violin(draw_quantiles = c(0.5)) +
    labs(title="Distribution of Re(Lambda2) per Iteration",
         x = "Iteration", y = "Re(Lambda2)", subtitle = title_str)
  print(eig_violin2)
  
  eig_scatter1_proj <- ggplot(full_eig_tibble, aes(x=x,y=y,color=Re(lambda1P))) +
    geom_point() +
    scale_color_viridis()  +
    labs(title="1st Projected Eigenvalue", x = "x", y = "y", subtitle = title_str) +
    facet_wrap(~iteration)
  print(eig_scatter1_proj)
  
  eig_scatter1_proj_bool <- ggplot(full_eig_tibble, aes(x=x,y=y,color=(Re(lambda1P)<0))) +
    geom_point() +
    labs(title="1st Projected Eigenvalue Sign", x = "x", y = "y", subtitle = title_str) +
    facet_wrap(~iteration)
  print(eig_scatter1_proj_bool)
  
  eig_violin1_proj <- ggplot(full_eig_tibble, aes(x=as.factor(iteration),y=Re(lambda1P))) +
    geom_violin(draw_quantiles = c(0.5)) +
    labs(title="Distribution of Re(Lambda1 Projected) per Iteration",
         x = "Iteration", y = "Re(Lambda1P)", subtitle = title_str)
  print(eig_violin1_proj)
  
  full_field_tibble <- do.call(rbind, lapply(spline_result_list, function(x){x$field}))
  
  rescale <- 0.005
  sample_tibble <- tibble(x = gd_points[,1], y = gd_points[,2], u = NA, v = NA)
  gradient_tibble <- tibble(x = rep(grid_matrix[,1],length(plot_indices)), y = rep(grid_matrix[,2],length(plot_indices)), 
                            u = rescale*full_field_tibble[,1], v =rescale*full_field_tibble[,2],
                            Batch = as.factor(full_field_tibble[,3]))
  
  field_plot <- ggplot(gradient_tibble, aes(x = x, y = y, u = u, v = v)) +
    geom_quiver(aes(color = Batch),alpha=.75) +
    geom_point(data = gd_points, aes(x = x, y = y, u = NA, v = NA), color = "#FDB515") +
    labs(title = "Gradient Field Evolution",subtitle = title_str)
  
  print(field_plot)

  return()
}

####################
# Helper Functions #
####################

get_gd_delta_tube <- function(data, radius = 0.2, sample_frac = 1){
  delta_ring_samples <- data[sort(sample(nrow(data),size=floor(nrow(data)*sample_frac),replace=FALSE)),]
  delta_ring_sideinfo <- t(apply(delta_ring_samples,1,get_normal_vectors, radius = radius))
  side_information_pos <- rbind(delta_ring_sideinfo[,c(1,2)],delta_ring_sideinfo[,c(3,4)])
  side_information_pos <- cbind(side_information_pos,rbind(data[,c(3,4)],data[,c(3,4)]))
  colnames(side_information_pos) <- c("x","y","f_x","f_y")
  return(side_information_pos)
}

get_jacobian_eig <- function(sample_index, basis_list, fdx, fdy){
  # TODO: replace with Jacobian helper function
  x_y_matrix <- outer(basis_list$x0[sample_index,],basis_list$y0[sample_index,])
  dx_y_matrix <- outer(basis_list$x1[sample_index,],basis_list$y0[sample_index,])
  x_dy_matrix <- outer(basis_list$x0[sample_index,], basis_list$y1[sample_index,])
  x_y_vec <- c(x_y_matrix)
  dx_y_vec <- c(dx_y_matrix)
  x_dy_vec <- c(x_dy_matrix)
  x_coefs <- c(fdx$coefs)
  y_coefs <- c(fdy$coefs)
  
  # get eigen-stuff from the original Jacobian
  jacobian <- matrix(c(x_coefs%*%dx_y_vec,y_coefs%*%dx_y_vec,x_coefs%*%x_dy_vec,y_coefs%*%x_dy_vec),nrow=2)
  jacobian_eigen <- eigen(jacobian)
  jacobian_eig_order <- order(Re(jacobian_eigen$values),decreasing=T)
  jacobian_eig_odered_values = jacobian_eigen$values[jacobian_eig_order]
  
  return(jacobian_eig_odered_values)
}

get_projected_jacobian_eig <- function(sample_index, basis_list, fdx, fdy, fx_eval, fy_eval){
  # TODO: replace with Jacobian helper function
  x_y_matrix <- outer(basis_list$x0[sample_index,],basis_list$y0[sample_index,])
  dx_y_matrix <- outer(basis_list$x1[sample_index,],basis_list$y0[sample_index,])
  x_dy_matrix <- outer(basis_list$x0[sample_index,], basis_list$y1[sample_index,])
  x_y_vec <- c(x_y_matrix)
  dx_y_vec <- c(dx_y_matrix)
  x_dy_vec <- c(x_dy_matrix)
  x_coefs <- c(fdx$coefs)
  y_coefs <- c(fdy$coefs)
  
  # get eigen-stuff from the original Jacobian
  jacobian <- matrix(c(x_coefs%*%dx_y_vec,y_coefs%*%dx_y_vec,x_coefs%*%x_dy_vec,y_coefs%*%x_dy_vec),nrow=2)
  jacobian_eigen <- eigen(jacobian)
  jacobian_eig_order <- order(Re(jacobian_eigen$values),decreasing=T)
  jacobian_eig_odered_values = jacobian_eigen$values[jacobian_eig_order]
  
  estimated_sample_gradient <- c(fx_eval[sample_index],fy_eval[sample_index])
  # TODO: check discrepancy is tolerable at initialization
  
  gradient_norm <- norm(estimated_sample_gradient, type="2")
  projection_matrix <- outer(estimated_sample_gradient,estimated_sample_gradient) / gradient_norm
  jacobian_residuals <- (diag(2) -  projection_matrix) %*% jacobian
  jacobian_residuals_eigen <- eigen(jacobian_residuals)
  jacobian_residuals_eig_order <- order(Re(jacobian_residuals_eigen$values),decreasing=T)
  jacobian_residuals_eig_ordered_values = jacobian_residuals_eigen$values[jacobian_residuals_eig_order]
  jacobian_residuals_eig_ordered_vectors = jacobian_residuals_eigen$vectors[,jacobian_residuals_eig_order]
  
  
  return(jacobian_residuals_eig_ordered_values)
}