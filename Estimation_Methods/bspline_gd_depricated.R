source(here::here("estimation_methods","bspline.R"))
library(MASS)

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

  batch_params <- gd_params$batching # extract batching parameters and build batching schedule
  if (batch_params$method == "random"){
    batch_program <- t(apply(matrix(1:nrow(gd_points), nrow= batch_params$num_updates*batch_params$cycle_size, ncol=nrow(gd_points), byrow=TRUE), 
                           1, sample, batch_params$batch_size))
    if (batch_params$batch_size == 1){ batch_program <- t(batch_program) } # transpose in this special case
    batch_results <- matrix(-5, nrow = nrow(batch_program), ncol = 2*batch_params$batch_size)
    #as.list(numeric(batch_params$num_updates*batch_params$cycle_size*batch_params$batch_size))
    #dim(batch_results) <- c(batch_params$num_updates*batch_params$cycle_size,batch_params$batch_size)

  } else if (batch_params$method == "cyclic"){
    batch_program <- matrix(1:nrow(gd_points), ncol=batch_params$batch_size) # Some points get run twice!
    batch_program <- rep(1,batch_params$num_updates) %x% batch_program
    
    batch_results <- matrix(-5, nrow = nrow(batch_program), ncol = 2*batch_params$batch_size)

  } else {stop("Batching Method Not Implemented")}
  
  gd_updated_splines <- batch_jacobian_gd(gd_params$eig_rank, init_spline_fit$fdx, init_spline_fit$fdy,
                                            gd_points, batch_params, batch_program, batch_results, eta=gd_params$eta)
  #gd_updated_splines <- random_point_jacobian_gd(init_spline_fit$fdx, init_spline_fit$fdy, gd_points,
  #                                               max_steps=gd_params$num_updates,eta=gd_params$eta)
  
  gd_spline_grid_x <- eval.bifd(x_grid,y_grid,gd_updated_splines$fdx)
  gd_spline_grid_y <- eval.bifd(x_grid,y_grid,gd_updated_splines$fdy)
  
  # return as |Grid| x 2 matrix
  spline_field <- matrix(c(c(gd_spline_grid_x),c(gd_spline_grid_y)),ncol=2)
  
  analyze_gd_convergence(gd_points, batch_params, batch_program, gd_updated_splines$batch_results, batch_params$num_updates)

  return(list(field = spline_field, fdx = gd_updated_splines$fdx, fdy = gd_updated_splines$fdy, 
              gd_log = gd_updated_splines$batch_results))
}

batch_jacobian_gd <- function(eig_rank, init_fdx, init_fdy, gd_points, batch_params, batch_program, batch_results, eta){
  # evalutate basis function and derivatives at knots
  # assumes identical basis functions for both directions, hence repeated use of fdx
  xbasis.eval = eval.basis(gd_points[,1],init_fdx$sbasis)
  ybasis.eval = eval.basis(gd_points[,2],init_fdx$tbasis)
  xbasis.deriv = eval.basis(gd_points[,1],init_fdx$sbasis, Lfdobj = 1)
  ybasis.deriv = eval.basis(gd_points[,2],init_fdx$tbasis, Lfdobj = 1) 
  basis_list = list(x0 = xbasis.eval, y0 = ybasis.eval, x1 = xbasis.deriv, y1 = ybasis.deriv)
  gd_results <- list(fdx = init_fdx, fdy = init_fdy, batch_results = batch_results)
  
  if (batch_params$algorithm == "Vanilla"){
    gd_algorithm <- run_gd
  } else if (batch_params$algorithm == "Projection"){
    gd_algorithm <- run_projection_gd
  } else{stop("GD algorithm not specified")}
  
  for (i in 1:nrow(batch_program)){
    batch_indices <- batch_program[i,] # only need the indices
    gd_results <- gd_algorithm(eig_rank, gd_results$fdx, gd_results$fdy, basis_list, batch_indices, i, gd_results$batch_results, eta)
  }
  #gd_results$batch_results <- cbind(gd_results$batch_results, rep(1:batch_params$num_updates, each = nrow(gd_results$batch_results)/batch_params$num_updates))
  
  return(gd_results)
}

run_gd <- function(eig_rank, init_fdx, init_fdy, basis_list, batch_indices, result_index, results_array, eta){
  fdx_updates <- matrix(0, nrow = nrow(init_fdx$coefs)*nrow(init_fdy$coefs), ncol = length(batch_indices))
  fdy_updates <- matrix(0, nrow = nrow(init_fdx$coefs)*nrow(init_fdy$coefs), ncol = length(batch_indices))
  
  for(i in 1:length(batch_indices)){
    sample_index <- batch_indices[i]
    dx_y = c(outer(basis_list$x1[sample_index,],basis_list$y0[sample_index,]))
    x_dy = c(outer(basis_list$x0[sample_index,], basis_list$y1[sample_index,]))
    x_coefs <- c(init_fdx$coefs)
    y_coefs <- c(init_fdy$coefs)
    
    jacobian <- matrix(c(x_coefs%*%dx_y,y_coefs%*%dx_y,x_coefs%*%x_dy,y_coefs%*%x_dy),nrow=2)
    jacobian_eigen <- eigen(jacobian)
    eig_order <- order(Re(jacobian_eigen$values),decreasing=T)
    ordered_values = jacobian_eigen$values[eig_order]
    ordered_vectors = jacobian_eigen$vectors[,eig_order]
    
    max_eig_outer <- outer(ordered_vectors[,eig_rank],ordered_vectors[,eig_rank])
    lambda_step.x <- Re(max_eig_outer[1,1]) * dx_y + Re(max_eig_outer[1,2]) * x_dy
    lambda_step.y <- Re(max_eig_outer[2,1]) * dx_y + Re(max_eig_outer[2,2]) * x_dy
    fdx_updates[,i] <- lambda_step.x
    fdy_updates[,i] <- lambda_step.y
    
    results_array[result_index,(2*i-1):(2*i)] = ordered_values # could be expanded to include eg eig vectors
  }

  init_fdx$coefs <- matrix(init_fdx$coefs - eta*apply(fdx_updates,1,mean), nrow = nrow(init_fdx$coefs))
  init_fdy$coefs <- matrix(init_fdy$coefs - eta*apply(fdy_updates,1,mean), nrow = nrow(init_fdy$coefs))
  
  return(list(fdx = init_fdx, fdy = init_fdy, batch_results = results_array))
}

run_projection_gd <- function(eig_rank, init_fdx, init_fdy, basis_list, batch_indices, result_index, results_array, eta){
  # this algorithm calculates the eigenvalue penalty of the Jacobian of residuals from the largest eigenvector 
  fdx_updates <- matrix(0, nrow = nrow(init_fdx$coefs)*nrow(init_fdy$coefs), ncol = length(batch_indices))
  fdy_updates <- matrix(0, nrow = nrow(init_fdx$coefs)*nrow(init_fdy$coefs), ncol = length(batch_indices))
  
  for(b in 1:length(batch_indices)){
    sample_index <- batch_indices[b]
    # these are now matrices!
    dx_y_matrix = outer(basis_list$x1[sample_index,],basis_list$y0[sample_index,])
    x_dy_matrix = outer(basis_list$x0[sample_index,], basis_list$y1[sample_index,])
    dx_y_vec = c(dx_y_matrix)
    x_dy_vec = c(x_dy_matrix)
    x_coefs <- c(init_fdx$coefs)
    y_coefs <- c(init_fdy$coefs)
    
    # get eigen-stuff from the original Jacobian
    jacobian <- matrix(c(x_coefs%*%c(dx_y_matrix),y_coefs%*%c(dx_y_matrix),x_coefs%*%c(x_dy_matrix),y_coefs%*%c(x_dy_matrix)),nrow=2)
    jacobian_eigen <- eigen(jacobian)
    jacobian_eig_order <- order(Re(jacobian_eigen$values),decreasing=T)
    jacobian_eig_odered_values = jacobian_eigen$values[jacobian_eig_order]
    jacobian_eig_ordered_vectors = jacobian_eigen$vectors[,jacobian_eig_order]
    
    # then for the residualized Jacobian
    projection_matrix <- outer(jacobian_eig_ordered_vectors[,projection_rank],jacobian_eig_ordered_vectors[,projection_rank])
    jacobian_residuals <- (diag(2) -  projection_matrix) %*% jacobian
    jacobian_residuals_eigen <- eigen(jacobian_residuals)
    jacobian_residuals_eig_order <- order(Re(jacobian_residuals_eigen$values),decreasing=T)
    jacobian_residuals_eig_ordered_values = jacobian_residuals_eigen$values[jacobian_residuals_eig_order]
    jacobian_residuals_eig_ordered_vectors = jacobian_residuals_eigen$vectors[,jacobian_residuals_eig_order]
    
    # the hermitian transpose is required for Thm. 2 of Magnus (1985)
    jacobian_residuals_conjugate <- Conj(t(jacobian_residuals))
    jacobian_residuals_conjugate_eigen <- eigen(jacobian_residuals_conjugate)
    jacobian_residuals_conjugate_eig_order <- order(Re(jacobian_residuals_conjugate_eigen$values),decreasing=T)
    jacobian_residuals_conjugate_eig_ordered_values = jacobian_residuals_conjugate_eigen$values[jacobian_residuals_conjugate_eig_order]
    jacobian_residuals_conjugate_eig_ordered_vectors = jacobian_residuals_conjugate_eigen$vectors[,jacobian_residuals_conjugate_eig_order]
    
    # these lists will contain the matrix of partial derivatives wrt each spline coefficient per direction
    jacobian_partial_list_x <- list()
    jacobian_residuals_partial_list_x <- list()
    jacobian_partial_list_y <- list()
    jacobian_residuals_partial_list_y <- list()
    
    # get partials of eigen vectors wrt Jacobian
    for (i in 1:length(x_coefs)){
      dJ_dc_x = matrix(c(dx_y_vec[i],0,x_dy_vec[i],0),2,2)
      jacobian_partial_list_x[[i]] <- dJ_dc_x
      dJ_dc_y = matrix(c(0,dx_y_vec[i],0,x_dy_vec[i]),2,2)
      jacobian_partial_list_y[[i]] <- dJ_dc_y
    }
    
    dP_dJ_tensor <- list() # partial of projection matrix wrt Jacobian
    for (i in 1:2){
      for (j in 1:2){
        ij_matrix <- get_dP_fixed_dJ(i, j , projection_matrix, jacobian, 
                                     jacobian_eig_ordered_vectors, jacobian_eig_odered_values[projection_rank], projection_rank)
        dP_dJ_tensor[[(2*(i-1)+j)]] <- ij_matrix
      }
    }
    
    dP_dC_x_11 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[1]]), jacobian_partial_list_x)
    dP_dC_y_11 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[1]]), jacobian_partial_list_y)
    dP_dC_x_21 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[2]]), jacobian_partial_list_x)
    dP_dC_y_21 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[2]]), jacobian_partial_list_y)
    dP_dC_x_12 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[3]]), jacobian_partial_list_x)
    dP_dC_y_12 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[3]]), jacobian_partial_list_y)
    dP_dC_x_22 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[4]]), jacobian_partial_list_x)
    dP_dC_y_22 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[4]]), jacobian_partial_list_y)
    
    dP_dC_x <- Map(function(w,x,y,z) matrix(c(w,x,y,z),2,2),dP_dC_x_11,dP_dC_x_12,dP_dC_x_21,dP_dC_x_22)
    dP_dC_y <- Map(function(w,x,y,z) matrix(c(w,x,y,z),2,2),dP_dC_y_11,dP_dC_y_12,dP_dC_y_21,dP_dC_y_22)
    
    dPJ_dC_x <- Map(function(x,y) x %*% jacobian + projection_matrix %*% y, dP_dC_x, jacobian_partial_list_x)
    dPJ_dC_y <- Map(function(x,y) x %*% jacobian + projection_matrix %*% y, dP_dC_y, jacobian_partial_list_y)
    
    dMJ_dC_x <- Map(function(x,y) x - y, jacobian_partial_list_x, dPJ_dC_x)
    dMJ_dC_y <- Map(function(x,y) x - y, jacobian_partial_list_y, dPJ_dC_y)
    
    w_0 <- jacobian_residuals_conjugate_eig_ordered_vectors[,penalty_ranks]
    v_0 <- jacobian_residuals_eig_ordered_vectors[,penalty_ranks]
    dL_dC_x <-unlist(Map(function(x) (1/sqrt(w_0%*%v_0) * (t(w_0) %*% x %*% v_0) ),dMJ_dC_x))
    dL_dC_y <- unlist(Map(function(x) (1/sqrt(w_0%*%v_0) * (t(w_0) %*% x %*% v_0) ),dMJ_dC_y))
    
    fdx_updates[,b] <- Re(dL_dC_x)
    fdy_updates[,b] <- Re(dL_dC_y)
    
    results_array[result_index,(2*b-1):(2*b)] = jacobian_eig_odered_values # could be expanded to include eg eig vectors
  }
  
  init_fdx$coefs <- matrix(init_fdx$coefs - eta*apply(fdx_updates,1,mean), nrow = nrow(init_fdx$coefs))
  init_fdy$coefs <- matrix(init_fdy$coefs - eta*apply(fdy_updates,1,mean), nrow = nrow(init_fdy$coefs))
  
  return(list(fdx = init_fdx, fdy = init_fdy, batch_results = results_array))
}

# projects off eigenvectors of Jacobian, which doesn't make much sense...
run_eig_projection_gd <- function(eig_rank, init_fdx, init_fdy, basis_list, batch_indices, result_index, results_array, eta){
  # this algorithm calculates the eigenvalue penalty of the Jacobian of residuals from the largest eigenvector 
  projection_rank <- eig_rank[1] # direction to project off
  penalty_ranks <- eig_rank[2] # direction to penalize in residuals
  fdx_updates <- matrix(0, nrow = nrow(init_fdx$coefs)*nrow(init_fdy$coefs), ncol = length(batch_indices))
  fdy_updates <- matrix(0, nrow = nrow(init_fdx$coefs)*nrow(init_fdy$coefs), ncol = length(batch_indices))
  
  for(b in 1:length(batch_indices)){
    sample_index <- batch_indices[b]
    # these are now matrices!
    dx_y_matrix = outer(basis_list$x1[sample_index,],basis_list$y0[sample_index,])
    x_dy_matrix = outer(basis_list$x0[sample_index,], basis_list$y1[sample_index,])
    dx_y_vec = c(dx_y_matrix)
    x_dy_vec = c(x_dy_matrix)
    x_coefs <- c(init_fdx$coefs)
    y_coefs <- c(init_fdy$coefs)
    
    # get eigen-stuff from the original Jacobian
    jacobian <- matrix(c(x_coefs%*%c(dx_y_matrix),y_coefs%*%c(dx_y_matrix),x_coefs%*%c(x_dy_matrix),y_coefs%*%c(x_dy_matrix)),nrow=2)
    jacobian_eigen <- eigen(jacobian)
    jacobian_eig_order <- order(Re(jacobian_eigen$values),decreasing=T)
    jacobian_eig_odered_values = jacobian_eigen$values[jacobian_eig_order]
    jacobian_eig_ordered_vectors = jacobian_eigen$vectors[,jacobian_eig_order]
    
    # then for the residualized Jacobian
    projection_matrix <- outer(jacobian_eig_ordered_vectors[,projection_rank],jacobian_eig_ordered_vectors[,projection_rank])
    jacobian_residuals <- (diag(2) -  projection_matrix) %*% jacobian
    jacobian_residuals_eigen <- eigen(jacobian_residuals)
    jacobian_residuals_eig_order <- order(Re(jacobian_residuals_eigen$values),decreasing=T)
    jacobian_residuals_eig_ordered_values = jacobian_residuals_eigen$values[jacobian_residuals_eig_order]
    jacobian_residuals_eig_ordered_vectors = jacobian_residuals_eigen$vectors[,jacobian_residuals_eig_order]
    
    # the hermitian transpose is required for Thm. 2 of Magnus (1985)
    jacobian_residuals_conjugate <- Conj(t(jacobian_residuals))
    jacobian_residuals_conjugate_eigen <- eigen(jacobian_residuals_conjugate)
    jacobian_residuals_conjugate_eig_order <- order(Re(jacobian_residuals_conjugate_eigen$values),decreasing=T)
    jacobian_residuals_conjugate_eig_ordered_values = jacobian_residuals_conjugate_eigen$values[jacobian_residuals_conjugate_eig_order]
    jacobian_residuals_conjugate_eig_ordered_vectors = jacobian_residuals_conjugate_eigen$vectors[,jacobian_residuals_conjugate_eig_order]
    
    # these lists will contain the matrix of partial derivatives wrt each spline coefficient per direction
    jacobian_partial_list_x <- list()
    jacobian_residuals_partial_list_x <- list()
    jacobian_partial_list_y <- list()
    jacobian_residuals_partial_list_y <- list()
    
    # get partials of eigen vectors wrt Jacobian
    for (i in 1:length(x_coefs)){
      dJ_dc_x = matrix(c(dx_y_vec[i],0,x_dy_vec[i],0),2,2)
      jacobian_partial_list_x[[i]] <- dJ_dc_x
      dJ_dc_y = matrix(c(0,dx_y_vec[i],0,x_dy_vec[i]),2,2)
      jacobian_partial_list_y[[i]] <- dJ_dc_y
    }
    
    dP_dJ_tensor <- list() # partial of projection matrix wrt Jacobian
    for (i in 1:2){
      for (j in 1:2){
        ij_matrix <- get_dP_fixed_dJ(i, j , projection_matrix, jacobian, 
                                     jacobian_eig_ordered_vectors, jacobian_eig_odered_values[projection_rank], projection_rank)
        dP_dJ_tensor[[(2*(i-1)+j)]] <- ij_matrix
      }
    }
    
    dP_dC_x_11 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[1]]), jacobian_partial_list_x)
    dP_dC_y_11 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[1]]), jacobian_partial_list_y)
    dP_dC_x_21 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[2]]), jacobian_partial_list_x)
    dP_dC_y_21 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[2]]), jacobian_partial_list_y)
    dP_dC_x_12 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[3]]), jacobian_partial_list_x)
    dP_dC_y_12 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[3]]), jacobian_partial_list_y)
    dP_dC_x_22 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[4]]), jacobian_partial_list_x)
    dP_dC_y_22 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[4]]), jacobian_partial_list_y)
    
    dP_dC_x <- Map(function(w,x,y,z) matrix(c(w,x,y,z),2,2),dP_dC_x_11,dP_dC_x_12,dP_dC_x_21,dP_dC_x_22)
    dP_dC_y <- Map(function(w,x,y,z) matrix(c(w,x,y,z),2,2),dP_dC_y_11,dP_dC_y_12,dP_dC_y_21,dP_dC_y_22)
    
    dPJ_dC_x <- Map(function(x,y) x %*% jacobian + projection_matrix %*% y, dP_dC_x, jacobian_partial_list_x)
    dPJ_dC_y <- Map(function(x,y) x %*% jacobian + projection_matrix %*% y, dP_dC_y, jacobian_partial_list_y)
    
    dMJ_dC_x <- Map(function(x,y) x - y, jacobian_partial_list_x, dPJ_dC_x)
    dMJ_dC_y <- Map(function(x,y) x - y, jacobian_partial_list_y, dPJ_dC_y)
    
    w_0 <- jacobian_residuals_conjugate_eig_ordered_vectors[,penalty_ranks]
    v_0 <- jacobian_residuals_eig_ordered_vectors[,penalty_ranks]
    dL_dC_x <-unlist(Map(function(x) (1/sqrt(w_0%*%v_0) * (t(w_0) %*% x %*% v_0) ),dMJ_dC_x))
    dL_dC_y <- unlist(Map(function(x) (1/sqrt(w_0%*%v_0) * (t(w_0) %*% x %*% v_0) ),dMJ_dC_y))
    
    fdx_updates[,b] <- Re(dL_dC_x)
    fdy_updates[,b] <- Re(dL_dC_y)
    
    results_array[result_index,(2*b-1):(2*b)] = jacobian_eig_odered_values # could be expanded to include eg eig vectors
  }
  
  init_fdx$coefs <- matrix(init_fdx$coefs - eta*apply(fdx_updates,1,mean), nrow = nrow(init_fdx$coefs))
  init_fdy$coefs <- matrix(init_fdy$coefs - eta*apply(fdy_updates,1,mean), nrow = nrow(init_fdy$coefs))
  
  return(list(fdx = init_fdx, fdy = init_fdy, batch_results = results_array))
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

get_dP_fixed_dJ <- function(i, j, P, J, X, eig_val, projection_rank){
  # TODO: Make more algorithmic
  dP_fixed_dJ <- matrix(rep(-5,4),ncol=2) # init matrix to return 
  if ((i == 1) & (j == 1)){
    if(projection_rank == 1){
      dP_fixed_dX <- matrix(c(2*X[1,1],0,0,0),ncol=2) 
    } else {
      dP_fixed_dX <- matrix(c(0,0,2*X[1,2],0),ncol=2) 
    }
  } else if (i != j){
    if(projection_rank ==  1){
      dP_fixed_dX <- matrix(c(X[2,1],X[1,1],0,0),ncol=2) 
    } else {
      dP_fixed_dX <- matrix(c(0,0,X[2,2],X[1,2]),ncol=2) 
    }
  } else if ((i == 2) & (j == 2)){
    if(projection_rank == 1){
      dP_fixed_dX <- matrix(c(0,2*X[2,1],0,0),ncol=2) 
    } else {
      dP_fixed_dX <- matrix(c(0,0,0,2*X[2,2]),ncol=2) 
    }
  }
  
  pseudo_inverse <- MASS::ginv(eig_val*diag(2) - J)
  
  for (i in 1:2){
    for (j in 1:2){
      selection_matrix <- matrix(replace(numeric(4), 2*(i-1)+j, 1), ncol = 2)
      dX_dJij <- pseudo_inverse %*% selection_matrix %*% X
      dP_fixed_dJ[i, j] <- sum(diag(t(dP_fixed_dX) %*% dX_dJij))
    }
  }
  
  return(dP_fixed_dJ)
}


analyze_gd_convergence <- function(samples, batch_params, batch_program, gd_log, num_updates){
  sample_data <- samples[t(batch_program),] # transpose means samples are grouped by batch
  sample_data$cycle <- rep(1:num_updates, each = (ncol(batch_program)*nrow(batch_program)/num_updates))
  sample_data$batch <- rep(rep(1:(nrow(batch_program)/num_updates),each=ncol(batch_program)),num_updates)
  eig_tibble <- tibble(lambda1 = c(t(gd_log[,c(TRUE,FALSE)])),
                       lambda2 = c(t(gd_log[,c(FALSE,TRUE)])))
  results_tibble <- cbind(as_tibble(sample_data),eig_tibble)
  results_tibble <- mutate(results_tibble, neg_eig = (Re(lambda1) < 0))
  
  title_str <- paste0("Alg:",batch_params$algorithm," Batching:",batch_params$method,
                      " Iterations:",batch_params$num_updates," Batch Size:",batch_params$batch_size,
                      " Batches/Iteration:",batch_params$cycle_size)
  
  eig_violin <- ggplot(results_tibble, aes(x=as.factor(cycle),y=Re(lambda1))) +
                  geom_violin(draw_quantiles = c(0.5)) +
                  labs(title="Distribution of Re(Lambda1) per Iteration",
                       x = "Iteration", y = "Re(Lambda1)", subtitle = title_str)
                  #scale_color_viridis() 
  print(eig_violin)
  
  eig_scatter <- ggplot(results_tibble, aes(x=x,y=y,color=Re(lambda1))) +
    geom_point() +
    scale_color_viridis()  +
    facet_wrap(~cycle) + 
    labs(title="1st Eigenvalue per Iteration", x = "x", y = "y", subtitle = title_str)
  print(eig_scatter)
  
  eig_scatter <- ggplot(results_tibble, aes(x=x,y=y,color=(Re(lambda1)<0))) +
    geom_point() +
    #scale_color_viridis()  +
    facet_wrap(~cycle) + 
    labs(title="1st Eigenvalue Sign per Iteration", x = "x", y = "y", subtitle = title_str)
  print(eig_scatter)
  
  eig_scatter <- ggplot(results_tibble, aes(x=x,y=y,color=Re(lambda2))) +
    geom_point() +
    scale_color_viridis()  +
    facet_wrap(~cycle) + 
    labs(title="2nd Eigenvalue per Iteration", x = "x", y = "y", subtitle = title_str)
  print(eig_scatter)
  
  eig_scatter <- ggplot(results_tibble, aes(x=x,y=y,color=(Re(lambda2)<0))) +
    geom_point() +
    #scale_color_viridis()  +
    facet_wrap(~cycle) + 
    labs(title="2nd Eigenvalue Sign per Iteration", x = "x", y = "y", subtitle = title_str)
  print(eig_scatter)
  
  #batch_map <- as_tibble(rep(1,num_updates) %x% batch_program)
  #batch_map$cycle <- rep(1:num_updates, each = nrow(batch_program))
  #gd_tibble <- as_tibble(gd_log)
  #View(batch_map)
  
  # plot precision graph
  # rolling_frac <- zoo::rollsum(gd_updated_splines$is_neg,1000)
  # precision_plot <- ggplot(data=tibble(x=1000:(999+length(rolling_frac)),y=rolling_frac),aes(x=x,y=rolling_frac/1000)) +
  #   geom_line() + 
  #   labs(title="Precision of Gradient Descent",y="% of Last 1000 Gradient Descent Samples with Negative Maximum Real Eigenvalue",x="Iteration")
  # print(precision_plot)
  
  return(NA)
}

##############
# Deprecated #
##############

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
    max_eig_outer <- outer(ordered_vectors[,2],ordered_vectors[,2]) # Changed for lambda_2
    lambda_step.x <- Re(max_eig_outer[1,1]) * dx_y + Re(max_eig_outer[1,2]) * x_dy
    lambda_step.y <- Re(max_eig_outer[2,1]) * dx_y + Re(max_eig_outer[2,2]) * x_dy
    spline.fd_x$coefs <- matrix(spline.fd_x$coefs - eta*lambda_step.x, nrow = nrow(spline.fd_x$coefs))
    spline.fd_y$coefs <- matrix(spline.fd_y$coefs - eta*lambda_step.y, nrow = nrow(spline.fd_y$coefs))
  }
  
  updated_fd <- list(fdx = spline.fd_x, fdy = spline.fd_y, is_neg = is_neg)
  return(updated_fd)
}


##############
# Deprecated #
##############

# random_point_jacobian_gd <- function(init_spline.fd_x,init_spline.fd_y,gd_points,max_steps=5000,eta=0.002){
#   # assumes same basis for f_x and f_y at the moment
#   xbasis.eval = eval.basis(gd_points[,1],init_spline.fd_x$sbasis)
#   ybasis.eval = eval.basis(gd_points[,2],init_spline.fd_x$tbasis)
#   xbasis.deriv = eval.basis(gd_points[,1],init_spline.fd_x$sbasis, Lfdobj = 1) # get derivative of basis functions
#   ybasis.deriv = eval.basis(gd_points[,2],init_spline.fd_x$tbasis, Lfdobj = 1) 
#   
#   spline.fd_x = init_spline.fd_x
#   spline.fd_y = init_spline.fd_y
#   is_neg <- c()
#   for(i in 1:max_steps){
#     sample_index <- sample(nrow(gd_points),1)
#     linearization_point <- gd_points[sample_index,]
#     dx_y = c(outer(xbasis.deriv[sample_index,],ybasis.eval[sample_index,]))
#     x_dy = c(outer(xbasis.eval[sample_index,], ybasis.deriv[sample_index,]))
#     x_coefs <- c(spline.fd_x$coefs)
#     y_coefs <- c(spline.fd_y$coefs)
#     jacobian <- matrix(c(x_coefs%*%dx_y,y_coefs%*%dx_y,x_coefs%*%x_dy,y_coefs%*%x_dy),nrow=2)
#     jacobian_eigen <- eigen(jacobian)
#     eig_order <- order(Re(jacobian_eigen$values),decreasing=T)
#     ordered_values = jacobian_eigen$values[eig_order]
#     ifelse(Re(ordered_values[1])<0,is_neg <- append(is_neg,1),is_neg <- append(is_neg,0))
#     ordered_vectors = jacobian_eigen$vectors[,eig_order]
#     max_eig_outer <- outer(ordered_vectors[,2],ordered_vectors[,2]) # Changed for lambda_2
#     lambda_step.x <- Re(max_eig_outer[1,1]) * dx_y + Re(max_eig_outer[1,2]) * x_dy
#     lambda_step.y <- Re(max_eig_outer[2,1]) * dx_y + Re(max_eig_outer[2,2]) * x_dy
#     spline.fd_x$coefs <- matrix(spline.fd_x$coefs - eta*lambda_step.x, nrow = nrow(spline.fd_x$coefs))
#     spline.fd_y$coefs <- matrix(spline.fd_y$coefs - eta*lambda_step.y, nrow = nrow(spline.fd_y$coefs))
#   }
#   
#   updated_fd <- list(fdx = spline.fd_x, fdy = spline.fd_y, is_neg = is_neg)
#   return(updated_fd)
# }


# # projects off eigenvectors of Jacobian, which doesn't make much sense...
# run_eig_projection_gd <- function(eig_rank, init_fdx, init_fdy, basis_list, batch_indices, result_index, results_array, eta){
#   # this algorithm calculates the eigenvalue penalty of the Jacobian of residuals from the largest eigenvector 
#   projection_rank <- eig_rank[1] # direction to project off
#   penalty_ranks <- eig_rank[2] # direction to penalize in residuals
#   fdx_updates <- matrix(0, nrow = nrow(init_fdx$coefs)*nrow(init_fdy$coefs), ncol = length(batch_indices))
#   fdy_updates <- matrix(0, nrow = nrow(init_fdx$coefs)*nrow(init_fdy$coefs), ncol = length(batch_indices))
#   
#   for(b in 1:length(batch_indices)){
#     sample_index <- batch_indices[b]
#     # these are now matrices!
#     dx_y_matrix = outer(basis_list$x1[sample_index,],basis_list$y0[sample_index,])
#     x_dy_matrix = outer(basis_list$x0[sample_index,], basis_list$y1[sample_index,])
#     dx_y_vec = c(dx_y_matrix)
#     x_dy_vec = c(x_dy_matrix)
#     x_coefs <- c(init_fdx$coefs)
#     y_coefs <- c(init_fdy$coefs)
#     
#     # get eigen-stuff from the original Jacobian
#     jacobian <- matrix(c(x_coefs%*%c(dx_y_matrix),y_coefs%*%c(dx_y_matrix),x_coefs%*%c(x_dy_matrix),y_coefs%*%c(x_dy_matrix)),nrow=2)
#     jacobian_eigen <- eigen(jacobian)
#     jacobian_eig_order <- order(Re(jacobian_eigen$values),decreasing=T)
#     jacobian_eig_odered_values = jacobian_eigen$values[jacobian_eig_order]
#     jacobian_eig_ordered_vectors = jacobian_eigen$vectors[,jacobian_eig_order]
#     
#     # then for the residualized Jacobian
#     projection_matrix <- outer(jacobian_eig_ordered_vectors[,projection_rank],jacobian_eig_ordered_vectors[,projection_rank])
#     jacobian_residuals <- (diag(2) -  projection_matrix) %*% jacobian
#     jacobian_residuals_eigen <- eigen(jacobian_residuals)
#     jacobian_residuals_eig_order <- order(Re(jacobian_residuals_eigen$values),decreasing=T)
#     jacobian_residuals_eig_ordered_values = jacobian_residuals_eigen$values[jacobian_residuals_eig_order]
#     jacobian_residuals_eig_ordered_vectors = jacobian_residuals_eigen$vectors[,jacobian_residuals_eig_order]
#     
#     # the hermitian transpose is required for Thm. 2 of Magnus (1985)
#     jacobian_residuals_conjugate <- Conj(t(jacobian_residuals))
#     jacobian_residuals_conjugate_eigen <- eigen(jacobian_residuals_conjugate)
#     jacobian_residuals_conjugate_eig_order <- order(Re(jacobian_residuals_conjugate_eigen$values),decreasing=T)
#     jacobian_residuals_conjugate_eig_ordered_values = jacobian_residuals_conjugate_eigen$values[jacobian_residuals_conjugate_eig_order]
#     jacobian_residuals_conjugate_eig_ordered_vectors = jacobian_residuals_conjugate_eigen$vectors[,jacobian_residuals_conjugate_eig_order]
#     
#     # these lists will contain the matrix of partial derivatives wrt each spline coefficient per direction
#     jacobian_partial_list_x <- list()
#     jacobian_residuals_partial_list_x <- list()
#     jacobian_partial_list_y <- list()
#     jacobian_residuals_partial_list_y <- list()
#     
#     # get partials of eigen vectors wrt Jacobian
#     for (i in 1:length(x_coefs)){
#       dJ_dc_x = matrix(c(dx_y_vec[i],0,x_dy_vec[i],0),2,2)
#       jacobian_partial_list_x[[i]] <- dJ_dc_x
#       dJ_dc_y = matrix(c(0,dx_y_vec[i],0,x_dy_vec[i]),2,2)
#       jacobian_partial_list_y[[i]] <- dJ_dc_y
#     }
#     
#     dP_dJ_tensor <- list() # partial of projection matrix wrt Jacobian
#     for (i in 1:2){
#       for (j in 1:2){
#         ij_matrix <- get_dP_fixed_dJ(i, j , projection_matrix, jacobian, 
#                                      jacobian_eig_ordered_vectors, jacobian_eig_odered_values[projection_rank], projection_rank)
#         dP_dJ_tensor[[(2*(i-1)+j)]] <- ij_matrix
#       }
#     }
#     
#     dP_dC_x_11 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[1]]), jacobian_partial_list_x)
#     dP_dC_y_11 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[1]]), jacobian_partial_list_y)
#     dP_dC_x_21 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[2]]), jacobian_partial_list_x)
#     dP_dC_y_21 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[2]]), jacobian_partial_list_y)
#     dP_dC_x_12 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[3]]), jacobian_partial_list_x)
#     dP_dC_y_12 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[3]]), jacobian_partial_list_y)
#     dP_dC_x_22 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[4]]), jacobian_partial_list_x)
#     dP_dC_y_22 <- Map(function(x,y)  sum(diag(t(x)%*%y)), list(dP_dJ_tensor[[4]]), jacobian_partial_list_y)
#     
#     dP_dC_x <- Map(function(w,x,y,z) matrix(c(w,x,y,z),2,2),dP_dC_x_11,dP_dC_x_12,dP_dC_x_21,dP_dC_x_22)
#     dP_dC_y <- Map(function(w,x,y,z) matrix(c(w,x,y,z),2,2),dP_dC_y_11,dP_dC_y_12,dP_dC_y_21,dP_dC_y_22)
#     
#     dPJ_dC_x <- Map(function(x,y) x %*% jacobian + projection_matrix %*% y, dP_dC_x, jacobian_partial_list_x)
#     dPJ_dC_y <- Map(function(x,y) x %*% jacobian + projection_matrix %*% y, dP_dC_y, jacobian_partial_list_y)
#     
#     dMJ_dC_x <- Map(function(x,y) x - y, jacobian_partial_list_x, dPJ_dC_x)
#     dMJ_dC_y <- Map(function(x,y) x - y, jacobian_partial_list_y, dPJ_dC_y)
#     
#     w_0 <- jacobian_residuals_conjugate_eig_ordered_vectors[,penalty_ranks]
#     v_0 <- jacobian_residuals_eig_ordered_vectors[,penalty_ranks]
#     dL_dC_x <-unlist(Map(function(x) (1/sqrt(w_0%*%v_0) * (t(w_0) %*% x %*% v_0) ),dMJ_dC_x))
#     dL_dC_y <- unlist(Map(function(x) (1/sqrt(w_0%*%v_0) * (t(w_0) %*% x %*% v_0) ),dMJ_dC_y))
#     
#     fdx_updates[,b] <- Re(dL_dC_x)
#     fdy_updates[,b] <- Re(dL_dC_y)
#     
#     results_array[result_index,(2*b-1):(2*b)] = jacobian_eig_odered_values # could be expanded to include eg eig vectors
#   }
#   
#   init_fdx$coefs <- matrix(init_fdx$coefs - eta*apply(fdx_updates,1,mean), nrow = nrow(init_fdx$coefs))
#   init_fdy$coefs <- matrix(init_fdy$coefs - eta*apply(fdy_updates,1,mean), nrow = nrow(init_fdy$coefs))
#   
#   return(list(fdx = init_fdx, fdy = init_fdy, batch_results = results_array))
# }

#function(eig_rank, init_fdx, init_fdy, basis_list, batch_indices, result_index, results_array, eta){
# run_gd <- function(gd_data, eig_rank, init_fdx, init_fdy, batch_params, eta, batch_num){ 
#   fdx_updates <- matrix(0, nrow = nrow(init_fdx$coefs)*nrow(init_fdy$coefs), ncol = length(batch_indices))
#   fdy_updates <- matrix(0, nrow = nrow(init_fdx$coefs)*nrow(init_fdy$coefs), ncol = length(batch_indices))
#   
#   for(i in 1:length(batch_indices)){
#     sample_index <- batch_indices[i]
#     dx_y = c(outer(basis_list$x1[sample_index,],basis_list$y0[sample_index,]))
#     x_dy = c(outer(basis_list$x0[sample_index,], basis_list$y1[sample_index,]))
#     x_coefs <- c(init_fdx$coefs)
#     y_coefs <- c(init_fdy$coefs)
#     
#     jacobian <- matrix(c(x_coefs%*%dx_y,y_coefs%*%dx_y,x_coefs%*%x_dy,y_coefs%*%x_dy),nrow=2)
#     jacobian_eigen <- eigen(jacobian)
#     eig_order <- order(Re(jacobian_eigen$values),decreasing=T)
#     ordered_values = jacobian_eigen$values[eig_order]
#     ordered_vectors = jacobian_eigen$vectors[,eig_order]
#     
#     max_eig_outer <- outer(ordered_vectors[,eig_rank],ordered_vectors[,eig_rank])
#     lambda_step.x <- Re(max_eig_outer[1,1]) * dx_y + Re(max_eig_outer[1,2]) * x_dy
#     lambda_step.y <- Re(max_eig_outer[2,1]) * dx_y + Re(max_eig_outer[2,2]) * x_dy
#     fdx_updates[,i] <- lambda_step.x
#     fdy_updates[,i] <- lambda_step.y
#     
#     results_array[result_index,(2*i-1):(2*i)] = ordered_values # could be expanded to include eg eig vectors
#   }
# 
#   init_fdx$coefs <- matrix(init_fdx$coefs - eta*apply(fdx_updates,1,mean), nrow = nrow(init_fdx$coefs))
#   init_fdy$coefs <- matrix(init_fdy$coefs - eta*apply(fdy_updates,1,mean), nrow = nrow(init_fdy$coefs))
#   
#   return(list(fdx = init_fdx, fdy = init_fdy, batch_results = results_array))
# }
