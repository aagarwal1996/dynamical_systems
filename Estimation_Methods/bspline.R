library(fda)
library(deSolve)
library(Rcpp)
library(ggplot2)
library(superheat)
sourceCpp(here::here('Estimation_Methods/bspline.cpp'))

calculate_spline_gradient_field <- function(data, x_grid, y_grid, side_info = list(), norder = 6, nbasis = 12,
                                            penalty_order = 2, lambda = 1e-8, plot_penalty = FALSE){
  for (i in 1:length(side_info)){
    if (length(side_info)==0){ break }
    if (side_info[[i]]$name=="boundary_box"){
      reps = 1000
      x_min <- min(x_grid)
      x_max <- max(x_grid)
      y_min <- min(y_grid)
      y_max <- max(y_grid)
      
      h_grid <- seq(x_min,x_max,length.out=reps)
      v_grid <- seq(y_min,y_max,length.out=reps)
      
      boundary_tibble <- as.matrix(tibble(x = c(h_grid,h_grid,rep(x_min,reps),rep(x_max,reps)),
                                y= c(rep(y_min,reps),rep(y_max,reps),v_grid,v_grid),
                                f_x = c(rep(0,reps),rep(0,reps),rep(1,reps),rep(-1,reps)),
                                f_y = c(rep(1,reps),rep(-1,reps),rep(0,reps),rep(0,reps))))
      data <- rbind(data, boundary_tibble)
    } else if (side_info[[i]]$name == "stable_fp"){
      # Currently Broken
      fixed_point = side_info[[i]]$fp
      x_vals = fixed_point[1] + 0.2*cos(seq(0,2*pi,length.out=250))
      y_vals = fixed_point[2] + 0.2*sin(seq(0,2*pi,length.out=250))
      fp_tibble <- as.matrix(tibble(x = x_vals, y = y_vals,
                     f_x = x_vals - fixed_point[1], f_y = y_vals - fixed_point[2]))
      data <- rbind(data, fp_tibble)
    } else if (side_info[[i]]$name == "delta_ring"){
      # expects list with items (name, radius, sample_frac, vector_magnitude)
      # for now we sample data points
      delta_ring_samples <- data[sort(sample(nrow(data),size=floor(nrow(data)*side_info[[i]]$sample_frac),replace=FALSE)),]
      delta_ring_sideinfo <- t(apply(delta_ring_samples,1,get_normal_vectors, radius = side_info[[i]]$radius))
      side_information_pos <- rbind(delta_ring_sideinfo[,c(1,2)],delta_ring_sideinfo[,c(3,4)])
      colnames(side_information_pos) <- c("x","y")
      #colnames(delta_ring_sideinfo) <- c("out.root.x", "out.root.y", "in.root.x", "in.root.y")
      shifted_data <- rbind(delta_ring_samples[(side_info[[i]]$shift+1):nrow(delta_ring_samples),], 
                            delta_ring_samples[0:(side_info[[i]]$shift),])
      #View(delta_ring_samples)
      #View(shifted_data)
      #View(delta_ring_sideinfo)
      side_information_grad <- rbind(shifted_data[,c(1,2)] - delta_ring_sideinfo[,c(1,2)], shifted_data[,c(1,2)] - delta_ring_sideinfo[,c(3,4)])
      #View(side_information_grad)
      side_information_grad <- (sqrt(side_info[[i]]$si_magnitude)/sqrt(2))*t(apply(side_information_grad, 1, function(x) x/norm(x)))
      colnames(side_information_grad) = c("f_x","f_y")
      side_information <- cbind(side_information_pos,side_information_grad)
      #View(side_information)
      data <- rbind(data, side_information) # extend data
    } 
  }
  # evaluate b-spline basis functions at coordinates in data (N x 2 matrix)
  bspline_basis_fns <- generate_bspline_basis(data, x_grid, y_grid, norder = norder, 
                                              nbasis = nbasis, penalty_order = penalty_order)
  
  if (!(lambda==0)){
    bspline_fit_coeffs <- fit_bsplines_cpp(bspline_basis_fns$xbasis.eval, bspline_basis_fns$ybasis.eval,
                                    bspline_basis_fns$xpenalty, bspline_basis_fns$ypenalty,
                                    data$f_x, data$f_y, lambda)

  } else{
    bspline_fit_coeffs_x <- smooth.basis(argvals=1:n, y, fdParobj)$fd$coeff
  }
  
  # create bivariate functional data objects for our fit splines 
  spline.fd_x <- bifd(t(matrix(bspline_fit_coeffs[,1],nbasis,nbasis)), 
                      bspline_basis_fns$xbasis,  bspline_basis_fns$ybasis)
  spline.fd_y <- bifd(t(matrix(bspline_fit_coeffs[,2],nbasis,nbasis)),
                      bspline_basis_fns$xbasis,  bspline_basis_fns$ybasis)
  
  # evaluate over user-specified grid
  spline_grid_x <- eval.bifd(x_grid,y_grid,spline.fd_x)
  spline_grid_y <- eval.bifd(x_grid,y_grid,spline.fd_y)
  
  # return as |Grid| x 2 matrix
  spline_field <- matrix(c(c(spline_grid_x),c(spline_grid_y)),ncol=2)
  
  if (plot_penalty){
    plot_spline_penalty(x_grid,y_grid, spline.fd_x, spline.fd_y, penalty_order)
  }
  
  # TODO: Currently on plots first type of side-info
  if (length(side_info)>0){
    if (side_info[[1]]$name == "delta_ring"){
      plot_ring_df <- side_information
      if (length(side_info[[1]]$lc_params) > 1){
        if (side_info[[1]]$lc_params$system == "van_der_pol"){
          # true field only works for VdP at the moment
          true_ring_grad <- t(apply(side_information[,c("x","y")], 1, van_der_pol_gradient_helper, mu = side_info[[i]]$lc_params$params$mu))
          true_ring_grad <- cbind(side_information[,c("x","y")],true_ring_grad)
          colnames(true_ring_grad) <- c("x","y","f_x","f_y")
          
          # convert to tibble and add source label
          plot_ring_df <- as_tibble(plot_ring_df) %>% mutate(source = "SI")
          true_ring_grad <- as_tibble(true_ring_grad) %>% mutate(source = "Truth")
          plot_ring_tibble <- bind_rows(true_ring_grad, plot_ring_df) %>% mutate(source = as.factor(source))
        } else {
          plot_ring_tibble <- as_tibble(plot_ring_df) %>% mutate(source = "SI")
        }
      } else {
        plot_ring_tibble <- as_tibble(plot_ring_df) %>% mutate(source = "SI")
      }
      ring_plot <- ggplot(plot_ring_tibble, aes(x = x, y = y, u = f_x, v = f_y, color = source)) +
        geom_quiver(aes(x = x, y = y, u = f_x, v = f_y),alpha = 0.5) +
        geom_point(data = plot_ring_tibble, aes(x = x, y = y), color = "red")
      print(ring_plot)
    }
  }

  return(list(field = spline_field, fdx = spline.fd_x, fdy = spline.fd_y))
}

get_normal_vectors <- function(sample, radius){
  rotation_matrix <- matrix(c(0,-1,1,0),ncol=2)
  magnitude = sqrt(sample["f_x"]^2 + sample["f_y"]^2)
  point <- c(sample["x"], sample["y"])
  normal_vec_out <- (rotation_matrix %*% c(sample["f_x"], sample["f_y"]))/magnitude

  ring_sample <- c(point + radius*normal_vec_out, point - radius*normal_vec_out) 
  # normal out x, y; normal in x, y
  return(ring_sample)
}
  
plot_spline_penalty <- function(x_grid, y_grid, spline_bifd_x, spline_bifd_y, penalty_order = 2){
  # plots nth-derivative penalty at each point in eval grid
  # get penalty at each point in grid
  
  #bspline_penalty <- get_bspline_penalty_cpp(xbasis, ybasis, xpen, ypen, xdata, ydata)
  #print(bspline_penalty)
  # build plots
  
  # eval.basis(x_grid,spline_bifd_x$sbasis)
  sbasis_penalty_x <- getbasismatrix(x_grid,spline_bifd_x$sbasis, penalty_order)
  tbasis_penalty_x <- getbasismatrix(y_grid,spline_bifd_x$tbasis, penalty_order)
  grid_penalty_x <- abs(sbasis_penalty_x %*% spline_bifd_x$coef %*% t(tbasis_penalty_x))
  # prep for plotting
  grid_penalty_x <- t(grid_penalty_x)
  rownames(grid_penalty_x) <- round(y_grid,2)
  colnames(grid_penalty_x) <- round(x_grid,2)
  superheat(grid_penalty_x,
            bottom.label.text.angle = 270,
            title = "Second Derivative Penalty for x Over Grid")
  
  sbasis_penalty_y <- getbasismatrix(x_grid,spline_bifd_y$sbasis, penalty_order)
  tbasis_penalty_y <- getbasismatrix(y_grid,spline_bifd_y$tbasis, penalty_order)
  grid_penalty_y <- abs(sbasis_penalty_y %*% spline_bifd_y$coef %*% t(tbasis_penalty_y))
  grid_penalty_y <- t(grid_penalty_y)
  rownames(grid_penalty_y) <- round(y_grid,2)
  colnames(grid_penalty_y) <- round(x_grid,2)
  superheat(grid_penalty_y,
            bottom.label.text.angle = 270,
            title = "Second Derivative Penalty for y Over Grid")
  
  #eval_grid <- unname(as.matrix(expand.grid(x_grid,y_grid)))
  #long_penalty_x <- cbind(eval_grid, c(t(grid_penalty_x)))
  #colnames(long_penalty_x) <- c("x", "y", "penalty")
  #x_penalty_plot <- ggplot(data.frame(long_penalty_x), aes(x, y, fill=penalty)) + 
  #  geom_raster() +
  #  labs(title="Second Derivative Penalty for x Over Grid")
  #print(x_penalty_plot)
  #long_penalty_y <- cbind(eval_grid, c(t(grid_penalty_y)))
  #colnames(long_penalty_y) <- c("x", "y", "penalty")
  #y_penalty_plot <- ggplot(data.frame(long_penalty_y), aes(x, y, fill=penalty)) + 
  #  geom_raster() +
  #  labs(title = "Second Derivative Penalty for y Over Grid")
  #print(y_penalty_plot)
  
  return()
}

generate_bspline_basis <- function(data, x_grid, y_grid, norder = 4, nbasis = 12,
                                   penalty_order = 2){
  # number of breaks equals `nbasis` - `norder` + 2; 10 for default parameters (8 interior knots)
  xbasis = create.bspline.basis(rangeval=c(min(x_grid),max(x_grid)),norder=norder,nbasis=nbasis)
  ybasis = create.bspline.basis(rangeval=c(min(y_grid),max(y_grid)),norder=norder,nbasis=nbasis)
  xbasis.vals = eval.basis(data$x,xbasis)
  ybasis.vals = eval.basis(data$y,ybasis)
  #View(solve(t(xbasis.vals) %*% xbasis.vals) %*% t(xbasis.vals) %*% data$x)
  # get roughness penalty for each axis
  xPen = eval.penalty(xbasis, penalty_order)
  yPen = eval.penalty(ybasis, penalty_order)
  
  spline_fit_list <- list(xbasis = xbasis, ybasis = ybasis, 
                          xbasis.eval = xbasis.vals, ybasis.eval = ybasis.vals,
                          xpenalty = xPen, ypenalty = yPen)
  return(spline_fit_list)
}

bifd_spline_gradient <- function(t,x,params){
  sfd_x = params$fdx
  sfd_y = params$fdy

  # catch solution paths which leave splines' support
  x_limits <- sfd_x$sbasis$rangeval
  y_limits <- sfd_x$tbasis$rangeval

  if((x[1] < x_limits[1]) | (x[1] > x_limits[2])){
    return(list(as.vector(c(0,0))))
  }
  if((x[2] < y_limits[1]) | (x[2] > y_limits[2])){
    return(list(as.vector(c(0,0))))
  }
  
  # Helper function for `lsoda` which evaluates derivatives of a spline fit
  dx = eval.bifd(x[1],x[2],sfd_x)
  dy = eval.bifd(x[1],x[2],sfd_y)
  return(list(as.vector(c(dx,dy))))
}


### New work --------------
# TODO: Integrate side information

return_spline_gradient_fn <- function(data_object, estimator_object){

	# evaluate b-spline basis functions at coordinates in data (N x 2 matrix)
	bspline_basis_fns <- generate_bspline_basis(data, x_grid, y_grid, norder = norder, 
		nbasis = nbasis, penalty_order = penalty_order)

	if (!(lambda==0)){
		bspline_fit_coeffs <- fit_bsplines_cpp(bspline_basis_fns$xbasis.eval, bspline_basis_fns$ybasis.eval,
			bspline_basis_fns$xpenalty, bspline_basis_fns$ypenalty,
			data$f_x, data$f_y, lambda)
	} else{
		bspline_fit_coeffs_x <- smooth.basis(argvals=1:n, y, fdParobj)$fd$coeff
	}

	# create bivariate functional data objects for our fit splines 
	spline.fd_x <- bifd(t(matrix(bspline_fit_coeffs[,1],nbasis,nbasis)), 
						bspline_basis_fns$xbasis,  bspline_basis_fns$ybasis)
	spline.fd_y <- bifd(t(matrix(bspline_fit_coeffs[,2],nbasis,nbasis)),
						bspline_basis_fns$xbasis,  bspline_basis_fns$ybasis)
}
