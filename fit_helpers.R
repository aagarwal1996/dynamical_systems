library(Rcpp)
sourceCpp(here::here('Estimation_Methods','bspline.cpp'))
sourceCpp(here::here('Estimation_Methods','loess.cpp'))
sourceCpp(here::here('Estimation_Methods','nw_regression.cpp'))
source(here::here('Estimation_Methods','knn.R'))
library(MASS)

###################
# Vanilla Splines #
###################

fit_vanilla_spline_field <- function(data, side_info = list(), norder = 6, nbasis = 12,
						penalty_order = 2, lambda = 1e-8){
	for (i in 1:length(side_info)){
		if (length(side_info)==0){ break }
		if (side_info[[i]]$name=="boundary_box"){
			reps = 1000
			x_range <- max(data$x) - min(data$x)
			x_min <- min(data$x) - 0.2*x_range
			x_max <- max(data$x) + 0.2*x_range
			y_range <- max(data$y) - min(data$y)
			y_min <- min(data$y) - 0.2*y_range
			y_max <- max(data$y) + 0.2*y_range
			
			h_grid <- seq(x_min,x_max,length.out=reps)
			v_grid <- seq(y_min,y_max,length.out=reps)
			
			boundary_tibble <- as.matrix(tibble(x = c(h_grid,h_grid,rep(x_min,reps),rep(x_max,reps)),
												y= c(rep(y_min,reps),rep(y_max,reps),v_grid,v_grid),
												f_x = c(rep(0,reps),rep(0,reps),rep(1,reps),rep(-1,reps)),
												f_y = c(rep(1,reps),rep(-1,reps),rep(0,reps),rep(0,reps))))
			data <- rbind(data, boundary_tibble)
		} else if (side_info[[i]]$name == "delta_ring"){
			# expects list with items (name, radius, sample_frac, vector_magnitude)
			delta_ring_samples <- data[sort(sample(nrow(data),size=floor(nrow(data)*side_info[[i]]$sample_frac),replace=FALSE)),]
			delta_ring_sideinfo <- t(apply(delta_ring_samples,1,get_normal_vectors, radius = side_info[[i]]$radius))
			side_information_pos <- rbind(delta_ring_sideinfo[,c(1,2)],delta_ring_sideinfo[,c(3,4)])
			colnames(side_information_pos) <- c("x","y")
			
			shifted_data <- rbind(delta_ring_samples[(side_info[[i]]$shift+1):nrow(delta_ring_samples),], 
								  delta_ring_samples[0:(side_info[[i]]$shift),])
			side_information_grad <- rbind(shifted_data[,c(1,2)] - delta_ring_sideinfo[,c(1,2)], shifted_data[,c(1,2)] - delta_ring_sideinfo[,c(3,4)])
			side_information_grad <- (sqrt(side_info[[i]]$si_magnitude)/sqrt(2))*t(apply(side_information_grad, 1, function(x) x/norm(x)))
			colnames(side_information_grad) = c("f_x","f_y")
			side_information <- cbind(side_information_pos,side_information_grad)
			data <- rbind(data, side_information) # extend data
		} 
	}
	# evaluate b-spline basis functions at coordinates in data (N x 2 matrix)
	bspline_basis_fns <- generate_bspline_basis(data, norder = norder, 
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

	return(list(fdx = spline.fd_x, fdy = spline.fd_y))
}

generate_bspline_basis <- function(data, norder = 4, nbasis = 12,
								   penalty_order = 2){
	x_range <- max(data$x) - min(data$x)
	y_range <- max(data$y) - min(data$y)
	
	# number of breaks equals `nbasis` - `norder` + 2; 10 for default parameters (8 interior knots)
	xbasis = fda::create.bspline.basis(rangeval=c(min(data$x) - 0.4*x_range,max(data$x) + 0.4*x_range),norder=norder,nbasis=nbasis)
	ybasis = fda::create.bspline.basis(rangeval=c(min(data$y) - 0.4*y_range,max(data$y) + 0.4*y_range),norder=norder,nbasis=nbasis)
	xbasis.vals = fda::eval.basis(data$x,xbasis)
	ybasis.vals = fda::eval.basis(data$y,ybasis)
	
	# get roughness penalty for each axis
	xPen = fda::eval.penalty(xbasis, penalty_order)
	yPen = fda::eval.penalty(ybasis, penalty_order)
	
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

nw_get_grad <- function(t,y,params){
	# helper function to get nw gradient in deSolve::lsoda
	data <- params[[1]]
	bw_matrix <- params[[2]]
	grad_pred <- NW_regression_cpp(unname(matrix(y,nrow = 1)), as.matrix(data), bw_matrix)
	grad_pred.lsoda <- list(c(x = grad_pred[1,1], y = grad_pred[1,2])) # reformat as list with named vector
	return(grad_pred.lsoda)
}

loess_get_grad <- function(t,y,params){
	data <- params[[1]]
	h <- params[[2]]
	grad_pred <- get_loess_pred(unname(matrix(y,nrow = 1)), as.matrix(data), h)
	grad_pred.lsoda <- list(c(x = grad_pred[1,1], y = grad_pred[1,2])) # reformat as list with named vector
	return(grad_pred.lsoda)
}

knn_get_grad <- function(t,y,params){
	data <- params[[1]]
	k <- params[[2]]
	grad_pred <- eval_knn_fit(unname(matrix(y,nrow = 1)), as.matrix(data), k)
	grad_pred.lsoda <- list(c(x = grad_pred[1,1], y = grad_pred[1,2])) # reformat as list with named vector
	return(grad_pred.lsoda)
}