library(fda)
library(deSolve)

generate_bspline_basis <- function(data, x_grid, y_grid, norder = 4, nbasis = 12){
  xbasis = create.bspline.basis(range=c(min(x_grid),max(x_grid)),norder=norder,nbasis=nbasis)
  ybasis = create.bspline.basis(range=c(min(y_grid),max(y_grid)),norder=norder,nbasis=nbasis)
  xbasis.vals = eval.basis(data$x,xbasis)
  ybasis.vals = eval.basis(data$y,ybasis)
  
  spline_fit_list <- list(xbasis = xbasis, ybasis = ybasis, 
                          xbasis.eval = xbasis.vals, ybasis.eval = ybasis.vals)
  return(spline_fit_list)
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
  
  
  traj = lsoda(c(data$x[1],data$y[1]),seq(0,100,by=0.1),sVdP,0,sfd_x=sfd_x,sfd_y=sfd_y)
  
  output_list <- list(x_grad_bifd = sfd_x, x_grad_eval = smat_x, y_grad_bifd = sfd_y, y_grad_eval = smat_y, 
                      trajectory = traj, second.deriv_penalty = allPen)
  return(output_list)
}


sVdP = function(t,x,p,sfd_x,sfd_y){
  # Helper function for `lsoda` which evaluates derivatives of a spline fit
  dx = eval.bifd(x[1],x[2],sfd_x)
  dy = eval.bifd(x[1],x[2],sfd_y)
  return(list(as.vector(c(dx,dy))))
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
