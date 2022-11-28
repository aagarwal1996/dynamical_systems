## Proof of concept for derivative of jacobian wrt max eigenvalue

eig_gd <- function(J_init, steps = 1000, eta = 0.1){
  J_matrix = J_init
  jacobian_eig = eigen(J_matrix)
  eig_order <- order(Re(jacobian_eig$values),decreasing=T)
  ordered_values = jacobian_eig$values[eig_order]
  ordered_vectors = jacobian_eig$vectors[,eig_order]
  
  for(i in 1:steps){
    v_1 = ordered_vectors[,1]
    right_eig = solve(t(ordered_vectors))[,1]
    gradient = outer(v_1, right_eig)
    print(gradient)
    J_matrix = J_matrix - eta * gradient # + matrix(rnorm(4, sd = 0.1),nrow=2)
    
    jacobian_eig = eigen(J_matrix)
    eig_order <- order(Re(jacobian_eig$values),decreasing=T)
    ordered_values = jacobian_eig$values[eig_order]
    ordered_vectors = jacobian_eig$vectors[,eig_order]
    
    new_max_eig <- ordered_values[1]
    print(ordered_values)
    print(paste0("Iteration ", i, " Max Eig: ", new_max_eig))
    if(Re(new_max_eig) < 0){break}
  }
  
  return(J_matrix)
}

test_jacobian <- matrix(rnorm(4),nrow=2)
eigs <- eigen(test_jacobian)
eig_order <- order(Re(eigs$values))
ordered_values = eigs$values[eig_order]
ordered_vectors = eigs$vectors[,eig_order]
print(eig_order)
print(test_jacobian)
print(eigen(test_jacobian))
print(eigen(test_jacobian)$vectors %*% diag(eigen(test_jacobian)$values) %*% solve(eigen(test_jacobian)$vectors))

test_jacobian <- matrix(rnorm(4),nrow=2)
result <- eig_gd(test_jacobian)
eigen(result)