library(nabor) # for kNN

average_nn_results <- function(nn_indices, data.matrix){
  nearest_neighbors <- matrix(data.matrix[nn_indices,], ncol = 2)
  averaged_result <- c(apply(nearest_neighbors, 2, mean))
  return(averaged_result)
}
eval_knn_fit <- function(query_matrix, data_matrix, k){
  data_coords <- data_matrix[,c(1,2)]
  data_grad <- data_matrix[,c(3,4)]
  knn_result <- knn(data_coords, query_matrix, k)
  nearest_neighbors <- t(apply(knn_result$nn.idx, 1, average_nn_results, data.matrix = data_grad)) #data_grad[knn_result$nn.idx,]
  return(nearest_neighbors)
}