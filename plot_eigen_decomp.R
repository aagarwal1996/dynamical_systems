here::i_am("plot_eigen_decomp.R")

library(dplyr)
library(ggplot2)
library(ggquiver)

get_jacobian <- function(gradient_field_fn, vector){
	if(is.null(gradient_field_fn)){
		return(NULL)
	}
	# Helper function to get the jacobian matrix at a given state
	jacobian_fn <- gradient_field_fn$jacobian
	jacobian <- jacobian_fn(vector)
	return(jacobian)
}

plot_eigen_decomp <- function(gradient_field_fn, limits,
	vector_density = 10, sample_df = NULL, title = ""){

	# Set up a grid of vectors for `ggquiver`
	x_partition <- seq(limits$x_lower,limits$x_upper,length.out=vector_density)
	y_partition <- seq(limits$y_lower,limits$y_upper,length.out=vector_density)
	vector_locations <- expand.grid(x=x_partition, y=y_partition)
	
	# TODO: Extract the Jacobian matrix at each vector
	vector_jacobians <- apply(vector_locations,1,function(vec){
		return(list(location = vec, jacobian = get_jacobian(gradient_field_fn,vec)))
	})

	# Apply the eigen-decomposition
	vector_decomps <- lapply(vector_jacobians,function(vec){
		return(list(location = vec$location,
					eig = eigen(vec$jacobian)))
	})
	
	# Convert the eigenvectors into `ggquiver` form
	vector_dfs <- lapply(vector_decomps,function(vec){
		first_eig_vec_real <- unlist(c(vec$location,
			0.1*Re(vec$eig$values[1])*Re(vec$eig$vectors[,1])))
		names(first_eig_vec_real) <- c("x","y","u","v")
		second_eig_vec_real <- unlist(c(vec$location,
			1*Re(vec$eig$values[2])*Re(vec$eig$vectors[,2])))
		names(second_eig_vec_real) <- c("x","y","u","v")
		vec_df <- rbind(first_eig_vec_real,second_eig_vec_real)
		return(list(vec_df))
	})
	vector_dfs <- unlist(vector_dfs,recursive=FALSE)
	vector_field_df <- as.data.frame(do.call(rbind,vector_dfs))
	vector_field_df$eig_rank <- as.factor(ifelse(
		grepl("first", rownames(vector_field_df), fixed = TRUE),"Eig1","Eig2"))
	rownames(vector_field_df) <- NULL # Remove rownames
	
	vector_field <- ggplot2::ggplot(vector_field_df,aes(x=x,y=y,u=u,v=v)) +
		ggquiver::geom_quiver(aes(color=eig_rank)) +
		labs(title=title)
	if(!is.null(sample_df)){
		sample_df$u = 0
		sample_df$v = 0
		vector_field <- vector_field + 
			ggplot2::geom_point(data=sample_df,aes(x=x,y=y,u=u,v=v))
	}
	return(vector_field)
}

plot_eigen_decomp_delta <- function(gradient_field_fn_1, gradient_field_fn_2, limits,
							  vector_density = 10, sample_df = NULL, title = ""){
	
	# Set up a grid of vectors for `ggquiver`
	x_partition <- seq(limits$x_lower,limits$x_upper,length.out=vector_density)
	y_partition <- seq(limits$y_lower,limits$y_upper,length.out=vector_density)
	vector_locations <- expand.grid(x=x_partition, y=y_partition)
	
	fields <- list(gradient_field_fn_1,gradient_field_fn_2)
	vector_field_list <- lapply(fields, function(gradient_field_fn){
		
		vector_jacobians <- apply(vector_locations,1,function(vec){
			return(list(location = vec,
				jacobian = get_jacobian(gradient_field_fn,vec)))
		})

		# Apply the eigen-decomposition
		vector_decomps <- lapply(vector_jacobians,function(vec){
			return(list(location = vec$location, eig = eigen(vec$jacobian)))
		})

		# Convert the eigenvectors into `ggquiver` form
		vector_dfs <- lapply(vector_decomps,function(vec){
			first_eig_vec_real <- unlist(c(vec$location,
				Re(vec$eig$values[1])*Re(vec$eig$vectors[,1])))
			names(first_eig_vec_real) <- c("x","y","u","v")
			second_eig_vec_real <- unlist(c(vec$location,
				Re(vec$eig$values[2])*Re(vec$eig$vectors[,2])))
			names(second_eig_vec_real) <- c("x","y","u","v")
			vec_df <- rbind(first_eig_vec_real,second_eig_vec_real)
			return(list(vec_df))
		})
		
		vector_dfs <- unlist(vector_dfs,recursive=FALSE)
		vector_field_df <- as.data.frame(do.call(rbind,vector_dfs))
		vector_field_df$eig_rank <- as.factor(ifelse(
			grepl("first", rownames(vector_field_df), fixed = TRUE),"Eig1","Eig2"))
		rownames(vector_field_df) <- NULL # Remove rownames
		return(vector_field_df)
	})
	
	vector_field_df_delta <- vector_field_list[[1]]
	vector_field_df_delta$u = vector_field_list[[1]]$u - vector_field_list[[2]]$u
	vector_field_df_delta$v = vector_field_list[[1]]$v - vector_field_list[[2]]$v
	vector_field <- ggplot2::ggplot(vector_field_df_delta,aes(x=x,y=y,u=u,v=v)) +
		ggquiver::geom_quiver(aes(color=eig_rank)) +
		labs(title=title)
	if(!is.null(sample_df)){
		sample_df$u = 0
		sample_df$v = 0
		vector_field <- vector_field + 
			ggplot2::geom_point(data=sample_df,aes(x=x,y=y,u=u,v=v))
	}
	return(vector_field)
}

#small_sample_df <- data.frame(x=c(0,0.25,0.5),y=c(0,0.25,0.5),u=rep(0,3),v=rep(0,3))
#limits_list <- list(x_lower = -1, x_upper = 1, y_lower = -1, y_upper = 1)
#return_plot <- plot_eigen_decomp(NULL,limits_list,sample_df = small_sample_df)
#return_plot <- plot_eigen_decomp(vdp_trial_fit$replicates[[1]]$estimators[[3]],limits_list)
#return_plot