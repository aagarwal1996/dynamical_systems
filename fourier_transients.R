# setup
library(here)
library(fda)
library(deSolve)
library(tidyverse)
library(seewave)
library(superheat)
#library(signal)

here::here()
source(here::here("Data_Generation","data_generation.R"))
source(here::here("simulation_fns.R"))

##############################
# Fourier Analysis Functions #
##############################

get_coordinate_spectrum <- function(full_data,transient_cutoff=50,gradient=FALSE){
  # we consider the full data, the transient, and the LC
  
  x_pos <- ifelse(rep(gradient,nrow(full_data)),full_data[,"f_x"],full_data[,"x"])
  x_transient <- x_pos[1:transient_cutoff]
  x_lc <- x_pos[transient_cutoff:length(x_pos)]
  y_pos <- ifelse(rep(gradient,nrow(full_data)),full_data[,"f_y"],full_data[,"y"])
  y_transient <- y_pos[1:transient_cutoff]
  y_lc <- y_pos[transient_cutoff:length(y_pos)]
  
  # get the FFT of each
  x_pos_fft <- stats::fft(x_pos)/sqrt(length(x_pos))
  x_pos_label <- rep("x_full",length(x_pos_fft))
  x_transient_fft <- stats::fft(x_transient)/sqrt(length(x_transient)) # use level not power
  x_transient_label <- rep("x_transient",length(x_transient_fft))
  x_lc_fft <- stats::fft(x_lc)/sqrt(length(x_lc))
  x_lc_label <- rep("x_lc",length(x_lc_fft))

  y_pos_fft <- stats::fft(y_pos)/sqrt(length(y_pos))
  y_pos_label <- rep("y_full",length(y_pos_fft))
  y_transient_fft <- stats::fft(y_transient)/sqrt(length(y_transient))
  y_transient_label <- rep("y_transient",length(y_transient_fft))
  y_lc_fft <- stats::fft(y_lc)/sqrt(length(y_lc))
  y_lc_label <- rep("y_lc",length(y_lc_fft))
  
  spectra <- c(x_pos_fft,x_transient_fft,x_lc_fft,y_pos_fft,y_transient_fft,y_lc_fft)
  labels <- c(x_pos_label,x_transient_label,x_lc_label,y_pos_label,y_transient_label,y_lc_label)
  spectral_tibble <- tibble(spectra=spectra,splice=labels)
  spectral_tibble <- spectral_tibble %>% group_by(splice) %>% mutate(index = row_number() - 1)
  
  return(spectral_tibble)
}

get_trajectory_plot <- function(full_samples,transient_cutoff,gradient=FALSE){

  x_var <- ifelse(gradient,"f_x","x")
  y_var <- ifelse(gradient,"f_y","y")
  
  # add indicator if in transient cutoff range
  full_samples_toplot <- full_samples %>% mutate(transient = ifelse(row_number() <= transient_cutoff,"Transient","Limit Cycle"))
  
  full_traj <- ggplot2::ggplot(full_samples_toplot,aes(x=.data[[x_var]],y=.data[[y_var]])) +
    geom_point(color = "#FDB515", alpha = 0.7) +
    facet_wrap(~transient) +
    labs(x = x_var, y = y_var, title = paste0("Trajectory with Cutoff at ",transient_cutoff, " Samples"))
  
  return(full_traj)
}

get_axis_plots <- function(full_samples,transient_cutoff,gradient=FALSE){
  
  x_var <- ifelse(gradient,"f_x","x")
  y_var <- ifelse(gradient,"f_y","y")
  
  full_samples_toplot <- full_samples %>% 
    mutate(index = row_number()) %>% 
    select(!!x_var,!!y_var,index) %>%
    pivot_longer(-index,names_to="axis",values_to="position")
  
  axis_plots <- ggplot2::ggplot(full_samples_toplot,aes(x=index,y=position)) +
    geom_point(alpha=0.7) +
    geom_vline(xintercept = transient_cutoff, linetype = "longdash")+
    facet_wrap(~axis) +
    labs(x = "Sample Index", y = "Position", title = "Projection onto Each Axis")
    
  return(axis_plots)
}

spectro_superheat <- function (wave, f, tlab = "Time (s)", flab = "Frequency (kHz)", 
                                 alab = "Amplitude\n(dB)\n", ...){
  # modified version of seewave::ggspectro
  
  spectrogram <- seewave::spectro(wave, f = f, plot = FALSE, ...)
  frequency <- rep(spectrogram$freq, times = ncol(spectrogram$amp))
  time <- rep(spectrogram$time, each = nrow(spectrogram$amp))
  amplitude <- as.vector(spectrogram$amp)
  
  df <- data.frame(time, frequency, amplitude)
  stft_df <- df %>% pivot_wider(names_from=time,values_from = amplitude)
  stft_df <- stft_df %>% column_to_rownames(var="frequency")
  stft_matrix <- as.matrix(stft_df)
  
  stft_plot <- superheat::superheat(stft_matrix,
                row.title = "Frequency",
                column.title = "Time (s)",
                left.label.size = 0.2,
                bottom.label.text.angle = 90)
  
  return(stft_plot)
}

analyze_lc_cutoff <- function(full_samples,sample_rate,transient_cutoff,gradient=FALSE){
  spectral_tibble <- get_coordinate_spectrum(full_samples,transient_cutoff,gradient=gradient)
  
  spectra_lollipop <- ggplot2::ggplot(spectral_tibble,aes(x=index,y=Mod(spectra))) + geom_point(alpha=0.3,color="blue") + 
    geom_segment(aes(x=index, xend=index, y=0, yend=Mod(spectra)),alpha=0.3) + facet_wrap(~splice,scales="free_x") +
    labs(x = "Freq. Index", y = "Normalized Amplitude", title = "FFT Coefficient Amplitude")
  phase_lollipop <- ggplot2::ggplot(spectral_tibble,aes(x=index,y=Arg(spectra))) + geom_point(alpha=0.3,color="blue") + 
    geom_segment(aes(x=index, xend=index, y=0, yend=Arg(spectra)),alpha=0.3) + facet_wrap(~splice,scales="free_x") +
    labs(x = "Freq. Index", y = "Radians", title = "FFT Coefficient Phase")
  
  trajectories <- get_trajectory_plot(full_samples,transient_cutoff,gradient=gradient)
  projections <- get_axis_plots(full_samples,transient_cutoff,gradient=gradient)
  
  x_spectrogram <- spectro_superheat(ifelse(rep(gradient,nrow(full_samples)),full_samples[,"f_x"],full_samples[,"x"]),
                                      f=sample_rate,wl=20,ovl=75,print.plot=T)
  y_spectrogram <- spectro_superheat(ifelse(rep(gradient,nrow(full_samples)),full_samples[,"f_y"],full_samples[,"y"]),
                                      f=sample_rate,wl=20,ovl=75,print.plot=T)
  
  combined_plot <- ggpubr::ggarrange(
    trajectories,
    projections,
    ggarrange(spectra_lollipop, phase_lollipop, ncol = 2), 
    ggarrange(x_spectrogram$plot,y_spectrogram$plot,ncol=2,labels=c("x","y")),
    nrow = 4
  )
  
  return(combined_plot)
}

derivative_comparison <- function(full_samples,sample_rate,transient_cutoff){
  f_plot <- analyze_lc_cutoff(full_samples,sample_rate,transient_cutoff,F)
  df_plot <- analyze_lc_cutoff(full_samples,sample_rate,transient_cutoff,T)
  
  output_plot <- combined_plot <- ggpubr::ggarrange(
    f_plot,
    df_plot,
    ncol = 2,
    labels = c("f","df")
  ) 
  
  output_plot
}

########################x
# Van der Pol No Noise #
########################

# Generate Van der Pol Data with No Noise 
sample_deltaT <- 0.1 # inverse of sampling rate
vdp_data_noiseless <- list(name = "VDP", system = "van_der_pol", params = list(mu=3),
                              n = 1000, lc_tail_n = 700, sample_density = sample_deltaT, var_x = 0.05, var_y = 0.05,
                              x_grid_size = 36, y_grid_size = 36, extrapolation_size = 0.5, 
                              smoother = "bspline", data_seed = 2, noise_seed = 2)
experiment_list <- list(vdp_data_noiseless)
experiment_data <- generate_data_object_model(experiment_list)

# Visualize true and spline-estimated solution paths with transients
truth <- list(method = "truth",  params = list())
spline <- list(method = "spline",  params = list(lambda = 1e-2, norder = 6, nbasis = 20, side_info = list()))

experiment_estimators <- list(truth, spline)
experiment_results <- evaluate_gradient_methods(experiment_data, experiment_estimators)

plot_delta <- list(type = "field_delta_paths", experiments = list(c(data = 1, estimator = 1, ref = 1),c(data = 1, estimator = 2, ref = 1)))
visualize_results(experiment_results, list(plot_delta))

#spectral_tibble <- get_coordinate_spectrum(experiment_data[[1]]$limit_cycle_samples)
derivative_comparison(experiment_data[[1]]$limit_cycle_samples,1/sample_deltaT,600)
analyze_lc_cutoff(experiment_data[[1]]$limit_cycle_samples,55)
analyze_lc_cutoff(experiment_data[[1]]$limit_cycle_samples,70)
analyze_lc_cutoff(experiment_data[[1]]$limit_cycle_samples,100)
analyze_lc_cutoff(experiment_data[[1]]$limit_cycle_samples,200)