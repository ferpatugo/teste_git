## Libraries
#source('/Users/arthurvaz/Desktop/CEFETRJ - Mestrado/Anomalia/Anomalia/carregar_bases.R')


###LOAD AND INSTALL
{
  # Library to load and install if ins't already installed
  #install.packages('pacman')
  library(pacman)
  
  vetor_pacotes <- c("ggplot2", "plotly", "dplyr", "kableExtra", "gridExtra",
                     "TSPred", "dbscan", "fpc", "factoextra", "fitdistrplus",
                     "logspline", "caret","RColorBrewer","R.filesets","glue")
  
  pacman::p_load(char = vetor_pacotes)
  
  library(daltoolbox)
  library(harbinger)
}

# Base directory for the datasets
base_dir <- "/Users/arthurvaz/Desktop/CEFETRJ - Mestrado/Anomalia/Anomalia/Dataset"


generate_time_series <- function(n, anomalies_qty, variance_size, anomaly_intensity,noise_sd = 0.5,seasonality = c(1,-1),trend_increment = 0.1 ,stationary = TRUE, window_anomaly = TRUE) {
  # Generate the time index
  time_index <- 1:n
  
  if (stationary) {
    # Generate a random time series with constant variance
    time_series <- rnorm(n, mean = 0, sd = variance_size)
  } else {
    set.seed(42)
    time_series <- create_nonstationary_ts(n = n, 
                                           trend = list(initial = 0, increment = trend_increment), 
                                           seasonality = list(pattern = seasonality), 
                                           noise_sd = noise_sd)
  }
  anomaly_idx_list <- c()
  # Add anomalies
  for (i in 1:anomalies_qty) {
    # Randomly select the position of the anomaly
    anomaly_position <- sample(1:n, 1)
    # Randomly select the size of the anomaly
    if (window_anomaly) {
      anomaly_size <- sample(1:10, 1)
      anomaly_indices <- anomaly_position:(anomaly_position + anomaly_size - 1)
      # Store the indices of the anomaly
      anomaly_idx_list <- c(anomaly_idx_list,anomaly_indices)
    } else {
      anomaly_indices <- anomaly_position
      # Store the indices of the anomaly
      anomaly_idx_list <- c(anomaly_idx_list,anomaly_indices)
    }
    # Randomly select the direction of the anomaly
    anomaly_direction <- sample(c(-1, 1), 1)
    # Add the anomaly to the time series
    time_series[anomaly_indices] <- time_series[anomaly_indices] + anomaly_direction * anomaly_intensity
  }
  
  # Transform the time series and anomalies into a data frame
  time_series_df <- data.frame(time = time_index, value = time_series)
  time_series_df$event <- FALSE
  if (window_anomaly) {
    time_series_df[anomaly_idx_list, 'event'] <- TRUE
  } else {
    time_series_df[anomaly_idx_list, 'event'] <- TRUE
  }
  
  return(time_series_df)
}

create_nonstationary_ts <- function(n, trend = NULL, seasonality = NULL, noise_sd = 1) {
  # Generate trend component
  if (!is.null(trend)) {
    trend_values <- seq(trend$initial, length.out = n, by = trend$increment)
  } else {
    trend_values <- rep(0, n)
  }
  
  # Generate seasonality component
  if (!is.null(seasonality)) {
    seasonality_values <- rep(seasonality$pattern, length.out = n)
  } else {
    seasonality_values <- rep(0, n)
  }
  
  # Generate noise
  noise <- rnorm(n, sd = noise_sd)
  
  # Combine components
  ts_values <- trend_values + seasonality_values + noise
  
  # Create time series object
  ts(ts_values)
}

###################### CRIACAO DAS FUNCOES ###########################

kmeans_func <- function(time_series, window, n_cluster = 2) {
  # disable the warnings
  options(warn = -1)
  model <- hanct_kmeans(seq = window, n_cluster)
  fitted_model <- fit(model, time_series)
  detection <- detect(fitted_model, time_series)
  seq_window <- detection %>% filter(event == TRUE)  
  if (nrow(seq_window) == 0) {
    print('No anomalies detected')
    return(time_series)
  }
  seq_window$idx_final <- seq_window$idx+seq_window$seqlen-1
  #print(seq_window)
  seq_window[seq_window$idx_final > nrow(time_series),'idx_final'] <- nrow(time_series)
  time_series$predicted_event <- FALSE
  time_series$seq_start <- 0
  time_series$seq_len <- 0
  for (i in 1:nrow(seq_window)) {
    time_series$predicted_event[seq_window$idx[i]:seq_window$idx_final[i]] <- TRUE
    time_series$seq_start[seq_window$idx[i]] <- 1
    time_series$seq_len[seq_window$idx[i]] <- window
  }
  return(time_series)
}

generate_plot <- function(time_series,predicted_event,model) {
  # check if the predicted_event column exists
  if (!predicted_event %in% colnames(time_series)) {
    print("The column does not exist in the data frame")
    return()
  }
  
  # Calculate the current range of y-values
  y_range <- range(time_series$value)
  # Calculate the new limits for the y-axis
  new_y_lim <- c(y_range[1] - 0.5 * diff(y_range), y_range[2] + 0.5 * diff(y_range))
  # Plot the time series with the detected anomalies
  plot(time_series$value, type = "l", col = "blue", lwd = 2, xlab = "Time", ylab = "Value", ylim = new_y_lim)
  # Add a scatter plot with the detected anomalies
  points(which(time_series[,predicted_event]), time_series$value[which(time_series[,predicted_event])], col = "red", pch = 25, cex = 0.8)
  points(which(time_series$event), time_series$value[which(time_series$event)], col = "green", pch = 19, cex = 0.5)
  # Add a legend to the plot
  legend("bottomright", legend = c("Actual Anomalies", "Detected Anomalies"), col = c("green", "red"), pch = c(19, 25), cex = 0.5)
  # Add a title to the plot
  title(glue("Time Series with Anomalies - {model}"), cex.main = 0.7)
}

dbscan_func <- function(time_series, window,minPts) {
  
  # Extract time series data with given sliding window size
  data_sw <- ts_data(time_series[, 'value'], window)
  
  # Calculate kNN distances
  kNNdist <- kNNdist(data_sw, k = minPts)
  
  # Sort the distances in ascending order
  sorted_distances <- sort(kNNdist)
  
  # Calculate the percentage change of the sorted distances
  percentage_change <- c(0, diff(sorted_distances) / sorted_distances[-length(sorted_distances)] * 100)
  
  # Find the index of the maximum percentage change
  max_index <- which.max(percentage_change)
  eps <- sorted_distances[max_index]
  
  # Convert data into matrix for dbscan
  datamatrix <- matrix(data_sw, ncol = window)
  
  # Perform DBSCAN clustering
  dbscanResult <- dbscan::dbscan(datamatrix, eps = eps, minPts = minPts)
  
  # Create a dataframe from datamatrix with clusters
  datamatrix_df <- data.frame(datamatrix)
  datamatrix_df$clusters <- dbscanResult$cluster
  
  idx <- datamatrix_df[datamatrix_df$clusters == 0,]%>%rownames()%>%as.numeric()
  
  if (length(idx) == 0) {
    print('No anomalies detected')
    return(time_series)
  }
  
  idx_final <- idx+window-1
  idx_final[idx_final > nrow(time_series)] <- nrow(time_series)
  # print(idx)
  # print(idx_final)
  time_series$predicted_event2 <- FALSE
  time_series$seq_start2 <- 0
  time_series$seq_len2 <- 0
  for (i in 1:length(idx)) {
    time_series$predicted_event2[idx[i]:idx_final[i]] <- TRUE
    time_series$seq_start2[idx[i]] <- 1
    time_series$seq_len2[idx[i]] <- window
  }
  
  return(time_series)
}

autoencoder_func <- function(time_series, window, encode = 1) {
  model <- han_autoencoder(window,encode)
  fitted_model <- fit(model, time_series[,'value'])
  detection <- detect(fitted_model,  time_series[,'value'])
  seq_window <- detection %>% filter(event == TRUE)  
  if (nrow(seq_window) == 0) {
    print('No anomalies detected')
    return(time_series)
  }
  #print(seq_window)
  seq_window$idx_final <- seq_window$idx+window-1
  seq_window[seq_window$idx_final > nrow(time_series),'idx_final'] <- nrow(time_series)
  #print(seq_window)
  time_series$predicted_event3 <- FALSE
  time_series$seq_start3 <- 0
  time_series$seq_len3 <- 0
  for (i in 1:nrow(seq_window)) {
    time_series$predicted_event3[seq_window$idx[i]:seq_window$idx_final[i]] <- TRUE
    time_series$seq_start3[seq_window$idx[i]] <- 1
    time_series$seq_len3[seq_window$idx[i]] <- window
  }
  return(time_series)
}

prob_dist <-  function(data, limiar) {
  distribution_list <- c("gamma", "lnorm", "exp", "gamma", "unif", "logis", "norm")
  fits <- list()
  
  for (dist in distribution_list) {
    tryCatch({
      fits[[dist]] <- fitdist(data, dist)
    }, error = function(e) {
      cat("Error fitting distribution", dist, ": ", conditionMessage(e), "\n")
    })
  }
  
  gof_stats <- lapply(fits, gofstat)
  best_fit_idx <- which.min(sapply(gof_stats, function(x) x$bic))
  best_fit <- names(best_fit_idx)
  print(best_fit)
  if (is.null(best_fit)) {
    print('No anomalies detected with the given limiar')
    return(data)
  } else {
    best_fit <- strsplit(best_fit, "-")[[1]][3]
  }
  
  
  params <- as.list(fits[[best_fit_idx]]$estimate)
  is_anomaly <- switch(
    best_fit,
    lognormal = dlnorm(data, meanlog = params$meanlog, sdlog = params$sdlog) < limiar,
    exp = dexp(data, rate = params$rate) < limiar,
    gamma = dgamma(data, shape = params$shape, rate = params$rate) < limiar,
    unif = dunif(data, min = params$min, max = params$max) < limiar,
    logis = dlogis(data, location = params$location, scale = params$scale) < limiar,
    dnorm(data, mean = params$mean, sd = params$sd) < limiar
  )
  
  anomalias <- data[is_anomaly]
  
  # Check if there are anomalies
  if (length(anomalias) == 0) {
    print('No anomalies detected using the given limiar')
    return()
  }
  
  return(anomalias)
}

anomaly_window_normalization <- function(time_series, method = 'dbscan', limiar = 0.2) {
  # seleciona as colunas de acordo com o modelo escolhido
  switch(method,
         'kmeans' = {
           predict_col <- 'predicted_event'
           seq_start_col <- 'seq_start'
           seq_len_col <- 'seq_len'
         },
         'dbscan' = {
           predict_col <- 'predicted_event2'
           seq_start_col <- 'seq_start2'
           seq_len_col <- 'seq_len2'
         },
         'autoencoder' = {
           predict_col <- 'predicted_event3'
           seq_start_col <- 'seq_start3'
           seq_len_col <- 'seq_len3'
         }
  )
  # filtra o data frame com as colunas selecionadas
  data <- time_series[time_series[, predict_col], c(predict_col, 'value', seq_start_col, seq_len_col)]
  # a coluna seq_start é igual a 1 (boolean) e para ele pegar a real posicao precisa olhar o índice do data frame
  data$seq_inicio <- rownames(data) %>% as.numeric() * data[, seq_start_col]
  data$seq_fim <- data$seq_inicio + data[, seq_len_col] - 1
  data$seq_fim <- data$seq_fim * data[, seq_start_col]
  
  z_scores_idx <- data[data[, seq_start_col] == 1, c('seq_inicio', 'seq_fim')]
  z_scores_list <- lapply(seq_along(data$seq_inicio), function(i) {
    seq_ <- data$value[rownames(data) %in% seq(data$seq_inicio[i], data$seq_fim[i])]
    mean_val <- mean(seq_)
    sd_val <- sd(seq_)
    z_scores <- scale(seq_, center = mean_val, scale = sd_val)
    z_scores[, 1]
  })
  
  y_ <- unlist(z_scores_list)
  # drop the NA values
  y <- y_[!is.na(y_)]
  anomalias <- prob_dist(y, limiar)
  
  # check if the anomalias is NULL, then print the message "não existe anomalias"
  if (is.null(anomalias)) {
    print('Não existe anomalias')
    return()
  }
  
  data2 <- data.frame(seq_anomaly = y_, dist_anomaly = y_ %in% anomalias)
  
  idx_anomalies <- lapply(seq(length(z_scores_idx$seq_inicio)), function(i) {
    seq(z_scores_idx$seq_inicio[i], z_scores_idx$seq_fim[i])
  })
  
  # filtrar para as posicoes existentes da serie levando em conta o seu tamanho
  anomalia_list <- unlist(idx_anomalies)
  anomalia_list <- anomalia_list[anomalia_list<= nrow(time_series)]
  
  # retorna a sequencia original
  data2$real_Seq <- anomalia_list
  
  # check the real_Seq if the dist_anomaly is False THEN real_Seq is Null
  data2$real_Seq <- ifelse(data2$dist_anomaly, data2$real_Seq, NA)
  
  
  return(data2)
}

window_abroad_normalization <- function(method, time_series, threshold_window, more = 1) {
  
  switch(method,
         'kmeans' = {
           predict_col <- 'predicted_event'
           seq_start_col <- 'seq_start'
           seq_len_col <- 'seq_len'
         },
         'dbscan' = {
           predict_col <- 'predicted_event2'
           seq_start_col <- 'seq_start2'
           seq_len_col <- 'seq_len2'
         },
         'autoencoder' = {
           predict_col <- 'predicted_event3'
           seq_start_col <- 'seq_start3'
           seq_len_col <- 'seq_len3'
         }
  )
  # filtra as colunas de acordo com o modelo escolhido
  data_total <- time_series[time_series[, predict_col], c(predict_col,'time','value', seq_start_col, seq_len_col)]
  # filtra a primeira posicao da sequencia anomala
  data_total <- data_total[data_total[,seq_start_col]==1,]
  # adiciona as colunas de aumento de janela
  data_total <- cbind(data_total, generate_length_window(data_total,method,threshold_window,predict_col,seq_start_col,seq_len_col, more))
  return(data_total)
}

generate_length_window <- function(data_total,method,threshold_window,predict_col,seq_start_col,seq_len_col, more = 1) {
  seq_start_left <- list()
  seq_end_left <- list()
  seq_len_left <- list()
  limit_left <- list()
  seq_start_right <- list()
  seq_end_right <- list()
  seq_len_right <- list()
  limit_right <- list()
  
  #reset index
  rownames(data_total) <- 1:nrow(data_total)
  
  # loop interaction
  for (i in 1:nrow(data_total)) {
    # print(i)
    data <- data_total[i,]
    # print(data)
    # if more equal to 1 ceiling will be on the right, if 0 floor will be on the right
    # a coluna time é a sequencia original
    if (more == 1) {
      if(floor((threshold_window - data[,seq_len_col])/2)+data[,'time']+data[,seq_len_col] > nrow(time_series)) {
        seq_start_right[i] <- (data[,'time']+data[,seq_len_col])+1
        seq_end_right[i] <- nrow(time_series)
        seq_len_right[i] <- seq_end_right[[i]][[1]] - seq_start_right[[i]][[1]]
        limit_right[i] <- TRUE
      } else {
        seq_start_right[i] <- (data[,'time']+data[,seq_len_col])+1
        seq_end_right[i] <- floor((threshold_window - data[,seq_len_col])/2)+data[,'time']+data[,seq_len_col]+1
        seq_len_right[i] <- seq_end_right[[i]][[1]] - seq_start_right[[i]][[1]]
        limit_right[i] <- FALSE
      }
      
      if(data[,'time']-ceiling((threshold_window - data[,seq_len_col])/2) < 0) {
        seq_start_left[i] <- 1
        seq_end_left[i] <- data[,'time']-seq_start_left[[i]][[1]]
        seq_len_left[i] <- seq_end_left[[i]][[1]] - seq_start_left[[i]][[1]]
        limit_left[i] <- TRUE
      } else {
        seq_start_left[i] <- data[,'time']-ceiling((threshold_window - data[,seq_len_col])/2)-1
        seq_end_left[i] <- seq_start_left[[i]][[1]]+ceiling((threshold_window - data[,seq_len_col])/2)
        seq_len_left[i] <- seq_end_left[[i]][[1]] - seq_start_left[[i]][[1]]
        limit_left[i] <- FALSE
      }
      
      if (limit_right[[i]][[1]] == TRUE) {
        residual_size <- floor((threshold_window - data[,seq_len_col])/2) - seq_len_right[[i]][[1]]
        seq_start_left[i] <- seq_start_left[[i]][[1]] - residual_size
        seq_len_left[i] <- seq_len_left[[i]][[1]] + residual_size
      } else if (limit_right[[i]][[1]] == TRUE) {
        residual_size <- ceiling((threshold_window - data[,seq_len_col])/2) - seq_len_left[[i]][[1]]
        seq_end_right[i] <- seq_end_right[[i]][[1]] + residual_size
        seq_len_right[i] <- seq_len_right[[i]][[1]] + residual_size
      }
      
    }
    else {
      if(ceiling((threshold_window - data[,seq_len_col])/2)+data[,'time']+data[,seq_len_col] > nrow(time_series)) {
        seq_start_right[i] <- (data[,'time']+data[,seq_len_col])+1
        seq_end_right[i] <- nrow(time_series)
        seq_len_right[i] <- seq_end_right[[i]][[1]] - seq_start_right[[i]][[1]]
        limit_right[i] <- TRUE
      } else {
        seq_start_right[i] <- (data[,'time']+data[,seq_len_col])
        seq_end_right[i] <- ceiling((threshold_window - data[,seq_len_col])/2)+data[,'time']+data[,seq_len_col]
        seq_len_right[i] <- seq_end_right[[i]][[1]] - seq_start_right[[i]][[1]]
        limit_right[i] <- FALSE
      }
      
      if(data[,'time']-floor((threshold_window - data[,seq_len_col])/2) < 0) {
        seq_start_left[i] <- 1
        seq_end_left[i] <- data[,'time']-seq_start_left[[i]][[1]]
        seq_len_left[i] <- seq_end_left[[i]][[1]] - seq_start_left[[i]][[1]]
        limit_left[i] <- TRUE
      } else {
        seq_start_left[i] <- data[,'time']-floor((threshold_window - data[,seq_len_col])/2)-1
        seq_end_left[i] <- seq_start_left[[i]][[1]]+floor((threshold_window - data[,seq_len_col])/2)
        seq_len_left[i] <- seq_end_left[[i]][[1]] - seq_start_left[[i]][[1]]
        limit_left[i] <- FALSE
      }
      
      if (limit_right[[i]][[1]] == TRUE) {
        residual_size <- ceiling((threshold_window - data[,seq_len_col])/2) - seq_len_right[[i]][[1]]
        seq_start_left[i] <- seq_start_left[[i]][[1]] - residual_size
        seq_len_left[i] <- seq_len_left[[i]][[1]] + residual_size
      } else if (limit_right[[i]][[1]] == TRUE) {
        residual_size <- floor((threshold_window - data[,seq_len_col])/2) - seq_len_left[[i]][[1]]
        seq_end_right[i] <- seq_end_right[[i]][[1]] + residual_size
        seq_len_right[i] <- seq_len_right[[i]][[1]] + residual_size
      }
    }
  }
  switch(method,
         'kmeans' = {
           result <- data.frame('seq_start_left' = unlist(seq_start_left), 'seq_end_left' = unlist(seq_end_left), 
                                'seq_len_left' = unlist(seq_len_left), 'seq_start_right' = unlist(seq_start_right),
                                'seq_end_right' = unlist(seq_end_right), 'seq_len_right' = unlist(seq_len_right))
         },
         'dbscan' = {
           result <- data.frame('seq_start_left2' = unlist(seq_start_left), 'seq_end_left2' = unlist(seq_end_left), 
                                'seq_len_left2' = unlist(seq_len_left), 'seq_start_right2' = unlist(seq_start_right),
                                'seq_end_right2' = unlist(seq_end_right), 'seq_len_right2' = unlist(seq_len_right))
         },
         'autoencoder' = {
           result <- data.frame('seq_start_left3' = unlist(seq_start_left), 'seq_end_left3' = unlist(seq_end_left), 
                                'seq_len_left3' = unlist(seq_len_left), 'seq_start_right3' = unlist(seq_start_right),
                                'seq_end_right3' = unlist(seq_end_right), 'seq_len_right3' = unlist(seq_len_right))
         }
  )
  
  
  
  return(result)
}

wide_window_anomaly <- function(method, time_series, threshold_window, more = 1, limiar) {
  data <- window_abroad_normalization(method, time_series, threshold_window, more = 1)
  anomaly <- list()
  
  switch(method,
         'kmeans' = {
           seq_start_left <- data$seq_start_left
           seq_end_left <- data$seq_end_left
           seq_start_right <- data$seq_start_right
           seq_end_right <- data$seq_end_right
         },
         'dbscan' = {
           seq_start_left <- data$seq_start_left2
           seq_end_left <- data$seq_end_left2
           seq_start_right <- data$seq_start_right2
           seq_end_right <- data$seq_end_right2
         },
         'autoencoder' = {
           seq_start_left <- data$seq_start_left3
           seq_end_left <- data$seq_end_left3
           seq_start_right <- data$seq_start_right3
           seq_end_right <- data$seq_end_right3
         }
  )
  
  
  
  for (i in 1:nrow(data)) {
    data_ <- time_series[time_series$time%in%c(seq_start_left[i]:seq_end_right[i]),'value']
    anomaly[[i]] <- prob_dist(data_,limiar)
  }
  
  return(anomaly)
}


###################### ANOMALY WINDOW STACIONAY ######################

set.seed(12345)

time_series <- generate_time_series(n = 100, 
                                    anomalies_qty = 1, 
                                    variance_size = 0.4, 
                                    stationary = TRUE,
                                    anomaly_intensity = 4)

write.csv2(time_series,"serie1_estacionaria.csv",row.names = F,sep = ";")
time_series = read.csv2("serie1_estacionaria.csv",sep=";",h=T)

# plot
plot(time_series$value, type = "l", col = "blue", lwd = 2, xlab = "Time", ylab = "Value",main = "Série Estacionária sintética com uma observação anômala")

# Example KMEANS:

time_series_anomaly <- kmeans_func(time_series, 10,1)

generate_plot(time_series_anomaly,'predicted_event','kmeans')

# Example DBSCAN:

time_series_anomaly <- dbscan_func(time_series_anomaly, 10,2)

generate_plot(time_series_anomaly,'predicted_event2','dbscan')

# Example Autoencoder:

time_series_anomaly <- autoencoder_func(time_series_anomaly, 10,1)

generate_plot(time_series_anomaly,'predicted_event3','autoencoder')

###################### PROBABILISTIC ANOMALIES ######################

#----------- Example KMEANS:

#### METODO 1
data <- time_series_anomaly[time_series_anomaly$predicted_event,'value']
anomalias <- prob_dist(data, 0.25)
time_series_anomaly$predicted_event_dist <- FALSE
time_series_anomaly[time_series_anomaly$value%in%anomalias,'predicted_event_dist'] <- TRUE
generate_plot(time_series_anomaly,'predicted_event_dist','kmeans com filtro dist unica')

#### METODO 2
data <- anomaly_window_normalization(time_series_anomaly, 'kmeans', 0.25)
time_series_anomaly$predicted_event_dist <- rownames(time_series_anomaly) %in% data$real_Seq
generate_plot(time_series_anomaly,'predicted_event_dist','kmeans com filtro dist unica norm')

#### METODO 3
anomalias <- wide_window_anomaly('kmeans', time_series_anomaly, threshold_window = 10, more = 1, limiar = 0.25)
time_series_anomaly$predicted_event_dist <- FALSE
time_series_anomaly[time_series_anomaly$value %in% unlist(anomalias),'predicted_event_dist'] <- TRUE
generate_plot(time_series_anomaly,'predicted_event_dist','kmeans com janela extendida dist ind')

#----------- Example DBSCAN:

#### METODO 1
data <- time_series_anomaly[time_series_anomaly$predicted_event2,'value']
anomalias <- prob_dist(data, 0.25)
time_series_anomaly$predicted_event2_dist <- FALSE
time_series_anomaly[time_series_anomaly$value%in%anomalias,'predicted_event2_dist'] <- TRUE
generate_plot(time_series_anomaly,'predicted_event2_dist','dbscan com filtro dist unica')

#### METODO 2
data <- anomaly_window_normalization(time_series_anomaly, 'dbscan', 0.25)
time_series_anomaly$predicted_event2_dist <- rownames(time_series_anomaly) %in% data$real_Seq
generate_plot(time_series_anomaly,'predicted_event2_dist','dbscan com filtro dist unica norm')

#### METODO 3
anomalias <- wide_window_anomaly('dbscan', time_series_anomaly, threshold_window = 10, more = 1, limiar = 0.25)
time_series_anomaly$predicted_event2_dist <- FALSE
time_series_anomaly[time_series_anomaly$value %in% unlist(anomalias),'predicted_event2_dist'] <- TRUE
generate_plot(time_series_anomaly,'predicted_event2_dist','dbscan com janela extendida dist ind')

#----------- Example AUTOENCODER:

#### METODO 1
data <- time_series_anomaly[time_series_anomaly$predicted_event3,'value']
anomalias <- prob_dist(data, 0.25)
time_series_anomaly$predicted_event3_dist <- FALSE
time_series_anomaly[time_series_anomaly$value%in%anomalias,'predicted_event3_dist'] <- TRUE
generate_plot(time_series_anomaly,'predicted_event3_dist','autoencoder com filtro dist unica')

#### METODO 2
data <- anomaly_window_normalization(time_series_anomaly, 'autoencoder', 0.25)
time_series_anomaly$predicted_event3_dist <- rownames(time_series_anomaly) %in% data$real_Seq
generate_plot(time_series_anomaly,'predicted_event3_dist','autoencoder com filtro dist unica norm')

#### METODO 3
anomalias <- wide_window_anomaly('autoencoder', time_series_anomaly, threshold_window = 10, more = 1, limiar = 0.25)
time_series_anomaly$predicted_event3_dist <- FALSE
time_series_anomaly[time_series_anomaly$value %in% unlist(anomalias),'predicted_event3_dist'] <- TRUE
generate_plot(time_series_anomaly,'predicted_event3_dist','autoencoder com janela extendida dist ind')


###################### ANOMALY SEQ STACIONARY ######################

set.seed(56789)

time_series2 <- generate_time_series(n = 100, 
                                    anomalies_qty = 4, 
                                    variance_size = 0.4, 
                                    stationary = TRUE,
                                    anomaly_intensity = 4,
                                    window_anomaly = FALSE)

write.csv2(time_series2,"serie2_estacionaria_seqanomaly.csv",row.names = F,sep = ";")

# plot
plot(time_series2$value, type = "l", col = "blue", lwd = 2, xlab = "Time", ylab = "Value", main = "Série Estacionária com 4 observações anômalas")

# Example KMEANS:

time_series_anomaly2 <- kmeans_func(time_series2, 10)

generate_plot(time_series_anomaly2,'predicted_event','kmeans')

# Example DBSCAN:

time_series_anomaly2 <- dbscan_func(time_series_anomaly2, 10,2)

generate_plot(time_series_anomaly2,'predicted_event2','dbscan')

# Example AUTOENCODER:

time_series_anomaly2 <- autoencoder_func(time_series_anomaly2, 10, 1)

generate_plot(time_series_anomaly2,'predicted_event3','autoencoder')

###################### PROBABILISTIC ANOMALIES ######################

#----------- Example KMEANS:
#nao precisa rodar nenhum método, pois o k-means não detectou nenhuma anomalia

#### METODO 1 
data <- time_series[time_series$predicted_event,'value']
anomalias <- prob_dist(data, 0.25)
time_series$predicted_event_dist <- FALSE
time_series[time_series$value%in%anomalias,'predicted_event_dist'] <- TRUE
generate_plot(time_series,'predicted_event_dist','dbscan com filtro dist unica')

#### METODO 2
data <- anomaly_window_normalization(time_series, 'kmeans', 0.25)
time_series$predicted_event_dist <- rownames(time_series) %in% data$real_Seq
generate_plot(time_series,'predicted_event_dist','kmeans com filtro dist unica norm')

#### METODO 3
anomalias <- wide_window_anomaly('kmeans', time_series, threshold_window = 10, more = 1, limiar = 0.25)
time_series$predicted_event_dist <- FALSE
time_series[time_series$value %in% unlist(anomalias),'predicted_event_dist'] <- TRUE
generate_plot(time_series,'predicted_event_dist','kmeans com janela extendida dist ind')

#----------- Example DBSCAN:

#### METODO 1
data <- time_series2[time_series_anomaly2$predicted_event2,'value']
anomalias <- prob_dist(data, 0.25)
time_series_anomaly2$predicted_event2_dist <- FALSE
time_series_anomaly2[time_series_anomaly2$value%in%anomalias,'predicted_event2_dist'] <- TRUE
generate_plot(time_series_anomaly2,'predicted_event2_dist','dbscan com filtro dist unica')

#### METODO 2
data <- anomaly_window_normalization(time_series_anomaly2, 'dbscan', 0.25)
time_series_anomaly2$predicted_event2_dist <- rownames(time_series_anomaly2) %in% data$real_Seq
generate_plot(time_series_anomaly2,'predicted_event2_dist','dbscan com filtro dist unica norm')

#### METODO 3
anomalias <- wide_window_anomaly('dbscan', time_series_anomaly2, threshold_window = 10, more = 1, limiar = 0.25)
time_series_anomaly2$predicted_event2_dist <- FALSE
time_series_anomaly2[time_series_anomaly2$value %in% unlist(anomalias),'predicted_event2_dist'] <- TRUE
generate_plot(time_series_anomaly2,'predicted_event2_dist','dbscan com janela extendida')

#----------- Example AUTOENCODER:

#### METODO 1
data <- time_series_anomaly2[time_series_anomaly2$predicted_event3,'value']
anomalias <- prob_dist(data, 0.25)
time_series_anomaly2$predicted_event3_dist <- FALSE
time_series_anomaly2[time_series_anomaly2$value%in%anomalias,'predicted_event3_dist'] <- TRUE
generate_plot(time_series_anomaly2,'predicted_event3_dist','autoencoder com filtro dist unica')

#### METODO 2
data <- anomaly_window_normalization(time_series_anomaly2, 'autoencoder', 0.25)
time_series_anomaly2$predicted_event3_dist <- rownames(time_series_anomaly2) %in% data$real_Seq
generate_plot(time_series_anomaly2,'predicted_event3_dist','autoencoder com filtro dist unica norm')

#### METODO 3
anomalias <- wide_window_anomaly('autoencoder', time_series_anomaly2, threshold_window = 10, more = 1, limiar = 0.25)
time_series_anomaly2$predicted_event3_dist <- FALSE
time_series_anomaly2[time_series_anomaly2$value %in% unlist(anomalias),'predicted_event3_dist'] <- TRUE
generate_plot(time_series_anomaly2,'predicted_event3_dist','autoencoder com janela extendida dist ind')


###################### ANOMALY POINT NO STACIONARY ######################

set.seed(91011)

time_series3 <- generate_time_series(n = 100, 
                                    anomalies_qty = 1, 
                                    variance_size = 0.4, 
                                    stationary = FALSE,
                                    anomaly_intensity = 4,
                                    window_anomaly = FALSE)

write.csv2(time_series3,"serie3_naoestacionaria.csv",row.names = F,sep = ";")

# plot
plot(time_series3$value, type = "l", col = "blue", lwd = 2, xlab = "Time", ylab = "Value")

# Example KMEANS:

time_series_anomaly3 <- kmeans_func(time_series3, 10,1)

generate_plot(time_series_anomaly3,'predicted_event','kmeans')

# Example DBSCAN:

time_series_anomaly3 <- dbscan_func(time_series_anomaly3, 10,2)

generate_plot(time_series_anomaly3,'predicted_event2','dbscan')

# Example AUTOENCODER:

time_series_anomaly3 <- autoencoder_func(time_series_anomaly3, 10, 1)

generate_plot(time_series_anomaly3,'predicted_event3','autoencoder')

###################### PROBABILISTIC ANOMALIES ######################

#----------- Example KMEANS:
#Não precisa rodar, pois o k-means não pegou nenhuma anomalia

#### METODO 1
data <- time_series_anomaly3[time_series$predicted_event,'value']
anomalias <- prob_dist(data, 0.2)
time_series$predicted_event_dist <- FALSE
time_series[time_series$value%in%anomalias,'predicted_event_dist'] <- TRUE
generate_plot(time_series,'predicted_event_dist','dbscan com filtro dist unica')

#### METODO 2
data <- anomaly_window_normalization(time_series, 'kmeans', 0.2)
time_series$predicted_event_dist <- rownames(time_series) %in% data$real_Seq
generate_plot(time_series,'predicted_event_dist','kmeans com filtro dist unica norm')

#### METODO 3
anomalias <- wide_window_anomaly('kmeans', time_series, threshold_window = 30, more = 1, limiar = 0.2)
time_series$predicted_event_dist <- FALSE
time_series[time_series$value %in% unlist(anomalias),'predicted_event_dist'] <- TRUE
generate_plot(time_series,'predicted_event_dist','kmeans com janela extendida dist ind')

#----------- Example DBSCAN:

#### METODO 1
data <- time_series_anomaly3[time_series_anomaly3$predicted_event2,'value']
anomalias <- prob_dist(data, 0.25)
time_series_anomaly3$predicted_event2_dist <- FALSE
time_series_anomaly3[time_series_anomaly3$value%in%anomalias,'predicted_event2_dist'] <- TRUE
generate_plot(time_series_anomaly3,'predicted_event2_dist','dbscan com filtro dist unica')

#### METODO 2 #não filtrou anomalias nesse limiar
data <- anomaly_window_normalization(time_series_anomaly3, 'dbscan', 0.25)
time_series_anomaly3$predicted_event2_dist <- rownames(time_series_anomaly3) %in% data$real_Seq
generate_plot(time_series_anomaly3,'predicted_event2_dist','dbscan com filtro dist unica norm')

#### METODO 3
anomalias <- wide_window_anomaly('dbscan', time_series_anomaly3, threshold_window = 10, more = 1, limiar = 0.25)
time_series_anomaly3$predicted_event2_dist <- FALSE
time_series_anomaly3[time_series_anomaly3$value %in% unlist(anomalias),'predicted_event2_dist'] <- TRUE
generate_plot(time_series_anomaly3,'predicted_event2_dist','dbscan com janela extendida')

#----------- Example AUTOENCODER:

#### METODO 1
data <- time_series_anomaly3[time_series_anomaly3$predicted_event3,'value']
anomalias <- prob_dist(data, 0.25)
time_series_anomaly3$predicted_event3_dist <- FALSE
time_series_anomaly3[time_series_anomaly3$value%in%anomalias,'predicted_event3_dist'] <- TRUE
generate_plot(time_series_anomaly3,'predicted_event3_dist','autoencoder com filtro dist unica')

#### METODO 2 #não filtrou anomalias nesse limiar
data <- anomaly_window_normalization(time_series_anomaly3, 'autoencoder', 0.25)
time_series_anomaly3$predicted_event3_dist <- rownames(time_series_anomaly3) %in% data$real_Seq
generate_plot(time_series_anomaly3,'predicted_event3_dist','autoencoder com filtro dist unica norm')

#### METODO 3
anomalias <- wide_window_anomaly('autoencoder', time_series_anomaly3, threshold_window = 10, more = 1, limiar = 0.25)
time_series_anomaly3$predicted_event3_dist <- FALSE
time_series_anomaly3[time_series_anomaly3$value %in% unlist(anomalias),'predicted_event3_dist'] <- TRUE
generate_plot(time_series_anomaly3,'predicted_event3_dist','autoencoder com janela extendida dist ind')
