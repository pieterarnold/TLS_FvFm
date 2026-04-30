# Calculate the temperature thresholds to bound data within
calculate_temp_threshold <- function(response, temperature, threshold,
                                     temp_min = 15, temp_max = 80,
                                     use_interpolation = TRUE,
                                     allow_extrapolation = TRUE) {
  
  valid_idx <- !is.na(response) & !is.na(temperature)
  response <- response[valid_idx]
  temperature <- temperature[valid_idx]
  
  if(length(response) < 3) return(NA)
  if(var(response, na.rm = TRUE) < 1e-6) return(NA)
  
  # Try interpolation first
  if(use_interpolation) {
    ord <- order(response)
    
    # Check if threshold is within data range
    if(threshold >= min(response) && threshold <= max(response)) {
      # Interpolate within range
      temp_at_threshold <- approx(x = response[ord], 
                                  y = temperature[ord], 
                                  xout = threshold,
                                  rule = 2,
                                  ties = mean)$y
    } else if(allow_extrapolation) {
      # Use linear model for extrapolation beyond data range
      m <- lm(temperature ~ response)
      temp_at_threshold <- predict(m, newdata = data.frame(response = threshold))
    } else {
      return(NA)
    }
  } else {
    # GLM approach (existing code)
    tryCatch({
      m <- suppressWarnings(glm(response ~ temperature, family = quasibinomial))
      slope <- coef(m)[2]
      if(abs(slope) < 0.0001) return(NA)
      temp_at_threshold <- (threshold - coef(m)[1]) / slope
    }, error = function(e) return(NA))
  }
  
  # Validate result
  if(is.finite(temp_at_threshold) && 
     temp_at_threshold >= temp_min && 
     temp_at_threshold <= temp_max) {
    return(temp_at_threshold)
  } else {
    return(NA)
  }
}


# Calculate CTmax and z-value from TDT data
calculate_ctmax_z <- function(temps, times, unit = "minute") {
  
  # Convert times to correct units before taking log10
  if(unit == "minute") {
    times_converted <- times 
  } else if(unit == "hour") {
    times_converted <- times / 60
  } else {
    stop("unit must be 'minute' or 'hour'")
  }
  
  # log10 of the converted times
  log_times <- log10(times_converted)
  
  # Fit linear model
  mod <- lm(log_times ~ temps)
  
  # Extract CTmax and z at the appropriate time unit
  list(
    ctmax = -coef(mod)[1] / coef(mod)[2],
    z = 1 / coef(mod)[2],
    r_squared = summary(mod)$r.squared
  )
}

# Generate prediction grid based on available data
generate_pred_grid <- function(sp) {
  
  # Get unique combinations that exist in the data
  combinations <- sp %>%
    dplyr::select(Temp, Time) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Time, Temp)
  
  pred <- data.frame(
    temp = combinations$Temp,
    time = combinations$Time
  )
  
  # Indexing for each time level
  time_levels <- sort(unique(pred$time))
  wtime_index <- lapply(time_levels, function(t) which(pred$time == t))
  
  list(pred = pred, wtime_index = wtime_index, time_levels = time_levels)
}

# Bootstrap predictions

bootstrap_predictions <- function(sp, pred_grid, response_vars, rep, resample = TRUE, progress = TRUE, n_cores = 4) {
  n_conditions <- nrow(pred_grid$pred)
  
  # Initialize output for each response variable
  predictions <- list()
  for(resp_var in names(response_vars)) {
    predictions[[resp_var]] <- matrix(NA, n_conditions, rep)
  }
  
  # Determine number of iterations
  n_iter <- if(resample) rep else 1
  
  # Setup parallel backend
  if(n_cores > 1 && n_iter > 1) {
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl))
    
    # Export necessary objects to workers
    clusterExport(cl, c("sp", "pred_grid", "response_vars", "resample", "n_conditions"), 
                  envir = environment())
    
    if(progress) {
      cat(sprintf("Running %d iterations on %d cores\n", n_iter, n_cores))
    }
    
    # Parallel bootstrap with progress bar
    results <- pblapply(1:n_iter, function(j) {
      
      if(resample) {
        dat <- sp[sample(1:nrow(sp), nrow(sp) - round(nrow(sp)/10, 0), replace = TRUE), ]
        dat$logTime <- log10(dat$Time)
      } else {
        dat <- sp
        dat$logTime <- log10(dat$Time)
      }
      
      # Fit models for each response variable
      models <- list()
      for(resp_var in names(response_vars)) {
        var_name <- response_vars[[resp_var]]$name
        
        tryCatch({
          models[[resp_var]] <- suppressWarnings(glm(
            as.formula(paste(var_name, "~ logTime * Temp")),
            data = dat, 
            family = quasibinomial
          ))
        }, error = function(e) {
          models[[resp_var]] <- NULL
        })
      }
      
      # Generate predictions
      iter_predictions <- list()
      for(resp_var in names(response_vars)) {
        iter_predictions[[resp_var]] <- rep(NA, n_conditions)
        
        if(!is.null(models[[resp_var]])) {
          for(i in 1:n_conditions) {
            newdata <- data.frame(
              Temp = pred_grid$pred$temp[i],
              logTime = log10(pred_grid$pred$time[i])
            )
            pred <- predict(models[[resp_var]], newdata = newdata, type = "response")
            iter_predictions[[resp_var]][i] <- pmax(0, pmin(1, pred))
          }
        }
      }
      
      iter_predictions
      
    }, cl = cl)  # Pass cluster for parallel execution
    
    # Reorganize results into matrix format
    for(resp_var in names(response_vars)) {
      for(j in 1:n_iter) {
        predictions[[resp_var]][, j] <- results[[j]][[resp_var]]
      }
    }
    
  } else {
    # Sequential version with progress
    for(j in 1:n_iter) {
      
      if(resample) {
        dat <- sp[sample(1:nrow(sp), nrow(sp) - round(nrow(sp)/10, 0), replace = TRUE), ]
        dat$logTime <- log10(dat$Time)
      } else {
        dat <- sp
        dat$logTime <- log10(dat$Time)
      }
      
      models <- list()
      for(resp_var in names(response_vars)) {
        var_name <- response_vars[[resp_var]]$name
        
        tryCatch({
          models[[resp_var]] <- suppressWarnings(glm(
            as.formula(paste(var_name, "~ logTime * Temp")),
            data = dat, 
            family = quasibinomial
          ))
        }, error = function(e) {
          models[[resp_var]] <- NULL
        })
      }
      
      for(i in 1:n_conditions) {
        newdata <- data.frame(
          Temp = pred_grid$pred$temp[i],
          logTime = log10(pred_grid$pred$time[i])
        )
        
        for(resp_var in names(response_vars)) {
          if(!is.null(models[[resp_var]])) {
            pred <- predict(models[[resp_var]], newdata = newdata, type = "response")
            predictions[[resp_var]][i, j] <- pmax(0, pmin(1, pred))
          }
        }
      }
      
      if(progress && j %% 10 == 0) {
        cat(sprintf("\r%d/%d iterations                    ", j, n_iter))
        flush.console()
      }
    }
    
    if(progress) cat("\n")
  }
  
  # If no resampling, replicate the single result across all columns
  if(!resample) {
    for(resp_var in names(response_vars)) {
      predictions[[resp_var]] <- matrix(
        rep(predictions[[resp_var]][, 1], rep),
        nrow = n_conditions,
        ncol = rep
      )
    }
  }
  
  return(predictions)
}



calculate_all_thresholds <- function(predictions, pred_grid, response_vars, 
                                     rep, temp_min = 15, temp_max = 80) {
  all_results <- list()
  
  # Process each response variable
  for(resp_var in names(response_vars)) {
    response_data <- predictions[[resp_var]]
    thresh_set <- response_vars[[resp_var]]$thresholds
    
    results <- list()
    
    for(thresh_name in names(thresh_set)) {
      thresh_value <- thresh_set[[thresh_name]]
      
      n_times <- length(pred_grid$time_levels)
      result_matrix <- matrix(NA, n_times, rep)
      
      for(i in seq_along(pred_grid$time_levels)) {
        time_idx <- pred_grid$wtime_index[[i]]
        
        for(j in 1:rep) {
          result_matrix[i, j] <- tryCatch({
            calculate_temp_threshold(
              response = response_data[time_idx, j],
              temperature = pred_grid$pred$temp[time_idx],
              threshold = thresh_value,
              temp_min = temp_min,
              temp_max = temp_max,
              use_interpolation = TRUE
            )
          }, error = function(e) NA)
        }
      }
      
      # Calculate mean and SD
      result_df <- data.frame(
        Time = pred_grid$time_levels,
        Mean = apply(result_matrix, 1, function(x) {
          x_valid <- x[is.finite(x) & x >= temp_min & x <= temp_max]
          if(length(x_valid) >= max(3, rep * 0.3)) {
            mean(x_valid, na.rm = TRUE)
          } else {
            NA
          }
        }),
        SD = apply(result_matrix, 1, function(x) {
          x_valid <- x[is.finite(x) & x >= temp_min & x <= temp_max]
          if(length(x_valid) >= max(3, rep * 0.3)) {
            sd(x_valid, na.rm = TRUE)
          } else {
            NA
          }
        })
      )
      
      names(result_df)[2:3] <- paste0(thresh_name, c("", "_sd"))
      results[[thresh_name]] <- result_df
    }
    
    all_results[[resp_var]] <- results
  }
  
  return(all_results)
}

# Calculate CTmax/z for all thresholds (with diagnostics)
calculate_tls_metrics <- function(threshold_results,
                                  min_r_squared = 0.5,
                                  temp_range = c(15, 80),
                                  z_range = c(-25, 5),
                                  min_temp_span = 1,
                                  warn = TRUE, verbose = FALSE) {
  tls_results <- list()
  
  for(thresh_name in names(threshold_results)) {
    df <- threshold_results[[thresh_name]]
    temp_col <- paste0(thresh_name, "")
    
    # Check if we have valid data
    temps <- df[[temp_col]]
    times <- df$Time
    
    valid_data <- is.finite(temps) & is.finite(times) & temps < 1e10
    
    if(verbose) {
      cat(sprintf("\n%s: temps = %s\n", thresh_name, 
                  paste(round(temps, 1), collapse = ", ")))
      cat(sprintf("%s: times = %s\n", thresh_name, 
                  paste(times, collapse = ", ")))
    }
    
    if(verbose) {
      cat(sprintf("%s: %d/%d valid data points\n", 
                  thresh_name, sum(valid_data), length(valid_data)))
    }
    
    if(sum(valid_data) >= 3) {
      # 1-hour units
      metrics_1h <- calculate_ctmax_z(temps, times, unit = "hour")
      
      # 1-minute units
      metrics_1m <- calculate_ctmax_z(temps, times, unit = "minute")
      
      # Check validity
      issues <- c()
      
      if(metrics_1h$r_squared < min_r_squared) {
        issues <- c(issues, sprintf("Poor fit (R² = %.3f)", metrics_1h$r_squared))
      }
      
      if(!is.na(metrics_1h$ctmax) && 
         (metrics_1h$ctmax < temp_range[1] || metrics_1h$ctmax > temp_range[2])) {
        issues <- c(issues, sprintf("CTmax out of range (%.1f°C)", metrics_1h$ctmax))
      }
      
      if(!is.na(metrics_1h$z) && 
         (metrics_1h$z < z_range[1] || metrics_1h$z > z_range[2])) {
        issues <- c(issues, sprintf("z out of range (%.1f)", metrics_1h$z))
      }
      
      # Check for flat line (near-zero slope)
      temp_range_data <- max(temps[valid_data]) - min(temps[valid_data])
      if(temp_range_data < min_temp_span) {
        issues <- c(issues, "Near-zero slope")
      }
      
      valid_fit <- length(issues) == 0
      
      if(warn && !valid_fit) {
        warning(sprintf("Invalid TLS fit for %s: %s", 
                        thresh_name, paste(issues, collapse = "; ")))
      }
      
      tls_results[[thresh_name]] <- data.frame(
        ctmax_1h = ifelse(valid_fit, metrics_1h$ctmax, NA),
        z_1h = ifelse(valid_fit, metrics_1h$z, NA),
        r_squared_1h = metrics_1h$r_squared,
        ctmax_1m = ifelse(valid_fit, metrics_1m$ctmax, NA),
        z_1m = ifelse(valid_fit, metrics_1m$z, NA),
        r_squared_1m = metrics_1m$r_squared,
        valid_fit = valid_fit,
        issues = ifelse(length(issues) > 0, paste(issues, collapse = "; "), NA)
      )
      
    } else {
      tls_results[[thresh_name]] <- data.frame(
        ctmax_1h = NA, z_1h = NA, r_squared_1h = NA,
        ctmax_1m = NA, z_1m = NA, r_squared_1m = NA,
        valid_fit = FALSE,
        issues = "Insufficient data points"
      )
    }
  }
  
  return(tls_results)
}


generate_diagnostic_plots <- function(sp, predictions, pred_grid, 
                                      threshold_results, tls_results, 
                                      k, response_vars) {
  
  # Get list of response types
  resp_types <- names(response_vars)
  
  cat(sprintf("\nGenerating plots for %d response types: %s\n", 
              length(resp_types), paste(resp_types, collapse = ", ")))
  
  # Create plots for each response type
  all_plots <- list()
  
  for(resp_var in resp_types) {
    cat(sprintf("  Processing %s...\n", resp_var))
    
    resp_name <- response_vars[[resp_var]]$description
    var_name <- response_vars[[resp_var]]$name
    
    # Check if variable exists in data
    if(!var_name %in% names(sp)) {
      cat(sprintf("  WARNING: %s not found in data. Skipping.\n", var_name))
      next
    }
    
    # Prepare prediction data for this response
    pred_df <- data.frame(
      temp = pred_grid$pred$temp,
      time = pred_grid$pred$time,
      fvfm_mean = apply(predictions[[resp_var]], 1, mean, na.rm = TRUE),
      fvfm_sd = apply(predictions[[resp_var]], 1, sd, na.rm = TRUE)
    )
    
    # Plot 1: Raw data with GLM fits (linear time scale)
    p1 <- ggplot(sp, aes(x = 10^logTime, y = .data[[var_name]], 
                         colour = factor(Temp))) +
      ggtitle(paste(sp$Species_sample[1], "\n", "-", resp_name)) +
      geom_point(size = 2) +
      geom_smooth(method = "glm", formula = y ~ x, 
                  method.args = list(family = "quasibinomial"), se = FALSE) +
      labs(y = expression(Scaled~italic(F)[v]/italic(F)[m]), 
           x = "Time (min)", colour = "Temperature") +
      scale_y_continuous(limits = c(0, 1)) +
      scale_colour_viridis_d(option = "inferno", end = 0.85) +
      theme_classic()
    
    # Plot 2: Raw data with linear fits (log time scale)
    p2 <- ggplot(sp, aes(x = logTime, y = .data[[var_name]], 
                         colour = factor(Temp))) +
      ggtitle(paste(sp$Species_sample[1], "\n", "-", resp_name)) +
      geom_point(size = 2) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      scale_x_continuous(
        breaks = log10(c(1, 5, 15, 30, 60, 120, 240)),
        labels = c("1", "5", "15", "30", "60", "120", "240")
      ) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(y = expression(Scaled~italic(F)[v]/italic(F)[m]), 
           x = expression(log[10]~Time~(min)), colour = "Temperature") +
      scale_colour_viridis_d(option = "inferno", end = 0.85) +
      theme_classic()
    
    # Plot 3: Temperature response curves
    p3 <- ggplot(pred_df, aes(x = temp, y = fvfm_mean, 
                              colour = factor(time), fill = factor(time))) +
      ggtitle(paste(sp$Species_sample[1], "\n", "-", resp_name)) +
      geom_point() +
      geom_hline(yintercept = c(0.1, 0.5, 0.9), 
                 linetype = 2, alpha = c(0.2, 0.5, 0.8)) +
      geom_smooth(method = "glm", formula = y ~ x,
                  method.args = list(family = "quasibinomial"), 
                  se = TRUE, fullrange = TRUE) +
      labs(y = expression(Scaled~italic(F)[v]/italic(F)[m]),
           x = expression(Temperature~(degree*C))) +
      scale_x_continuous(limits = c(20, 70)) +
      scale_y_continuous(limits = c(0, 1)) +
      scale_colour_viridis_d(option = "mako", direction = -1, end = 0.9, 
                             "Duration (min)") +
      scale_fill_viridis_d(option = "mako", direction = -1, end = 0.9, 
                           "Duration (min)") +
      theme_bw()
    
    threshold_plots <- list()
    
    # Get threshold names from the response_vars config
    thresh_names <- names(response_vars[[resp_var]]$thresholds)
    
    for(thresh_name in thresh_names) {  # Changed from hardcoded c("T50", "T90", "T10")
      if(thresh_name %in% names(threshold_results[[resp_var]])) {
        df <- threshold_results[[resp_var]][[thresh_name]]
        temp_col <- paste0(thresh_name, "")
        sd_col <- paste0(thresh_name, "_sd")
        
        # Get TLS metrics if available
        label_text <- ""
        label_colour <- "blue"
        pred_data <- NULL
        ci_data <- NULL
        
        if(!is.null(tls_results[[resp_var]]) && 
           thresh_name %in% names(tls_results[[resp_var]])) {
          tls <- tls_results[[resp_var]][[thresh_name]]
          
          # Get valid data points
          temps <- df[[temp_col]]
          times <- df$Time
          valid_idx <- is.finite(temps) & is.finite(times) & temps < 1e10
          
          if(sum(valid_idx) >= 3) {
            # Fit the model used to calculate CTmax_1h and z_1h
            times_hours <- times[valid_idx] / 60
            log_times <- log10(times_hours)
            temps_clean <- temps[valid_idx]
            
            # Fit model: log10(time_hours) ~ temp
            mod <- lm(log_times ~ temps_clean)
            
            # Get coefficients
            intercept <- coef(mod)[1]
            slope <- coef(mod)[2]
            
            # Generate predictions across time range (1 to 260 minutes)
            time_seq_minutes <- 10^seq(log10(1), log10(260), length.out = 260)
            time_seq_hours <- time_seq_minutes / 60
            log_time_hours <- log10(time_seq_hours)
            
            # Calculate temperature: temp = (log_time - intercept) / slope
            temp_pred <- (log_time_hours - intercept) / slope
            
            # Calculate standard errors for confidence interval
            temp_range <- range(temps_clean)
            temp_for_se <- seq(min(temp_pred, temp_range[1], na.rm = TRUE), 
                               max(temp_pred, temp_range[2], na.rm = TRUE), 
                               length.out = 260)
            
            pred_se <- predict(mod, 
                               newdata = data.frame(temps_clean = temp_for_se),
                               se.fit = TRUE)
            
            # Convert SE back to temperature space (approximate)
            temp_lower <- (pred_se$fit - 1.96 * pred_se$se.fit - intercept) / slope
            temp_upper <- (pred_se$fit + 1.96 * pred_se$se.fit - intercept) / slope
            
            pred_data <- data.frame(
              log_time = log10(time_seq_minutes),
              temp = temp_pred
            )
            
            pred_for_ci <- predict(mod, 
                                   newdata = data.frame(temps_clean = temp_pred),
                                   se.fit = TRUE)
            
            se_temp <- pred_for_ci$se.fit / abs(slope)
            
            ci_data <- data.frame(
              log_time = log10(time_seq_minutes),
              temp_lower = temp_pred - 1.96 * se_temp,
              temp_upper = temp_pred + 1.96 * se_temp
            )
          }
          
          # Create label text
          if(!is.na(tls$valid_fit) && tls$valid_fit) {
            label_text <- sprintf("CTmax_1h = %.1f°C\nCTmax_1m = %.1f°C\nz = %.1f\nR² = %.2f", 
                                  tls$ctmax_1h, tls$ctmax_1m, tls$z_1h, tls$r_squared_1h)
            label_colour <- "blue"
          } else {
            label_text <- sprintf("INVALID FIT\nR² = %.2f", 
                                  tls$r_squared_1h)
            label_colour <- "red"
          }
        }
        
        p <- ggplot(df, aes(x = log10(Time), y = .data[[temp_col]])) +
          ggtitle(paste(sp$Species_sample[1], "\n", thresh_name, "-", resp_name))
        
        # Add CI ribbon if available
        if(!is.null(pred_data) && !is.null(ci_data)) {
          # Match CI data to prediction data by log_time
          pred_data_merged <- merge(pred_data, ci_data, by = "log_time", all.x = TRUE)
          
          p <- p + geom_ribbon(data = pred_data_merged,
                               aes(x = log_time, ymin = temp_lower, ymax = temp_upper),
                               alpha = 0.2, fill = label_colour,
                               inherit.aes = FALSE)
        }
        
        # Add fitted line
        if(!is.null(pred_data)) {
          p <- p + geom_line(data = pred_data, 
                             aes(x = log_time, y = temp),
                             colour = label_colour, linewidth = 1,
                             inherit.aes = FALSE)
        }
        
        # Add data points on top
        p <- p + geom_pointrange(aes(ymin = .data[[temp_col]] - .data[[sd_col]], 
                                     ymax = .data[[temp_col]] + .data[[sd_col]]), 
                                 size = 0.35) +
          expand_limits(x = log10(c(1, 260)), y = c(40, 60)) +
          scale_x_continuous(
            breaks = log10(c(1, 5, 15, 30, 60, 120, 240)),
            labels = c("1", "5", "15", "30", "60", "120", "240")
          ) +
          labs(x = "Duration (min)", 
               y = expression(Heat~tolerance~(degree*C))) +
          theme_classic()
        
        # Add label if it exists
        if(label_text != "") {
          p <- p + annotate("text", x = Inf, y = Inf, label = label_text, 
                            hjust = 1.05, vjust = 1.5, colour = label_colour, size = 3)
        }
        
        # Store the plot
        threshold_plots[[thresh_name]] <- p
        
      } 
    } 
    
    # Combine plots for this response type
    response_plots <- c(list(p1, p2, p3), threshold_plots)
    all_plots[[resp_var]] <- response_plots
    
    cat(sprintf("  Created %d plots for %s\n", length(response_plots), resp_var))
  } 
  
  cat(sprintf("Total response types with plots: %d\n", length(all_plots)))
  
  return(all_plots)
}



generate_comparison_plots <- function(sp, predictions, pred_grid, 
                                      threshold_results, tls_results, 
                                      k, response_vars) {
  
  # Combine predictions from all response types for comparison
  pred_combined <- bind_rows(lapply(names(response_vars), function(resp_var) {
    data.frame(
      temp = pred_grid$pred$temp,
      time = pred_grid$pred$time,
      fvfm_mean = apply(predictions[[resp_var]], 1, mean, na.rm = TRUE),
      fvfm_sd = apply(predictions[[resp_var]], 1, sd, na.rm = TRUE),
      response_type = response_vars[[resp_var]]$description
    )
  }))
  
  # Temperature response curves comparison
  p_comp <- ggplot(pred_combined, 
                   aes(x = temp, y = fvfm_mean, 
                       colour = factor(time), 
                       linetype = response_type)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ response_type, ncol = 1) +
    geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.5) +
    labs(y = expression(Scaled~italic(F)[v]/italic(F)[m]),
         x = expression(Temperature~(degree*C)),
         colour = "Duration (min)",
         linetype = "Recovery Type",
         title = paste("Sample", k, "- Recovery Comparison")) +
    scale_colour_viridis_d(option = "mako", direction = -1, end = 0.9) +
    theme_bw()
  
  # TDT comparison for T50
  tdt_combined <- bind_rows(lapply(names(response_vars), function(resp_var) {
    if("T50" %in% names(threshold_results[[resp_var]])) {
      df <- threshold_results[[resp_var]]$T50
      df$response_type <- response_vars[[resp_var]]$description
      return(df)
    }
    return(NULL)
  }))
  
  p_tdt <- ggplot(tdt_combined, aes(x = log10(Time), y = T50, 
                                    colour = response_type)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = T50 - T50_sd, ymax = T50 + T50_sd), 
                  width = 0.1, alpha = 0.5) +
    scale_x_continuous(
      breaks = log10(c(1, 5, 15, 30, 60, 120, 240)),
      labels = c("1", "5", "15", "30", "60", "120", "240")
    ) +
    labs(x = "Duration (min)",
         y = expression(T[50]~(degree*C)),
         colour = "Recovery Type",
         title = paste("Sample", k, "- T50 Comparison")) +
    theme_bw()
  
  # Combine comparison plots
  p_comp / p_tdt + plot_layout(heights = c(2, 1))
}



# Function to compile results for a specific response type
compile_results_by_response <- function(all_results, response_var, threshold_name) {
  
  result_list <- lapply(all_results, function(x) {
    if(!is.null(x$thresholds[[response_var]][[threshold_name]])) {
      df <- x$thresholds[[response_var]][[threshold_name]]
      df$unique_sample <- x$tls_metrics[[response_var]][[threshold_name]]$unique_sample[1]
      return(df)
    }
    return(NULL)
  })
  
  result_list <- result_list[!sapply(result_list, is.null)]
  do.call(rbind, result_list)
}

# Compile TLS metrics for each response type
compile_tls_by_response <- function(all_results, response_var, threshold_name) {
  
  tls_list <- lapply(all_results, function(x) {
    if(!is.null(x$tls_metrics[[response_var]][[threshold_name]])) {
      return(x$tls_metrics[[response_var]][[threshold_name]])
    }
    return(NULL)
  })
  
  tls_list <- tls_list[!sapply(tls_list, is.null)]
  do.call(rbind, tls_list)
}

# Function to flatten nested results_summary into a single dataframe
flatten_results_summary <- function(results_summary) {
  
  # Get all response variables
  response_vars <- names(results_summary)
  
  all_data <- list()
  
  for(resp_var in response_vars) {
    
    # Extract TLS metrics for all thresholds
    tls_data <- results_summary[[resp_var]]$tls
    thresholds <- names(tls_data)
    
    for(thresh in thresholds) {
      
      # Get the dataframe for this threshold
      df <- tls_data[[thresh]]
      
      # Add threshold column if not present
      if(!"threshold" %in% names(df)) {
        df$threshold <- thresh
      }
      
      # Add response_type column if not present
      if(!"response_type" %in% names(df)) {
        df$response_type <- resp_var
      }
      
      # Store in list
      key <- paste(resp_var, thresh, sep = "_")
      all_data[[key]] <- df
    }
  }
  
  # Combine all dataframes
  combined_df <- bind_rows(all_data)
  
  # Reorder columns for better readability
  col_order <- c("unique_sample", "response_type", "threshold", 
                 "ctmax_1h", "ctmax_1m", "z_1h", "z_1m",
                 "r_squared_1h", "r_squared_1m", 
                 "valid_fit", "issues")
  
  # Select columns that exist
  available_cols <- intersect(col_order, names(combined_df))
  other_cols <- setdiff(names(combined_df), available_cols)
  
  combined_df <- combined_df %>%
    select(all_of(c(available_cols, other_cols)))
  
  return(combined_df)
}

# Function to create summary statistics by response_type and threshold
summarize_results <- function(flattened_df) {
  
  summary_stats <- flattened_df %>%
    filter(valid_fit == TRUE) %>%
    group_by(response_type, threshold) %>%
    summarise(
      n_valid = n(),
      
      # CTmax_1h stats
      ctmax_1h_mean = mean(ctmax_1h, na.rm = TRUE),
      ctmax_1h_sd = sd(ctmax_1h, na.rm = TRUE),
      ctmax_1h_min = min(ctmax_1h, na.rm = TRUE),
      ctmax_1h_max = max(ctmax_1h, na.rm = TRUE),
      
      # CTmax_1m stats
      ctmax_1m_mean = mean(ctmax_1m, na.rm = TRUE),
      ctmax_1m_sd = sd(ctmax_1m, na.rm = TRUE),
      ctmax_1m_min = min(ctmax_1m, na.rm = TRUE),
      ctmax_1m_max = max(ctmax_1m, na.rm = TRUE),
      
      # z_1h stats
      z_1h_mean = mean(z_1h, na.rm = TRUE),
      z_1h_sd = sd(z_1h, na.rm = TRUE),
      z_1h_min = min(z_1h, na.rm = TRUE),
      z_1h_max = max(z_1h, na.rm = TRUE),
      
      # R² stats
      r_squared_1h_mean = mean(r_squared_1h, na.rm = TRUE),
      r_squared_1h_sd = sd(r_squared_1h, na.rm = TRUE),
      r_squared_1h_min = min(r_squared_1h, na.rm = TRUE),
      r_squared_1h_max = max(r_squared_1h, na.rm = TRUE),
      
      .groups = "drop"
    )
  
  return(summary_stats)
}

# Function to compare valid vs invalid fits
compare_valid_invalid <- function(flattened_df) {
  
  validity_summary <- flattened_df %>%
    group_by(response_type, threshold) %>%
    summarise(
      n_total = n(),
      n_valid = sum(valid_fit == TRUE, na.rm = TRUE),
      n_invalid = sum(valid_fit == FALSE, na.rm = TRUE),
      pct_valid = round(100 * n_valid / n_total, 1),
      .groups = "drop"
    )
  
  return(validity_summary)
}

# Function to extract threshold data (T50, T90, T10 values across time)
flatten_threshold_data <- function(results_summary) {
  
  response_vars <- names(results_summary)
  
  all_threshold_data <- list()
  
  for(resp_var in response_vars) {
    
    # Extract threshold data
    threshold_data <- results_summary[[resp_var]]$thresholds
    thresholds <- names(threshold_data)
    
    for(thresh in thresholds) {
      
      df <- threshold_data[[thresh]]
      
      # Add metadata columns
      df$response_type <- resp_var
      df$threshold <- thresh
      
      # Rename the threshold column to "temperature" for consistency
      # Updated regex to match T followed by digits and optional decimal
      thresh_col <- names(df)[grepl("^T[0-9.]+$", names(df))]  # ← Changed here
      if(length(thresh_col) > 0) {
        names(df)[names(df) == thresh_col] <- "temperature"
        
        # Also rename the SD column
        sd_col <- paste0(thresh_col, "_sd")
        if(sd_col %in% names(df)) {
          names(df)[names(df) == sd_col] <- "temperature_sd"
        }
      }
      
      key <- paste(resp_var, thresh, sep = "_")
      all_threshold_data[[key]] <- df
    }
  }
  
  # Combine all dataframes
  combined_threshold_df <- bind_rows(all_threshold_data)
  
  # Reorder columns
  col_order <- c("unique_sample", "response_type", "threshold", 
                 "Time", "temperature", "temperature_sd")
  
  available_cols <- intersect(col_order, names(combined_threshold_df))
  other_cols <- setdiff(names(combined_threshold_df), available_cols)
  
  combined_threshold_df <- combined_threshold_df %>%
    select(all_of(c(available_cols, other_cols)))
  
  return(combined_threshold_df)
}

# Colourblind friendly colour palette options
cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbp3 <- c("#56B4E9", "#D55E00", "#CC79A7", "#0072B2", 
          "#F0E442", "#E69F00", "#009E73", "#000000")
