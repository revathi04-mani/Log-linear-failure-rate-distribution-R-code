library(parallel)

qil_ll<- function(n, alpha, beta, alpha1) {
  u <- runif(n)
  if (any(u <= 0 | u >= 1)) stop("u must be between 0 and 1")
  if (alpha <= 0 || beta <= 0 || alpha1 <= 0) stop("Parameters must be positive")
  q_vals<- (1/beta) * (sqrt(alpha^2 - 2 * beta * log((1 - alpha1^(1-u))/(1 - alpha1))) - alpha)
  return(q_vals)
}

nloglik_il_ll <- function(params, data) {
  alpha <- params[1]
  beta <- params[2]
  alpha1<- params[3]
  n <- length(data)
  
  if (alpha <= 0 || beta <= 0 || alpha1 <= 0 ) {
    return(Inf)
  }
  
  ll <- n * log((alpha1 - 1) / log(alpha1)) + 
    sum(log(alpha + beta * data)) - 
    sum(alpha * data + (beta/2) * data^2) - 
    sum(log(1 - (1 - alpha1) * exp(-(alpha * data + (beta/2) * data^2))))
  
  return(-ll)   
}

fit_il_ll <- function(data, init_values = c(0.1, 0.1, 0.1)) {
  opt_result <- try(optim(
    init_values, 
    nloglik_il_ll, 
    data = data,
    method = "L-BFGS-B",
    lower = c(0.01, 0.01, 0.01),
    upper = c(5, 5, 5),
    control = list(maxit = 1000)
  ), silent = TRUE)
  
  if (!inherits(opt_result, "try-error") && opt_result$convergence == 0) {
    return(list(
      estimates = opt_result$par,
      success = TRUE
    ))
  } else {
    init_values2 <- c(0.2, 0.2, 0.2)
    opt_result2 <- try(optim(
      init_values2, 
      nloglik_il_ll, 
      data = data,
      method = "L-BFGS-B",
      lower = c(0.01, 0.01, 0.01),
      upper = c(5, 5, 5),
      control = list(maxit = 1000)
    ), silent = TRUE)
    
    if (!inherits(opt_result2, "try-error") && opt_result2$convergence == 0) {
      return(list(
        estimates = opt_result2$par,
        success = TRUE
      ))
    } else {
      return(list(
        estimates = c(NA, NA, NA),
        success = FALSE
      ))
    }
  }
}

rel_function <- function(x, alpha, beta, alpha1) {
  if (any(x < 0)) stop("x must be >=0")
  if (alpha <= 0 || beta <= 0 || alpha1 <= 0) stop("parameters must be positive")
  log(1 - (1 - alpha1) * exp(-(alpha * x + (beta/2) * x^2))) / log(alpha1)
}

haz_function <- function(x, alpha, beta, alpha1) {
  if (any(x < 0))  stop("x must be >=0")
  if (alpha <= 0 || beta <= 0 || alpha1 <= 0) stop("parameters must be positive")
  num <- (alpha1-1) * (alpha + beta * x) * exp(-(alpha * x + (beta/2) * x^2))
  denom <- (1 - (1 - alpha1) * exp(-(alpha * x + (beta/2) * x^2))) * 
    log(1 - (1 - alpha1) * exp(-(alpha * x + (beta/2) * x^2)))
  return(num / denom)
}

run_single_replication <- function(args) {
  set.seed(args$seed)
  n <- args$n
  true_alpha <- args$true_alpha
  true_beta <- args$true_beta
  true_alpha1 <- args$true_alpha1
  x_values <- args$x_values
  
  sample_data <- qil_ll(n, true_alpha, true_beta, true_alpha1)
  
  fit_result <- fit_il_ll(sample_data)
  
  if (fit_result$success) {
    est_alpha <- fit_result$estimates[1]
    est_beta <- fit_result$estimates[2]
    est_alpha1 <- fit_result$estimates[3]
    
    
    result <- list(
      n = n,
      alpha_est = est_alpha,
      beta_est = est_beta,
      alpha1_est = est_alpha1,
      success = TRUE
    )
    
    
    for (x in x_values) {
      rel_true <- rel_function(x, true_alpha, true_beta, true_alpha1)
      rel_est <- rel_function(x, est_alpha, est_beta, est_alpha1)
      haz_true <- haz_function(x, true_alpha, true_beta, true_alpha1)
      haz_est <- haz_function(x, est_alpha, est_beta, est_alpha1)
      
      result[[paste0("rel_true_", x)]] <- rel_true
      result[[paste0("rel_est_", x)]] <- rel_est
      result[[paste0("haz_true_", x)]] <- haz_true
      result[[paste0("haz_est_", x)]] <- haz_est
    }
    
    return(result)
  } else {
    result <- list(
      n = n,
      alpha_est = NA,
      beta_est = NA,
      alpha1_est = NA,
      success = FALSE
    )
    
    for (x in x_values) {
      result[[paste0("rel_true_", x)]] <- NA
      result[[paste0("rel_est_", x)]] <- NA
      result[[paste0("haz_true_", x)]] <- NA
      result[[paste0("haz_est_", x)]] <- NA
    }
    
    return(result)
  }
}

# Function to calculate confidence intervals
calculate_ci <- function(values, confidence_level = 0.95) {
  if (all(is.na(values)) || length(values) == 0) {
    return(list(lower = NA, upper = NA))
  }
  
  alpha_level <- 1 - confidence_level
  lower_percentile <- alpha_level / 2
  upper_percentile <- 1 - alpha_level / 2
  
  values_clean <- values[!is.na(values)]
  
  if (length(values_clean) == 0) {
    return(list(lower = NA, upper = NA))
  }
  
  ci_lower <- quantile(values_clean, lower_percentile, na.rm = TRUE)
  ci_upper <- quantile(values_clean, upper_percentile, na.rm = TRUE)
  
  return(list(lower = ci_lower, upper = ci_upper))
}

run_simulation <- function(true_alpha, true_beta, true_alpha1,
                           sample_sizes = c(50, 100, 200, 350),
                           replications = 1000,
                           x_values = c(2.5),
                           num_cores = NULL,
                           confidence_level = 0.95){
  
  cat("Running Monte Carlo simulation with parameters:\n")
  cat("alpha =", true_alpha, ", beta =", true_beta, ", alpha1 =", true_alpha1,  "\n")
  cat("Replications:", replications, "\n")
  cat("Confidence Level:", confidence_level, "\n")
  
  if (is.null(num_cores)) {
    num_cores <- max(1, detectCores() - 1)
  }
  cat("Using", num_cores, "cores\n")
  cl <- makeCluster(num_cores)
  clusterExport(cl, c("qil_ll", "nloglik_il_ll", "fit_il_ll", "run_single_replication", 
                      "rel_function", "haz_function"), 
                envir = environment())
  
  
  args_list <- list()
  counter <- 1
  for (n in sample_sizes) {
    for (i in 1:replications) {
      args_list[[counter]] <- list(
        n = n,
        true_alpha = true_alpha,
        true_beta = true_beta,
        true_alpha1= true_alpha1,
        x_values = x_values,
        seed = i * 1000 + counter 
      )
      counter <- counter + 1
    }
  }
  
  
  cat("Starting simulations...\n")
  start_time <- Sys.time()
  
  results_list <- parLapply(cl, args_list, run_single_replication)
  
  
  stopCluster(cl)
  
  end_time <- Sys.time()
  cat("Simulations completed in", difftime(end_time, start_time, units = "mins"), "minutes\n")
  
  all_cols <- unique(unlist(lapply(results_list, names)))
  results_df <- data.frame(matrix(NA, nrow = length(results_list), ncol = length(all_cols)))
  names(results_df) <- all_cols
  
  for (i in 1:length(results_list)) {
    result <- results_list[[i]]
    for (col in names(result)) {
      results_df[i, col] <- result[[col]]
    }
  }
  
  stats_table <- data.frame()
  
  for (n in sample_sizes) {
    n_results <- results_df[results_df$n == n & results_df$success == TRUE, ]
    n_count <- nrow(n_results)
    
    if (n_count > 0) {
      # Alpha parameter
      alpha_mean <- mean(n_results$alpha_est, na.rm = TRUE)
      alpha_bias <- alpha_mean - true_alpha
      alpha_mse <- mean((n_results$alpha_est - true_alpha)^2, na.rm = TRUE)
      alpha_ci <- calculate_ci(n_results$alpha_est, confidence_level)
      
      # Beta parameter
      beta_mean <- mean(n_results$beta_est, na.rm = TRUE)
      beta_bias <- beta_mean - true_beta
      beta_mse <- mean((n_results$beta_est - true_beta)^2, na.rm = TRUE)
      beta_ci <- calculate_ci(n_results$beta_est, confidence_level)
      
      # Alpha1 parameter
      alpha1_mean <- mean(n_results$alpha1_est, na.rm = TRUE)
      alpha1_bias <- alpha1_mean - true_alpha1
      alpha1_mse <- mean((n_results$alpha1_est - true_alpha1)^2, na.rm = TRUE)
      alpha1_ci <- calculate_ci(n_results$alpha1_est, confidence_level)
      
      
      param_rows <- data.frame(
        Type = rep("Parameter", 3),
        Parameter = c("alpha", "beta", "alpha1"),
        n = rep(n, 3),
        x = rep(NA, 3),
        Bias = c(alpha_bias, beta_bias, alpha1_bias),
        MSE = c(alpha_mse, beta_mse, alpha1_mse),
        Estimate = c(alpha_mean, beta_mean, alpha1_mean),
        TrueValue = c(true_alpha, true_beta, true_alpha1),
        CI_Lower = c(alpha_ci$lower, beta_ci$lower, alpha1_ci$lower),
        CI_Upper = c(alpha_ci$upper, beta_ci$upper, alpha1_ci$upper),
        CI_Width = c(alpha_ci$upper - alpha_ci$lower, 
                     beta_ci$upper - beta_ci$lower, 
                     alpha1_ci$upper - alpha1_ci$lower),
        Coverage = c(
          ifelse(!is.na(alpha_ci$lower) && !is.na(alpha_ci$upper), 
                 true_alpha >= alpha_ci$lower && true_alpha <= alpha_ci$upper, NA),
          ifelse(!is.na(beta_ci$lower) && !is.na(beta_ci$upper), 
                 true_beta >= beta_ci$lower && true_beta <= beta_ci$upper, NA),
          ifelse(!is.na(alpha1_ci$lower) && !is.na(alpha1_ci$upper), 
                 true_alpha1 >= alpha1_ci$lower && true_alpha1 <= alpha1_ci$upper, NA)
        )
      )
      
      stats_table <- rbind(stats_table, param_rows)
      
      for (x in x_values) {
        rel_true_col <- paste0("rel_true_", x)
        rel_est_col <- paste0("rel_est_", x)
        haz_true_col <- paste0("haz_true_", x)
        haz_est_col <- paste0("haz_est_", x)
        
        
        rel_true <- mean(n_results[[rel_true_col]], na.rm = TRUE)
        rel_mean <- mean(n_results[[rel_est_col]], na.rm = TRUE)
        rel_bias <- rel_mean - rel_true
        rel_mse <- mean((n_results[[rel_est_col]] - n_results[[rel_true_col]])^2, na.rm = TRUE)
        rel_ci <- calculate_ci(n_results[[rel_est_col]], confidence_level)
        
        
        haz_true <- mean(n_results[[haz_true_col]], na.rm = TRUE)
        haz_mean <- mean(n_results[[haz_est_col]], na.rm = TRUE)
        haz_bias <- haz_mean - haz_true
        haz_mse <- mean((n_results[[haz_est_col]] - n_results[[haz_true_col]])^2, na.rm = TRUE)
        haz_ci <- calculate_ci(n_results[[haz_est_col]], confidence_level)
        
        func_rows <- data.frame(
          Type = c("Reliability", "Hazard"),
          Parameter = c("R(x)", "h(x)"),
          n = c(n, n),
          x = c(x, x),
          Bias = c(rel_bias, haz_bias),
          MSE = c(rel_mse, haz_mse),
          Estimate = c(rel_mean, haz_mean),
          TrueValue = c(rel_true, haz_true),
          CI_Lower = c(rel_ci$lower, haz_ci$lower),
          CI_Upper = c(rel_ci$upper, haz_ci$upper),
          CI_Width = c(rel_ci$upper - rel_ci$lower, haz_ci$upper - haz_ci$lower),
          Coverage = c(
            ifelse(!is.na(rel_ci$lower) && !is.na(rel_ci$upper), 
                  rel_true >= rel_ci$lower && rel_true <= rel_ci$upper, NA),
            ifelse(!is.na(haz_ci$lower) && !is.na(haz_ci$upper), 
                   haz_true >= haz_ci$lower && haz_true <= haz_ci$upper, NA)
          )
        )
        
        stats_table <- rbind(stats_table, func_rows)
      }
    } else {
      param_rows <- data.frame(
        Type = rep("Parameter", 3),
        Parameter = c("alpha", "beta", "alpha1"),
        n = rep(n, 3),
        x = rep(NA, 3),
        Bias = rep(NA, 3),
        MSE = rep(NA, 3),
        Estimate = rep(NA, 3),
        TrueValue = c(true_alpha, true_beta, true_alpha1),
        CI_Lower = rep(NA, 3),
        CI_Upper = rep(NA, 3),
        CI_Width = rep(NA, 3),
        Coverage = rep(NA, 3)
      )
      
      stats_table <- rbind(stats_table, param_rows)
      
      for (x in x_values) {
        func_rows <- data.frame(
          Type = c("Reliability", "Hazard"),
          Parameter = c("R(x)", "h(x)"),
          n = c(n, n),
          x = c(x, x),
          Bias = c(NA, NA),
          MSE = c(NA, NA),
          Estimate = c(NA, NA),
          TrueValue = c(
            rel_function(x, true_alpha, true_beta, true_alpha1),
            haz_function(x, true_alpha, true_beta, true_alpha1)
          ),
          CI_Lower = c(NA, NA),
          CI_Upper = c(NA, NA),
          CI_Width = c(NA, NA),
          Coverage = c(NA, NA)
        )
        
        stats_table <- rbind(stats_table, func_rows)
      }
    }
  }
  
  return(stats_table)
}


format_table <- function(results) {
  
  formatted <- results
  formatted$Bias <- sprintf("%.4f", formatted$Bias)
  formatted$MSE <- sprintf("%.4f", formatted$MSE)
  formatted$Estimate <- sprintf("%.4f", formatted$Estimate)
  formatted$TrueValue <- sprintf("%.4f", formatted$TrueValue)
  formatted$CI_Lower <- sprintf("%.4f", formatted$CI_Lower)
  formatted$CI_Upper <- sprintf("%.4f", formatted$CI_Upper)
  formatted$CI_Width <- sprintf("%.4f", formatted$CI_Width)
  formatted$Coverage <- ifelse(is.na(formatted$Coverage), "NA", 
                               ifelse(formatted$Coverage, "Yes", "No"))
  
  # Add confidence interval as a combined column
  formatted$CI_Interval <- paste0("(", formatted$CI_Lower, ", ", formatted$CI_Upper, ")")
  
  
  type_order <- c("Parameter", "Reliability", "Hazard")
  param_order <- c("alpha", "beta", "alpha1", "R(x)", "h(x)")
  
  formatted$Type <- factor(formatted$Type, levels = type_order)
  formatted$Parameter <- factor(formatted$Parameter, levels = param_order)
  
  
  formatted <- formatted[order(formatted$n, formatted$Type, formatted$Parameter, formatted$x), ]
  
  return(formatted)
}

# Function to create a summary table with confidence intervals
create_summary_table <- function(results) {
  summary_cols <- c("Type", "Parameter", "n", "x", "Bias", "MSE", "Estimate", 
                    "TrueValue", "CI_Interval", "CI_Width", "Coverage")
  return(results[, summary_cols])
}

run_and_print_results <- function(true_alpha, true_beta, true_alpha1, 
                                  sample_sizes = c(50, 100, 200, 350),
                                  replications = 1000,
                                  x_values = c(2.5),
                                  confidence_level = 0.95) {
  
  results <- run_simulation(
    true_alpha = true_alpha,
    true_beta = true_beta,
    true_alpha1= true_alpha1,
    sample_sizes = sample_sizes,
    replications = replications,
    x_values = x_values,
    confidence_level = confidence_level
  )
  
  
  formatted_table <- format_table(results)
  summary_table <- create_summary_table(formatted_table)
  
  cat("\n=== SIMULATION RESULTS WITH CONFIDENCE INTERVALS ===\n")
  print(summary_table, row.names = FALSE)
  
  # Print coverage summary
  cat("\n=== COVERAGE SUMMARY ===\n")
  for (param in c("alpha", "beta", "alpha1")) {
    param_results <- results[results$Parameter == param & results$Type == "Parameter", ]
    coverage_rate <- mean(param_results$Coverage, na.rm = TRUE)
    cat(sprintf("%s Coverage Rate: %.3f\n", param, coverage_rate))
  }
  
  return(list(full_results = formatted_table, summary = summary_table))
}

# case I
#results <- run_and_print_results(
 # true_alpha = 0.4,
  #true_beta = 0.5, 
 # true_alpha1 = 0.2,
  #sample_sizes = c(50, 100, 200, 350),
  #replications = 1000,
 # x_values = c(2.5),
 #x_values = c(3),
 # confidence_level = 0.95
#)

# case II
results <- run_and_print_results(
  true_alpha = 0.3,
  true_beta = 0.8, 
  true_alpha1 = 0.1,
  sample_sizes = c(50, 100, 200, 350),
  replications = 1000,
   #x_values = c(2.5),
  x_values = c(3),
  confidence_level = 0.95
)
