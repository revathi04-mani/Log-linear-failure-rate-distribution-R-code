#----------------------------------------------------------------------
code for dataset-I of TTT plot, PP-plot, QQ-plot of Log-LFR distribution
#-----------------------------------------------------------------------
rm(list = ls(all = TRUE))
library(nleqslv)
#data set I(failure time datasets)
data<- c(0.040, 1.866, 2.385, 3.443, 0.301, 1.876, 2.481, 3.467, 0.309, 1.899,
2.610, 3.478, 0.557, 1.911, 2.625, 3.578, 0.943, 1.912, 2.632, 3.595, 1.070, 1.914, 2.646, 3.699,
1.124, 1.981, 2.661, 3.779, 1.248, 2.010, 2.688, 3.924, 1.281, 2.038, 2.823, 4.035, 1.281, 2.085,
2.890, 4.121, 1.303, 2.089, 2.902, 4.167, 1.432, 2.097, 2.934, 4.240, 1.480, 2.135, 2.962, 4.255,
1.505, 2.154, 2.964, 4.278, 1.506, 2.190, 3.000, 4.305, 1.568, 2.194, 3.103, 4.376, 1.615, 2.223,
3.114, 4.449, 1.619, 2.224, 3.117, 4.485, 1.652, 2.229, 3.166, 4.570, 1.652, 2.300, 3.344, 4.602,
1.757, 2.324, 3.376, 4.663)

n <- length(data)
# Define the PDF of the Log-lfr distribution
pdf_llfr <- function(x, alpha, beta, p) {
  if (alpha <= 0 || beta <=0 || p<=0 ) return(rep(1e-10, length(x)))
  pdf1<- (p-1)*(alpha + (beta*x))* exp(-((alpha*x) + (beta/2)* x^2))
  pdf2 <-(1-(1-p)*exp(-((alpha*x) + (beta/2) *x^2)))
  pdf<- pdf1/pdf2
  pdf[is.na(pdf) | is.infinite(pdf) | pdf <= 0] <- 1e-10  # Ensure support is x >= 0
  return(pdf)
}

# Define the CDF of the Log-lfr distribution
cdf_llfr<- function(x, alpha, beta, p) {
  if (alpha <= 0 || beta <=0 || p <=0) return(rep(1e-10, length(x)))
  #cdf <- 1 - exp(-((alpha*x) +((beta/2) *x^2)))
  cdf<- 1 - (log(1-(1-p)*exp(-((alpha*x) + (beta/2)* x^2)))/log(p))
  cdf[x < 0] <- 0 # Ensure support is x >= 0
  return(cdf)
}

# Negative log-likelihood function for parameter estimation
neg_log_likelihood <- function(params, data) {
  alpha <- params[1]
  beta<- params[2]
  p<- params[3]
  if (is.na(alpha)  ||  is.na(beta) || is.na(p) || alpha <= 0 || beta<=0 || p<=0) return(Inf)
  pdf_values <- sapply(data, pdf_llfr, alpha=alpha, beta=beta, p=p)
  if (any(pdf_values <= 0)) return(Inf)
  return(-sum(log(pdf_values)))
}


neg_log_likelihood_safe <- function(params, data) {
  result <- neg_log_likelihood(params, data)
  if (!is.finite(result)) {
    return(1e-10)  # Return a large but finite value
  }
  return(result)
}

# Fit the distribution using optim for MLE
fit_llfr<- function(data) {
  init_params <-  c(1, 1, 1)
  fit <- optim(init_params, neg_log_likelihood, data=data, method= "L-BFGS-B", 
               lower=c(0.01, 0.01, 0.01), upper=c(3.5, 3.5, 3.5), control = list(maxit = 1000, pgtol = 1e-8))
  return(fit)
}

# Kolmogorov-Smirnov test for goodness-of-fit
ks_test <- function(data, alpha, beta, p){
  empirical_cdf <- ecdf(data)
  theoretical_cdf <- function(x) sapply(x, function(xi) cdf_llfr(xi, alpha, beta, p))
  ks <- ks.test(data, theoretical_cdf)
  return(ks)
}

# Compute Information Criteria
compute_criteria <- function(log_likelihood, k, n) {
  aic <- 2 * k - 2 * log_likelihood
  bic <- k * log(n) - 2 * log_likelihood
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  hqic <- 2 * k * log(log(n)) - 2 * log_likelihood
  return(list(AIC = aic, BIC = bic, AICc = aicc, HQIC = hqic))
}

# Fit the model
fit <- fit_llfr(data)
alpha_hat <- fit$par[1]
beta_hat <- fit$par[2]
p_hat <- fit$par[3]
log_likelihood <- -fit$value

# Compute goodness-of-fit statistics and criteria
criteria <- compute_criteria(log_likelihood, k = 3, n = n)
ks_result <- ks_test(data, alpha_hat, beta_hat, p_hat)

# Sorting the data for plotting
sorted_data <- sort(data)
empcdf <- seq(1:n) / n

# Histogram with PDF overlay
hist(data, breaks = 20, freq = FALSE, col = "grey", border = "black", 
     xlab = "Data", ylab = "Density", main = "Histogram with Fitted PDF")
curve(pdf_llfr(x, alpha_hat, beta_hat, p_hat), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = c("Fitted PDF"), col = c("red"), lwd = c(2))

# P-P Plot
fitted_cdf <- sapply(sorted_data, function(x) cdf_llfr(x, alpha_hat, beta_hat, p_hat))
plot(fitted_cdf, empcdf, xlim = c(0, 1), ylim = c(0, 1), col = "red", 
     xlab = "Theoretical probabilities", main="P-P Plot", ylab = "Empirical probabilities")
abline(0, 1, col = "black", lwd = 2)

# Q-Q Plot
fitted_quantiles <- numeric(n)
for (i in 1:n) {
  quantile_eq <- function(xi) cdf_llfr(xi, alpha_hat, beta_hat, p_hat) - i / n
  result <- tryCatch({
    nleqslv(1, quantile_eq)$x
  }, error = function(e) {
    warning("Quantile calculation failed for i = ", i, ": ", e$message)
    return(NA)
  })
  fitted_quantiles[i] <- result
}
# Remove NA values if any
valid_indices <- !is.na(fitted_quantiles)
fitted_quantiles <- fitted_quantiles[valid_indices]
sorted_data <- sorted_data[valid_indices]

# Plot Q-Q Plot
plot(fitted_quantiles, sorted_data, xlim = c(0, max(data)), ylim = c(0, max(data)), 
     col = "red", xlab = "Theoretical Quantiles", main="Q-Q Plot", ylab = "Empirical Quantiles")

abline(0, 1, col = "black", lwd = 2)

# Total Time on Test (TTT) Plot
t_sorted <- sort(data)
n <- length(t_sorted)
cum_sum <- cumsum(t_sorted)
t_sum <- sum(t_sorted)
t_cum <- cum_sum / t_sum
p_cum <- (1:n) / n
plot(p_cum, t_cum, type = "l", col = "black", lwd = 2,
     xlab = "Proportion of Failures (i/n)",
     #ylab = "Cumulative TTT (Σt_i / Σt)",
     #ylab = "G(i/n)",
     main = "Total Time on Test (TTT) Plot)",
     ylab = "G(i/n)")
abline(0, 1, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("TTT Curve", "Reference Line"),
       col = c("black", "red"), lwd = 2, lty = c(1, 2))

# Print Results
cat("Estimated Parameter:\n")
cat("alpha:", alpha_hat, "\n")
cat("beta:", beta_hat, "\n")
cat("p:", p_hat, "\n")
cat("Negative Log-Likelihood:", -log_likelihood, "\n")
cat("AIC:", criteria$AIC, "\n")
cat("BIC:", criteria$BIC, "\n")
cat("AICc:", criteria$AICc, "\n")
cat("HQIC:", criteria$HQIC, "\n")
cat("KS Test Statistic:", ks_result$statistic, "\n")
cat("KS Test P-Value:", ks_result$p.value, "\n")
#-----------------------------------------------------------------------------------------------
#code for fitted Log-LFR distribution with other competing distribution models using Dataset-I
#-----------------------------------------------------------------------------------------------
 # Required libraries
library(stats)
library(ggplot2)


#dataset I(aircraft_failuretime datasets)
aircraft_data<- c(0.040, 1.866, 2.385, 3.443, 0.301, 1.876, 2.481, 3.467, 0.309, 1.899,
       2.610, 3.478, 0.557, 1.911, 2.625, 3.578, 0.943, 1.912, 2.632, 3.595, 1.070, 1.914, 2.646, 3.699,
      1.124, 1.981, 2.661, 3.779, 1.248, 2.010, 2.688, 3.924, 1.281, 2.038, 2.823, 4.035, 1.281, 2.085,
     2.890, 4.121, 1.303, 2.089, 2.902, 4.167, 1.432, 2.097, 2.934, 4.240, 1.480, 2.135, 2.962, 4.255,
    1.505, 2.154, 2.964, 4.278, 1.506, 2.190, 3.000, 4.305, 1.568, 2.194, 3.103, 4.376, 1.615, 2.223,
   3.114, 4.449, 1.619, 2.224, 3.117, 4.485, 1.652, 2.229, 3.166, 4.570, 1.652, 2.300, 3.344, 4.602,
  1.757, 2.324, 3.376, 4.663)


# Parameter estimates (dataset I)
params <- list(
LLFR = list(alpha =0.0306, beta =0.3326, p=3.5),
LFR= list(alpha=0.0107, beta =0.25),
RAY = list(alpha = 0.75),
EXP = list(lambda = 0.3910),
TLFR=list(beta=0.0767, theta=0.3037, lambda=-0.7675)
)

#  PDF function(LLFR)
pdf_llfr <- function(x, alpha, beta, p) {
  pdf1<- (p-1)*(alpha + (beta*x))* exp(-((alpha*x) + (beta/2)* x^2))
  pdf2 <-(1-(1-p)*exp(-((alpha*x) + (beta/2) *x^2)))
  pdf<- pdf1/pdf2
  return(pdf)
}

#  PDF function (LFR)
pdf_lfr <- function(x, alpha, beta) {
  pdf4 <- (alpha+(beta*x))*exp(-((alpha*x) + (beta/2)*x^2))
  return(pdf4)
}

#  PDF function(RAY)
pdf_ray <- function(x, alpha) {
  pdf5<-(x/alpha)*exp(-(x^2 /(2*alpha)))
  return(pdf5)
}

#  PDF function (EXP)
pdf_exp <- function(x,lambda) {
  pdf6<- lambda*exp(-(lambda*x))
  return(pdf6)
}
# PDF function (TLFR)
pdf_tlfr<- function(x, beta, theta, lambda){
  z <- beta*x+(theta/2)*x^2
  pdf7<- (beta+theta*x) *exp(-z)*(1-lambda+2*lambda*exp(-z))
  return(pdf7)
}

# CDF functions 
#  CDF (LLFR)
cdf_llfr <- function(x, alpha, beta, p) {
  cdf1<-  1 - (log(1-(1-p)*exp(-((alpha*x) + (beta/2)* x^2)))/log(p))
  return(cdf1)
}

#  CDF(LFR)
cdf_lfr<- function(x,alpha, beta) {
  cdf4<- 1- exp(-(alpha*x + (beta/2)*x^2))
  return(cdf4)
}

# CDF (RAY)
cdf_ray<- function(x, alpha) {
  cdf5<- 1-exp(-(x^2/(2*alpha)))
  return(cdf5)
}

#  CDF (EXP)
cdf_exp <- function(x,lambda) {
  cdf6<- 1-exp(-(lambda*x))
  return(cdf6)
}
#CDF (TLFR)Transmuted LFR
cdf_tlfr<- function(x, beta, theta, lambda){
  z <-  beta*x+(theta/2)*x^2
  cdf7 <-(1-exp(-z))*(1+lambda*exp(-z))
  return(cdf7)
}

# Create a sequence for plotting
x_seq <- seq(0, 5, by = 0.1)

# Calculate PDFs
pdf_values <- list(
  LLFR = sapply(x_seq, function(x) pdf_llfr(x, params$LLFR$alpha, params$LLFR$beta, 
                                            params$LLFR$p)),
  
  LFR= sapply(x_seq, function(x) pdf_lfr(x, params$LFR$alpha, params$LFR$beta)), 
  
  RAY= sapply(x_seq, function(x) pdf_ray(x, params$RAY$alpha)),
  EXP= sapply(x_seq, function(x) pdf_exp(x, params$EXP$lambda)),
  TLFR = sapply(x_seq, function(x) pdf_tlfr(x,params$TLFR$beta, 
                                            params$TLFR$theta, params$TLFR$lambda))
  
  )

# Calculate CDFs
cdf_values <- list(
  LLFR = sapply(x_seq, function(x) cdf_llfr(x, params$LLFR$alpha, params$LLFR$beta, 
                                            params$LLFR$p)),
  
  LFR = sapply(x_seq, function(x) cdf_lfr(x, params$LFR$alpha, params$LFR$beta)), 
  
  RAY= sapply(x_seq, function(x) cdf_ray(x, params$RAY$alpha)),
  EXP = sapply(x_seq, function(x) cdf_exp(x, params$EXP$lambda)),
  TLFR = sapply(x_seq, function(x) cdf_tlfr(x,params$TLFR$beta, 
                                            params$TLFR$theta, params$TLFR$lambda))
  
  
)


# Create data frames for plotting
pdf_df <- data.frame(
  x = rep(x_seq, 5),
  density = c(pdf_values$LLFR, pdf_values$LFR, 
              pdf_values$RAY, pdf_values$EXP, pdf_values$TLFR),
  Distribution = factor(rep(c("LLFR", "LFR", "RAY", "EXP", "TLFR"), 
                            each = length(x_seq)))
)

cdf_df <- data.frame(
  x = rep(x_seq, 5),
  probability = c(cdf_values$LLFR, cdf_values$LFR,
                  cdf_values$RAY, cdf_values$EXP, cdf_values$TLFR),
  Distribution = factor(rep(c("LLFR", "LFR", "RAY", "EXP", "TLFR"), 
                            each = length(x_seq)))
)

# Create colors for the plots
line_colors <- c("red", "green","blue", "magenta", "black")
line_types  <- c("twodash", "dashed", "dotdash", "longdash", "dotted")
#  Fitted PDFs
pdf_plot <- ggplot() +
  geom_histogram(data = data.frame(x = aircraft_data), aes(x = x, y = ..density..), 
                 bins = 15, color = "black", fill = "grey", alpha = 0.5) +
  geom_line(data = pdf_df, aes(x = x, y = density, color = Distribution, linetype=Distribution), size = 1) +
  scale_color_manual(values = line_colors) +  scale_color_manual(values = line_colors) +
  scale_linetype_manual(values = line_types) +
  labs(title = " ",
       x = "x",
       y = "probability density function") +
  theme_minimal() +
  xlim(0, 5.5) +
  theme(legend.position = "top")

# Create empirical CDF function
ecdf_values <- ecdf(aircraft_data)(x_seq)
ecdf_df <- data.frame(x = x_seq, probability = ecdf_values, Distribution = "Empirical CDF")

# Fitted CDFs
cdf_plot <- ggplot() +
  geom_step(data = ecdf_df, aes(x = x, y = probability), color = "black") +
  geom_line(data = cdf_df, aes(x = x, y = probability, color = Distribution), size = 1) +
  scale_color_manual(values = c(line_colors)) +
  labs(title = " ",
       x = "x",
       y = "CDF") +
  theme_minimal() +
  xlim(0, 5.5) +
  ylim(0, 1) +
  theme(legend.position = "right")


# Display plots

print(pdf_plot)
print(cdf_plot)

# To save the plots as shown in the paper
ggsave("Figure4_PDF.png", pdf_plot, width = 10, height = 6)
ggsave("Figure5_CDF.png", cdf_plot, width = 10, height = 6)
               
#----------------------------------------------------------------------
code for dataset-II of TTT plot, PP-plot, QQ-plot of Log-LFR distribution
#-----------------------------------------------------------------------
               
rm(list = ls(all = TRUE))
library(nleqslv)
#data set II(service time datasets)
data<- c(0.046, 1.436, 2.592, 0.140, 1.492, 2.600, 0.150, 1.580, 2.670, 0.248, 1.719,
         2.717, 0.280, 1.794, 2.819, 0.313, 1.915, 2.820, 0.389, 1.920, 2.878, 0.487, 1.963, 2.950, 0.622, 1.978, 3.003, 0.900,
         2.053, 3.102, 0.952, 2.065, 3.304, 0.996, 2.117, 3.483,1.003, 2.137, 3.500, 1.010, 2.141, 3.622, 1.085, 2.163, 3.665,
         1.092, 2.183, 3.695, 1.152, 2.240, 4.015, 1.183, 2.341, 4.628, 1.244, 2.435, 4.806, 1.249, 2.464, 4.881,1.262, 2.543,
         5.140)
n <- length(data)
# Define the PDF of the Log-lfr distribution
pdf_llfr <- function(x, alpha, beta, p) {
  if (alpha <= 0 || beta <=0 || p<=0 ) return(rep(1e-10, length(x)))
  pdf1<- (p-1)*(alpha + (beta*x))* exp(-((alpha*x) + (beta/2)* x^2))
  pdf2 <-(1-(1-p)*exp(-((alpha*x) + (beta/2) *x^2)))
  pdf<- pdf1/pdf2
  pdf[is.na(pdf) | is.infinite(pdf) | pdf <= 0] <- 1e-10  # Ensure support is x >= 0
  return(pdf)
}

# Define the CDF of the Log-lfr distribution
cdf_llfr<- function(x, alpha, beta, p) {
  if (alpha <= 0 || beta <=0 || p <=0) return(rep(1e-10, length(x)))
  #cdf <- 1 - exp(-((alpha*x) +((beta/2) *x^2)))
  cdf<- 1 - (log(1-(1-p)*exp(-((alpha*x) + (beta/2)* x^2)))/log(p))
  cdf[x < 0] <- 0 # Ensure support is x >= 0
  return(cdf)
}

# Negative log-likelihood function for parameter estimation
neg_log_likelihood <- function(params, data) {
  alpha <- params[1]
  beta<- params[2]
  p<- params[3]
  if (is.na(alpha)  ||  is.na(beta) || is.na(p) || alpha <= 0 || beta<=0 || p<=0) return(Inf)
  pdf_values <- sapply(data, pdf_llfr, alpha=alpha, beta=beta, p=p)
  if (any(pdf_values <= 0)) return(Inf)
  return(-sum(log(pdf_values)))
}


neg_log_likelihood_safe <- function(params, data) {
  result <- neg_log_likelihood(params, data)
  if (!is.finite(result)) {
    return(1e-10)  # Return a large but finite value
  }
  return(result)
}

# Fit the distribution using optim for MLE
fit_llfr<- function(data) {
  init_params <-  c(1, 1, 1)
  fit <- optim(init_params, neg_log_likelihood, data=data, method= "L-BFGS-B", 
               lower=c(0.01, 0.01, 0.01), upper=c(3.5, 3.5, 3.5), control = list(maxit = 1000, pgtol = 1e-8))
  return(fit)
}

# Kolmogorov-Smirnov test for goodness-of-fit
ks_test <- function(data, alpha, beta, p){
  empirical_cdf <- ecdf(data)
  theoretical_cdf <- function(x) sapply(x, function(xi) cdf_llfr(xi, alpha, beta, p))
  ks <- ks.test(data, theoretical_cdf)
  return(ks)
}

# Compute Information Criteria
compute_criteria <- function(log_likelihood, k, n) {
  aic <- 2 * k - 2 * log_likelihood
  bic <- k * log(n) - 2 * log_likelihood
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  hqic <- 2 * k * log(log(n)) - 2 * log_likelihood
  return(list(AIC = aic, BIC = bic, AICc = aicc, HQIC = hqic))
}

# Fit the model
fit <- fit_llfr(data)
alpha_hat <- fit$par[1]
beta_hat <- fit$par[2]
p_hat <- fit$par[3]
log_likelihood <- -fit$value

# Compute goodness-of-fit statistics and criteria
criteria <- compute_criteria(log_likelihood, k = 3, n = n)
ks_result <- ks_test(data, alpha_hat, beta_hat, p_hat)

# Sorting the data for plotting
sorted_data <- sort(data)
empcdf <- seq(1:n) / n

# Histogram with PDF overlay
hist(data, breaks = 20, freq = FALSE, col = "grey", border = "black", 
     xlab = "Data", ylab = "Density", main = "Histogram with Fitted PDF")
curve(pdf_llfr(x, alpha_hat, beta_hat, p_hat), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = c("Fitted PDF"), col = c("red"), lwd = c(2))

# P-P Plot
fitted_cdf <- sapply(sorted_data, function(x) cdf_llfr(x, alpha_hat, beta_hat, p_hat))
plot(fitted_cdf, empcdf, xlim = c(0, 1), ylim = c(0, 1), col = "red", 
     xlab = "Theoretical probabilities", main="P-P Plot", ylab = "Empirical probabilities")
abline(0, 1, col = "black", lwd = 2)

# Q-Q Plot
fitted_quantiles <- numeric(n)
for (i in 1:n) {
  quantile_eq <- function(xi) cdf_llfr(xi, alpha_hat, beta_hat, p_hat) - i / n
  result <- tryCatch({
    nleqslv(1, quantile_eq)$x
  }, error = function(e) {
    warning("Quantile calculation failed for i = ", i, ": ", e$message)
    return(NA)
  })
  fitted_quantiles[i] <- result
}
# Remove NA values if any
valid_indices <- !is.na(fitted_quantiles)
fitted_quantiles <- fitted_quantiles[valid_indices]
sorted_data <- sorted_data[valid_indices]

# Plot Q-Q Plot
plot(fitted_quantiles, sorted_data, xlim = c(0, max(data)), ylim = c(0, max(data)), 
     col = "red", xlab = "Theoretical Quantiles", main="Q-Q Plot", ylab = "Empirical Quantiles")
    
abline(0, 1, col = "black", lwd = 2)

# Total Time on Test (TTT) Plot
t_sorted <- sort(data)
n <- length(t_sorted)
cum_sum <- cumsum(t_sorted)
t_sum <- sum(t_sorted)
t_cum <- cum_sum / t_sum
p_cum <- (1:n) / n
plot(p_cum, t_cum, type = "l", col = "black", lwd = 2,
     xlab = "Proportion of Failures (i/n)",
     #ylab = "Cumulative TTT (Σt_i / Σt)",
     #ylab = "G(i/n)",
     main = "Total Time on Test (TTT) Plot)",
     ylab = "G(i/n)")
abline(0, 1, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("TTT Curve", "Reference Line"),
       col = c("black", "red"), lwd = 2, lty = c(1, 2))

# Print Results
cat("Estimated Parameter:\n")
cat("alpha:", alpha_hat, "\n")
cat("beta:", beta_hat, "\n")
cat("p:", p_hat, "\n")
cat("Negative Log-Likelihood:", -log_likelihood, "\n")
cat("AIC:", criteria$AIC, "\n")
cat("BIC:", criteria$BIC, "\n")
cat("AICc:", criteria$AICc, "\n")
cat("HQIC:", criteria$HQIC, "\n")
cat("KS Test Statistic:", ks_result$statistic, "\n")
cat("KS Test P-Value:", ks_result$p.value, "\n")
#-----------------------------------------------------------------------------------------------
#code for fitted Log-LFR distribution with other competing distribution models using Dataset-II
#-----------------------------------------------------------------------------------------------
# Required libraries
library(stats)
library(ggplot2)
#dataset II(aircraft_servicetime datasets)
aircraft_data<- c(0.046, 1.436, 2.592, 0.140, 1.492, 2.600, 0.150, 1.580, 2.670, 0.248, 1.719,
         2.717, 0.280, 1.794, 2.819, 0.313, 1.915, 2.820, 0.389, 1.920, 2.878, 0.487, 1.963, 2.950, 0.622, 1.978, 3.003, 0.900,
         2.053, 3.102, 0.952, 2.065, 3.304, 0.996, 2.117, 3.483,1.003, 2.137, 3.500, 1.010, 2.141, 3.622, 1.085, 2.163, 3.665,
         1.092, 2.183, 3.695, 1.152, 2.240, 4.015, 1.183, 2.341, 4.628, 1.244, 2.435, 4.806, 1.249, 2.464, 4.881,1.262, 2.543,
         5.140)
 #Parameter estimates (dataset II)
params <- list(
  LLFR = list(alpha =0.2782, beta =0.2653, p=3.5),
  LFR= list(alpha=0.1316, beta =0.2469),
  RAY = list(alpha = 0.75),
  EXP = list(lambda =0.4795),
  TLFR=list(beta=0.2275, theta=0.2435, lambda=-0.3849)
)
#  PDF function(LLFR)
pdf_llfr <- function(x, alpha, beta, p) {
  pdf1<- (p-1)*(alpha + (beta*x))* exp(-((alpha*x) + (beta/2)* x^2))
  pdf2 <-(1-(1-p)*exp(-((alpha*x) + (beta/2) *x^2)))
  pdf<- pdf1/pdf2
  return(pdf)
}

#  PDF function (LFR)
pdf_lfr <- function(x, alpha, beta) {
  pdf4 <- (alpha+(beta*x))*exp(-((alpha*x) + (beta/2)*x^2))
  return(pdf4)
}

#  PDF function(RAY)
pdf_ray <- function(x, alpha) {
  pdf5<-(x/alpha)*exp(-(x^2 /(2*alpha)))
  return(pdf5)
}

#  PDF function (EXP)
pdf_exp <- function(x,lambda) {
  pdf6<- lambda*exp(-(lambda*x))
  return(pdf6)
}
# PDF function (TLFR)
pdf_tlfr<- function(x, beta, theta, lambda){
  z <- beta*x+(theta/2)*x^2
  pdf7<- (beta+theta*x) *exp(-z)*(1-lambda+2*lambda*exp(-z))
  return(pdf7)
}
# CDF functions 
#  CDF (LLFR)
cdf_llfr <- function(x, alpha, beta, p) {
  cdf1<-  1 - (log(1-(1-p)*exp(-((alpha*x) + (beta/2)* x^2)))/log(p))
  return(cdf1)
}
#  CDF(LFR)
cdf_lfr<- function(x,alpha, beta) {
  cdf4<- 1- exp(-(alpha*x + (beta/2)*x^2))
  return(cdf4)
}
# CDF (RAY)
cdf_ray<- function(x, alpha) {
  cdf5<- 1-exp(-(x^2/(2*alpha)))
  return(cdf5)
}
#  CDF (EXP)
cdf_exp <- function(x,lambda) {
  cdf6<- 1-exp(-(lambda*x))
  return(cdf6)
}
#CDF (TLFR)Transmuted LFR
cdf_tlfr<- function(x, beta, theta, lambda){
  z <-  beta*x+(theta/2)*x^2
  cdf7 <-(1-exp(-z))*(1+lambda*exp(-z))
  return(cdf7)
}
# Create a sequence for plotting
x_seq <- seq(0, 5, by = 0.1)

# Calculate PDFs
pdf_values <- list(
  LLFR = sapply(x_seq, function(x) pdf_llfr(x, params$LLFR$alpha, params$LLFR$beta, 
                                            params$LLFR$p)),

  LFR= sapply(x_seq, function(x) pdf_lfr(x, params$LFR$alpha, params$LFR$beta)), 
  
  RAY= sapply(x_seq, function(x) pdf_ray(x, params$RAY$alpha)),
  EXP= sapply(x_seq, function(x) pdf_exp(x, params$EXP$lambda)),
  TLFR = sapply(x_seq, function(x) pdf_tlfr(x,params$TLFR$beta, 
                                            params$TLFR$theta, params$TLFR$lambda))
  
)

# Calculate CDFs
cdf_values <- list(
  LLFR = sapply(x_seq, function(x) cdf_llfr(x, params$LLFR$alpha, params$LLFR$beta, 
                                            params$LLFR$p)),
   
 LFR = sapply(x_seq, function(x) cdf_lfr(x, params$LFR$alpha, params$LFR$beta)), 
  
  RAY= sapply(x_seq, function(x) cdf_ray(x, params$RAY$alpha)),
  EXP = sapply(x_seq, function(x) cdf_exp(x, params$EXP$lambda)),
 TLFR = sapply(x_seq, function(x) cdf_tlfr(x,params$TLFR$beta, 
                                           params$TLFR$theta, params$TLFR$lambda)) 
)
# Create data frames for plotting
pdf_df <- data.frame(
  x = rep(x_seq, 5),
  density = c(pdf_values$LLFR, pdf_values$LFR, 
              pdf_values$RAY, pdf_values$EXP, pdf_values$TLFR),
  Distribution = factor(rep(c("LLFR", "LFR", "RAY", "EXP", "TLFR"), 
                            each = length(x_seq)))
)

cdf_df <- data.frame(
  x = rep(x_seq, 5),
  probability = c(cdf_values$LLFR, cdf_values$LFR,
                  cdf_values$RAY, cdf_values$EXP, cdf_values$TLFR),
  Distribution = factor(rep(c("LLFR", "LFR", "RAY", "EXP", "TLFR"), 
                            each = length(x_seq)))
)

# Create colors for the plots
line_colors <- c("red", "green","blue", "magenta", "black")
line_types  <- c("twodash", "dashed", "dotdash", "longdash", "dotted")
#  Fitted PDFs
pdf_plot <- ggplot() +
  geom_histogram(data = data.frame(x = aircraft_data), aes(x = x, y = ..density..), 
                 bins = 15, color = "black", fill = "grey", alpha = 0.5) +
  geom_line(data = pdf_df, aes(x = x, y = density, color = Distribution, linetype=Distribution), size = 1) +
  scale_color_manual(values = line_colors) +  scale_color_manual(values = line_colors) +
  scale_linetype_manual(values = line_types) +
  labs(title = " ",
       x = "x",
       y = "probability density function") +
  theme_minimal() +
  xlim(0, 5.5) +
  theme(legend.position = "top")

# Create empirical CDF function
ecdf_values <- ecdf(aircraft_data)(x_seq)
ecdf_df <- data.frame(x = x_seq, probability = ecdf_values, Distribution = "Empirical CDF")

# Fitted CDFs
cdf_plot <- ggplot() +
  geom_step(data = ecdf_df, aes(x = x, y = probability), color = "black") +
  geom_line(data = cdf_df, aes(x = x, y = probability, color = Distribution), size = 1) +
  scale_color_manual(values = c(line_colors)) +
  labs(title = " ",
       x = "x",
       y = "CDF") +
  theme_minimal() +
  xlim(0, 5.5) +
  ylim(0, 1) +
  theme(legend.position = "right")
# Display plots
print(pdf_plot)
print(cdf_plot)
# To save the plots as shown in the paper
ggsave("Figure4_PDF.png", pdf_plot, width = 10, height = 6)
ggsave("Figure5_CDF.png", cdf_plot, width = 10, height = 6)

















