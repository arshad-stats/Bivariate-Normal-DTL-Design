# ==============================================================================
# Title: Estimation Following Surrogate Based Selection in Multi-Arm 
#        Drop-the-Losers Designs
# Description: Simulation study evaluating Bias and MSE of estimators (Naive, 
#              Plug-in, UMVCUE) for the primary endpoint across varying 
#              surrogate effect separations (Delta mu_1).
# ==============================================================================

library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------------------------
# CORE SIMULATION FUNCTION
# ------------------------------------------------------------------------------
simulate_one <- function(n, m, mu_mat, k, sigma, Nsim = 20000) {
  set.seed(pi)
  # Initialize storage
  sum_bias <- rep(0, 4)
  sum_mse  <- rep(0, 4)
  sum_scaled_mse <- rep(0, 4)
  
  # Constant degrees of freedom and factor for UMVCUE
  nu <- k * n + m - k - 1
  factor_const <- sqrt(n / (m * (n + m)))
  
  for (i in 1:Nsim) {
    # 1. First Stage Sampling
    xj_mean <- matrix(0, nrow = k, ncol = 2)
    Ssum_x <- matrix(0, 2, 2)
    
    for (j in 1:k) {
      xj <- mvrnorm(n, mu_mat[j, ], sigma)
      xj_mean[j, ] <- colMeans(xj)
      # Sum of squares: (n-1) * var(xj)
      Ssum_x <- Ssum_x + (t(xj) %*% xj - n * tcrossprod(xj_mean[j, ]))
    }
    
    # 2. Selection (Based on the first component/surrogate)
    idx_sorted <- order(xj_mean[, 1], decreasing = TRUE)
    q <- idx_sorted[1]
    q_prime <- idx_sorted[2]
    
    mu_2q <- mu_mat[q, 2]
    xbar_1_qprime <- xj_mean[q_prime, 1]
    
    # 3. Second Stage Sampling (from the selected treatment arm)
    y <- mvrnorm(m, mu_mat[q, ], sigma)
    y_mean <- colMeans(y)
    Ssum_y <- t(y) %*% y - m * tcrossprod(y_mean)
    
    # 4. Sufficient Statistics & Pooling
    s_within <- Ssum_x + Ssum_y
    z_vec <- n * xj_mean[q, ] + m * y_mean
    z_over <- z_vec / (n + m)
    
    delta_diff <- z_over - y_mean
    s_pooled <- s_within + ((m * (n + m)) / n) * tcrossprod(delta_diff)
    
    s11 <- s_pooled[1, 1]
    s12 <- s_pooled[1, 2]
    
    # 5. Estimators
    # Naive 1 (Maximum Likelihood Estimator)
    delta_naive1 <- z_over[2]
    
    # Naive 2 & 3 (Plugin Adjustments)
    arg1 <- sqrt((n * (n + m)) / m) * (z_over[1] - xbar_1_qprime) / sqrt(s_within[1, 1] / nu)
    arg2 <- sqrt((n * (n + m)) / m) * (z_over[1] - xbar_1_qprime) / sqrt(s11 / (nu + 1))
    
    delta_naive2 <- z_over[2] - factor_const * ((s_within[1, 2] * dnorm(arg1)) / (sqrt(s_within[1, 1]) * pnorm(arg1)))
    delta_naive3 <- z_over[2] - factor_const * ((s12 * dnorm(arg2)) / (sqrt(s11) * pnorm(arg2)))
    
    # UMVCUE (delta_U)
    u_star <- sqrt(n * (n + m) / m) * (z_over[1] - xbar_1_qprime)
    min_val <- min(1, u_star / sqrt(s11))
    
    pow_term <- (1 - min_val^2)^(nu / 2)
    c_par <- nu / 2
    F_arg <- (1 + min_val) / 2
    
    F_val <- pbeta(F_arg, c_par, c_par)
    denom <- (2^(nu - 1)) * nu * sqrt(s11) * beta(c_par, c_par) * F_val
    
    delta_U <- z_over[2] - (factor_const * s12 * pow_term / denom)
    
    # 6. Accumulate results
    errors <- c(delta_naive1 - mu_2q, 
                delta_naive2 - mu_2q, 
                delta_naive3 - mu_2q, 
                delta_U - mu_2q)
    
    sum_bias <- sum_bias + errors
    sum_mse  <- sum_mse + errors^2
    sum_scaled_mse <- sum_scaled_mse + (errors^2 / sigma[2, 2])
  }
  
  return(data.frame(
    n = n, m = m,
    mu1_1 = mu_mat[1,1],
    estimator = c("naive1", "naive2", "naive3", "UMVCUE"),
    mean_bias = sum_bias / Nsim,
    mean_mse = sum_mse / Nsim,
    mean_scaled_mse = sum_scaled_mse / Nsim,
    stringsAsFactors = FALSE
  ))
}

# ------------------------------------------------------------------------------
# EXPERIMENT CONFIGURATION
# ------------------------------------------------------------------------------
# Fixed values
k <- 2
n <- 10
m <- 5
mu_mat <- matrix(0, nrow = k, ncol = 2)

fixed_s12 <- 0.1  # Set covariance sigma12 as a constant
mu1_seq <- seq(0, 2, by = 0.2) # Increased density for smoother line plots

results <- list()

# ------------------------------------------------------------------------------
# EXECUTION LOOP
# ------------------------------------------------------------------------------
for (i in seq_along(mu1_seq)) {
  mu1_val <- mu1_seq[i]
  mu1 <- c(mu1_val, 0)
  mu_mat[1, ] <- mu1
  
  sigma <- matrix(c(2, fixed_s12,
                    fixed_s12, 2), 2, 2, byrow = TRUE)
  
  cat(sprintf("Running for Delta mu1 = %.1f (fixed sigma12 = %.1f)\n", mu1_val, fixed_s12))
  
  # Nsim = 50000 provides highly stable results
  results[[i]] <- simulate_one(n = n, m = m, k = k, mu_mat = mu_mat, sigma = sigma, Nsim = 50000)
}

results_df <- do.call(rbind, results)

# Refactor estimators for consistent plotting
results_df$estimator <- factor(
  results_df$estimator,
  levels = c("naive1", "naive2", "naive3", "UMVCUE")
)

# ------------------------------------------------------------------------------
# VISUALIZATION
# ------------------------------------------------------------------------------

# Define common styling attributes to avoid repetition
plot_colors <- c("naive1" = "black", "naive2" = "green", "naive3" = "red", "UMVCUE" = "skyblue")
plot_labels <- c("naive1" = expression(delta[N]), "naive2" = expression(delta[P[1]]),
                 "naive3" = expression(delta[P[2]]), "UMVCUE" = expression(delta[U]))
plot_linetypes <- c("naive1" = "solid", "naive2" = "longdash", "naive3" = "dotted", "UMVCUE" = "dotdash")
plot_shapes <- c("naive1" = 16, "naive2" = 15, "naive3" = 17, "UMVCUE" = 3)

base_theme <- theme_bw(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    legend.background = element_rect(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.key.width = unit(1.6, "cm")
  )

legend_guide <- guides(
  color = guide_legend(
    override.aes = list(linetype = plot_linetypes, shape = plot_shapes, size = 1.2)
  ),
  linetype = "none",
  shape = "none"
)

# 1. Plot Bias
p_bias <- results_df %>%
  ggplot(aes(x = mu1_1, y = mean_bias, color = estimator, linetype = estimator, shape = estimator)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.6, stroke = 1.2) +
  scale_color_manual(values = plot_colors, labels = plot_labels) +
  scale_linetype_manual(values = plot_linetypes) +
  scale_shape_manual(values = plot_shapes) +
  labs(x = expression(Delta * mu[1]), y = "Bias", color = "Estimators") +
  legend_guide +
  base_theme

print(p_bias)

# 2. Plot MSE
p_mse <- results_df %>%
  ggplot(aes(x = mu1_1, y = mean_mse, color = estimator, linetype = estimator, shape = estimator)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.6, stroke = 1.2) +
  scale_color_manual(values = plot_colors, labels = plot_labels) +
  scale_linetype_manual(values = plot_linetypes) +
  scale_shape_manual(values = plot_shapes) +
  labs(x = expression(Delta * mu[1]), y = "MSE", color = "Estimators") +
  legend_guide +
  base_theme

print(p_mse)