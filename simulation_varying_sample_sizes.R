# ==============================================================================
# Title: Estimation Following Surrogate Based Selection in Multi-Arm 
#        Drop-the-Losers Designs
# Description: Simulation study for estimating the primary endpoint after 
#              surrogate-based selection in a bivariate normal framework.
#              Calculates and plots Bias and MSE for Naive, Plug-in, and UMVCUE 
#              estimators across varying second-stage sample sizes (m).
# ==============================================================================

library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------------------------
# CORE SIMULATION FUNCTION
# ------------------------------------------------------------------------------
simulate_one <- function(n, m, mu_mat, k, sigma, Nsim = 20000) {
  set.seed(pi) # Set standard integer seed for reproducibility
  
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
    
    # 2. Selection (Based on the first component / surrogate endpoint)
    idx_sorted <- order(xj_mean[, 1], decreasing = TRUE)
    q <- idx_sorted[1]
    q_prime <- idx_sorted[2]
    
    mu_2q <- mu_mat[q, 2]
    xbar_1_qprime <- xj_mean[q_prime, 1]
    
    # 3. Second Stage Sampling (from the selected arm)
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
    
    # Naive 2 & 3 (Plug-in Adjustments)
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
k <- 2
N <- 20
m_seq <- seq(2, 18, 2)

# Parameter Setup mapping to Figure 1 
mu_mat <- matrix(0, nrow = k, ncol = 2)
sigma <- matrix(c(2, 1, 1, 2), 2, 2, byrow = TRUE)

results <- list()

# Execute Simulation
for (i in seq_along(m_seq)) {
  m <- m_seq[i]
  n <- 20 - m
  cat("Running for m =", m, " \n")
  
  # Nsim = 50000 provides highly stable results 
  results[[i]] <- simulate_one(n = n, m = m, k = k, mu_mat = mu_mat, sigma = sigma, Nsim = 50000)
}

results_df <- do.call(rbind, results)
results_df$estimator <- factor(results_df$estimator, levels = c("naive1", "naive2", "naive3", "UMVCUE"))

# ------------------------------------------------------------------------------
# VISUALIZATION 
# ------------------------------------------------------------------------------

# General plotting attributes for reuse
plot_colors <- c(naive1 = "black", naive2 = "green", naive3 = "red", UMVCUE = "skyblue")
plot_labels <- c(naive1 = expression(delta[N]), naive2 = expression(delta[P[1]]), 
                 naive3 = expression(delta[P[2]]), UMVCUE = expression(delta[U]))
plot_linetypes <- c(naive1 = "solid", naive2 = "longdash", naive3 = "dotted", UMVCUE = "dotdash")
plot_shapes <- c(naive1 = 16, naive2 = 15, naive3 = 17, UMVCUE = 3)

# 1. Plot Bias
p_bias <- results_df %>%
  ggplot(aes(x = m, y = mean_bias, color = estimator, linetype = estimator, shape = estimator)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.6, stroke = 1.2) +
  scale_color_manual(values = plot_colors, labels = plot_labels) +
  scale_linetype_manual(values = plot_linetypes) +
  scale_shape_manual(values = plot_shapes) +
  labs(x = "m", y = "Bias", color = "Estimators") +
  guides(
    color = guide_legend(override.aes = list(linetype = plot_linetypes, shape = plot_shapes, size = 1.2)),
    linetype = "none",
    shape = "none"
  ) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    legend.background = element_rect(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.key.width = unit(1.6, "cm")
  )

print(p_bias)

# 2. Plot MSE
p_mse <- results_df %>%
  ggplot(aes(x = m, y = mean_mse, color = estimator, linetype = estimator, shape = estimator)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.6, stroke = 1.2) +
  scale_color_manual(values = plot_colors, labels = plot_labels) +
  scale_linetype_manual(values = plot_linetypes) +
  scale_shape_manual(values = plot_shapes) +
  labs(x = "m", y = "MSE", color = "Estimators") +
  guides(
    color = guide_legend(override.aes = list(linetype = plot_linetypes, shape = plot_shapes, size = 1.2)),
    linetype = "none",
    shape = "none"
  ) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    legend.background = element_rect(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.key.width = unit(1.6, "cm")
  )

print(p_mse)
