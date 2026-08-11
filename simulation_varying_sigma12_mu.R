# ==============================================================================
# Title: Estimation Following Surrogate Based Selection in Multi-Arm 
#        Drop-the-Losers Designs
# Description: Simulation study generating contour plots to evaluate the Bias 
#              and MSE of estimators (Naive, Plug-in, UMVCUE) for the primary 
#              endpoint. This script evaluates a k=10 arm design across varying 
#              covariance (sigma_12) and surrogate effect separations (mu_1).
# ==============================================================================

# --- 1. Libraries ---
library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales) 

# --- 2. Simulation Function ---
simulate_one <- function(n, m, mu_mat, k, sigma, Nsim = 20000) {
  set.seed(pi)
  sum_bias <- rep(0, 4)
  sum_mse  <- rep(0, 4)
  sum_scaled_mse <- rep(0, 4)
  nu <- k * n + m - k - 1
  factor_const <- sqrt(n / (m * (n + m)))
  
  for (i in 1:Nsim) {
    xj_mean <- matrix(0, nrow = k, ncol = 2)
    Ssum_x <- matrix(0, 2, 2)
    
    for (j in 1:k) {
      xj <- mvrnorm(n, mu_mat[j, ], sigma)
      xj_mean[j, ] <- colMeans(xj)
      Ssum_x <- Ssum_x + (t(xj) %*% xj - n * tcrossprod(xj_mean[j, ]))
    }
    idx_sorted <- order(xj_mean[, 1], decreasing = TRUE)
    q <- idx_sorted[1]
    q_prime <- idx_sorted[2]
    
    mu_2q <- mu_mat[q, 2]
    xbar_1_qprime <- xj_mean[q_prime, 1]
    y <- mvrnorm(m, mu_mat[q, ], sigma)
    y_mean <- colMeans(y)
    Ssum_y <- t(y) %*% y - m * tcrossprod(y_mean)
    s_within <- Ssum_x + Ssum_y
    z_vec <- n * xj_mean[q, ] + m * y_mean
    z_over <- z_vec / (n + m)
    
    delta_diff <- z_over - y_mean
    s_pooled <- s_within + ((m * (n + m)) / n) * tcrossprod(delta_diff)
    
    s11 <- s_pooled[1, 1]
    s12 <- s_pooled[1, 2]
    
    # 5. Estimators
    # Naive 1
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
    mu1_1 = mu_mat[1,1], sigma12 = sigma[1,2],
    estimator = c("naive1", "naive2", "naive3", "UMVCUE"),
    mean_bias = sum_bias / Nsim,
    mean_mse = sum_mse / Nsim,
    mean_scaled_mse = sum_scaled_mse / Nsim,
    stringsAsFactors = FALSE
  ))
}

# --- 3. Configuration & Execution ---
n <- 10
m <- 10
k <- 10
# Increased density for smoother contour rendering
mu1_seq <- seq(-2, 2, by = 0.25) 
sigma12_seq <- seq(0.5, 2, by = 0.25)

results_list <- list()

for (mu1_val in mu1_seq) {
  mu_mat <- matrix(0, nrow = k, ncol = 2)
  mu_mat[1, ] <- c(mu1_val, mu1_val) 
  
  for (s12 in sigma12_seq) {
    sigma <- matrix(c(2, s12, 
                      s12, 2), nrow = 2, byrow = TRUE)
    
    cat(sprintf("Processing: mu1 = %5.2f | sigma12 = %5.2f\n", mu1_val, s12))
    
    # Lowered Nsim to balance the denser grid execution time
    sim_res <- simulate_one(n = n, m = m, k = k, mu_mat = mu_mat, sigma = sigma, Nsim = 10000)
    results_list[[length(results_list) + 1]] <- sim_res
  }
}

results_df <- bind_rows(results_list)

# --- 4. Plotting Setup & Function ---
est_labs_char <- c(
  naive1 = "delta[N]",
  naive2 = "delta[P[1]]",
  naive3 = "delta[P[2]]",
  UMVCUE = "delta[U]"
)

plot_simulation <- function(data, metric_name, title_label) {
  plot_data <- data %>%
    rename(mu1 = mu1_1, value = !!sym(metric_name))
  
  # Lock factor levels to maintain uniform facet order
  plot_data$estimator <- factor(plot_data$estimator, levels = c("naive1", "naive2", "naive3", "UMVCUE"))
  
  ggplot(plot_data, aes(x = mu1, y = sigma12, z = value)) +
    geom_contour_filled(bins = 12, show.legend = FALSE) + 
    geom_contour(color = "white", alpha = 0.3, linewidth = 0.2) +
    
    # Transparent dummy layer to force a continuous legend
    geom_point(aes(color = value), alpha = 0, size = 0.1) +
    
    facet_wrap(~ estimator, labeller = as_labeller(est_labs_char, label_parsed)) +
    
    scale_fill_viridis_d() +
    scale_color_viridis_c(
      name = title_label,
      breaks = scales::pretty_breaks(n = 8),
      guide = guide_colorbar(
        barheight = unit(15, "lines"), 
        barwidth = unit(1.2, "lines"),  
        frame.colour = "black",
        frame.linewidth = 0.5
      )
    ) +
    coord_cartesian(expand = FALSE) +
    labs(
      x = expression(mu[11] == mu[21]),
      y = expression(sigma[12])
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold", margin = margin(b = 10)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
      strip.text = element_text(face = "bold"),
      panel.spacing.x = unit(1.5, "lines") 
    )
}

# --- 5. Generate and Display Plots ---

# Bias 
plot_bias_all <- plot_simulation(results_df, "mean_bias", "Bias Values")
suppressWarnings(print(plot_bias_all))

# MSE 
plot_mse_all <- plot_simulation(results_df, "mean_mse", "MSE Values")
suppressWarnings(print(plot_mse_all))