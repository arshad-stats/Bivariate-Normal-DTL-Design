# ==============================================================================
# Title: Real Data Application - Prostate Cancer Dataset
# Description: Applies the two-stage drop-the-losers design estimation 
#              methodology to a clinical dataset. Evaluates bivariate normality, 
#              covariance homogeneity, and calculates primary endpoint 
#              estimators (Naive, Plug-in, UMVCUE) after surrogate selection.
# ==============================================================================

# URL for the raw prostate cancer data file
data_url <- "https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data"
prostate_data <- read.table(data_url, header = TRUE, sep = "\t")
prostate_data <- prostate_data[, -1]

# Group 1: Gleason Score equals 6
data_group_1 <- prostate_data[prostate_data$gleason == 6, ]

X1 <- data_group_1$lweight
Y1 <- data_group_1$lcavol

# Group 2: Gleason Score more than 6 (Select rows where 'gleason' is greater than 6)
data_group_2 <- prostate_data[prostate_data$gleason > 6, ]

X2 <- data_group_2$lweight
Y2 <- data_group_2$lcavol

# --- 3. Bivariate Normality Test (Royston Test) ---
library(MVN)
data_X1Y1 <- data.frame(lweight = X1, lcavol = Y1)
royston_test_1 <- mvn(data = data_X1Y1, mvn_test = "royston")

# Combine X2 and Y2 into a single data frame for the Royston test
data_X2Y2 <- data.frame(lweight = X2, lcavol = Y2)
royston_test_2 <- mvn(data = data_X2Y2, mvn_test = "royston")

cat("\n--- Bivariate Normality Test Results (Royston Test) ---\n")

print(royston_test_1$multivariate_normality)
print(royston_test_2$multivariate_normality)


# --- 4. Test for Equality of Covariance Matrices (Box's M Test) ---
library(biotools)

combined_data <- data.frame(
  lweight = c(X1, X2), 
  lcavol = c(Y1, Y2), 
  Group = factor(c(rep("Gleason_6", length(X1)), rep("Gleason_GT_6", length(X2))))
)

box_m_test <- boxM(combined_data[, c("lweight", "lcavol")], combined_data$Group)
print(box_m_test)




#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
# Considering n1 = n2 = 30

data_url <- "https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data"
prostate_data <- read.table(data_url, header = TRUE, sep = "\t")
prostate_data <- prostate_data[, -1]

# Group 1: Gleason Score equals 6
data_group_1 <- prostate_data[prostate_data$gleason == 6, ]
X1 <- data_group_1$lweight[1:30]
Y1 <- data_group_1$lcavol[1:30]

data_group_2 <- prostate_data[prostate_data$gleason > 6, ]
X2 <- data_group_2$lweight[1:30]
Y2 <- data_group_2$lcavol[1:30]

library(MVN)
data_X1Y1 <- data.frame(lweight = X1, lcavol = Y1)
royston_test_1 <- mvn(data = data_X1Y1, mvn_test = "royston")

data_X2Y2 <- data.frame(lweight = X2, lcavol = Y2)
royston_test_2 <- mvn(data = data_X2Y2, mvn_test = "royston")

cat("\n--- Bivariate Normality Test Results (Royston Test) ---\n")

print(royston_test_1$multivariate_normality)
print(royston_test_2$multivariate_normality)


# --- 4. Test for Equality of Covariance Matrices (Box's M Test) ---
library(biotools)

combined_data <- data.frame(
  lweight = c(X1, X2), 
  lcavol = c(Y1, Y2), 
  Group = factor(c(rep("Gleason_6", length(X1)), rep("Gleason_GT_6", length(X2))))
)

box_m_test <- boxM(combined_data[, c("lweight", "lcavol")], combined_data$Group)
print(box_m_test)
#----------------------------------
#----------------------------------
#----------------------------------

Estimates <- function(x1_mat,x2_mat, q, y_mat) {
  set.seed(pi)
  n = length(x1_mat[,1])
  m = length(y_mat[,1])
  k=2
  nu <- k * n + m - k - 1
  factor_const <- sqrt(n / (m * (n + m)))
  x1_mean = colMeans(x1_mat)
  x2_mean = colMeans(x2_mat)
  y_mean = colMeans(y_mat)
  
  Ssum_x <- (t(x1_mat) %*% x1_mat - n * tcrossprod(x1_mean))
  + (t(x2_mat) %*% x2_mat - n * tcrossprod(x2_mean))
  Ssum_y <- t(y_mat) %*% y_mat - m * tcrossprod(y_mean)
  s_within <- Ssum_x + Ssum_y
  
  
  if(q ==1){
    z_vec = n*x1_mean + m*y_mean
    xbar_1_qprime <- x2_mean[1]
  }else{
    z_vec = n*x2_mean + m*y_mean
    xbar_1_qprime <- x1_mean[1]
  }
  
  z_over <- z_vec / (n + m)
  
  delta_diff <- z_over - y_mean
  s_pooled <- s_within + ((m * (n + m)) / n) * tcrossprod(delta_diff)
  
  s11 <- s_pooled[1, 1]
  s12 <- s_pooled[1, 2]
  
  # Naive 1
  delta_naive1 <- z_over[2]
  
  # Naive 2 & 3 (Plugin Adjustments)
  arg1 <- sqrt((n * (n + m)) / m) * (z_over[1] - xbar_1_qprime) / sqrt(s_within[1, 1] / nu)
  arg2 <- sqrt((n * (n + m)) / m) * (z_over[1] - xbar_1_qprime) / sqrt(s11 / (nu + 1))
  delta_naive2 <- z_over[2] - factor_const * ((s_within[1, 2] * dnorm(arg1)) / (sqrt(s_within[1, 1]) * pnorm(arg1)))
  delta_naive3 <- z_over[2] - factor_const * ((s12 * dnorm(arg2)) / (sqrt(s11) * pnorm(arg2)))
  
  # UMVCUE (delta_U)
  u_star <- sqrt(n * (n + m) / m) * (z_over[1] - xbar_1_qprime)
  min_val <-min(1, u_star / sqrt(s11))
  
  pow_term <- (1 - min_val^2)^(nu / 2)
  c_par <- nu / 2
  F_arg <- (1 + min_val) / 2
  
  F_val <- pbeta(F_arg, c_par, c_par)
  denom <- (2^(nu - 1)) * nu * sqrt(s11) * beta(c_par, c_par) * F_val
  delta_U <- z_over[2] - (factor_const * s12 * pow_term / denom)
  return(c(delta_naive1, delta_naive2, delta_naive3, delta_U))
}


n=5; m=25
x1_mat <- cbind(X1[1:n],Y1[1:n])   # treatment arm 1
x2_mat <- cbind(X2[1:n],Y2[1:n])   # treatment arm 2

x1_mean <- colMeans(x1_mat)
x2_mean <- colMeans(x2_mat)


if (x1_mean[1] > x2_mean[1]) {
  q <- 1
  y_mat <- cbind(X1[(n+1):(n+m)],Y1[(n+1):(n+m)])
} else {
  q <- 2
  y_mat <- cbind(X2[(n+1):(n+m)],Y2[(n+1):(n+m)])
}
Estimates(x1_mat = x1_mat, x2_mat = x2_mat, q = q, y_mat = y_mat)