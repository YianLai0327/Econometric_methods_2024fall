#simulate the plug-in estimator for the mean of IID random variables

# Load necessary library

library(MASS)  # For generating random samples from t-distribution
library(ggplot2)  # For plotting
#simulate the finite sample distribution of the statics: sample mean
#test for n = 10, 50, 100
# Experiment parameters
r <- 1000
n_values <- c(10, 50, 100, 10000) #n=10000 is a extra simulation for the large sample size
mode <- c("y_bar", "y_bar_scaled")
random_variable_mode <- c("t", "n")

for (m in mode) {
  for (rv in random_variable_mode) {
    for (n in n_values) {
      # Initialize vectors for storing the sample means
      Y_bar <- rep(0, r)
      
      for (i in 1:r) {
        if (rv == "t") {
          x <- rt(n, df = 2)
        } else {
          x <- rnorm(n)
        }
        Y_bar[i] <- mean(x)
      }
      
      # Calculate scaled sample means
      if (m == "y_bar_scaled") {
        Y_bar <- Y_bar * sqrt(n)
      }
      # Plot for Y_bar_scaled
      png(paste0("HW6/outputs/", n, "_", rv, "_", m, ".png"))

      # Set the range for x-axis
      x_range <- seq(-6, 6, length.out = 300)

      # Calculate the density for Y_bar_sqrt_n within the range [-6, 6]
      density_Y_bar <- density(Y_bar, from = -6, to = 6)
      
      # Calculate the normal distribution over the same range
      normal_pdf <- dnorm(x_range)
      y_max <- max(max(normal_pdf), max(density_Y_bar$y))

      # Plot the density of Y_bar_sqrt_n
      plot(density_Y_bar$x, density_Y_bar$y, type = "l", col = "red", lwd = 2, xlab = "x", ylab = "Density", 
      main = paste("Density of ", m, "with distribution=", rv, "n =", n), xlim = c(-6, 6), ylim = c(0, y_max*1.2))

      # Overlay the normal distribution
      lines(x_range, normal_pdf, col = "blue", lwd = 2)

      # Add legend
      legend("topright", legend = c("Y_bar * sqrt(n)", "N(0, 1)"), col = c("red", "blue"), lwd = 2)

      dev.off()
    }
  }
}

dt <- read.csv("./HW6/Equity_Premium.csv")
time <- dt$Time
y <- dt$y
ones <- rep(1, length(y))
x_dfy <- dt$x_dfy
x_infl <- dt$x_infl
x_svar <- dt$x_svar
x_tms <- dt$x_tms
x_tbl <- dt$x_tbl
x_dfr <- dt$x_dfr
x_dp <- dt$x_dp
x_ltr <- dt$x_ltr
x_ep <- dt$x_ep
x_bmr <- dt$x_bmr
x_ntis <- dt$x_ntis

x <- cbind(ones, x_dfy, x_infl, x_svar, x_tms, x_tbl, x_dfr, x_dp, x_ltr, x_ep, x_bmr, x_ntis)

#Do the Wald's test for each beta_hat_j with null hypothesis beta_hat_j = 0

fit <- lm(y ~ x - 1)
beta_hat <- coef(fit)
vcov_matrix <- vcov(fit)

wald_test_results <- data.frame(
  beta = beta_hat,
  wald_statistic = beta_hat^2 / diag(vcov_matrix),
  p_value = 1 - pchisq(beta_hat^2 / diag(vcov_matrix), df = 1)
)

#reject the null hypothesis if p_value < alpha
alpha <- 0.05
wald_test_results$reject_null <- wald_test_results$p_value < alpha

print(wald_test_results)


cat("\n=====================================================\n")
#conducted joint wald's test with null hypothesis beta_hat_1 = 0 and beta_hat_2 + beta_hat_3 = 0
# Calculate the Wald statistic
R <- matrix(c(1, 0, 0, 0, 1, 1), nrow = 2, ncol = 3)
r <- c(0, 0)
R_beta <- R %*% beta_hat[c(1, 2, 3)]
R_vcov <- R %*% vcov_matrix[c(1, 2, 3), c(1, 2, 3)] %*% t(R)
wald_statistic <- t(R_beta) %*% solve(R_vcov) %*% (R_beta - r)

p_value = 1 - pchisq(wald_statistic, df = 2)

cat(paste("Wald statistic:", wald_statistic, "\n"))
cat(paste("P-value:", p_value, "\n"))
cat(paste("Reject null hypothesis if p_value < alpha: ", p_value < alpha, "\n"))
cat("=====================================================\n")