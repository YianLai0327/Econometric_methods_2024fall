getwd()
dt <- read.csv("./HW3/Equity_Premium.csv")
time <- dt$Time
y <- dt$y
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

x <- cbind(1, x_dfy, x_infl, x_svar, x_tms, x_tbl, x_dfr, x_dp, x_ltr)
X_transpose_X <- t(x) %*% x
beta_hat <- solve(X_transpose_X) %*% t(x) %*% y
beta_df <- data.frame(
  Estimate = beta_hat
)

cat("=====================================================\n")
cat("beta_hat\n")
print(beta_df)

e_hat <- y - x %*% beta_hat
s_square_hat <- t(e_hat) %*% e_hat / (nrow(x) - ncol(x))


s_beta_hat <- sqrt(diag(s_square_hat[1] * solve(X_transpose_X)))

s_beta_df <- data.frame(
  Estimate = s_beta_hat
)

cat("=====================================================\n")
cat("s_beta_hat:\n")
print(s_beta_df)

e_hat_square_xi_xi_transpose <- matrix(0, ncol(x), ncol(x))

for(i in 1:nrow(x)) {
    x_i <- x[i, ]
    e_hat_square_xi_xi_transpose <- e_hat_square_xi_xi_transpose + e_hat[i]^2 * x_i %*% t(x_i)
}

hetero_var_beta_hat <- solve(X_transpose_X) %*% e_hat_square_xi_xi_transpose %*% solve(X_transpose_X)
s_beta_hat_hetero <- sqrt(diag(hetero_var_beta_hat))

s_beta_hat_hetero_df <- data.frame(
  Estimate = s_beta_hat_hetero
)
cat("=====================================================\n")
cat("s_beta_hat_hetero:\n")
print(s_beta_hat_hetero_df)

t_stats <- beta_hat / s_beta_hat

t_stats_df <- data.frame(
  t_statistic = t_stats
)

cat("=====================================================\n")
cat("t_stats:\n")
print(t_stats_df)

cat("=====================================================\n")

n <- nrow(x)
p <- ncol(x)

alpha_1 <- qt(1 - 0.01 / 2, df = n - p)
alpha_5 <- qt(1 - 0.05 / 2, df = n - p)
alpha_10 <- qt(1 - 0.10 / 2, df = n - p)

reject_1 <- abs(t_stats) > alpha_1
reject_5 <- abs(t_stats) > alpha_5
reject_10 <- abs(t_stats) > alpha_10

cat("alpha rejection 1% level:", alpha_1, "\n")
cat("alpha rejection 5% level:", alpha_5, "\n")
cat("alpha rejection 10% level:", alpha_10, "\n")

reject_1_df <- data.frame(
  Reject_H0 = reject_1
)

reject_5_df <- data.frame(
  Reject_H0 = reject_5
)

reject_10_df <- data.frame(
  Reject_H0 = reject_10
)

cat("=====================================================\n")
cat("Reject at 1% level:\n")
print(reject_1_df)

cat("=====================================================\n")
cat("Reject at 5% level:\n")
print(reject_5_df)

cat("=====================================================\n")
cat("Reject at 10% level:\n")
print(reject_10_df)
cat("=====================================================\n")

sigma_square_hat <- (1/n) * sum(e_hat^2)

e_hat_standardized <- e_hat / sqrt(sigma_square_hat)  # Standardized residuals

skewness <- (1/n) * sum(e_hat_standardized^3)
kurtosis <- (1/n) * sum(e_hat_standardized^4)

# Compute the Jarque-Bera test statistic
JB_stat <- n * (skewness^2 / 6 + (kurtosis - 3)^2 / 24)

cat("Skewness:", skewness, "\n")
cat("Kurtosis:", kurtosis, "\n")
cat("Jarque-Bera statistic:", JB_stat, "\n\n")

alpha_1 <- qchisq(0.99, df = 2)
alpha_5 <- qchisq(0.95, df = 2)
alpha_10 <- qchisq(0.90, df = 2)

cat("Critical value at 1% level:", alpha_1, "\n")
cat("Critical value at 5% level:", alpha_5, "\n")
cat("Critical value at 10% level:", alpha_10, "\n")

reject_1 <- JB_stat > alpha_1
reject_5 <- JB_stat > alpha_5
reject_10 <- JB_stat > alpha_10

cat("\nReject at 1% level:", reject_1, "\n")
cat("Reject at 5% level:", reject_5, "\n")
cat("Reject at 10% level:", reject_10, "\n")
cat("=====================================================\n")

library(ggplot2)
density_resid <- density(e_hat_standardized)
x_vals <- seq(min(e_hat_standardized), max(e_hat_standardized), length.out = 100)
normal_pdf <- dnorm(x_vals)

# Save the plot
png("./HW5/Kernel_Density_vs_Normal.png", width = 800, height = 600)
plot(density_resid, main="Kernel Density vs Normal(0,1)", xlab="Standardized Residuals")
lines(x_vals, normal_pdf, col="red")
legend("topright", legend=c("Kernel Density", "N(0,1) PDF"), col=c("black", "red"), lty=1)
dev.off()