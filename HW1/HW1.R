dt <- read.csv("Equity_Premium.csv")
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

#Q1. Plot the time series and the histograms (with the estimated normal density functions) of y and the x's.
library(ggplot2)
library(tidyr)
library(purrr)

# Function to create time series plot
# plot_time_series <- function(data, var_name) {
#   ggplot(data, aes(x = Time, y = .data[[var_name]])) +
#     geom_line() +
#     labs(title = paste("Time Series of", var_name), y = var_name)
# }

# # Function to create histogram with normal density
# plot_histogram <- function(data, var_name) {
#   ggplot(data, aes(x = .data[[var_name]])) +
#     geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
#     stat_function(fun = dnorm, args = list(mean = mean(data[[var_name]]), sd = sd(data[[var_name]])), color = "red") +
#     labs(title = paste("Histogram of", var_name), x = var_name)
# }

# # List of variables to plot
# vars <- c("y", "x_dfy", "x_infl", "x_svar", "x_tms", "x_tbl", "x_dfr", "x_dp", "x_ltr", "x_ep", "x_bmr", "x_ntis")

# # Create and save plots
# walk(vars, ~{
#   time_series_plot <- plot_time_series(dt, .x)
#   histogram_plot <- plot_histogram(dt, .x)
  
#   ggsave(paste0("./plots/time_series_", .x, ".png"), time_series_plot, width = 10, height = 6)
#   ggsave(paste0("./plots/histogram_", .x, ".png"), histogram_plot, width = 10, height = 6)
# })

#Q2. Given the same size n and the X defined in this running example, show: 
#Q2.1. trace(X(X'X)^-1X') = ?
X <- cbind(x_dfy, x_infl, x_svar, x_tms, x_tbl, x_dfr, x_dp, x_ltr, x_ep, x_bmr, x_ntis)

# Calculate X'X
X_transpose_X <- t(X) %*% X

# Calculate (X'X)^-1
X_transpose_X_inv <- solve(X_transpose_X)

# Calculate X(X'X)^-1X'
hat_matrix <- X %*% X_transpose_X_inv %*% t(X)

# Calculate trace(X(X'X)^-1X')
trace_result_Q2_1 <- sum(diag(hat_matrix))
cat("Q2-1: The result of trace(X(X'X)^-1X') is:", trace_result_Q2_1, "\n")

#Q2.2. trace(I_n - X(X'X)^-1X') = ?
identity_matrix <- diag(nrow(X))
trace_result_Q2_2 <- sum(diag(identity_matrix - hat_matrix))
cat("Q2-2: The result of trace(I_n - X(X'X)^-1X') is:", trace_result_Q2_2, "\n")

#Q3. Following #2, let λ_j be an eigenvalue of X'X, for j = 1, 2, ..., k. Show the "scree plot" of the eigenvalues. (The horizontal axis is j, and the vertical axis is λ_j).
eigenvalues <- eigen(X_transpose_X)$values

# Create a data frame for plotting
scree_data <- data.frame(
  j = 1:length(eigenvalues),
  lambda = eigenvalues
)

# Create the scree plot
scree_plot <- ggplot(scree_data, aes(x = j, y = lambda)) +
  geom_line() +
  geom_point() +
  geom_text(aes(label = round(lambda, 2)), vjust = -0.5, hjust = 0.5) +
  labs(title = "Scree Plot of Eigenvalues",
       x = "j",
       y = "λ_j (Eigenvalue)") +
  scale_x_continuous(breaks = 1:11) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  )

# Save the plot
ggsave("scree_plots/scree_plot_q3.png", scree_plot, width = 12, height = 8, dpi = 300)

#Q4. Following #3, let X̄ be an n × k matrix. The j-th column of X̄ is defined by "standardizing" the j-th column of X, with the sample mean 0 and the sample variance 1, for each j. Compare the scree plot of the eigenvalues of X̄'X̄ with that of X'X.

# Standardize each column of X
X_standardized <- scale(X)

# Calculate X̄'X̄
X_bar_transpose_X_bar <- t(X_standardized) %*% X_standardized

# Calculate eigenvalues of X̄'X̄
eigenvalues_X_bar <- eigen(X_bar_transpose_X_bar)$values

# Create a data frame for plotting
scree_data_X_bar <- data.frame(
  j = 1:length(eigenvalues_X_bar),
  lambda = eigenvalues_X_bar
)

# Create the scree plot for X̄'X̄
scree_plot_X_bar <- ggplot(scree_data_X_bar, aes(x = j, y = lambda)) +   
  geom_line() +
  geom_point() +
  geom_text(aes(label = round(lambda, 2)), vjust = -0.5, hjust = 0.5) +
  labs(title = "Scree Plot of Eigenvalues",
       x = "j",
       y = "λ_j (Eigenvalue)") +
  scale_x_continuous(breaks = 1:11) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  )

# Save the plot
#ggsave("scree_plots/scree_plot_q4.png", scree_plot_X_bar, width = 12, height = 8, dpi = 300)

#Q5. Following #4, compute A = (X̄'X̄)⁻¹ using the spectral decomposition and verify that AA⁻¹ = Iₖ.

# Spectral decomposition of X̄'X̄
eigen_decomp <- eigen(X_bar_transpose_X_bar)
eigenvalues <- eigen_decomp$values
eigenvectors <- eigen_decomp$vectors

# Compute A using spectral decomposition
A <- eigenvectors %*% diag(1/eigenvalues) %*% t(eigenvectors)

# Verify that AA⁻¹ = Iₖ
Verify_matrix <- A %*% solve(A)
Verify_matrix <- as.matrix(Verify_matrix)
colnames(Verify_matrix) <- NULL
rownames(Verify_matrix) <- NULL

cat("Q5. Verify that AA⁻¹ equals to Iₖ by calculating the maximum element-wise difference:\n")
cat("Maximum element-wise difference:", max_diff)
is_effectively_identity <- max_diff < 1e-12
if (!is_effectively_identity) {
  cat(".\n=> The result is not close to the identity matrix within the specified tolerance 10^-12.\n")
} else {
  cat(".\n=> The result is effectively the identity matrix within the specified tolerance 10^-12.\n")
}

#Q6. Following #5, let y be an n × 1 vector. The i-th element of y corresponds to the i-th observation of y, for i = 1, 2, ..., n.
#Consider the linear equation:y = X̄b̂
#where b̂ is a k × 1 vector. Show b̂ = ?

#Estimate b by projection matrix
b <- A %*% t(X_standardized) %*% y
#Print the result
print("Q6: Result of b is:")
print(b)