set.seed(12345)

n_values <- c(50, 100, 200, 500)
k <- 5
R <- matrix(c(1, 0, 0, 0, 0, 
              1, 1, 0, 0, 0), nrow = 2, byrow = TRUE)
theta_0 <- c(1, 2)
alpha <- 0.05
replications <- 1000

library(MASS)

generate_data <- function(n, k) {
    data <- mvrnorm(n, mu = rep(0, k + 1), Sigma = diag(k + 1))
    #the first k columns are the X's, the last column is the e
    X <- data[, 1:k]
    e <- data[, k + 1]
    list(X = X, e = e)
}

compute_test_statistic <- function(X, Y, R, theta_0) {
  beta_hat <- solve(t(X) %*% X) %*% t(X) %*% Y
  sigma_hat <- sum((Y - X %*% beta_hat)^2) / (nrow(X) - ncol(X))
  V_theta_hat <- sigma_hat * R %*% solve(t(X) %*% X) %*% t(R)

  T <- t(R %*% beta_hat - theta_0) %*% solve(V_theta_hat) %*% (R %*% beta_hat - theta_0)
  as.numeric(T)
}

critical_value <- qchisq(1 - alpha, df = 2)

simulation <- function(n, R, theta_0, beta, critical_value, replications) {
  rejections <- 0
  for (i in 1:replications) {
    data <- generate_data(n, k)
    X <- data$X
    e <- data$e
    Y <- X %*% beta + e
    T <- compute_test_statistic(X, Y, R, theta_0)
    if (T > critical_value) {
      rejections <- rejections + 1
    }
  }
  rejections / replications
}

# Run the simulation for each n to get the empirical size
beta_alt_2a <- rep(1, k)
size_results <- sapply(n_values, function(n) simulation(n, R, theta_0, beta_alt_2a, critical_value, replications))
results_2a <- data.frame(n = n_values, size = size_results)
cat("Empirical size for each n:\n")
print(results_2a)

# Run the simulation for each n with beta = (1, 2, 3, ..., k)
beta_alt_2b <- 1:k
power_results_2b <- sapply(n_values, function(n) simulation(n, R, theta_0, beta_alt_2b, critical_value, replications))
cat("Empirical power for each n with beta = (1, 2, 3, ..., k):\n")
print(data.frame(n = n_values, power = power_results_2b))

# Range of h values for Problem 2c
h_values <- 1:10

# Run the simulation for each n and each h in 2c
power_results_2c <- sapply(n_values, function(n) {
  sapply(h_values, function(h) {
    beta <- rep(1 + h / sqrt(n), k)  # Construct beta for each h and n
    simulation(n, R, theta_0, beta, critical_value, replications)
  })
})

power_results_2c <- t(power_results_2c)

# Convert power results to a more readable format
colnames(power_results_2c) <- paste("h =", h_values)
rownames(power_results_2c) <- paste("n =", n_values)
cat("Empirical power for each n and h with beta = (1 + n^(-1/2) * h, 2, 3, ..., k):\n")
print(data.frame(power_results_2c))

