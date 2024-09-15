# Load necessary library
library(MASS)  # For generating random samples from t-distribution

# Define the function for Y = 1 / (1 + X^4)
Y_func <- function(X) {
  return(1 / (1 + X^4))
}

# Generate random samples for X from t-distribution with 3 degrees of freedom
set.seed(0)
num_samples <- 100000
X_samples <- rt(num_samples, df = 3)

# Calculate the corresponding Y values
Y_samples <- Y_func(X_samples)

# Estimate the necessary expectations
E_X <- mean(X_samples)
E_Y <- mean(Y_samples)
E_X2 <- mean(X_samples^2)
E_XY <- mean(X_samples * Y_samples)

# Calculate the linear regression coefficient b
b <- (E_XY - E_X * E_Y) / (E_X2 - (E_X)^2)

# Print the result
print(paste("The estimated value of b is:", b))