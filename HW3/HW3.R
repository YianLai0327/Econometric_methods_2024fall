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

x <- cbind(1, x_dfy, x_infl, x_svar, x_tms, x_tbl, x_dfr)
#Y = X'beta + epsilon
#please compute LS estimater of beta by using the formula beta_hat = (X'X)^-1X'y

# Calculate X'X
X_transpose_X <- t(x) %*% x

beta_hat <- solve(X_transpose_X) %*% t(x) %*% y
print("beta_hat")
print(beta_hat)

#please compute the LS estimator of beta by FWL theorem
x1 <- cbind(1, x_dfy, x_infl, x_svar)
x2 <- cbind(x_tms, x_tbl, x_dfr)

m1 <- diag(nrow(x1)) - x1 %*% solve(t(x1) %*% x1) %*% t(x1)
m2 <- diag(nrow(x2)) - x2 %*% solve(t(x2) %*% x2) %*% t(x2)

beta_hat1 <- solve(t(m2 %*% x1) %*% m2 %*% x1) %*% t(m2 %*% x1) %*% m2 %*% y
beta_hat2 <- solve(t(m1 %*% x2) %*% m1 %*% x2) %*% t(m1 %*% x2) %*% m1 %*% y

beta_hat_fwl <- matrix(t(c(beta_hat1, beta_hat2)))
print("beta_hat_fwl")
print(beta_hat_fwl)

#compute R square of a set of variables
compute_R2 <- function(x, y) {
  # Calculate beta_hat
  beta_hat <- solve(t(x) %*% x) %*% t(x) %*% y
  
  # Calculate y_hat
  y_hat <- x %*% beta_hat
  
  # Calculate R^2
  R2 <- 1 - sum((y - y_hat)^2) / sum((y - mean(y))^2)
  
  return(R2)
}

rsquare_list <- c()

for(i in 1:7) {
    x_subset <- x[, 1:i]
    R2 <- compute_R2(x_subset, y)
    rsquare_list <- c(rsquare_list, R2)
}

png("./HW3/plot.png")
plot(1:7, rsquare_list, type = "o", xlab = "j", ylab = "R^2", main = "R^2 vs. set of j")
dev.off()
#save the plot
