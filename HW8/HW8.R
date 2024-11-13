# Load necessary libraries
library(combinat)  # To generate all combinations of predictors

# Load data
dt <- read.csv("Equity_Premium.csv")
y <- dt$y
ones <- rep(1, length(y))
x_dfy <- dt$x_dfy
x_infl <- dt$x_infl
x_svar <- dt$x_svar
x_tms <- dt$x_tms
x_tbl <- dt$x_tbl

x_dfy_squared <- x_dfy^2
x_infl_squared <- x_infl^2
x_svar_squared <- x_svar^2
x_tms_squared <- x_tms^2
x_tbl_squared <- x_tbl^2

# Create data frame with all variables
data <- data.frame(y = y, ones = ones, x_dfy = x_dfy, x_infl = x_infl,
                   x_svar = x_svar, x_tms = x_tms, x_tbl = x_tbl,
                   x_dfy_squared = x_dfy_squared,
                   x_infl_squared = x_infl_squared,
                   x_svar_squared = x_svar_squared,
                   x_tms_squared = x_tms_squared,
                   x_tbl_squared = x_tbl_squared)

# Define the function to collect criteria
collect_criterian <- function(y, x) {
    n <- length(y)
    k <- ncol(x)
    x_matrix <- as.matrix(x)
    beta_hat <- solve(t(x_matrix) %*% x_matrix) %*% t(x_matrix) %*% y
    y_hat <- x_matrix %*% beta_hat
    residuals <- y - y_hat
    rss <- sum(residuals^2)
    tss <- sum((y - mean(y))^2)
    sigma_hat <- rss / (n - k)
    r2 <- 1 - rss / tss
    adjusted_r2 <- 1 - ((1 - r2) * (n - 1) / (n - k))
    aic <- n * log(rss / n) + 2 * k
    bic <- n * log(rss / n) + k * log(n)
    cp <- rss + 2 * k * sigma_hat
    h_ii <- rowSums((x_matrix %*% solve(t(x_matrix) %*% x_matrix)) * x_matrix)
    loo_cv <- mean((residuals / (1 - h_ii))^2)
    
    model_results <- list(r2 = r2, adjusted_r2 = adjusted_r2, aic = aic,
                          bic = bic, cp = cp, loo_cv = loo_cv)
    return(model_results)
}

# Get all combinations of predictors
predictor <- c("ones", "x_dfy", "x_infl", "x_svar", "x_tms", "x_tbl",
               "x_dfy_squared", "x_infl_squared", "x_svar_squared",
               "x_tms_squared", "x_tbl_squared")

predictor_combinations <- unlist(
    lapply(1:length(predictor), function(i) combn(predictor, i, simplify = FALSE)),
    recursive = FALSE)

# Initialize lists
models <- vector("list", length(predictor_combinations))
model_predictors <- vector("list", length(predictor_combinations))

# Build models and collect criteria
for (i in seq_along(predictor_combinations)) {
    predictors <- predictor_combinations[[i]]
    x <- data[, predictors]
    result <- collect_criterian(y, x)
    if (is.null(result)) {
        # Skip singular models
        next
    }
    models[[i]] <- result
    model_predictors[[i]] <- predictors
}

# Remove NULL entries from models and model_predictors
valid_indices <- which(!sapply(models, is.null))
models <- models[valid_indices]
model_predictors <- model_predictors[valid_indices]

# Extract criteria values
r2_values <- sapply(models, function(model) model$r2)
adjusted_r2_values <- sapply(models, function(model) model$adjusted_r2)
aic_values <- sapply(models, function(model) model$aic)
bic_values <- sapply(models, function(model) model$bic)
cp_values <- sapply(models, function(model) model$cp)
loo_cv_values <- sapply(models, function(model) model$loo_cv)

#Plot the selection model statistics of all models by 折線圖
#define the output figure's height and width
png("HW8/selection_model_statistics.png", width = 900, height = 600)
par(mfrow = c(3, 2))
plot(1:length(model_predictors), r2_values, type = "l", xlab = "Model Index", ylab = "R^2", main = "R^2")
plot(1:length(model_predictors), adjusted_r2_values, type = "l", xlab = "Model Index", ylab = "Adjusted R^2", main = "Adjusted R^2")
plot(1:length(model_predictors), aic_values, type = "l", xlab = "Model Index", ylab = "AIC", main = "AIC")
plot(1:length(model_predictors), bic_values, type = "l", xlab = "Model Index", ylab = "BIC", main = "BIC")
plot(1:length(model_predictors), cp_values, type = "l", xlab = "Model Index", ylab = "Cp", main = "Cp")
plot(1:length(model_predictors), loo_cv_values, type = "l", xlab = "Model Index", ylab = "LOO CV", main = "LOO CV")
dev.off()


# Find best models
best_indices <- list(
    r2 = which.max(r2_values),
    adjusted_r2 = which.max(adjusted_r2_values),
    aic = which.min(aic_values),
    bic = which.min(bic_values),
    cp = which.min(cp_values),
    loo_cv = which.min(loo_cv_values)
)

# Print best models and predictors
for (criterion in names(best_indices)) {
    index <- best_indices[[criterion]]
    cat(sprintf("\nBest Model Based on %s:\n", toupper(criterion)))
    cat("Predictors:\n")
    print(model_predictors[[index]])
    # Print the reference criterion value
    cat(sprintf("%s: %.4f\n", criterion, models[[index]][[criterion]]))

}

