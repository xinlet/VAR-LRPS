library(MASS)
library(vars)
library(glmnet)
library(parallel)

var_fit <- function(ts, lag) { 
  if (is.null(colnames(ts))) {
    colnames(ts) <- paste0("Y", seq_len(ncol(ts)))
  }
  fitvar <- VAR(ts, p = lag, type = "none")  
  Bhat <- Bcoef(fitvar)  
  return(list(coef = Bhat, ts = ts))  
}

LRfun <- function(y_train, k = ncol(y_train)) {
  pc <- prcomp(y_train, center = FALSE, scale. = FALSE)
  y_train <- t(t(pc$x[, 1:k] %*% t(pc$rotation[, 1:k])))
  return(y_train)
}

LRPS_varfit <- function(ts, lag, k){
  n <- nrow(ts)
  Y <- ts[-(1:lag), ]  
  X <- ts[1:(n-lag),]
  Ys <- LRfun(Y, k)
  coef <- solve(t(X) %*% X) %*% t(X) %*% Ys
  return(list(coef = coef, Ys = Ys, ts = ts))
}

svar_fit <- function(ts, lag, lambda){
  n <- nrow(ts)
  Y <- ts[-(1:lag), ]  
  X <- ts[1:(n-lag),]
  p <- ncol(Y)
  B_sparse <- matrix(0, ncol(X), p)
  
  for (i in 1:p) {
    fit <- glmnet(X, Y[, i], alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE)
    B_sparse[, i] <- as.vector(coef(fit))[-1]
  }
  return(list(coef = B_sparse, ts = ts))  
}

compute_mse <- function(B_true, B_pred) {
  return(mean((B_true - B_pred)^2))
}

cv_time_series <- function(ts, lag, k_values, init_train_size, model_type = c("LRPS", "VAR")) {
  model_type <- match.arg(model_type)
  
  ts <- as.matrix(ts)
  ts <- apply(ts, 2, as.numeric)
  
  n <- nrow(ts)
  p <- ncol(ts)
  
  num_steps <- n - init_train_size - 1
  mse_k <- matrix(NA, nrow = num_steps, ncol = length(k_values))
  
  for (t in (init_train_size + 1):(n - 1)) {
    X_train <- ts[1:t, , drop = FALSE]
    Y_train <- ts[-(1:lag), ][1:t, , drop = FALSE]
    
    X_test <- ts[t + 1, , drop = FALSE]
    Y_test <- ts[t + 1, , drop = FALSE]
    
    if (model_type == "VAR") {
      model <- var_fit(X_train, lag)
      B_hat <- model$coef
      Y_pred <- X_test %*% B_hat
      mse_k[t - init_train_size, ] <- rep(mean((Y_test - Y_pred)^2), length(k_values))
    } else {
      for (j in seq_along(k_values)) {
        k <- k_values[j]
        model <- LRPS_varfit(X_train, lag, k)
        B_hat <- model$coef
        Y_pred <- X_test %*% B_hat
        mse_k[t - init_train_size, j] <- mean((Y_test - Y_pred)^2)
      }
    }
  }
  
  # For VAR: no k selection, just one MSPE series
  if (model_type == "VAR") {
    mspe <- rowMeans(mse_k, na.rm = TRUE)
    return(list(best_k = NA, mspe = mspe))
  }
  
  # For LRPS: average MSPE per k, then choose best
  avg_mse_per_k <- colMeans(mse_k, na.rm = TRUE)
  best_k <- k_values[which.min(avg_mse_per_k)]
  mspe <- mse_k[, which.min(avg_mse_per_k)]
  
  return(list(best_k = best_k, mspe = mspe, mse_k_matrix = mse_k, avg_mse_per_k = avg_mse_per_k))
}

library(data.table)
determine_best_k_parallel <- function(train_files, lag, k_values, init_train_size, n_cores = detectCores() - 1) {
  # Each file returns the best k (mode) from its rolling CV
  best_k_lrps_list <- mclapply(train_files, function(file) {
    ts_data <- fread(file)
    ts_data <- scale(as.matrix(ts_data))
    
    result <- cv_time_series(ts_data, lag, k_values, init_train_size, model_type = "LRPS")
    best_k <- result$best_k  # scalar (most frequent for that file)
    return(best_k)
  }, mc.cores = n_cores)
  
  all_k_lrps <- unlist(best_k_lrps_list)
  final_k_lrps <- as.numeric(names(which.max(table(all_k_lrps))))
  
  return(list(k_lrps = final_k_lrps, all_k_lrps = all_k_lrps))
}

data_dir <- "your directory"  
files <- list.files(data_dir, full.names = TRUE, pattern = "*.csv")

train_files <- files[1:10]  
lag <- 1
k_values <- 1:10  # Possible k values
init_train_size <- 150  # Initial training size

# Determine best k values
best_k_results <- determine_best_k_parallel(train_files, lag, k_values, init_train_size)
print(best_k_results)

library(data.table)
library(glmnet)
library(zoo)


# Function to fit sparse VAR(1) using Lasso (column-wise regression)
sparse_var_fit <- function(X, Y, lambda) {
  p <- ncol(Y)
  B_sparse <- matrix(0, ncol(X), p)
  
  for (j in 1:p) {
    fit <- glmnet(X, Y[, j], alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE)
    B_sparse[, j] <- as.vector(coef(fit))[-1]  # exclude intercept
  }
  return(B_sparse)
}

# Cross-validation to choose lambda
cv_sparse_var <- function(ts_data, lag, lambda_seq, init_train_size) {
  n <- nrow(ts_data)
  p <- ncol(ts_data)
  mspe_matrix <- matrix(NA, nrow = length(lambda_seq), ncol = n - init_train_size - 1)
  
  for (l in seq_along(lambda_seq)) {
    lambda <- lambda_seq[l]
    
    for (t in (init_train_size + 1):(n - 1)) {
      X_train <- ts_data[1:t, , drop = FALSE]
      Y_train <- ts_data[-(1:lag), ][1:t, , drop = FALSE]
      
      X_test <- ts_data[t + 1, , drop = FALSE]
      Y_test <- ts_data[t + 1, , drop = FALSE]
      
      B_hat <- sparse_var_fit(X_train, Y_train, lambda)
      Y_pred <- X_test %*% B_hat
      mspe_matrix[l, t - init_train_size] <- mean((Y_test - Y_pred)^2)
    }
  }
  
  avg_mspe <- rowMeans(mspe_matrix, na.rm = TRUE)
  best_lambda <- lambda_seq[which.min(avg_mspe)]
  return(list(best_lambda = best_lambda, mspe_path = avg_mspe))
}


lambda_seq <- 10^seq(-4, 0, length.out = 20)  # e.g., range of lambda values

files <- list.files(data_dir, full.names = TRUE, pattern = "*.csv")

train_files <- files[1:10]


init_train_size <- 150
lag <- 1

all_best_lambda <- c()

for (file in train_files) {
  ts_data <- fread(file) |> as.matrix() |> scale()
  if (any(is.na(ts_data))) ts_data <- na.aggregate(ts_data)
  lambda_cv_result <- cv_sparse_var(ts_data, lag, lambda_seq, init_train_size)
  all_best_lambda <- c(all_best_lambda, lambda_cv_result$best_lambda)
}

# Choose the most frequently selected lambda (or mean/median)
final_lambda <- median(all_best_lambda)



#parallel version of cv lambda====
library(data.table)
library(glmnet)
library(zoo)
library(doParallel)
library(foreach)

# Set up parallel backend
n_cores <- parallel::detectCores() - 1  # use all but one core
cl <- makeCluster(n_cores)
registerDoParallel(cl)

lambda_seq <- 10^seq(-4, 0, length.out = 20)
files <- list.files(data_dir, full.names = TRUE, pattern = "*.csv")
train_files <- files[1:10]
init_train_size <- 150
lag <- 1

sparse_var_fit <- function(X, Y, lambda) {
  p <- ncol(Y)
  B_sparse <- matrix(0, ncol(X), p)
  
  for (j in 1:p) {
    fit <- glmnet(X, Y[, j], alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE)
    B_sparse[, j] <- as.vector(coef(fit))[-1]
  }
  return(B_sparse)
}

cv_sparse_var <- function(ts_data, lag, lambda_seq, init_train_size) {
  n <- nrow(ts_data)
  p <- ncol(ts_data)
  mspe_matrix <- matrix(NA, nrow = length(lambda_seq), ncol = n - init_train_size - 1)
  
  for (l in seq_along(lambda_seq)) {
    lambda <- lambda_seq[l]
    
    for (t in (init_train_size + 1):(n - 1)) {
      X_train <- ts_data[1:t, , drop = FALSE]
      Y_train <- ts_data[-(1:lag), ][1:t, , drop = FALSE]
      
      X_test <- ts_data[t + 1, , drop = FALSE]
      Y_test <- ts_data[t + 1, , drop = FALSE]
      
      B_hat <- sparse_var_fit(X_train, Y_train, lambda)
      Y_pred <- X_test %*% B_hat
      mspe_matrix[l, t - init_train_size] <- mean((Y_test - Y_pred)^2)
    }
  }
  
  avg_mspe <- rowMeans(mspe_matrix, na.rm = TRUE)
  best_lambda <- lambda_seq[which.min(avg_mspe)]
  return(list(best_lambda = best_lambda, mspe_path = avg_mspe))
}

# Run in parallel over files
all_best_lambda <- foreach(file = train_files, .combine = c, .packages = c("data.table", "glmnet", "zoo")) %dopar% {
  ts_data <- fread(file) |> as.matrix() |> scale()
  if (any(is.na(ts_data))) ts_data <- na.aggregate(ts_data)
  lambda_cv_result <- cv_sparse_var(ts_data, lag, lambda_seq, init_train_size)
  lambda_cv_result$best_lambda
}

# Stop the cluster
stopCluster(cl)

# Summary statistic for selected lambdas
final_lambda <- median(all_best_lambda)
final_lambda







#fmri====
test_files <- setdiff(files, train_files)

save_matrix <- function(mat, file, out_dir, method_name) {
  # Compose output file path
  out_file <- file.path(out_dir, paste0(method_name, "_", basename(file)))
  write.csv(mat, file = out_file, row.names = FALSE)
}

best_k <- 1 #10 for axcpt, 10 for cuedts
lambda <- 1
lag <- 1
output_dir <- "your directory"
dir.create(output_dir, showWarnings = FALSE)

library(data.table)


for (test_file in test_files) {
  cat("Processing:", basename(test_file), "\n")
  ts_data <- fread(test_file)
  ts_data <- as.matrix(ts_data)
  ts_data <- apply(ts_data, 2, as.numeric)
  ts_data <- scale(ts_data)
  ts_data <- na.omit(ts_data)  # Remove any rows with NAs
  
  if (nrow(ts_data) <= lag) {
    warning(paste("Skipping:", basename(test_file), "- Not enough data after NA removal"))
    next
  }
  
  # Fit VAR
  var_result <- var_fit(ts_data, lag)
  save_matrix(var_result$coef, test_file, output_dir, "VAR1")
  
  # Fit LRPS (using best k)
  lrps_result <- LRPS_varfit(ts_data, lag, best_k)
  save_matrix(lrps_result$coef, test_file, output_dir, "LRPS")
  
  # Fit SVAR (sparse VAR)
  svar_result <- svar_fit(ts_data, lag, lambda)
  save_matrix(svar_result$coef, test_file, output_dir, "SVAR")
}
#===================================
#Until here, the code is used for generating VAR matrix for each method
#===================================

#===================================
# The following code is used for generating predictive performance for each method
#===================================

#use more metric====
library(data.table)
library(parallel)
library(zoo)  # for na.aggregate

evaluate_prediction_fixed_k <- function(test_files, lag, k_lrps, init_train_size) {
  lambda_svar <- 0.0546
  
  process_file <- function(file) {
    ts_data <- fread(file)
    ts_data <- as.matrix(ts_data)
    ts_data <- apply(ts_data, 2, as.numeric)
    ts_data <- scale(ts_data)
    
    if (any(is.na(ts_data))) {
      message("Handling missing values in file: ", file)
      ts_data <- na.aggregate(ts_data, FUN = mean)
    }
    
    n <- nrow(ts_data)
    rmse_var <- c(); mae_var <- c()
    rmse_lrps <- c(); mae_lrps <- c()
    rmse_svar <- c(); mae_svar <- c()
    
    for (t in (init_train_size + 1):(n - 1)) {
      X_train <- ts_data[1:t, , drop = FALSE]
      Y_train <- ts_data[-(1:lag), ][1:t, , drop = FALSE]
      
      X_test <- ts_data[t + 1, , drop = FALSE]
      Y_test <- ts_data[t + 1, , drop = FALSE]
      
      # VAR
      var_model <- var_fit(X_train, lag)
      if (!is.null(var_model)) {
        B_var <- var_model$coef
        Y_pred_var <- X_test %*% B_var
        rmse_var <- c(rmse_var, sqrt(mean((Y_test - Y_pred_var)^2)))
        mae_var <- c(mae_var, mean(abs(Y_test - Y_pred_var)))
      }
      
      # LRPS
      lrps_model <- LRPS_varfit(X_train, lag, k_lrps)
      if (!is.null(lrps_model)) {
        B_lrps <- lrps_model$coef
        Y_pred_lrps <- X_test %*% B_lrps
        rmse_lrps <- c(rmse_lrps, sqrt(mean((Y_test - Y_pred_lrps)^2)))
        mae_lrps <- c(mae_lrps, mean(abs(Y_test - Y_pred_lrps)))
      }
      
      # Sparse VAR
      svar_model <- svar_fit(X_train, lag, lambda = lambda_svar)
      if (!is.null(svar_model)) {
        B_svar <- svar_model$coef
        Y_pred_svar <- X_test %*% B_svar
        rmse_svar <- c(rmse_svar, sqrt(mean((Y_test - Y_pred_svar)^2)))
        mae_svar <- c(mae_svar, mean(abs(Y_test - Y_pred_svar)))
      }
    }
    
    return(list(
      rmse_var = rmse_var, mae_var = mae_var,
      rmse_lrps = rmse_lrps, mae_lrps = mae_lrps,
      rmse_svar = rmse_svar, mae_svar = mae_svar
    ))
  }
  
  # Start timing
  start_time <- Sys.time()
  
  # Parallel processing
  results <- mclapply(test_files, process_file, mc.cores = detectCores() - 1)
  
  # End timing
  end_time <- Sys.time()
  total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Combine results
  extract_metric <- function(res_list, name) {
    unlist(lapply(res_list, function(x) x[[name]]))
  }
  
  rmse_var_all <- extract_metric(results, "rmse_var")
  mae_var_all <- extract_metric(results, "mae_var")
  rmse_lrps_all <- extract_metric(results, "rmse_lrps")
  mae_lrps_all <- extract_metric(results, "mae_lrps")
  rmse_svar_all <- extract_metric(results, "rmse_svar")
  mae_svar_all <- extract_metric(results, "mae_svar")
  
  return(list(
    avg_rmse_var = mean(rmse_var_all, na.rm = TRUE),
    avg_mae_var = mean(mae_var_all, na.rm = TRUE),
    sd_rmse_var = sd(rmse_var_all, na.rm = TRUE),
    sd_mae_var = sd(mae_var_all, na.rm = TRUE),
    
    avg_rmse_lrps = mean(rmse_lrps_all, na.rm = TRUE),
    avg_mae_lrps = mean(mae_lrps_all, na.rm = TRUE),
    sd_rmse_lrps = sd(rmse_lrps_all, na.rm = TRUE),
    sd_mae_lrps = sd(mae_lrps_all, na.rm = TRUE),
    
    avg_rmse_svar = mean(rmse_svar_all, na.rm = TRUE),
    avg_mae_svar = mean(mae_svar_all, na.rm = TRUE),
    sd_rmse_svar = sd(rmse_svar_all, na.rm = TRUE),
    sd_mae_svar = sd(mae_svar_all, na.rm = TRUE),
    
    total_time_seconds = total_time
  ))
}

run_fixed_k_prediction <- function(data_dir, lag, k_lrps, init_train_size = 600) {
  files <- list.files(data_dir, full.names = TRUE, pattern = "*.csv")
  
  train_files <- files[1:10]
  test_files <- setdiff(files, train_files)
  
  # Perform rolling prediction on test files with best k values
  results <- evaluate_prediction_fixed_k(test_files, lag, k_lrps, init_train_size)
  
  return(results)
}
k_lrps <- 1  # Use pre-determined best k for LRPS
lambda_svar = 0.054
init_train_size <- 150
prediction_results <- run_fixed_k_prediction(data_dir, lag, k_lrps, init_train_size)


# Print results
print(prediction_results)




#==================================
#The following code is used for generating CV plot for each method
#==================================
                  
#mspe vs k during cv====
cv_time_series <- function(ts, lag, k_values, init_train_size, model_type = c("LRPS", "RRR")) {
  model_type <- match.arg(model_type)
  
  ts <- as.matrix(ts)  # Convert dataframe to numeric matrix
  ts <- apply(ts, 2, as.numeric)  # Ensure all columns are numeric
  
  n <- nrow(ts)
  p <- ncol(ts)
  
  mse_k <- matrix(NA, nrow = n - init_train_size - 1, ncol = length(k_values))  # Store MSE for each k
  best_k_choices <- rep(NA, n - init_train_size - 1)  # Track best k choices
  
  for (t in (init_train_size + 1):(n - 1)) {
    X_train <- ts[1:t, , drop = FALSE]
    Y_train <- ts[-(1:lag), ][1:t, , drop = FALSE]
    
    X_test <- ts[t + 1, , drop = FALSE]
    Y_test <- ts[t + 1, , drop = FALSE]
    
    for (j in seq_along(k_values)) {
      k <- k_values[j]
      if (model_type == "LRPS") {
        model <- LRPS_varfit(X_train, lag, k)
      } else if (model_type == "RRR") {
        model <- RRR_varfit(X_train, lag, k)
      }
      
      B_hat <- model$coef
      Y_pred <- X_test %*% B_hat
      mse_k[t - init_train_size, j] <- mean((Y_test - Y_pred)^2)
    }
    
    best_k_choices[t - init_train_size] <- k_values[which.min(mse_k[t - init_train_size, ])]
  }
  
  mspe_per_k <- colMeans(mse_k, na.rm = TRUE)  # NEW: MSPE averaged over time for each k
  
  return(list(
    best_k_choices = best_k_choices,
    mspe = apply(mse_k, 1, min, na.rm = TRUE),
    mspe_per_k = mspe_per_k,  # NEW
    k_values = k_values       # NEW
  ))
}
plot_mspe_vs_k <- function(mspe_per_k, k_values, model_name = "LRPS") {
  plot(k_values, mspe_per_k, type = "b", pch = 19,
       col = "darkred", xlab = "k", ylab = "Average cvMSPE")
  grid()
}
data_dir <- "your directory"   
files <- list.files(data_dir, full.names = TRUE, pattern = "*.csv")

train_files <- files[1:10]  
lag <- 1
k_values <- 1:10  # Possible k values
init_train_size <- 150  # Initial training size

ts_data <- fread(train_files[1])  # Example file
ts_data <- as.matrix(ts_data)
ts_data <- scale(apply(ts_data, 2, as.numeric))

result <- cv_time_series(ts_data, lag, k_values, init_train_size, model_type = "LRPS")
plot_mspe_vs_k(result$mspe_per_k, k_values, model_name = "LRPS")



#mspe vs lambda during cv====
library(data.table)
library(glmnet)
library(zoo)
library(doParallel)
library(foreach)

# Function to fit sparse VAR using Lasso 
sparse_var_fit <- function(X, Y, lambda) {
  p <- ncol(Y)
  B_sparse <- matrix(0, ncol(X), p)
  
  for (j in 1:p) {
    fit <- glmnet(X, Y[, j], alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE)
    B_sparse[, j] <- as.vector(coef(fit))[-1]  # Exclude intercept
  }
  return(B_sparse)
}

# Parallelised CV function for Sparse VAR
cv_sparse_var_parallel <- function(ts_data, lag, lambda_seq, init_train_size) {
  n <- nrow(ts_data)
  p <- ncol(ts_data)
  
  mspe_list <- foreach(lambda = lambda_seq,
                       .combine = rbind,
                       .packages = c("glmnet"),
                       .export = c("sparse_var_fit")) %dopar% {
                         mspe_each_time <- rep(NA, n - init_train_size - 1)
                         
                         for (t in (init_train_size + 1):(n - 1)) {
                           X_train <- ts_data[1:t, , drop = FALSE]
                           Y_train <- ts_data[-(1:lag), ][1:t, , drop = FALSE]
                           
                           X_test <- ts_data[t + 1, , drop = FALSE]
                           Y_test <- ts_data[t + 1, , drop = FALSE]
                           
                           B_hat <- sparse_var_fit(X_train, Y_train, lambda)
                           Y_pred <- X_test %*% B_hat
                           mspe_each_time[t - init_train_size] <- mean((Y_test - Y_pred)^2)
                         }
                         
                         mspe_each_time
                       }
  
  avg_mspe <- rowMeans(mspe_list, na.rm = TRUE)
  best_lambda <- lambda_seq[which.min(avg_mspe)]
  
  return(list(best_lambda = best_lambda, mspe_path = avg_mspe))
}

# Plot function: lambda vs MSPE
plot_lambda_vs_mspe <- function(lambda_seq, mspe_path, log_scale = TRUE) {
  if (log_scale) {
    plot(log10(lambda_seq), mspe_path, type = "b", pch = 19, col = "darkblue",
         xlab = expression(log[10](lambda)), ylab = "Average cvMSPE"
         #main = "MSPE vs. Lambda (Sparse VAR)"
    )
  } else {
    plot(lambda_seq, mspe_path, type = "b", pch = 19, col = "darkblue",
         xlab = expression(lambda), ylab = "Average cvMSPE"
         # main = "MSPE vs. Lambda (Sparse VAR)"
    )
  }
  grid()
}

# === SETUP PARALLEL BACKEND ===
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# === EXAMPLE USAGE ===
data_dir <- "your directory"   
files <- list.files(data_dir, full.names = TRUE, pattern = "*.csv")
train_files <- files[1:10]

lambda_seq <- 10^seq(-4, 0, length.out = 20)
lag <- 1
init_train_size <- 150

# Use first subject as example
ts_data <- fread(train_files[1]) |> as.matrix() |> scale()
if (any(is.na(ts_data))) ts_data <- na.aggregate(ts_data)

# Run parallel CV
lambda_cv_result <- cv_sparse_var_parallel(ts_data, lag, lambda_seq, init_train_size)

# Plot
plot_lambda_vs_mspe(lambda_seq, lambda_cv_result$mspe_path)

# === CLEANUP PARALLEL BACKEND ===
stopCluster(cl)


