library(MASS)
library(vars)
library(tsDyn)

##flat eigenvalues====
singleSim <- function(n = 50, p = 2, rho = 0.2, lag = 1) {
  Sigma_E <- diag(p)
  
  # off if 0.2 lrps is good
  # off if 0.01 lrps maybe bad
  
  off_diag <- 0.01
  B <- matrix(off_diag, nrow = p, ncol = p)
  diag(B) <- rho 
  
  # Ensure stationarity by normalizing if necessary
  max_eigen <- max(abs(eigen(B)$values))
  if (max_eigen >= 1) {
    B <- B / (max_eigen + 0.1)  # Scale down to keep largest eigenvalue < 1
  }
  
  ts <- VAR.sim(B=B, n=n, lag = 1, include="none", varcov = Sigma_E)
  
  Y <- ts[-(1:lag), ]  
  X <- ts[1:(n-lag),]
  
  signal_frobenius <- norm(X %*% t(B), "F")  
  noise_frobenius <- norm(Sigma_E, "F")  
  snr <- signal_frobenius / noise_frobenius 
  
  return(list(ts = ts, B = B, Y = Y, X = X, snr = snr))
}


##not flat eigenvalues====
singleSim <- function(n = 50, p = 2, rho = 0.2, lag = 1) {
  Sigma_E <- diag(p)
  
  buildS <- function(n, rho) {
    idi <- idj <- 1:n
    tmp <- rho^(abs(outer(idi, idi, "-")))
    return(tmp)
  }
  
  B <- buildS(p, rho) 
  
  # Ensure stationarity by normalizing if necessary
  max_eigen <- max(abs(eigen(B)$values))
  if (max_eigen >= 1) {
    B <- B / (max_eigen + 0.1)  # Scale down to keep largest eigenvalue < 1
  }
  
  ts <- VAR.sim(B=B, n=n, lag = 1, include="none", varcov = Sigma_E)
  
  Y <- ts[-(1:lag), ]  
  X <- ts[1:(n-lag),]
  
  signal_frobenius <- norm(X %*% t(B), "F")  
  noise_frobenius <- norm(Sigma_E, "F")  
  snr <- signal_frobenius / noise_frobenius 
  
  return(list(ts = ts, B = B, Y = Y, X = X, snr = snr))
}

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

# Function to compute Mean Squared Error
compute_mse <- function(B_true, B_pred) {
  return(mean((B_true - B_pred)^2))
}

##Cross validation====
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

fixed_k_prediction <- function(ts, lag, k, init_train_size, model_type = c("LRPS", "VAR")) {
  model_type <- match.arg(model_type)
  
  ts <- as.matrix(ts)
  ts <- apply(ts, 2, as.numeric)
  
  n <- nrow(ts)
  mspe <- rep(NA, n - init_train_size - 1)
  
  for (t in (init_train_size + 1):(n - 1)) {
    X_train <- ts[1:t, , drop = FALSE]
    Y_train <- ts[-(1:lag), ][1:t, , drop = FALSE]
    
    X_test <- ts[t + 1, , drop = FALSE]
    Y_test <- ts[t + 1, , drop = FALSE]
    
    if (model_type == "VAR") {
      model <- var_fit(X_train, lag)
      Y_pred <- X_test %*% model$coef
    } else {
      model <- LRPS_varfit(X_train, lag, k)
      Y_pred <- X_test %*% model$coef
    }
    mspe[t - init_train_size] <- mean((Y_test - Y_pred)^2)
  }
  
  return(mspe)
}
library(foreach)
library(doParallel)

n_cores <- parallel::detectCores() - 1  # use one less than total
cl <- makeCluster(n_cores)
registerDoParallel(cl)

run_global_k_experiment <- function(n, p, rho_values, lag, k_values, init_train_size, reps = 10) {
  all_results <- list()
  
  for (rho in rho_values) {
    # Phase 1: parallel CV to choose best k per subject
    selected_ks <- foreach(r = 1:reps, .combine = c,
                           .packages = c("vars", "MASS", "tsDyn"),
                           .export = c("singleSim", "cv_time_series", "var_fit", "LRPS_varfit", "LRfun")) %dopar% {
                             sim_train <- singleSim(n = n, p = p, rho = rho, lag = lag)
                             cv_res <- cv_time_series(sim_train$ts, lag, k_values, init_train_size, model_type = "LRPS")
                             cv_res$best_k
                           }
    
    global_k <- as.numeric(names(which.max(table(selected_ks))))
    
    # Phase 2: parallel evaluation using fixed global_k
    eval_results <- foreach(r = 1:reps, .combine = rbind,
                            .packages = c("vars", "MASS", "tsDyn"),
                            .export = c("singleSim", "fixed_k_prediction", "var_fit", "LRPS_varfit", "LRfun")) %dopar% {
                              sim_eval <- singleSim(n = n, p = p, rho = rho, lag = lag)
                              mspe_lrps <- fixed_k_prediction(sim_eval$ts, lag, global_k, init_train_size, model_type = "LRPS")
                              mspe_var  <- fixed_k_prediction(sim_eval$ts, lag, NULL, init_train_size, model_type = "VAR")
                              cbind(mspe_lrps, mspe_var)
                            }
    
    
    all_results[[paste0("rho_", rho)]] <- list(
      global_k = global_k,
      selected_ks = selected_ks,
      mspe_lrps = as.vector(eval_results[, 1]),
      mspe_var  = as.vector(eval_results[, 2])
    )
  }
  
  return(all_results)
}


plot_k_histogram_global <- function(results, rho_values, p) {
  par(mfrow = c(1, length(rho_values)), mar = c(4, 4, 3, 1))
  
  max_y <- max(sapply(rho_values, function(rho) {
    table(factor(results[[paste0("rho_", rho)]]$selected_ks, levels = 1:(p - 1))) |> max()
  }))
  
  for (rho in rho_values) {
    k_vals <- results[[paste0("rho_", rho)]]$selected_ks
    freq <- table(factor(k_vals, levels = 1:(p - 1)))
    global_k <- results[[paste0("rho_", rho)]]$global_k
    
    barplot(freq, col = "blue", border = "black",
            names.arg = 1:(p - 1), xlab = "Chosen k",
            ylab = if (rho == rho_values[1]) "Frequency" else "",
            main = bquote(rho == .(rho) ~ "(" ~ hat(k) == .(global_k) * ")"),
            ylim = c(0, max_y + 2))
  }
}

plot_mspe_density_multi <- function(results, rho_values) {
  par(mfrow = c(1, length(rho_values)), mar = c(4, 4, 3, 1))
  
  all_x <- all_y <- c()
  densities <- list()
  
  for (rho in rho_values) {
    res <- results[[paste0("rho_", rho)]]
    d_lrps <- density(res$mspe_lrps, na.rm = TRUE)
    d_var  <- density(res$mspe_var, na.rm = TRUE)
    
    densities[[as.character(rho)]] <- list(lrps = d_lrps, var = d_var)
    all_x <- c(all_x, d_lrps$x, d_var$x)
    all_y <- c(all_y, d_lrps$y, d_var$y)
  }
  
  xlim <- range(all_x)
  ylim <- range(all_y)
  
  for (i in seq_along(rho_values)) {
    rho <- rho_values[i]
    d_lrps <- densities[[as.character(rho)]]$lrps
    d_var  <- densities[[as.character(rho)]]$var
    
    plot(d_lrps, col = "blue", lwd = 2, xlim = xlim, ylim = ylim,
         main = bquote(rho == .(rho)), xlab = "MSPE",
         ylab = if (i == 1) "Density" else "")
    lines(d_var, col = "red", lwd = 2, lty = 2)
    
    legend("topright", legend = c("LRPS", "VAR"),
           col = c("blue", "red"), lwd = 2, lty = c(1, 2))
  }
}
set.seed(123)

rho_values <- c(0.1, 0.5, 0.9)
k_values <- 1:29
n <- 100
p <- 30
init_train_size <- 50
lag <- 1
reps <- 50  # number of subjects

results_global <- run_global_k_experiment(n, p, rho_values, lag, k_values, init_train_size, reps)
stopCluster(cl)

# Plots
plot_k_histogram_global(results_global, rho_values, p)
plot_mspe_density_multi(results_global, rho_values)

pdf("cv_rank_0.01.pdf", width = 7, height = 3)
plot_k_histogram_global(results_global, rho_values, p)
dev.off()


pdf("mspe_0.01.pdf", width = 7, height = 3)
plot_mspe_density_multi(results_global, rho_values)
dev.off()