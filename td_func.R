
computeStat <- function(res, start_date = ymDate(1987,3)) {
  stat_rt <- map(names(res), function(var) {
    var_tsl <- sub("balance\\.composite", "balance", var)
    rr <- res[[var]]
    rr <- rr[as.Date(names(rr)) >= start_date]
    rt <- ts(map_dbl(rr, ~tail(.x$ytd,1)), start = dateToTime(names(rr)[1]), frequency = 12)
    rt <- window(rt, start=dateToTime(addMonths(startDate(rt), 0)))
    ym <- tsl[[var_tsl]]
    res <- list(
      rt = rt,
      rt_rmse_m23=rmse(rt, ym, c(2,3,5,6,8,9,11,12)),
      rt_rmse_m2=rmse(rt, ym, c(2,5,8,11)),
      rt_rmse_m3=rmse(rt, ym, c(3,6,9,12)),
      post_rmse_m23 = sqrt(mean(map_dbl(rr, ~mse(.x$ytd, ym, c(2,3,5,6,8,9,11,12))))),
      post_rmse_m2 = sqrt(mean(map_dbl(rr, ~mse(.x$ytd, ym, c(2,5,8,11))))),
      post_rmse_m3 = sqrt(mean(map_dbl(rr, ~mse(.x$ytd, ym, c(3,6,9,12))))))
    res$rt_rrmse_m2_vs_m3 <- res$rt_rmse_m2 / res$rt_rmse_m3
    res$post_rrmse_m2_vs_m3 <- res$post_rmse_m2 / res$post_rmse_m3
    res
  })
  names(stat_rt) <- names(res)
  stat_rt
}

lassoBIC <- function(y, X, penalty_factors = rep(1, ncol(X)))
{
  glm_result <- glmnet(X, y, alpha = 1, penalty.factor = penalty_factors)
  pred_result <- predict(glm_result, X)
  n <- length(y)
  lambda_bics <- map_dbl(seq_along(glm_result$lambda), function(j) {
    log(mean((y - pred_result[,j])^2)) + glm_result$df[j] * log(n)/n
  })
  best_index <- which.min(lambda_bics)
  min_bic <- lambda_bics[best_index]
  best_lambda <- glm_result$lambda[best_index]
  list(beta=as.vector(coef(glm_result, s = best_lambda)), min_bic=min_bic)
}

lassoCV <- function(y, tsl)
{
  X <- do.call(cbind, tsl)
  res <- cv.glmnet(x = X, y = as.vector(y), alpha=1, nfolds = min(floor(length(y) / 8), 10), pmax = length(y)-1)
  as.vector(coef(res, s = res$lambda.min))
}

estimAR1 <- function(y, X, beta)
{
  r <- y[-1] - X[-1,] %*% beta[-1] - beta[1]
  r_l1 <- y[-length(y)] - X[-nrow(X),] %*% beta[-1] - beta[1]
  lm_res <- lm(r ~ r_l1)
  max(0, lm_res$coefficients[2])
}

tempDisagg <- function(y, tsl, method = "chow-lin-maxlog")
{
  if(length(tsl)==0) {
    res <- tempdisagg::td(formula = y ~ 1, method = method, conversion="first")
  } else {
    X <- do.call(cbind, tsl)
    res <- tempdisagg::td(formula = y ~ X, method = method, conversion="first")
  }
  predict(res)
}

em <- function(y, tsl, nr_comp = 1)
{
  sd_y <- sd(y, na.rm = T)
  mean_y <- mean(y, na.rm = T)
  
  X <- do.call(cbind, c(list(y), tsl))
  X <- scale(X, T, T)
  A <- is.na(X)
  X_orig <- X
  
  # EM algorithm
  X[A] <- 0
  has_converged = F
  iter <- 0
  while(!has_converged && iter < 1000) {
    X <- scale(X, T, T)
    res_prcomp <- prcomp(X, center = F, scale = F)
    L <- res_prcomp$rotation[,1:nr_comp]
    S <- X %*% L
    I <- S %*% t(L)
    X <- X_orig
    X[A] <- I[A]
    has_converged <- all(abs(X[A] - I[A]) < 1e-6)
    iter <- iter+1
  }
  
  ts(X[,1] * sd_y + mean_y, start = start(y), frequency = 12)
}

chowlin <- function(y, tsl, beta, rho, chunk = c(1, 0, 0))
{
  rho <- rho^(1/3)
  # Create matrix from indicator series.
  X <- do.call(cbind, tsl)
  
  # Compute distribution matrix.
  zeros <- rep(0, 3)
  nl <- length(y)
  n <- 3*nl
  C <- t(sapply(1:nl, function(i) c(rep(zeros, i-1), chunk, rep(zeros, nl-i))))
  row <- sapply(0:(n-1), function(e) rho^e)
  S <- sapply(n:1, function(i) c(rep(NA, n-i),row[1:i]))
  A <- is.na(S)
  S[A] <- t(S)[A]
  S <- S / (1+rho^2)
  D <- S %*% t(C) %*% solve(C %*% S %*% t(C))
  
  # Do temp disagg.
  p <- X %*% beta[-1] + beta[1]
  u <- y - C %*% p
  ts(as.vector(p + D %*% u), start=start(tsl[[1]]), frequency = 12)
}

lassoAR1 <- function(y, tsl)
{
  X <- do.call(cbind, tsl)
  
  # Init with ridge regression with lambda ~ 0.
  glm_result <- glmnet(X, y, alpha = 0, lambda = 0)
  beta <- as.vector(coef(glm_result, s = 0))
  beta_zero <- abs(beta[-1]) < 1e-4
  beta[c(F, beta_zero)] <- 0
  
  # Estimate rho.
  rho <- estimAR1(y, X, beta)

  iter_result <- list()
  beta_diff <- 1
  iter <- 0
  while(beta_diff > 1e-4 && iter < 100 && sum(!beta_zero) > 1) {
    
    # Do Lasso with estimated rho.
    y_diff <- y[-1] - rho * y[-length(y)]
    X_diff <- X[-1,] - rho * X[-nrow(X),]
    
    # Remove columns with zero beta.
    X_diff <- X_diff[,!beta_zero]
    
    penalty_factors <- 1 / abs(beta[-1][!beta_zero])
    # Lasso with intercept.
    lasso_res <- lassoBIC(y_diff, X_diff, penalty_factors)
    
    # Save previous beta/rho.
    prev_beta <- beta
    prev_rho <- rho
    
    beta[c(T, !beta_zero)] <- lasso_res$beta
    
    # Correct intercept.
    beta[1] <- beta[1] / (1 - rho)
    
    beta_zero <- abs(beta[-1]) < 1e-4
    beta[c(F, beta_zero)] <- 0
    
    # Estimate rho.
    rho <- estimAR1(y, X, beta)
    
    beta_diff <- sum(abs(beta[-1] - prev_beta[-1])) + abs(rho - prev_rho)
    
    iter <- iter+1
    iter_result[[iter]] <- list(rho = rho, beta = beta, min_bic = lasso_res$min_bic)
  }
  
  if(iter == 100) {
    c(iter_result[[which.min(map_dbl(iter_result, ~.x$min_bic))]], list(iter=iter))
  } else {
    c(iter_result[[iter]], list(iter=iter))
  }
}
