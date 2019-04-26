
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
  res <- cv.glmnet(x = X, y = as.vector(y), alpha=1)
  as.vector(coef(res, s = res$lambda.min))
}

estimAR1 <- function(y, X, beta)
{
  r <- y[-1] - X[-1,] %*% beta[-1] - beta[1]
  r_l1 <- y[-length(y)] - X[-nrow(X),] %*% beta[-1] - beta[1]
  lm_res <- lm(r ~ r_l1)
  lm_res$coefficients[2]
}

tempDisagg <- function(y, tsl, method = "chow-lin-maxlog")
{
  X <- do.call(cbind, tsl)
  res <- tempdisagg::td(formula = y ~ X, method = method, conversion="first")
  predict(res)
}

em <- function(y, tsl, nr_comp = 1)
{
  sd_y <- sd(y)
  mean_y <- mean(y)
  
  X <- do.call(cbind, c(list(toMonthlySeries(y)), tsl))
  X <- apply(X, 2, stdize)
  
  # EM algorithm
  A <- is.na(X)
  X[A] <- 0
  has_converged = F
  iter <- 0
  while(!has_converged && iter < 1000) {
    res_prcomp <- prcomp(X, center = F, scale = F)
    L <- res_prcomp$rotation[,1:nr_comp]
    S <- X %*% L
    I <- S %*% t(L)
    has_converged <- all(abs(X[A] - I[A]) < 1e-6)
    X[A] <- I[A]
    iter <- iter+1
  }
  
  ts(X[,1] * sd_y + mean_y, start = start(y), frequency = 12)
}

chowlin <- function(y, tsl, beta, rho, chunk = c(1, 0, 0))
{
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
