# ==========================================
# functions_nongroup.R
# Unit-level Indicator Models for Multi-type SAE
# ==========================================
library(BayesLogit)
library(Matrix)
library(stats)

# Core precision-based MVN sampler using Cholesky decomposition
rmvn_prec <- function(b, Q) {
  Q <- (Q + t(Q)) / 2 
  L <- chol(Q)
  y <- forwardsolve(t(L), as.vector(b))
  z <- stats::rnorm(length(y))
  as.vector(backsolve(L, y + z))
}

# --- Univariate Gaussian Pseudo-Likelihood Model ---
unis_gaus <- function(X, Y, S, sig2b = 1000, wgt = NULL, 
                      predX, predS, nburn = 1000, nsim = 5000) {
  N <- length(Y); p <- ncol(X); r <- ncol(S); npred <- nrow(predX)
  if (is.null(wgt)) wgt <- rep(1, N)
  
  # Pre-compute fixed components
  X_w_X <- crossprod(X, wgt * X)
  S_w_S <- crossprod(S, wgt * S)
  X_w_S <- crossprod(X, wgt * S)
  S_w_X <- t(X_w_S)
  tX_wY <- crossprod(X, wgt * Y)
  tS_wY <- crossprod(S, wgt * Y)
  w_sum <- sum(wgt)
  Ip <- diag(p); Ir <- diag(r)

  # Initialization
  Beta <- rep(0, p); U <- rep(0, r); Mu <- rep(0, N); sig2 <- 1; Sigma2_u <- 1
  Beta_chain <- matrix(0, nrow = p, ncol = nsim)
  U_chain <- matrix(0, nrow = r, ncol = nsim)
  sig2_chain <- numeric(nsim); preds_chain <- matrix(0, nrow = npred, ncol = nsim)

  for (iter in 1:(nburn + nsim)) {
    # 1. Update sig2
    sig2 <- 1 / stats::rgamma(1, shape = 0.1 + w_sum / 2, 
                              rate = 0.1 + 0.5 * sum(wgt * (Y - Mu)^2))
    
    # 2. Update Beta
    Q_Beta <- X_w_X / sig2 + Ip / sig2b
    b_Beta <- (tX_wY - X_w_S %*% U) / sig2
    Beta   <- rmvn_prec(b_Beta, Q_Beta)
    
    # 3. Update Area Effect U & its variance Sigma2_u
    Sigma2_u <- 1 / stats::rgamma(1, shape = 0.1 + r / 2, rate = 0.1 + 0.5 * sum(U^2))
    Q_U <- S_w_S / sig2 + Ir / Sigma2_u
    b_U <- (tS_wY - S_w_X %*% Beta) / sig2
    U   <- rmvn_prec(b_U, Q_U)
    U   <- U - mean(U) # Identifiability constraint
    
    Mu <- as.vector(X %*% Beta + S %*% U)
    
    if (iter > nburn) {
      pos <- iter - nburn
      Beta_chain[, pos] <- Beta; U_chain[, pos] <- U; sig2_chain[pos] <- sig2
      preds_chain[, pos] <- predX %*% Beta + predS %*% U
    }
  }
  list(Beta = Beta_chain, U = U_chain, sig2 = sig2_chain, Preds = preds_chain)
}

# --- Univariate Binomial Pseudo-Likelihood Model (Pólya-Gamma) ---
unis_bios <- function(X, Y, S, sig2b = 1000, wgt = NULL, n_trials = NULL, 
                      predX, predS, nburn = 1000, nsim = 5000) {
  N <- length(Y); p <- ncol(X); r <- ncol(S); npred <- nrow(predX)
  if (is.null(wgt)) wgt <- rep(1, N)
  if (is.null(n_trials)) n_trials <- rep(1, N)
  
  kappa <- wgt * (Y - n_trials / 2)
  tX_k <- crossprod(X, kappa); tS_k <- crossprod(S, kappa)
  Ip <- diag(p); Ir <- diag(r)

  Beta <- rep(0, p); U <- rep(0, r); Mu <- rep(0, N); Sigma2_u <- 1
  Beta_chain <- matrix(0, nrow = p, ncol = nsim)
  preds_chain <- matrix(0, nrow = npred, ncol = nsim)

  for (iter in 1:(nburn + nsim)) {
    # 1. Update PG latent variables
    omega <- BayesLogit::rpg.gamma(N, wgt * n_trials, Mu)
    
    X_w_X <- crossprod(X, omega * X)
    S_w_S <- crossprod(S, omega * S)
    S_w_X <- crossprod(S, omega * X)

    # 2. Update Beta
    Q_Beta <- X_w_X + Ip / sig2b
    b_Beta <- tX_k - t(S_w_X) %*% U
    Beta   <- rmvn_prec(b_Beta, Q_Beta)

    # 3. Update Area Effect U & Variance
    Sigma2_u <- 1 / stats::rgamma(1, shape = 0.1 + r / 2, rate = 0.1 + 0.5 * sum(U^2))
    Q_U <- S_w_S + Ir / Sigma2_u
    b_U <- tS_k - S_w_X %*% Beta
    U   <- rmvn_prec(b_U, Q_U)
    U   <- U - mean(U)
    
    Mu <- as.vector(X %*% Beta + S %*% U)
    
    if (iter > nburn) {
      pos <- iter - nburn
      Beta_chain[, pos] <- Beta
      preds_chain[, pos] <- stats::plogis(predX %*% Beta + predS %*% U)
    }
  }
  list(Beta = Beta_chain, Preds = preds_chain)
}

# --- Multi-type Unit-Level Joint Model ---
MTSM_br <- function(X_1, X_2, Z_1, Z_2, S, sig2b = 1000, wgt = NULL, 
                    predX, predS, nburn = 1000, nsim = 5000) {
  N <- nrow(X_1); p_1 <- ncol(X_1); p_2 <- ncol(X_2); r <- ncol(S); npred <- nrow(predX)
  if (is.null(wgt)) wgt <- rep(1, N)

  # Gaussian pre-computes
  w_sum <- sum(wgt)
  X1_w_X1 <- crossprod(X_1, wgt * X_1); S_w_S <- crossprod(S, wgt * S)
  X1_w_S <- crossprod(X_1, wgt * S); S_w_X1 <- t(X1_w_S)
  tX1_wZ1 <- crossprod(X_1, wgt * Z_1); tS_wZ1 <- crossprod(S, wgt * Z_1)
  
  # Binomial pre-computes (Bernoulli case implicitly n_trials=1)
  kappa <- wgt * (Z_2 - 0.5)
  tX2_k <- crossprod(X_2, kappa); tS_k <- crossprod(S, kappa)

  Ip1 <- diag(p_1); Ip2 <- diag(p_2); Ir <- diag(r)

  # Init parameters
  Beta_1 <- rep(0, p_1); Beta_2 <- rep(0, p_2)
  eta <- rep(0, r); zeta <- rep(0, r); tau_1 <- 1
  Mu_1 <- rep(0, N); Mu_2 <- rep(0, N); sig2 <- 1
  
  preds_gaus_chain <- matrix(0, nrow = npred, ncol = nsim)
  preds_bios_chain <- matrix(0, nrow = npred, ncol = nsim)
  sig2_chain <- numeric(nsim)

  for (iter in 1:(nburn + nsim)) {
    # 1. Update Gaussian Residual Variance
    sig2 <- 1 / stats::rgamma(1, shape = 0.1 + w_sum / 2, 
                              rate = 0.1 + 0.5 * sum(wgt * (Z_1 - Mu_1)^2))

    # 2. PG Data Augmentation for Binomial
    omega <- BayesLogit::rpg.gamma(N, wgt, Mu_2)
    OM <- crossprod(S, omega * S)
    X2_w_X2 <- crossprod(X_2, omega * X_2)
    S_w_X2 <- crossprod(S, omega * X_2)

    # 3. Shared Area Effect (eta)
    sig2e <- 1 / stats::rgamma(1, shape = 0.1 + r / 2, rate = 0.1 + 0.5 * sum(eta^2))
    Q_eta <- tau_1^2 * (S_w_S / sig2) + OM + Ir / sig2e
    b_eta <- tau_1 * (tS_wZ1 - S_w_X1 %*% Beta_1) / sig2 + (tS_k - S_w_X2 %*% Beta_2 - OM %*% zeta)
    eta <- rmvn_prec(b_eta, Q_eta)
    eta <- eta - mean(eta)

    # 4. Gaussian Fixed Effects (Beta_1)
    Q_Beta_1 <- (X1_w_X1 / sig2) + Ip1 / sig2b
    b_Beta_1 <- (tX1_wZ1 - tau_1 * X1_w_S %*% eta) / sig2
    Beta_1 <- rmvn_prec(b_Beta_1, Q_Beta_1)

    # 5. Binomial Specific Area Effect (zeta)
    sig2zeta <- 1 / stats::rgamma(1, shape = 0.1 + r / 2, rate = 0.1 + 0.5 * sum(zeta^2))
    Q_zeta <- OM + Ir / sig2zeta
    b_zeta <- tS_k - S_w_X2 %*% Beta_2 - OM %*% eta
    zeta <- rmvn_prec(b_zeta, Q_zeta)
    zeta <- zeta - mean(zeta)

    # 6. Binomial Fixed Effects (Beta_2)
    Q_Beta_2 <- X2_w_X2 + Ip2 / sig2b
    b_Beta_2 <- tX2_k - t(S_w_X2) %*% (eta + zeta)
    Beta_2 <- rmvn_prec(b_Beta_2, Q_Beta_2)

    # 7. Scale Parameter (tau_1)
    M_sig2 <- S_w_S / sig2
    var_tau_1 <- 1 / (sum(eta * (M_sig2 %*% eta)) + 1 / 10) # sig2t = 10 prior
    mean_tau_1 <- var_tau_1 * sum(eta * ((tS_wZ1 - S_w_X1 %*% Beta_1) / sig2))
    tau_1 <- stats::rnorm(1, mean_tau_1, sqrt(var_tau_1))

    # Update Linear Predictors
    Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta)
    Mu_2 <- as.vector(X_2 %*% Beta_2 + S %*% eta + S %*% zeta)

    if (iter > nburn) {
      pos <- iter - nburn
      sig2_chain[pos] <- sig2
      preds_gaus_chain[, pos] <- predX %*% Beta_1 + tau_1 * predS %*% eta
      preds_bios_chain[, pos] <- stats::plogis(predX %*% Beta_2 + predS %*% eta + predS %*% zeta)
    }
  }
  list(Preds_Gaus = preds_gaus_chain, Preds_Bios = preds_bios_chain, sig2 = sig2_chain)
}

