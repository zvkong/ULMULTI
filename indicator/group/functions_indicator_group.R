# ==========================================
# functions_group.R
# Grouped Exact-Collapse Models for Multi-type SAE
# ==========================================
# (Make sure to source functions_nongroup.R first to get rmvn_prec and unis_gaus)

# --- Grouped Univariate Binomial Model ---
unis_bios_grouped <- function(X_grp, y_sum, n_sum, S_grp, sig2b = 1000, 
                              predX, predS, nburn = 1000, nsim = 5000) {
  G <- nrow(X_grp); p <- ncol(X_grp); r <- ncol(S_grp); npred <- nrow(predX)

  kappa <- y_sum - n_sum / 2
  tX_k <- crossprod(X_grp, kappa); tS_k <- crossprod(S_grp, kappa)
  Ip <- diag(p); Ir <- diag(r)

  Mu <- rep(0, G); Beta <- rep(0, p); U <- rep(0, r)
  preds_chain <- matrix(0, nrow = npred, ncol = nsim)

  for (iter in 1:(nburn + nsim)) {
    # PG Latent Variable - dimension G instead of N!
    omega <- BayesLogit::rpg.gamma(G, n_sum, Mu) 
    
    X_w_X <- crossprod(X_grp, omega * X_grp)
    S_w_S <- crossprod(S_grp, omega * S_grp)
    S_w_X <- crossprod(S_grp, omega * X_grp)

    Q_Beta <- X_w_X + Ip / sig2b
    b_Beta <- tX_k - t(S_w_X) %*% U
    Beta <- rmvn_prec(b_Beta, Q_Beta)

    Sigma2_u <- 1 / stats::rgamma(1, shape = 0.1 + r / 2, rate = 0.1 + 0.5 * sum(U^2))
    Q_U <- S_w_S + Ir / Sigma2_u
    b_U <- tS_k - S_w_X %*% Beta
    U <- rmvn_prec(b_U, Q_U)
    U <- U - mean(U)

    Mu <- as.vector(X_grp %*% Beta + S_grp %*% U)

    if (iter > nburn) {
      preds_chain[, iter - nburn] <- stats::plogis(predX %*% Beta + predS %*% U)
    }
  }
  list(Preds = preds_chain)
}

# --- Multi-type Joint Model (Gaussian Unit-level + Binomial Grouped) ---
MTSM_br_grouped <- function(X_1, Z_1, S_1, wgt_1, 
                            X_grp, y_sum, n_sum, S_grp, 
                            sig2b = 1000, predX, predS, nburn = 1000, nsim = 5000) {
  
  N1 <- nrow(X_1); G <- nrow(X_grp); p_1 <- ncol(X_1); p_2 <- ncol(X_grp); r <- ncol(S_1)
  npred <- nrow(predX)

  # Gaussian (Unit-level) Pre-computes
  w_sum <- sum(wgt_1)
  X1_w_X1 <- crossprod(X_1, wgt_1 * X_1); S1_w_S1 <- crossprod(S_1, wgt_1 * S_1)
  X1_w_S1 <- crossprod(X_1, wgt_1 * S_1); S1_w_X1 <- t(X1_w_S1)
  tX1_wZ1 <- crossprod(X_1, wgt_1 * Z_1); tS1_wZ1 <- crossprod(S_1, wgt_1 * Z_1)

  # Binomial (Grouped) Pre-computes
  kappa <- y_sum - n_sum / 2
  tX2_k <- crossprod(X_grp, kappa); tS2_k <- crossprod(S_grp, kappa)

  Ip1 <- diag(p_1); Ip2 <- diag(p_2); Ir <- diag(r)

  Beta_1 <- rep(0, p_1); Beta_2 <- rep(0, p_2)
  eta <- rep(0, r); zeta <- rep(0, r); tau_1 <- 1
  Mu_1 <- rep(0, N1); Mu_2 <- rep(0, G); sig2 <- 1
  
  preds_gaus_chain <- matrix(0, nrow = npred, ncol = nsim)
  preds_bios_chain <- matrix(0, nrow = npred, ncol = nsim)
  sig2_chain <- numeric(nsim)

  for (iter in 1:(nburn + nsim)) {
    # 1. Gaussian Residual Variance (Unit Level)
    sig2 <- 1 / stats::rgamma(1, shape = 0.1 + w_sum / 2, rate = 0.1 + 0.5 * sum(wgt_1 * (Z_1 - Mu_1)^2))

    # 2. PG Latent Variable (Grouped Level)
    omega <- BayesLogit::rpg.gamma(G, n_sum, Mu_2)
    OM <- crossprod(S_grp, omega * S_grp)
    X2_w_X2 <- crossprod(X_grp, omega * X_grp)
    S2_w_X2 <- crossprod(S_grp, omega * X_grp)

    # 3. Shared Area Effect
    sig2e <- 1 / stats::rgamma(1, shape = 0.1 + r / 2, rate = 0.1 + 0.5 * sum(eta^2))
    Q_eta <- (tau_1^2) * (S1_w_S1 / sig2) + OM + Ir / sig2e
    b_eta <- tau_1 * (tS1_wZ1 - S1_w_X1 %*% Beta_1) / sig2 + (tS2_k - S2_w_X2 %*% Beta_2 - OM %*% zeta)
    eta <- rmvn_prec(b_eta, Q_eta)
    eta <- eta - mean(eta) 

    # 4. Gaussian Fixed Effects
    Q_Beta_1 <- X1_w_X1 / sig2 + Ip1 / sig2b
    b_Beta_1 <- (tX1_wZ1 - tau_1 * X1_w_S1 %*% eta) / sig2
    Beta_1 <- rmvn_prec(b_Beta_1, Q_Beta_1)

    # 5. Binomial Specific Area Effect (zeta)
    sig2zeta <- 1 / stats::rgamma(1, shape = 0.1 + r / 2, rate = 0.1 + 0.5 * sum(zeta^2))
    Q_zeta <- OM + Ir / sig2zeta
    b_zeta <- tS2_k - S2_w_X2 %*% Beta_2 - OM %*% eta
    zeta <- rmvn_prec(b_zeta, Q_zeta)
    zeta <- zeta - mean(zeta)

    # 6. Binomial Fixed Effects
    Q_Beta_2 <- X2_w_X2 + Ip2 / sig2b
    b_Beta_2 <- tX2_k - t(S2_w_X2) %*% (eta + zeta)
    Beta_2 <- rmvn_prec(b_Beta_2, Q_Beta_2)

    # 7. Scale Parameter
    M_sig2 <- S1_w_S1 / sig2
    var_tau_1 <- 1 / (sum(eta * (M_sig2 %*% eta)) + 1 / 10)
    mean_tau_1 <- var_tau_1 * sum(eta * ((tS1_wZ1 - S1_w_X1 %*% Beta_1) / sig2))
    tau_1 <- stats::rnorm(1, mean_tau_1, sqrt(var_tau_1))

    Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S_1 %*% eta)
    Mu_2 <- as.vector(X_grp %*% Beta_2 + S_grp %*% eta + S_grp %*% zeta)

    if (iter > nburn) {
      pos <- iter - nburn
      sig2_chain[pos] <- sig2
      preds_gaus_chain[, pos] <- predX %*% Beta_1 + tau_1 * predS %*% eta
      preds_bios_chain[, pos] <- stats::plogis(predX %*% Beta_2 + predS %*% eta + predS %*% zeta)
    }
  }
  list(Preds_Gaus = preds_gaus_chain, Preds_Bios = preds_bios_chain, sig2 = sig2_chain)
}

