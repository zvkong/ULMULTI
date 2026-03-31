suppressPackageStartupMessages({
  library(Matrix)
  library(MASS)
  library(BayesLogit)
})

rmvn <- function(n, mu, Sigma) {
  p <- nrow(Sigma)
  L <- Matrix::chol(Sigma)
  z <- Matrix::Matrix(stats::rnorm(n * p), p, n)
  t(mu + L %*% z)
}

interval_score <- function(lower, upper, truth, alpha = 0.05) {
  width <- upper - lower
  penalty_low <- (2 / alpha) * pmax(0, lower - truth)
  penalty_high <- (2 / alpha) * pmax(0, truth - upper)
  width + penalty_low + penalty_high
}

gaus_post <- function(preds, sig2chain, true_mean, region, popsize) {
  region <- as.character(region)
  regions <- unique(region)
  C <- length(regions)
  nsim <- ncol(preds)

  stopifnot(length(sig2chain) == nsim)

  idx_by_reg <- split(seq_along(region), region)
  N_by_reg <- sapply(idx_by_reg[regions], function(id) sum(popsize[id]))

  post_mean <- matrix(NA_real_, nrow = C, ncol = nsim, dimnames = list(regions, NULL))
  post <- matrix(NA_real_, nrow = C, ncol = nsim, dimnames = list(regions, NULL))

  for (r in seq_len(nsim)) {
    theta_j <- preds[, r]
    sig2_r <- sig2chain[r]

    for (c in seq_along(regions)) {
      ids <- idx_by_reg[[regions[c]]]
      Nk <- N_by_reg[c]

      mu_mean <- sum(popsize[ids] * theta_j[ids]) / Nk
      post_mean[c, r] <- mu_mean
      post[c, r] <- stats::rnorm(1, mean = mu_mean, sd = sqrt(sig2_r / Nk))
    }
  }

  if (!is.null(names(true_mean))) {
    true_mean <- true_mean[regions]
  }

  est <- rowMeans(post)
  sigma2 <- apply(post, 1, stats::var)
  lb <- apply(post, 1, stats::quantile, probs = 0.025)
  ub <- apply(post, 1, stats::quantile, probs = 0.975)
  cr <- mean(lb <= true_mean & true_mean <= ub)

  est_param <- rowMeans(post_mean)
  var_param <- apply(post_mean, 1, stats::var)
  lb_param <- apply(post_mean, 1, stats::quantile, probs = 0.025)
  ub_param <- apply(post_mean, 1, stats::quantile, probs = 0.975)
  cr_param <- mean(lb_param <= true_mean & true_mean <= ub_param)

  list(
    est = est,
    lb = lb,
    ub = ub,
    cr = cr,
    sigma2 = sigma2,
    post = post,
    est_nonoise = est_param,
    lb_nonoise = lb_param,
    ub_nonoise = ub_param,
    cr_nonoise = cr_param,
    sigma2_nonoise = var_param,
    post_nonoise = post_mean
  )
}

bios_post <- function(preds, true_mean, region, popsize) {
  region <- as.character(region)
  regions <- unique(region)
  C <- length(regions)
  R <- ncol(preds)

  idx_by_reg <- split(seq_along(region), region)
  N_by_reg <- sapply(idx_by_reg[regions], function(id) sum(popsize[id]))

  post_mean <- matrix(NA_real_, nrow = C, ncol = R, dimnames = list(regions, NULL))
  post <- matrix(NA_real_, nrow = C, ncol = R, dimnames = list(regions, NULL))

  N_int <- as.integer(round(popsize))

  for (r in seq_len(R)) {
    p_j <- preds[, r]

    for (c in seq_along(regions)) {
      ids <- idx_by_reg[[regions[c]]]
      post_mean[c, r] <- sum(popsize[ids] * p_j[ids]) / N_by_reg[c]
    }

    s_j <- stats::rbinom(n = length(N_int), size = N_int, prob = p_j)

    for (c in seq_along(regions)) {
      ids <- idx_by_reg[[regions[c]]]
      post[c, r] <- sum(s_j[ids]) / N_by_reg[c]
    }
  }

  if (!is.null(names(true_mean))) {
    true_mean <- true_mean[regions]
  }

  est <- rowMeans(post)
  sigma2 <- apply(post, 1, stats::var)
  lb <- apply(post, 1, stats::quantile, probs = 0.025)
  ub <- apply(post, 1, stats::quantile, probs = 0.975)
  cr <- mean(lb <= true_mean & true_mean <= ub)

  est_param <- rowMeans(post_mean)
  var_param <- apply(post_mean, 1, stats::var)
  lb_param <- apply(post_mean, 1, stats::quantile, probs = 0.025)
  ub_param <- apply(post_mean, 1, stats::quantile, probs = 0.975)
  cr_param <- mean(lb_param <= true_mean & true_mean <= ub_param)

  list(
    est = est,
    lb = lb,
    ub = ub,
    cr = cr,
    sigma2 = sigma2,
    post = post,
    est_nonoise = est_param,
    lb_nonoise = lb_param,
    ub_nonoise = ub_param,
    cr_nonoise = cr_param,
    sigma2_nonoise = var_param,
    post_nonoise = post_mean
  )
}

unis_bios <- function(X, Y, S,
                      sig2b = 1000, wgt = NULL, n = NULL,
                      predX, predS,
                      nburn = 1000, nsim = 5000, nthin = 1,
                      a = 0.1, b = 0.1) {
  N <- length(Y)
  p <- ncol(X)
  r <- ncol(S)
  npred <- nrow(predX)

  if (is.null(wgt)) wgt <- rep(1, N)
  if (is.null(n)) n <- rep(1, N)

  k <- wgt * (Y - n / 2)
  Ip <- Matrix::Diagonal(p)
  Ir <- Matrix::Diagonal(r)

  Mu <- rep(1, N)
  Beta <- rep(1, p)
  Sigma2_u <- 1
  U <- rep(1, r)

  n_keep <- nsim / nthin
  Beta.chain <- array(0, dim = c(p, n_keep))
  Sigma2_u.chain <- numeric(n_keep)
  U.chain <- array(0, dim = c(r, n_keep))
  Mu.chain <- array(0, dim = c(N, n_keep))
  preds.chain <- array(0, dim = c(npred, n_keep))
  logit_bios.chain <- array(0, dim = c(npred, n_keep))

  pb <- utils::txtProgressBar(min = 0, max = nsim + nburn, style = 3)

  for (index in seq_len(nsim + nburn)) {
    omega <- BayesLogit::rpg.gamma(N, wgt * n, Mu)
    Omega <- Matrix::Diagonal(x = omega)
    SO <- Matrix::t(S) %*% Omega
    OM <- SO %*% S
    gamma <- k / omega

    var.Beta <- solve(Matrix::t(X) %*% Omega %*% X + Ip / sig2b)
    mean.Beta <- var.Beta %*% (Matrix::t(X) %*% Omega %*% (gamma - S %*% U))
    Beta <- as.vector(MASS::mvrnorm(1, mu = as.vector(mean.Beta), Sigma = as.matrix(var.Beta)))

    Sigma2_u <- 1 / stats::rgamma(
      1,
      shape = a + r / 2,
      rate = b + 0.5 * drop(Matrix::t(U) %*% U)
    )

    var.U <- solve(OM + Ir / Sigma2_u)
    mean.U <- var.U %*% SO %*% (gamma - X %*% Beta)
    U <- as.vector(rmvn(1, as.vector(mean.U), var.U))

    Mu <- as.vector(X %*% Beta + S %*% U)
    logit_bios <- as.vector(predX %*% Beta + predS %*% U)
    preds <- stats::plogis(logit_bios)

    utils::setTxtProgressBar(pb, index)

    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) / nthin
      Beta.chain[, pos] <- Beta
      U.chain[, pos] <- U
      Sigma2_u.chain[pos] <- Sigma2_u
      Mu.chain[, pos] <- Mu
      preds.chain[, pos] <- preds
      logit_bios.chain[, pos] <- logit_bios
    }
  }

  close(pb)

  list(
    Beta.chain = Beta.chain,
    U.chain = U.chain,
    Sigma2_u.chain = Sigma2_u.chain,
    Mu.chain = Mu.chain,
    Preds = preds.chain,
    logit_bios.chain = logit_bios.chain
  )
}

unis_gaus <- function(X, Y, S,
                      sig2b = 1000, wgt = NULL, n = NULL,
                      predX, predS,
                      nburn = 1000, nsim = 1000, nthin = 1,
                      a = 0.1, b = 0.1,
                      a_eps = 0.1, b_eps = 0.1) {
  N <- length(Y)
  p <- ncol(X)
  r <- ncol(S)
  npred <- nrow(predX)

  if (is.null(wgt)) wgt <- rep(1, N)
  if (is.null(n)) n <- rep(1, N)

  w <- sum(wgt)
  Wgt <- Matrix::Diagonal(x = wgt)
  Ip <- Matrix::Diagonal(p)
  Ir <- Matrix::Diagonal(r)

  Beta <- rep(1, p)
  Sigma2_u <- 1
  U <- rep(0, r)
  Mu <- rep(0, N)
  sig2 <- 1

  n_keep <- nsim / nthin
  Beta.chain <- array(0, dim = c(p, n_keep))
  Sigma2_u.chain <- numeric(n_keep)
  U.chain <- array(0, dim = c(r, n_keep))
  Mu.chain <- array(0, dim = c(N, n_keep))
  sig2.chain <- numeric(n_keep)
  preds.chain <- array(0, dim = c(npred, n_keep))

  pb <- utils::txtProgressBar(min = 0, max = nsim + nburn, style = 3)
  tX <- Matrix::t(X)

  for (index in seq_len(nsim + nburn)) {
    sig2 <- 1 / stats::rgamma(
      1,
      shape = a_eps + w / 2,
      rate = b_eps + 0.5 * drop(Matrix::t(Y - Mu) %*% Wgt %*% (Y - Mu))
    )

    d <- Wgt / sig2
    SD <- Matrix::t(S) %*% d
    M <- SD %*% S
    XD <- tX %*% d %*% X

    var.Beta <- solve(XD + Ip / sig2b)
    mean.Beta <- var.Beta %*% tX %*% d %*% (Y - S %*% U)
    Beta <- as.vector(MASS::mvrnorm(1, mu = as.vector(mean.Beta), Sigma = as.matrix(var.Beta)))

    Sigma2_u <- 1 / stats::rgamma(
      1,
      shape = a + r / 2,
      rate = b + 0.5 * drop(Matrix::t(U) %*% U)
    )

    var.U <- solve(M + Ir / Sigma2_u)
    mean.U <- var.U %*% (SD %*% (Y - X %*% Beta))
    U <- as.vector(rmvn(1, as.vector(mean.U), var.U))

    Mu <- as.vector(X %*% Beta + S %*% U)
    preds <- as.vector(predX %*% Beta + predS %*% U)

    utils::setTxtProgressBar(pb, index)

    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) / nthin
      Beta.chain[, pos] <- Beta
      U.chain[, pos] <- U
      Sigma2_u.chain[pos] <- Sigma2_u
      Mu.chain[, pos] <- Mu
      sig2.chain[pos] <- sig2
      preds.chain[, pos] <- preds
    }
  }

  close(pb)

  list(
    Beta.chain = Beta.chain,
    U.chain = U.chain,
    Sigma2_u.chain = Sigma2_u.chain,
    Mu.chain = Mu.chain,
    sig2.chain = sig2.chain,
    Preds = preds.chain
  )
}

MTSM_br <- function(X_1, X_2, Z_1, Z_2, S,
                    sig2b = 1000, wgt = NULL, n = NULL,
                    predX, predS, n_preds = NULL,
                    nburn = 1000, nsim = 5000, nthin = 1,
                    sig2t = 10, sig2e = 10, tau_1_init = 1, tau_2_init = 1,
                    a_eps = 0.1, b_eps = 0.1,
                    aeta = 0.1, beta = 0.1,
                    alambda = 0.1, blambda = 0.1) {
  N <- nrow(X_1)
  r <- ncol(S)
  npred <- nrow(predX)

  if (is.null(wgt)) wgt <- rep(1, N)
  if (is.null(n)) n <- rep(1, N)

  w <- sum(wgt)
  Wgt <- Matrix::Diagonal(x = wgt)
  p_1 <- ncol(X_1)
  p_2 <- ncol(X_2)
  Ip1 <- Matrix::Diagonal(p_1)
  Ip2 <- Matrix::Diagonal(p_2)
  Ir <- Matrix::Diagonal(r)
  t_X_1 <- Matrix::t(X_1)
  t_X_2 <- Matrix::t(X_2)
  k <- wgt * (Z_2 - n / 2)

  tau_1 <- tau_1_init
  tau_2 <- tau_2_init
  Beta_1 <- rep(1, p_1)
  Beta_2 <- rep(1, p_2)
  lambda <- rep(0, r)
  eta <- rep(1, r)
  Mu_1 <- rep(1, N)
  Mu_2 <- rep(1, N)
  sig2 <- 1
  sig2l <- 1

  n_keep <- nsim / nthin
  tau_1.chain <- numeric(n_keep)
  tau_2.chain <- numeric(n_keep)
  Beta_1.chain <- array(0, dim = c(p_1, n_keep))
  Beta_2.chain <- array(0, dim = c(p_2, n_keep))
  Sigma2_lambda.chain <- numeric(n_keep)
  Sigma2_eta.chain <- numeric(n_keep)
  sig2.chain <- numeric(n_keep)
  lambda.chain <- array(0, dim = c(r, n_keep))
  eta.chain <- array(0, dim = c(r, n_keep))
  Mu_1.chain <- array(0, dim = c(N, n_keep))
  Mu_2.chain <- array(0, dim = c(N, n_keep))
  preds_gaus.chain <- array(0, dim = c(npred, n_keep))
  preds_bios.chain <- array(0, dim = c(npred, n_keep))
  logit_bios.chain <- array(0, dim = c(npred, n_keep))

  pb <- utils::txtProgressBar(min = 0, max = nsim + nburn, style = 3)

  for (index in seq_len(nsim + nburn)) {
    sig2 <- 1 / stats::rgamma(
      1,
      shape = a_eps + w / 2,
      rate = b_eps + 0.5 * drop(Matrix::t(Z_1 - Mu_1) %*% Wgt %*% (Z_1 - Mu_1))
    )

    d <- Wgt / sig2
    SD <- Matrix::t(S) %*% d
    M <- SD %*% S
    XD <- t_X_1 %*% d %*% X_1

    omega <- BayesLogit::rpg.gamma(N, wgt * n, Mu_2)
    Omega <- Matrix::Diagonal(x = omega)
    SO <- Matrix::t(S) %*% Omega
    OM <- SO %*% S
    gamma <- k / omega

    sig2e <- 1 / stats::rgamma(
      1,
      shape = aeta + r / 2,
      rate = beta + 0.5 * drop(Matrix::t(eta) %*% eta)
    )

    var.eta <- solve(
      tau_1 * M * tau_1 + Ir / sig2e + tau_2 * OM * tau_2,
      sparse = TRUE
    )

    mean.eta <- var.eta %*% (
      tau_1 * SD %*% (Z_1 - X_1 %*% Beta_1) +
        tau_2 * SO %*% (gamma - X_2 %*% Beta_2 - S %*% lambda)
    )

    eta <- as.vector(rmvn(1, as.vector(mean.eta), var.eta))

    var.Beta_1 <- solve(XD + Ip1 / sig2b)
    mean.Beta_1 <- var.Beta_1 %*% t_X_1 %*% d %*% (Z_1 - tau_1 * S %*% eta)
    Beta_1 <- as.vector(MASS::mvrnorm(1, mu = as.vector(mean.Beta_1), Sigma = as.matrix(var.Beta_1)))

    sig2l <- 1 / stats::rgamma(
      1,
      shape = alambda + r / 2,
      rate = blambda + 0.5 * drop(Matrix::t(lambda) %*% lambda)
    )

    var.lambda <- solve(OM + Ir / sig2l, sparse = TRUE)
    mean.lambda <- var.lambda %*% (SO %*% (gamma - X_2 %*% Beta_2 - tau_2 * S %*% eta))
    lambda <- as.vector(rmvn(1, as.vector(mean.lambda), var.lambda))
    lambda <- lambda - mean(lambda)

    var.Beta_2 <- solve(t_X_2 %*% Omega %*% X_2 + Ip2 / sig2b)
    mean.Beta_2 <- var.Beta_2 %*% t_X_2 %*% Omega %*% (gamma - tau_2 * S %*% eta - S %*% lambda)
    Beta_2 <- as.vector(MASS::mvrnorm(1, mu = as.vector(mean.Beta_2), Sigma = as.matrix(var.Beta_2)))

    var_tau_2 <- as.numeric(solve(drop(Matrix::t(eta) %*% OM %*% eta) + 1 / sig2t))
    mean_tau_2 <- as.numeric(
      var_tau_2 * drop(Matrix::t(eta) %*% SO %*% (gamma - X_2 %*% Beta_2 - S %*% lambda))
    )
    tau_2 <- stats::rnorm(1, mean_tau_2, sqrt(var_tau_2))

    Mu_1 <- as.vector(X_1 %*% Beta_1 + tau_1 * S %*% eta)
    Mu_2 <- as.vector(X_2 %*% Beta_2 + tau_2 * S %*% eta + S %*% lambda)

    preds_gaus <- as.vector(predX %*% Beta_1 + tau_1 * predS %*% eta)
    logit_bios <- as.vector(predX %*% Beta_2 + tau_2 * predS %*% eta + predS %*% lambda)
    preds_bios <- stats::plogis(logit_bios)

    utils::setTxtProgressBar(pb, index)

    if (index > nburn && (index - nburn) %% nthin == 0) {
      pos <- (index - nburn) / nthin
      tau_1.chain[pos] <- tau_1
      tau_2.chain[pos] <- tau_2
      Beta_1.chain[, pos] <- Beta_1
      Beta_2.chain[, pos] <- Beta_2
      lambda.chain[, pos] <- lambda
      eta.chain[, pos] <- eta
      sig2.chain[pos] <- sig2
      Sigma2_lambda.chain[pos] <- sig2l
      Sigma2_eta.chain[pos] <- sig2e
      Mu_1.chain[, pos] <- Mu_1
      Mu_2.chain[, pos] <- Mu_2
      preds_gaus.chain[, pos] <- preds_gaus
      preds_bios.chain[, pos] <- preds_bios
      logit_bios.chain[, pos] <- logit_bios
    }
  }

  close(pb)

  list(
    Beta_1.chain = Beta_1.chain,
    Beta_2.chain = Beta_2.chain,
    lambda.chain = lambda.chain,
    eta.chain = eta.chain,
    Sigma2_lambda.chain = Sigma2_lambda.chain,
    Sigma2_eta.chain = Sigma2_eta.chain,
    sig2.chain = sig2.chain,
    Mu_1.chain = Mu_1.chain,
    Mu_2.chain = Mu_2.chain,
    preds_gaus.chain = preds_gaus.chain,
    preds_bios.chain = preds_bios.chain,
    logit_bios.chain = logit_bios.chain,
    tau_1.chain = tau_1.chain,
    tau_2.chain = tau_2.chain
  )
}
