suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(sampling)
  library(ggplot2)
  library(patchwork)
  library(tibble)
})
source("packages.R")
source("functions_indicator_group_revised.R")

pums23 <- readr::read_csv("IL21.csv", show_col_types = FALSE)

pums1 <- pums23 |>
  dplyr::select(PUMA, PWGTP, SEX, RACWHT, POVPIP, PINCP, SCHL) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    PUMA    = sprintf("%05d", as.integer(PUMA)),
    PWGTP   = as.numeric(PWGTP),
    POVPIP  = as.numeric(POVPIP),
    SEX     = factor(SEX),
    RACE    = factor(RACWHT),
    LOGINCO = log(as.numeric(PINCP)),
    INCOME  = as.numeric(PINCP),
    POV     = dplyr::case_when(
      POVPIP <= 100 ~ 1,
      POVPIP > 100  ~ 0
    ),
    BACH = dplyr::case_when(
      SCHL >= 21 ~ 1,
      SCHL < 21  ~ 0
    ) |> factor()
  )

pums1 <- pums1[is.finite(pums1$LOGINCO), ]

pums <- pums1 |>
  dplyr::select(PUMA, PWGTP, SEX, RACE, POV, LOGINCO, INCOME, BACH) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    POV = as.numeric(POV),
    EDU = as.numeric(BACH) - 1
  )

pums$INCO <- pums$LOGINCO
pums$INCO <- (pums$INCO - min(pums$INCO)) / (max(pums$INCO) - min(pums$INCO))

puma_sf <- tigris::pumas(state = "IL", cb = TRUE, year = 2020, class = "sf")

if ("PUMACE20" %in% names(puma_sf)) {
  puma_sf$PUMA <- sprintf("%05d", as.integer(puma_sf$PUMACE20))
} else if ("PUMACE10" %in% names(puma_sf)) {
  puma_sf$PUMA <- sprintf("%05d", as.integer(puma_sf$PUMACE10))
} else if ("GEOID20" %in% names(puma_sf)) {
  geo_tmp <- as.character(puma_sf$GEOID20)
  puma_sf$PUMA <- substr(geo_tmp, nchar(geo_tmp) - 4, nchar(geo_tmp))
} else if ("GEOID" %in% names(puma_sf)) {
  geo_tmp <- as.character(puma_sf$GEOID)
  puma_sf$PUMA <- substr(geo_tmp, nchar(geo_tmp) - 4, nchar(geo_tmp))
} else {
  stop("Cannot identify the PUMA code column in puma_sf.")
}

puma_sf <- puma_sf |>
  dplyr::filter(PUMA %in% unique(pums$PUMA)) |>
  dplyr::arrange(PUMA)

area_levels <- puma_sf$PUMA

pums <- pums |>
  dplyr::filter(PUMA %in% area_levels) |>
  dplyr::arrange(factor(PUMA, levels = area_levels), SEX, BACH)

truth <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(
    INCO = mean(INCO),
    POV  = mean(POV),
    .groups = "drop"
  ) |>
  dplyr::arrange(factor(PUMA, levels = area_levels))

pcells <- pums |>
  dplyr::group_by(PUMA, SEX, BACH) |>
  dplyr::summarise(popsize = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(factor(PUMA, levels = area_levels), SEX, BACH)

pcells_region <- dplyr::tibble(PUMA = area_levels) |>
  dplyr::left_join(
    pums |>
      dplyr::group_by(PUMA) |>
      dplyr::summarise(popsize = dplyr::n(), .groups = "drop"),
    by = "PUMA"
  )

nb <- spdep::poly2nb(puma_sf, queen = TRUE)
A <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
A <- (A + t(A)) / 2

eig <- eigen(A, symmetric = TRUE)
pos_idx <- which(eig$values > 1e-8)

basis_mat <- eig$vectors[, pos_idx, drop = FALSE]
rownames(basis_mat) <- area_levels
colnames(basis_mat) <- paste0("bf", seq_len(ncol(basis_mat)))

predX <- model.matrix(~ SEX + BACH - 1, data = pcells)
predPsi <- basis_mat[match(pcells$PUMA, rownames(basis_mat)), , drop = FALSE]

n_sim <- 100
nsim  <- 500
nburn <- 500
nthin <- 1

n_area <- length(area_levels)

ubios_pre_b   <- array(NA_real_, dim = c(n_area, n_sim))
mbios_pre_b   <- array(NA_real_, dim = c(n_area, n_sim))
ugaus_pre_b   <- array(NA_real_, dim = c(n_area, n_sim))
mgaus_pre_b   <- array(NA_real_, dim = c(n_area, n_sim))
dgaus_pre     <- array(NA_real_, dim = c(n_area, n_sim))
dbio_pre      <- array(NA_real_, dim = c(n_area, n_sim))

ubios_qual_b  <- array(NA_real_, dim = c(n_area, 2, n_sim))
mbios_qual_b  <- array(NA_real_, dim = c(n_area, 2, n_sim))
ugaus_qual_b  <- array(NA_real_, dim = c(n_area, 2, n_sim))
mgaus_qual_b  <- array(NA_real_, dim = c(n_area, 2, n_sim))

cor_set <- numeric(n_sim)

for (k in seq_len(n_sim)) {
  set.seed(k)
  ss <- 1000

  prob <- sampling::inclusionprobabilities(
    pums$PWGTP * (1 + 5 * (pums$POV == 1)),
    ss
  )

  Ind <- sampling::UPsystematic(prob)

  samp <- pums[as.logical(Ind), ]
  samp$P <- prob[as.logical(Ind)]
  samp$W <- 1 / samp$P
  samp$scaledWGT <- samp$W * nrow(samp) / sum(samp$W)

  compare_df <- dplyr::tibble(PUMA = area_levels) |>
    dplyr::left_join(
      samp |>
        dplyr::group_by(PUMA) |>
        dplyr::summarise(
          unweighted_means     = mean(INCO),
          unweighted_means_POV = mean(POV),
          weighted_means       = stats::weighted.mean(INCO, W),
          weighted_means_POV   = stats::weighted.mean(POV, W),
          .groups = "drop"
        ),
      by = "PUMA"
    )

  cor_set[k] <- stats::cor(samp$POV, samp$INCO)

  modwgt <- samp$scaledWGT
  modX   <- model.matrix(~ SEX + BACH - 1, data = samp)
  modPsi <- basis_mat[match(samp$PUMA, rownames(basis_mat)), , drop = FALSE]
  modY   <- samp$INCO
  modZ   <- samp$POV

  group_binom <- samp |>
    dplyr::group_by(PUMA, SEX, BACH) |>
    dplyr::summarise(
      y_sum = sum(scaledWGT * POV),
      n_sum = sum(scaledWGT),
      .groups = "drop"
    ) |>
    dplyr::arrange(factor(PUMA, levels = area_levels), SEX, BACH)

  groupX   <- model.matrix(~ SEX + BACH - 1, data = group_binom)
  groupPsi <- basis_mat[match(group_binom$PUMA, rownames(basis_mat)), , drop = FALSE]

  unis_wage_b <- unis_gaus(
    X      = modX,
    Y      = modY,
    S      = modPsi,
    sig2b  = 1000,
    wgt    = modwgt,
    n      = NULL,
    predX  = predX,
    predS  = predPsi,
    nburn  = nburn,
    nsim   = nsim,
    nthin  = nthin,
    a      = 0.5,
    b      = 0.5,
    a_eps  = 0.1,
    b_eps  = 0.1
  )

  unis_pov_b <- unis_bios_grouped(
    X      = groupX,
    y_sum  = group_binom$y_sum,
    n_sum  = group_binom$n_sum,
    S      = groupPsi,
    sig2b  = 1000,
    predX  = predX,
    predS  = predPsi,
    nburn  = nburn,
    nsim   = nsim,
    nthin  = nthin,
    a      = 0.1,
    b      = 0.1
  )

  mult_pov_b <- MTSM_br_grouped(
    X_1      = modX,
    Z_1      = modY,
    S_1      = modPsi,
    wgt_1    = modwgt,
    X_2      = groupX,
    y_sum    = group_binom$y_sum,
    n_sum    = group_binom$n_sum,
    S_2      = groupPsi,
    sig2b    = 1000,
    predX    = predX,
    predS    = predPsi,
    n_preds  = NULL,
    nburn    = nburn,
    nsim     = nsim,
    nthin    = nthin,
    sig2t    = 5,
    tau_1    = 1,
    a_eps    = 0.1,
    b_eps    = 0.1,
    aeta     = 0.1,
    beta     = 0.1,
    alambda  = 2,
    blambda  = 1
  )

  results_ug_b <- gaus_post(
    preds      = unis_wage_b$Preds,
    sig2chain  = unis_wage_b$sig2.chain,
    true_mean  = truth$INCO,
    region     = pcells$PUMA,
    popsize    = pcells$popsize
  )

  results_mg_b <- gaus_post(
    preds      = mult_pov_b$preds_gaus.chain,
    sig2chain  = mult_pov_b$sig2.chain,
    true_mean  = truth$INCO,
    region     = pcells$PUMA,
    popsize    = pcells$popsize
  )

  results_ub_b <- bios_post(
    preds      = unis_pov_b$Preds,
    true_mean  = truth$POV,
    region     = pcells$PUMA,
    popsize    = pcells$popsize
  )

  results_mb_b <- bios_post(
    preds      = mult_pov_b$preds_bios.chain,
    true_mean  = truth$POV,
    region     = pcells$PUMA,
    popsize    = pcells$popsize
  )

  ubios_pre_b[, k] <- results_ub_b$est
  mbios_pre_b[, k] <- results_mb_b$est
  ugaus_pre_b[, k] <- results_ug_b$est
  mgaus_pre_b[, k] <- results_mg_b$est
  dgaus_pre[, k]   <- compare_df$weighted_means
  dbio_pre[, k]    <- compare_df$weighted_means_POV

  ubios_qual_b[, 1, k] <- results_ub_b$lb
  ubios_qual_b[, 2, k] <- results_ub_b$ub
  mbios_qual_b[, 1, k] <- results_mb_b$lb
  mbios_qual_b[, 2, k] <- results_mb_b$ub
  ugaus_qual_b[, 1, k] <- results_ug_b$lb
  ugaus_qual_b[, 2, k] <- results_ug_b$ub
  mgaus_qual_b[, 1, k] <- results_mg_b$lb
  mgaus_qual_b[, 2, k] <- results_mg_b$ub

  cat("Finished", k, "simulation dataset\n")
}

cr_ub_b <- numeric(n_sim)
cr_mb_b <- numeric(n_sim)
cr_ug_b <- numeric(n_sim)
cr_mg_b <- numeric(n_sim)

is_ub_b <- numeric(n_sim)
is_mb_b <- numeric(n_sim)
is_ug_b <- numeric(n_sim)
is_mg_b <- numeric(n_sim)

mse_ub_b <- numeric(n_sim)
mse_mb_b <- numeric(n_sim)
mse_ug_b <- numeric(n_sim)
mse_mg_b <- numeric(n_sim)
mse_dg   <- numeric(n_sim)
mse_db   <- numeric(n_sim)

for (i in seq_len(n_sim)) {
  cr_ub_b[i] <- mean(
    ubios_qual_b[, 1, i] < truth$POV &
      truth$POV < ubios_qual_b[, 2, i]
  )
  cr_mb_b[i] <- mean(
    mbios_qual_b[, 1, i] < truth$POV &
      truth$POV < mbios_qual_b[, 2, i]
  )
  cr_ug_b[i] <- mean(
    ugaus_qual_b[, 1, i] < truth$INCO &
      truth$INCO < ugaus_qual_b[, 2, i]
  )
  cr_mg_b[i] <- mean(
    mgaus_qual_b[, 1, i] < truth$INCO &
      truth$INCO < mgaus_qual_b[, 2, i]
  )

  is_ub_b[i] <- mean(
    interval_score(
      ubios_qual_b[, 1, i],
      ubios_qual_b[, 2, i],
      truth$POV
    )
  )
  is_mb_b[i] <- mean(
    interval_score(
      mbios_qual_b[, 1, i],
      mbios_qual_b[, 2, i],
      truth$POV
    )
  )
  is_ug_b[i] <- mean(
    interval_score(
      ugaus_qual_b[, 1, i],
      ugaus_qual_b[, 2, i],
      truth$INCO
    )
  )
  is_mg_b[i] <- mean(
    interval_score(
      mgaus_qual_b[, 1, i],
      mgaus_qual_b[, 2, i],
      truth$INCO
    )
  )

  mse_ub_b[i] <- mean((ubios_pre_b[, i] - truth$POV)^2, na.rm = TRUE)
  mse_mb_b[i] <- mean((mbios_pre_b[, i] - truth$POV)^2, na.rm = TRUE)
  mse_ug_b[i] <- mean((ugaus_pre_b[, i] - truth$INCO)^2, na.rm = TRUE)
  mse_mg_b[i] <- mean((mgaus_pre_b[, i] - truth$INCO)^2, na.rm = TRUE)
  mse_dg[i]   <- mean((dgaus_pre[, i] - truth$INCO)^2, na.rm = TRUE)
  mse_db[i]   <- mean((dbio_pre[, i] - truth$POV)^2, na.rm = TRUE)
}

mse_ratio_mg_b <- mse_mg_b / mse_ug_b
mse_ratio_mb_b <- mse_mb_b / mse_ub_b
mse_ratio_ug_b <- mse_ug_b / mse_dg
mse_ratio_ub_b <- mse_ub_b / mse_db

is_ratio_mg_b <- is_mg_b / is_ug_b
is_ratio_mb_b <- is_mb_b / is_ub_b

df_mse <- dplyr::bind_rows(
  data.frame(Ratio = mse_ratio_mg_b, Comparison = "Multi / Univariate", Type = "Gaussian"),
  data.frame(Ratio = mse_ratio_mb_b, Comparison = "Multi / Univariate", Type = "Binomial"),
  data.frame(Ratio = mse_ratio_ug_b, Comparison = "Univariate / Direct", Type = "Gaussian"),
  data.frame(Ratio = mse_ratio_ub_b, Comparison = "Univariate / Direct", Type = "Binomial")
)

df_mse$Comparison <- factor(
  df_mse$Comparison,
  levels = c("Multi / Univariate", "Univariate / Direct")
)

df_is <- dplyr::bind_rows(
  data.frame(Ratio = is_ratio_mg_b, Comparison = "Multi / Univariate", Type = "Gaussian"),
  data.frame(Ratio = is_ratio_mb_b, Comparison = "Multi / Univariate", Type = "Binomial")
)

df_is$Comparison <- factor(
  df_is$Comparison,
  levels = c("Multi / Univariate")
)

df_cr_plot <- dplyr::bind_rows(
  data.frame(CR = cr_mg_b, Method = "Multi-type", Type = "Gaussian"),
  data.frame(CR = cr_mb_b, Method = "Multi-type", Type = "Binomial"),
  data.frame(CR = cr_ug_b, Method = "Univariate", Type = "Gaussian"),
  data.frame(CR = cr_ub_b, Method = "Univariate", Type = "Binomial")
)

df_cr_plot$Method <- factor(
  df_cr_plot$Method,
  levels = c("Multi-type", "Univariate")
)

p_mse <- ggplot(df_mse, aes(x = Comparison, y = Ratio, fill = Type)) +
  geom_hline(yintercept = 1, color = "#E63946", linetype = "dashed", linewidth = 1) +
  geom_boxplot(
    alpha = 0.7,
    width = 0.6,
    outlier.shape = NA,
    position = position_dodge(width = 0.75),
    fatten = 1
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    color = "grey50",
    shape = 16,
    alpha = 0.6,
    size = 0.8
  ) +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Binomial" = "#F4A261")) +
  labs(title = "MSE Ratios", y = "Ratio to Baseline", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 12)
  )

p_is <- ggplot(df_is, aes(x = Comparison, y = Ratio, fill = Type)) +
  geom_hline(yintercept = 1, color = "#E63946", linetype = "dashed", linewidth = 1) +
  geom_boxplot(
    alpha = 0.7,
    width = 0.6,
    outlier.shape = NA,
    position = position_dodge(width = 0.75),
    fatten = 1
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    color = "grey50",
    shape = 16,
    alpha = 0.6,
    size = 0.8
  ) +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Binomial" = "#F4A261")) +
  labs(title = "Interval Score Ratios", y = "Ratio to Baseline", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 12)
  )

p_cr <- ggplot(df_cr_plot, aes(x = Method, y = CR, fill = Type)) +
  geom_hline(yintercept = 0.95, color = "#E63946", linetype = "dashed", linewidth = 1) +
  geom_boxplot(
    alpha = 0.7,
    width = 0.6,
    outlier.shape = NA,
    position = position_dodge(width = 0.75),
    fatten = 1
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    color = "grey50",
    shape = 16,
    alpha = 0.6,
    size = 0.8
  ) +
  coord_cartesian(ylim = c(0.8, 1.0)) +
  scale_y_continuous(breaks = seq(0.8, 1.0, by = 0.05)) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Binomial" = "#F4A261")) +
  labs(title = "Empirical Coverage Rates", y = "Coverage Rate", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(face = "bold", size = 12)
  )

final_plot <- p_mse / (p_is + p_cr) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(final_plot)
save.image("empirical_basis_group.RData")