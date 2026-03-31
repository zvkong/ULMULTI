suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(sampling)
  library(ggplot2)
  library(patchwork)
  library(tibble)
})

source("functions_indicator_nongroup_revised.R")

pums21 <- readr::read_csv("IL21.csv", show_col_types = FALSE)

pums <- pums21 |>
  dplyr::select(PUMA, PWGTP, SEX, POVPIP, PINCP, SCHL) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    PWGTP = as.numeric(PWGTP),
    POVPIP = as.numeric(POVPIP),
    SEX = factor(SEX),
    LOGINCO = log(as.numeric(PINCP)),
    POV = dplyr::if_else(POVPIP <= 100, 1, 0),
    BACH = factor(dplyr::if_else(SCHL >= 21, 1, 0))
  ) |>
  dplyr::filter(is.finite(LOGINCO)) |>
  dplyr::mutate(
    INCO = (LOGINCO - min(LOGINCO)) / (max(LOGINCO) - min(LOGINCO))
  )

truth <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(
    INCO = mean(INCO),
    POV = mean(POV),
    .groups = "drop"
  ) |>
  dplyr::arrange(PUMA)

area_levels <- truth$PUMA

pcells <- pums |>
  dplyr::group_by(PUMA, SEX, BACH) |>
  dplyr::summarise(popsize = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(factor(PUMA, levels = area_levels), SEX, BACH)

predX <- model.matrix(~ SEX + BACH - 1, data = pcells)
predPsi <- model.matrix(~ factor(PUMA, levels = area_levels) - 1, data = pcells)

n_sim <- 100
nsim <- 500
nburn <- 500
nthin <- 1
sample_size <- 1000

n_area <- nrow(truth)

ubios_pre <- array(NA_real_, dim = c(n_area, n_sim))
mbios_pre <- array(NA_real_, dim = c(n_area, n_sim))
ugaus_pre <- array(NA_real_, dim = c(n_area, n_sim))
mgaus_pre <- array(NA_real_, dim = c(n_area, n_sim))
dgaus_pre <- array(NA_real_, dim = c(n_area, n_sim))
dbio_pre <- array(NA_real_, dim = c(n_area, n_sim))

ubios_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))
mbios_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))
ugaus_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))
mgaus_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))

cor_set <- numeric(n_sim)

for (k in seq_len(n_sim)) {
  set.seed(k)

  prob <- sampling::inclusionprobabilities(
    pums$PWGTP * (1 + 5 * (pums$POV == 1)),
    sample_size
  )

  Ind <- sampling::UPsystematic(prob)
  samp <- pums[as.logical(Ind), , drop = FALSE]
  samp$P <- prob[as.logical(Ind)]
  samp$W <- 1 / samp$P
  samp$scaledWGT <- samp$W * nrow(samp) / sum(samp$W)

  compare_df <- samp |>
    dplyr::group_by(PUMA) |>
    dplyr::summarise(
      weighted_means = stats::weighted.mean(INCO, W),
      weighted_means_POV = stats::weighted.mean(POV, W),
      .groups = "drop"
    )

  dgaus_pre[, k] <- compare_df$weighted_means[match(area_levels, compare_df$PUMA)]
  dbio_pre[, k] <- compare_df$weighted_means_POV[match(area_levels, compare_df$PUMA)]

  cor_set[k] <- stats::cor(samp$POV, samp$INCO)

  modwgt <- samp$scaledWGT
  modX <- model.matrix(~ SEX + BACH - 1, data = samp)
  modPsi <- model.matrix(~ factor(PUMA, levels = area_levels) - 1, data = samp)
  modY <- samp$INCO
  modZ <- samp$POV

  fit_ug <- unis_gaus(
    X = modX,
    Y = modY,
    S = modPsi,
    sig2b = 1000,
    wgt = modwgt,
    n = NULL,
    predX = predX,
    predS = predPsi,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    a = 0.5,
    b = 0.5,
    a_eps = 0.1,
    b_eps = 0.1
  )

  fit_ub <- unis_bios(
    X = modX,
    Y = modZ,
    S = modPsi,
    sig2b = 1000,
    wgt = modwgt,
    n = NULL,
    predX = predX,
    predS = predPsi,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    a = 0.1,
    b = 0.1
  )

  fit_mt <- MTSM_br(
    X_1 = modX,
    X_2 = modX,
    Z_1 = modY,
    Z_2 = modZ,
    S = modPsi,
    sig2b = 1000,
    wgt = modwgt,
    n = NULL,
    predX = predX,
    predS = predPsi,
    n_preds = NULL,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    sig2t = 5,
    sig2e = 10,
    tau_1_init = 1,
    tau_2_init = -0.1,
    a_eps = 0.1,
    b_eps = 0.1,
    aeta = 0.1,
    beta = 0.1,
    alambda = 2,
    blambda = 1
  )

  res_ug <- gaus_post(
    preds = fit_ug$Preds,
    sig2chain = fit_ug$sig2.chain,
    true_mean = truth$INCO,
    region = pcells$PUMA,
    popsize = pcells$popsize
  )

  res_mg <- gaus_post(
    preds = fit_mt$preds_gaus.chain,
    sig2chain = fit_mt$sig2.chain,
    true_mean = truth$INCO,
    region = pcells$PUMA,
    popsize = pcells$popsize
  )

  res_ub <- bios_post(
    preds = fit_ub$Preds,
    true_mean = truth$POV,
    region = pcells$PUMA,
    popsize = pcells$popsize
  )

  res_mb <- bios_post(
    preds = fit_mt$preds_bios.chain,
    true_mean = truth$POV,
    region = pcells$PUMA,
    popsize = pcells$popsize
  )

  ugaus_pre[, k] <- res_ug$est
  mgaus_pre[, k] <- res_mg$est
  ubios_pre[, k] <- res_ub$est
  mbios_pre[, k] <- res_mb$est

  ugaus_qual[, 1, k] <- res_ug$lb
  ugaus_qual[, 2, k] <- res_ug$ub
  mgaus_qual[, 1, k] <- res_mg$lb
  mgaus_qual[, 2, k] <- res_mg$ub
  ubios_qual[, 1, k] <- res_ub$lb
  ubios_qual[, 2, k] <- res_ub$ub
  mbios_qual[, 1, k] <- res_mb$lb
  mbios_qual[, 2, k] <- res_mb$ub

  cat("Finished", k, "simulation dataset\n")
}

cr_ub <- cr_mb <- cr_ug <- cr_mg <- numeric(n_sim)
is_ub <- is_mb <- is_ug <- is_mg <- numeric(n_sim)
mse_ub <- mse_mb <- mse_ug <- mse_mg <- mse_dg <- mse_db <- numeric(n_sim)

for (i in seq_len(n_sim)) {
  cr_ub[i] <- mean(ubios_qual[, 1, i] < truth$POV & truth$POV < ubios_qual[, 2, i], na.rm = TRUE)
  cr_mb[i] <- mean(mbios_qual[, 1, i] < truth$POV & truth$POV < mbios_qual[, 2, i], na.rm = TRUE)
  cr_ug[i] <- mean(ugaus_qual[, 1, i] < truth$INCO & truth$INCO < ugaus_qual[, 2, i], na.rm = TRUE)
  cr_mg[i] <- mean(mgaus_qual[, 1, i] < truth$INCO & truth$INCO < mgaus_qual[, 2, i], na.rm = TRUE)

  is_ub[i] <- mean(interval_score(ubios_qual[, 1, i], ubios_qual[, 2, i], truth$POV), na.rm = TRUE)
  is_mb[i] <- mean(interval_score(mbios_qual[, 1, i], mbios_qual[, 2, i], truth$POV), na.rm = TRUE)
  is_ug[i] <- mean(interval_score(ugaus_qual[, 1, i], ugaus_qual[, 2, i], truth$INCO), na.rm = TRUE)
  is_mg[i] <- mean(interval_score(mgaus_qual[, 1, i], mgaus_qual[, 2, i], truth$INCO), na.rm = TRUE)

  mse_ub[i] <- mean((ubios_pre[, i] - truth$POV)^2, na.rm = TRUE)
  mse_mb[i] <- mean((mbios_pre[, i] - truth$POV)^2, na.rm = TRUE)
  mse_ug[i] <- mean((ugaus_pre[, i] - truth$INCO)^2, na.rm = TRUE)
  mse_mg[i] <- mean((mgaus_pre[, i] - truth$INCO)^2, na.rm = TRUE)
  mse_db[i] <- mean((dbio_pre[, i] - truth$POV)^2, na.rm = TRUE)
  mse_dg[i] <- mean((dgaus_pre[, i] - truth$INCO)^2, na.rm = TRUE)
}

summary_table <- tibble::tibble(
  `Response / Model` = c(
    "Bernoulli, HT",
    "Bernoulli, Univariate",
    "Bernoulli, Multi-type",
    "Gaussian, HT",
    "Gaussian, Univariate",
    "Gaussian, Multi-type"
  ),
  MSE = c(
    mean(mse_db, na.rm = TRUE),
    mean(mse_ub, na.rm = TRUE),
    mean(mse_mb, na.rm = TRUE),
    mean(mse_dg, na.rm = TRUE),
    mean(mse_ug, na.rm = TRUE),
    mean(mse_mg, na.rm = TRUE)
  ),
  Ratio_to_direct = c(
    1,
    mean(mse_ub / mse_db, na.rm = TRUE),
    mean(mse_mb / mse_db, na.rm = TRUE),
    1,
    mean(mse_ug / mse_dg, na.rm = TRUE),
    mean(mse_mg / mse_dg, na.rm = TRUE)
  ),
  Mean_IS = c(
    NA_real_,
    mean(is_ub, na.rm = TRUE),
    mean(is_mb, na.rm = TRUE),
    NA_real_,
    mean(is_ug, na.rm = TRUE),
    mean(is_mg, na.rm = TRUE)
  )
)

mse_ratio_mg <- mse_mg / mse_ug
mse_ratio_mb <- mse_mb / mse_ub
mse_ratio_ug <- mse_ug / mse_dg
mse_ratio_ub <- mse_ub / mse_db

is_ratio_mg <- is_mg / is_ug
is_ratio_mb <- is_mb / is_ub

df_mse <- dplyr::bind_rows(
  data.frame(Ratio = mse_ratio_mg, Comparison = "Multi / Univariate", Type = "Gaussian"),
  data.frame(Ratio = mse_ratio_mb, Comparison = "Multi / Univariate", Type = "Binomial"),
  data.frame(Ratio = mse_ratio_ug, Comparison = "Univariate / Direct", Type = "Gaussian"),
  data.frame(Ratio = mse_ratio_ub, Comparison = "Univariate / Direct", Type = "Binomial")
)
df_mse$Comparison <- factor(df_mse$Comparison, levels = c("Multi / Univariate", "Univariate / Direct"))

df_is <- dplyr::bind_rows(
  data.frame(Ratio = is_ratio_mg, Comparison = "Multi / Univariate", Type = "Gaussian"),
  data.frame(Ratio = is_ratio_mb, Comparison = "Multi / Univariate", Type = "Binomial")
)
df_is$Comparison <- factor(df_is$Comparison, levels = c("Multi / Univariate"))

df_cr_plot <- dplyr::bind_rows(
  data.frame(CR = cr_mg, Method = "Multi-type", Type = "Gaussian"),
  data.frame(CR = cr_mb, Method = "Multi-type", Type = "Binomial"),
  data.frame(CR = cr_ug, Method = "Univariate", Type = "Gaussian"),
  data.frame(CR = cr_ub, Method = "Univariate", Type = "Binomial")
)
df_cr_plot$Method <- factor(df_cr_plot$Method, levels = c("Multi-type", "Univariate"))

p_mse <- ggplot(df_mse, aes(x = Comparison, y = Ratio, fill = Type)) +
  geom_hline(yintercept = 1, color = "#E63946", linetype = "dashed", linewidth = 1) +
  geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA, position = position_dodge(width = 0.75), fatten = 1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), color = "grey50", shape = 16, alpha = 0.6, size = 0.8) +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Binomial" = "#F4A261")) +
  labs(title = "MSE Ratios", y = "Ratio to Baseline", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(face = "bold", size = 12))

p_is <- ggplot(df_is, aes(x = Comparison, y = Ratio, fill = Type)) +
  geom_hline(yintercept = 1, color = "#E63946", linetype = "dashed", linewidth = 1) +
  geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA, position = position_dodge(width = 0.75), fatten = 1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), color = "grey50", shape = 16, alpha = 0.6, size = 0.8) +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2)) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Binomial" = "#F4A261")) +
  labs(title = "Interval Score Ratios", y = "Ratio to Baseline", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(face = "bold", size = 12))

p_cr <- ggplot(df_cr_plot, aes(x = Method, y = CR, fill = Type)) +
  geom_hline(yintercept = 0.95, color = "#E63946", linetype = "dashed", linewidth = 1) +
  geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA, position = position_dodge(width = 0.75), fatten = 1) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), color = "grey50", shape = 16, alpha = 0.6, size = 0.8) +
  coord_cartesian(ylim = c(0.8, 1.0)) +
  scale_y_continuous(breaks = seq(0.8, 1.0, by = 0.05)) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Binomial" = "#F4A261")) +
  labs(title = "Empirical Coverage Rates", y = "Coverage Rate", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(face = "bold", size = 12))

final_plot <- p_mse / (p_is + p_cr) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(summary_table)
print(final_plot)

saveRDS(
  list(
    summary_table = summary_table,
    mse = list(db = mse_db, ub = mse_ub, mb = mse_mb, dg = mse_dg, ug = mse_ug, mg = mse_mg),
    is = list(ub = is_ub, mb = is_mb, ug = is_ug, mg = is_mg),
    cr = list(ub = cr_ub, mb = cr_mb, ug = cr_ug, mg = cr_mg),
    cor_set = cor_set
  ),
  "empirical_indicator_nongroup_results.rds"
)

ggplot2::ggsave("empirical_indicator_nongroup_plot.png", final_plot, width = 12, height = 6, dpi = 300)
save.image("empirical_indicator_nongroup_results.RData")
