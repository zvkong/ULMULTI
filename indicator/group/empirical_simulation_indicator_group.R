library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(sampling)
library(BayesLogit)
library(Matrix)

source("indicator/group/functions_indicator_group.R")

# pums21 <- readr::read_csv("IL21.csv", show_col_types = FALSE)
pums_vars <- c(
  "PINCP",
  "PWGTP",     # person weight
  "PUMA",      # PUMA code (2010 PUMA for 2021 ACS)
  "AGEP",
  "SEX",
  "RAC1P",
  "SCHL",
  "POVPIP"     # poverty ratio
)

pums21 <- tidycensus::get_pums(
  variables = pums_vars,
  state = "IL",
  year = 2021,
  survey = "acs1",
  recode = TRUE
)
pums <- pums21 |>
  dplyr::select(PUMA, PWGTP, SEX, POVPIP, PINCP, SCHL) |>
  tidyr::drop_na() |>
  dplyr::mutate(
    PUMA = as.character(PUMA),
    PWGTP = as.numeric(PWGTP),
    POVPIP = as.numeric(POVPIP),
    PINCP = as.numeric(PINCP),
    SEX = factor(SEX),
    LOGINCO = log(PINCP),
    POV = dplyr::if_else(POVPIP < 100, 1, 0),
    BACH = factor(dplyr::if_else(as.numeric(SCHL) >= 21, 1, 0))
  ) |>
  dplyr::filter(is.finite(LOGINCO)) |>
  dplyr::select(PUMA, PWGTP, SEX, BACH, POV, LOGINCO) |>
  dplyr::arrange(PUMA, SEX, BACH)

pums$INCO <- (pums$LOGINCO - min(pums$LOGINCO)) /
  (max(pums$LOGINCO) - min(pums$LOGINCO))

truth <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(
    INCO = mean(INCO),
    POV = mean(POV),
    .groups = "drop"
  ) |>
  dplyr::arrange(PUMA)

puma_levels <- truth$PUMA
sex_levels <- levels(pums$SEX)
bach_levels <- levels(pums$BACH)

pcells <- pums |>
  dplyr::mutate(
    PUMA_F = factor(PUMA, levels = puma_levels),
    SEX = factor(SEX, levels = sex_levels),
    BACH = factor(BACH, levels = bach_levels)
  ) |>
  dplyr::group_by(PUMA, PUMA_F, SEX, BACH) |>
  dplyr::summarise(popsize = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(PUMA, SEX, BACH)

area_sizes <- pums |>
  dplyr::group_by(PUMA) |>
  dplyr::summarise(N_area = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(match(PUMA, puma_levels))

predX <- model.matrix(~ SEX + BACH - 1, data = pcells)
predPhi <- model.matrix(~ PUMA_F - 1, data = pcells)

n_sim <- 10
nsim <- 1000
nburn <- 1000
nthin <- 1
n_area <- nrow(truth)

ubios_pre <- array(NA_real_, dim = c(n_area, n_sim))
mbios_pre <- array(NA_real_, dim = c(n_area, n_sim))
ugaus_pre <- array(NA_real_, dim = c(n_area, n_sim))
mgaus_pre <- array(NA_real_, dim = c(n_area, n_sim))
dbio_pre <- array(NA_real_, dim = c(n_area, n_sim))
dgaus_pre <- array(NA_real_, dim = c(n_area, n_sim))

ubios_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))
mbios_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))
ugaus_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))
mgaus_qual <- array(NA_real_, dim = c(n_area, 2, n_sim))

for (k in seq_len(n_sim)) {
  set.seed(k)

  prob <- sampling::inclusionprobabilities(
    pums$PWGTP * (1 + 5 * (pums$POV == 1)),
    1000
  )

  Ind <- sampling::UPsystematic(prob)

  samp <- pums[as.logical(Ind), ] |>
    dplyr::mutate(
      PUMA_F = factor(PUMA, levels = puma_levels),
      SEX = factor(SEX, levels = sex_levels),
      BACH = factor(BACH, levels = bach_levels)
    )

  samp$pi <- prob[as.logical(Ind)]
  samp$W <- 1 / samp$pi
  samp$scaledWGT <- samp$W * nrow(samp) / sum(samp$W)

  direct_df <- samp |>
    dplyr::group_by(PUMA) |>
    dplyr::summarise(
      ht_inco_total = sum(W * INCO),
      ht_pov_total = sum(W * POV),
      .groups = "drop"
    ) |>
    dplyr::right_join(area_sizes, by = "PUMA") |>
    dplyr::mutate(
      ht_inco_total = tidyr::replace_na(ht_inco_total, 0),
      ht_pov_total = tidyr::replace_na(ht_pov_total, 0),
      ht_inco = ht_inco_total / N_area,
      ht_pov = ht_pov_total / N_area
    ) |>
    dplyr::arrange(match(PUMA, puma_levels))

  modX1 <- model.matrix(~ SEX + BACH - 1, data = samp)
  modPhi1 <- model.matrix(~ PUMA_F - 1, data = samp)
  modY <- samp$INCO
  modwgt <- samp$scaledWGT

  grp_dat <- samp |>
    dplyr::group_by(PUMA, PUMA_F, SEX, BACH) |>
    dplyr::summarise(
      y_sum = sum(scaledWGT * POV),
      n_sum = sum(scaledWGT),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      SEX = factor(SEX, levels = sex_levels),
      BACH = factor(BACH, levels = bach_levels),
      PUMA_F = factor(PUMA, levels = puma_levels)
    ) |>
    dplyr::arrange(PUMA, SEX, BACH)

  modX2 <- model.matrix(~ SEX + BACH - 1, data = grp_dat)
  modPhi2 <- model.matrix(~ PUMA_F - 1, data = grp_dat)

  fit_ug <- fit_uni_gaussian_indicator(
    X = modX1,
    Y = modY,
    Phi = modPhi1,
    wgt = modwgt,
    predX = predX,
    predPhi = predPhi,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    sig2b = 1000,
    a_u = 0.5,
    b_u = 0.5,
    a_eps = 0.1,
    b_eps = 0.1
  )

  fit_ub <- fit_uni_binomial_indicator_group(
    X = modX2,
    y_sum = grp_dat$y_sum,
    n_sum = grp_dat$n_sum,
    Phi = modPhi2,
    predX = predX,
    predPhi = predPhi,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    sig2b = 1000,
    a_u = 0.1,
    b_u = 0.1
  )

  fit_mt <- fit_multi_indicator_group(
    X_1 = modX1,
    Z_1 = modY,
    Phi_1 = modPhi1,
    wgt_1 = modwgt,
    X_2 = modX2,
    y_sum = grp_dat$y_sum,
    n_sum = grp_dat$n_sum,
    Phi_2 = modPhi2,
    predX = predX,
    predPhi = predPhi,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    sig2b = 1000,
    sig2t = 10,
    a_eps = 0.1,
    b_eps = 0.1,
    a_eta = 0.1,
    b_eta = 0.1,
    a_zeta = 0.1,
    b_zeta = 0.1
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
  dgaus_pre[, k] <- direct_df$ht_inco
  dbio_pre[, k] <- direct_df$ht_pov

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

cr_ub <- numeric(n_sim)
cr_mb <- numeric(n_sim)
cr_ug <- numeric(n_sim)
cr_mg <- numeric(n_sim)

is_ub <- numeric(n_sim)
is_mb <- numeric(n_sim)
is_ug <- numeric(n_sim)
is_mg <- numeric(n_sim)

mse_ub <- numeric(n_sim)
mse_mb <- numeric(n_sim)
mse_ug <- numeric(n_sim)
mse_mg <- numeric(n_sim)
mse_db <- numeric(n_sim)
mse_dg <- numeric(n_sim)

for (i in seq_len(n_sim)) {
  cr_ub[i] <- mean(ubios_qual[, 1, i] <= truth$POV & truth$POV <= ubios_qual[, 2, i])
  cr_mb[i] <- mean(mbios_qual[, 1, i] <= truth$POV & truth$POV <= mbios_qual[, 2, i])
  cr_ug[i] <- mean(ugaus_qual[, 1, i] <= truth$INCO & truth$INCO <= ugaus_qual[, 2, i])
  cr_mg[i] <- mean(mgaus_qual[, 1, i] <= truth$INCO & truth$INCO <= mgaus_qual[, 2, i])

  is_ub[i] <- mean(interval_score(ubios_qual[, 1, i], ubios_qual[, 2, i], truth$POV))
  is_mb[i] <- mean(interval_score(mbios_qual[, 1, i], mbios_qual[, 2, i], truth$POV))
  is_ug[i] <- mean(interval_score(ugaus_qual[, 1, i], ugaus_qual[, 2, i], truth$INCO))
  is_mg[i] <- mean(interval_score(mgaus_qual[, 1, i], mgaus_qual[, 2, i], truth$INCO))

  mse_ub[i] <- mean((ubios_pre[, i] - truth$POV)^2)
  mse_mb[i] <- mean((mbios_pre[, i] - truth$POV)^2)
  mse_ug[i] <- mean((ugaus_pre[, i] - truth$INCO)^2)
  mse_mg[i] <- mean((mgaus_pre[, i] - truth$INCO)^2)
  mse_db[i] <- mean((dbio_pre[, i] - truth$POV)^2)
  mse_dg[i] <- mean((dgaus_pre[, i] - truth$INCO)^2)
}

save.image("empirical_simulation_indicator_group.RData")

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

df_mse$Comparison <- factor(
  df_mse$Comparison,
  levels = c("Multi / Univariate", "Univariate / Direct")
)

df_is <- dplyr::bind_rows(
  data.frame(Ratio = is_ratio_mg, Comparison = "Multi / Univariate", Type = "Gaussian"),
  data.frame(Ratio = is_ratio_mb, Comparison = "Multi / Univariate", Type = "Binomial")
)

df_is$Comparison <- factor(
  df_is$Comparison,
  levels = c("Multi / Univariate")
)

df_cr <- dplyr::bind_rows(
  data.frame(CR = cr_mg, Method = "Multi-type", Type = "Gaussian"),
  data.frame(CR = cr_mb, Method = "Multi-type", Type = "Binomial"),
  data.frame(CR = cr_ug, Method = "Univariate", Type = "Gaussian"),
  data.frame(CR = cr_ub, Method = "Univariate", Type = "Binomial")
)

df_cr$Method <- factor(df_cr$Method, levels = c("Multi-type", "Univariate"))

p_mse <- ggplot(df_mse, aes(x = Comparison, y = Ratio, fill = Type)) +
  geom_hline(yintercept = 1, color = "#E63946", linetype = "dashed", linewidth = 1) +
  geom_boxplot(
    alpha = 0.7,
    width = 0.6,
    outlier.shape = NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    color = "grey50",
    alpha = 0.6,
    size = 0.8
  ) +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Binomial" = "#F4A261")) +
  labs(title = "MSE Ratios", y = "Ratio to Baseline", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p_is <- ggplot(df_is, aes(x = Comparison, y = Ratio, fill = Type)) +
  geom_hline(yintercept = 1, color = "#E63946", linetype = "dashed", linewidth = 1) +
  geom_boxplot(
    alpha = 0.7,
    width = 0.6,
    outlier.shape = NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    color = "grey50",
    alpha = 0.6,
    size = 0.8
  ) +
  coord_cartesian(ylim = c(0, 1.1)) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Binomial" = "#F4A261")) +
  labs(title = "Interval Score Ratios", y = "Ratio to Baseline", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p_cr <- ggplot(df_cr, aes(x = Method, y = CR, fill = Type)) +
  geom_hline(yintercept = 0.95, color = "#E63946", linetype = "dashed", linewidth = 1) +
  geom_boxplot(
    alpha = 0.7,
    width = 0.6,
    outlier.shape = NA,
    position = position_dodge(width = 0.75)
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    color = "grey50",
    alpha = 0.6,
    size = 0.8
  ) +
  coord_cartesian(ylim = c(0.8, 1.0)) +
  scale_fill_manual(values = c("Gaussian" = "#A8DADC", "Binomial" = "#F4A261")) +
  labs(title = "Empirical Coverage Rates", y = "Coverage Rate", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

final_plot <- p_mse / (p_is + p_cr) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(final_plot)

ggsave("Combined_Metrics_Plot_indicator_group.png", final_plot, width = 12, height = 6, dpi = 300)