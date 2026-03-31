# ==============================================================================
# run_simulation.R
# End-to-End Empirical Simulation for Multi-type Unit-level SAE
# ==============================================================================

# 1. 载入依赖与自定义函数
# ==============================================================================
source("indicator/nongroup/functions_indicator_nongroup.R")
source("indicator/group/functions_indicator_group.R")
source("evaluation_tools.R")

library(tidycensus)
library(dplyr)
library(tidyr)
library(sampling)
library(ggplot2)
library(patchwork)

# 2. 数据获取与清洗 (Data Loading & Preprocessing)
# ==============================================================================
message(">>> Loading and preprocessing ACS PUMS data...")

pums_vars <- c("PINCP", "PWGTP", "PUMA", "AGEP", "SEX", "RAC1P", "SCHL", "POVPIP")

# 获取 2021 年伊利诺伊州微观数据
pums_raw <- tidycensus::get_pums(
  variables = pums_vars,
  state = "IL",
  year = 2021,
  survey = "acs1",
  recode = TRUE
)

# 数据清洗与变量转换
pums_clean <- pums_raw %>%
  dplyr::select(PUMA, PWGTP, SEX, POVPIP, PINCP, SCHL) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(
    PWGTP = as.numeric(PWGTP),
    POVPIP = as.numeric(POVPIP),
    SEX = factor(SEX),
    LOGINCO = log(as.numeric(PINCP)),
    INCOME = as.numeric(PINCP),
    # 构建二值化贫困指标
    POV = dplyr::case_when(POVPIP <= 100 ~ 1, POVPIP > 100 ~ 0),
    # 构建教育程度指标 (Bachelor or higher)
    BACH = factor(dplyr::case_when(as.numeric(SCHL) >= 21 ~ 1, as.numeric(SCHL) < 21 ~ 0))
  ) %>%
  dplyr::filter(LOGINCO != -Inf) %>%
  dplyr::select(PUMA, PWGTP, SEX, POV, LOGINCO, INCOME, BACH) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(
    POV = as.numeric(POV),
    EDU = as.numeric(BACH) - 1
  ) %>%
  dplyr::arrange(PUMA, SEX, BACH)

# 高斯响应变量 (INCO) 的 Min-Max 缩放
min_log <- min(pums_clean$LOGINCO)
max_log <- max(pums_clean$LOGINCO)
pums_clean$INCO <- (pums_clean$LOGINCO - min_log) / (max_log - min_log)

# 3. 构建有限总体与后分层结构 (Finite Population & Poststratification)
# ==============================================================================
message(">>> Constructing finite population and poststratification cells...")

# 计算每个 PUMA 的真实均值 (Target Truth)
truth <- pums_clean %>%
  dplyr::group_by(PUMA) %>%
  dplyr::summarise(
    INCO = mean(INCO),
    POV = mean(POV),
    .groups = "drop"
  ) %>%
  dplyr::arrange(PUMA)

# 构建后分层单元格 (Poststratification Cells)
pcells <- pums_clean %>%
  dplyr::group_by(PUMA, SEX, BACH) %>%
  dplyr::summarise(popsize = dplyr::n(), .groups = "drop")

# 提取用于预测的设计矩阵
predX <- model.matrix(~ SEX + BACH - 1, data = pcells)
puma_levels <- levels(factor(pcells$PUMA))
predPsi <- model.matrix(~ factor(pcells$PUMA, levels = puma_levels) - 1)

# 4. 模拟循环设置 (Simulation Setup)
# ==============================================================================
n_sim <- 5      # 总的 Empirical Datasets 数量
ss <- 1000      # 样本量
nburn <- 1000   # MCMC Burn-in
nsim_mcmc <- 1000 # MCMC 迭代次数

# 存储指标的初始化... (与之前代码保持一致)
mse_dg <- mse_db <- numeric(n_sim)
mse_ug_gp <- mse_mg_gp <- mse_ub_gp <- mse_mb_gp <- numeric(n_sim)
is_ug_gp <- is_mg_gp <- is_ub_gp <- is_mb_gp <- numeric(n_sim)
cr_ug_gp <- cr_mg_gp <- cr_ub_gp <- cr_mb_gp <- numeric(n_sim)

# ==============================================================================
# 5. 核心模拟循环 (Main Simulation Loop with Dual-level Progress)
# ==============================================================================
# 记录整体开始时间
start_time_total <- Sys.time()

for (k in 1:n_sim) {
  # >>> [进度追踪 1：整体 Dataset 进度] <<<
  cat(sprintf("\n================================================================\n"))
  cat(sprintf(">>> OVERALL PROGRESS: Dataset %d / %d (%.1f%%) \n", k, n_sim, (k-1)/n_sim*100))
  cat(sprintf(">>> Current Time: %s \n", format(Sys.time(), "%H:%M:%S")))
  cat(sprintf("================================================================\n"))
  
  set.seed(k)
  
  # --- A. PPS 抽样 & B. 提取直接估计量 & C. 构建模型输入矩阵 ---
  # (这部分代码与上一版完全一致，此处省略以节省空间，直接复制你现有的 A, B, C 部分即可)
  # ... [PPS Sampling & Matrix Construction Code] ...
  
  # 为了示例完整性，这里假定 modX, modY, modPsi 等变量已准备好

  # --- D. 运行 Grouped 模型并显示内部进度 ---
  cat("\n[1/3] Fitting Univariate Gaussian Model...\n")
  fit_uni_g_ng <- unis_gaus(X = modX, Y = modY, S = modPsi, wgt = modwgt, 
                            predX = predX, predS = predPsi, nburn = nburn, nsim = nsim_mcmc)

  cat("\n[2/3] Fitting Univariate Binomial Model (Grouped)...\n")
  fit_uni_b_gp <- unis_bios_grouped(X_grp = modX_grp, y_sum = samp_grp$y_sum, n_sum = samp_grp$n_sum, 
                                    S_grp = modPsi_grp, predX = predX, predS = predPsi, 
                                    nburn = nburn, nsim = nsim_mcmc)

  cat("\n[3/3] Fitting Multi-type Joint Model (Grouped)...\n")
  fit_multi_gp <- MTSM_br_grouped(
    X_1 = modX, Z_1 = modY, S_1 = modPsi, wgt_1 = modwgt,
    X_grp = modX_grp, y_sum = samp_grp$y_sum, n_sum = samp_grp$n_sum, S_grp = modPsi_grp,
    predX = predX, predS = predPsi, nburn = nburn, nsim = nsim_mcmc
  )

  # --- E. Grouped 后分层计算 ---
  cat("\n--> Poststratifying and computing metrics...\n")
  res_ug_gp <- gaus_post(fit_uni_g_ng$Preds, fit_uni_g_ng$sig2, truth$INCO, pcells$PUMA, pcells$popsize) 
  res_ub_gp <- bios_post(fit_uni_b_gp$Preds, truth$POV, pcells$PUMA, pcells$popsize)
  res_mg_gp <- gaus_post(fit_multi_gp$Preds_Gaus, fit_multi_gp$sig2, truth$INCO, pcells$PUMA, pcells$popsize)
  res_mb_gp <- bios_post(fit_multi_gp$Preds_Bios, truth$POV, pcells$PUMA, pcells$popsize)

  # --- F. 指标计算与记录 ---
  # ... [指标计算代码与上一版完全一致] ...
}

# 结束时打印总耗时
end_time_total <- Sys.time()
cat(sprintf("\n================================================================\n"))
cat(">>> SIMULATION COMPLETE!\n")
cat(">>> Total Run Time: ", format(round(end_time_total - start_time_total, 2)), "\n")
cat(sprintf("================================================================\n"))

# 6. 数据可视化 (Visualization - using Grouped results)
# ==============================================================================
message("\n>>> Generating Final Evaluation Plots...")

# 整理绘图用的 DataFrames
df_mse <- data.frame(
  Ratio = c(mse_mg_gp / mse_ug_gp, mse_mb_gp / mse_ub_gp, mse_ug_gp / mse_dg, mse_ub_gp / mse_db),
  Comparison = factor(rep(c("Multi / Univariate", "Univariate / Direct"), each = 2 * n_sim), 
                      levels = c("Multi / Univariate", "Univariate / Direct")),
  Type = rep(c("Gaussian", "Binomial", "Gaussian", "Binomial"), each = n_sim)
)

df_is <- data.frame(
  Ratio = c(is_mg_gp / is_ug_gp, is_mb_gp / is_ub_gp),
  Comparison = factor("Multi / Univariate", levels = c("Multi / Univariate")),
  Type = rep(c("Gaussian", "Binomial"), each = n_sim)
)

df_cr_plot <- data.frame(
  CR = c(cr_mg_gp, cr_mb_gp, cr_ug_gp, cr_ub_gp),
  Method = factor(rep(c("Multi-type", "Univariate"), each = 2 * n_sim), 
                  levels = c("Multi-type", "Univariate")),
  Type = rep(c("Gaussian", "Binomial", "Gaussian", "Binomial"), each = n_sim)
)

# 调用你在 evaluation_tools.R 中的作图函数
final_plot <- generate_plots(df_mse, df_is, df_cr_plot)

# 打印图像
print(final_plot)

# 保存图像 (取消注释以保存)
# ggsave("Simulation_Results_Grouped.png", final_plot, width = 12, height = 6, dpi = 300)
message(">>> Pipeline execution complete.")