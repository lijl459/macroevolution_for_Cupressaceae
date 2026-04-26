# ============================================================
# 多元回归分析 (OLS)
# —— 全局 + 生物地理两大板块 + 各生物地理分区 + 统一作图
# ============================================================

# ---------- 加载依赖 ----------
library(readr)
library(dplyr)
library(car)
library(MASS)
library(broom)
library(ggplot2)
library(caret)
library(tidyr)
library(stringr)
library(forcats)
library(scales)
library(ggrepel)
library(corrplot)
library(grid)

# ---------- 参数 ----------
INPUT_CSV <- "grid_100km_env_values_185sp.csv"
OUT_GLOBAL <- "ols_global_coefficients_185sp.csv"
OUT_BIOGEO <- "ols_biogeo_coefficients_185sp.csv"
OUT_REGION <- "ols_region_coefficients_185sp.csv"
OUT_PLOT   <- "ols_all_regions_unified_predictors_185sp.jpg"
OUT_PDF    <- "ols_all_regions_unified_predictors_185sp.pdf"

VIF_CUTOFF <- 5
COR_CUTOFF <- 0.7

# ============================================================
# 通用：变量分类函数（气候 / 地形 / 土壤 / 地表覆盖 / 其他）
# —— 全局、Laurasia/Gondwana、region 以及画图统一使用
# ============================================================
term_to_group <- function(x) {
  case_when(
    # ---- 气候类细分 ----
    str_detect(x, "^BIO([1-9]|1[0-1])$") ~ "Temperature",        # BIO1–BIO11
    str_detect(x, "^BIO(1[2-9])$") |
      str_detect(x, "(?i)^AI$")              ~ "Precipitation",   # BIO12–BIO19 + AI
    
    # ---- 其他类型 ----
    str_detect(x, "(?i)^(ELEV|elevation|slope|topo)")        ~ "Topography",
    str_detect(x, "^LC_|(?i)landcover|consensus")            ~ "Landcover",
    str_detect(x, "(?i)(SAND|CLAY|SILT|CEC|PH|ORG|N|SOIL)")  ~ "Soil",
    TRUE ~ "Other"
  )
}

# ============================================================
# 1️⃣ 数据读取与准备
# ============================================================
dat <- read_csv(INPUT_CSV, show_col_types = FALSE)
cat("✅ 数据读取成功，行数:", nrow(dat), "列数:", ncol(dat), "\n")

# 保证响应变量存在且无 NA
dat <- dat %>% filter(!is.na(species_richness))

# 对物种丰富度取 log
dat$log_richness <- log(dat$species_richness + 1)

# 排除不作为自变量的列
exclude_vars <- c(
  "region", "grid_id", "lon", "lat",
  "species_richness", "mean_div_time", "log_richness",
  "PET", 
  "SAND", "SILT", "CLAY", "ORG_CARBON", "PH_WATER",
  "BULK", "TOTAL_N", "CEC_SOIL", "SWR", "LEC_COND", "ELEC_COND"
)

# 初始自变量集合
predictors <- setdiff(names(dat), exclude_vars)

# 仅保留数值且有方差的列（响应变量 log_richness 也会被保留）
dat <- dat[, sapply(dat, function(x) is.numeric(x) && var(x, na.rm = TRUE) > 0)]

# ============================================================
# 2️⃣ 全局 OLS：先对土地利用变量做 PCA，再去相关 + VIF + stepAIC
# ============================================================

# -------------------------------
# PCA：土地利用变量分析与可视化
# -------------------------------
lc_vars <- grep("^LC_", names(dat), value = TRUE)

if (length(lc_vars) > 1) {
  cat("▶ 对土地利用变量执行 PCA...\n")
  
  lc_mat <- dat[, lc_vars, drop = FALSE]
  
  # 删除零方差变量
  lc_mat <- lc_mat[, sapply(lc_mat, function(x) var(x, na.rm = TRUE) > 0), drop = FALSE]
  
  # 仅用完整行做 PCA
  lc_complete <- complete.cases(lc_mat)
  lc_pca <- prcomp(lc_mat[lc_complete, , drop = FALSE], scale. = TRUE)
  
  # ---- 计算方差解释率 ----
  var_expl <- summary(lc_pca)$importance["Proportion of Variance", 1:2] * 100
  
  # 为所有栅格构建 PC1/PC2 分数（缺失值行设为 NA）
  scores <- matrix(NA, nrow = nrow(dat), ncol = 2)
  scores[lc_complete, ] <- lc_pca$x[, 1:2]
  
  dat$LC_PC1 <- scores[, 1]
  dat$LC_PC2 <- scores[, 2]
  
  cat("✅ 提取 LC_PC1, LC_PC2 (方差解释率: ",
      round(var_expl[1], 1), "%, ",
      round(var_expl[2], 1), "% )\n", sep = "")
  
  # 用 PC 替换原始 LC 变量
  predictors <- c(setdiff(predictors, lc_vars), "LC_PC1", "LC_PC2")
  
  # 导出 PCA 变量贡献
  eig <- lc_pca$sdev[1:2]
  loadings <- lc_pca$rotation[, 1:2]
  corvar <- sweep(loadings, 2, eig, `*`)
  
  df_var <- data.frame(
    var = rownames(corvar),
    PC1 = corvar[, 1],
    PC2 = corvar[, 2]
  )
  df_var$Cos2 <- df_var$PC1^2 + df_var$PC2^2
  
  write.csv(df_var, "PCA_LC_variable_contributions.csv", row.names = FALSE)
  
  # 3. 绘制 PCA 相关圆图
  vexp1 <- round(var_expl[1], 1)
  vexp2 <- round(var_expl[2], 1)
  
  circle <- data.frame(
    x = cos(seq(0, 2*pi, length.out = 200)),
    y = sin(seq(0, 2*pi, length.out = 200))
  )
  
  p_pca <- ggplot() +
    geom_path(data = circle, aes(x = x, y = y),
              linewidth = 0.6, linetype = 2) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
    geom_segment(
      data = df_var,
      aes(x = 0, y = 0, xend = PC1, yend = PC2),
      arrow = arrow(length = unit(0.18, "cm")),
      linewidth = 0.6
    ) +
    geom_text_repel(
      data = df_var,
      aes(x = PC1, y = PC2, label = var),
      max.overlaps = Inf, size = 3.8
    ) +
    coord_equal(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05), expand = TRUE) +
    labs(
      title = "Principal Component Analysis (PCA) of Land-Use Variables",
      subtitle = sprintf("PC1: %.1f%% variance; PC2: %.1f%% variance", vexp1, vexp2),
      x = "PC1",
      y = "PC2"
    ) +
    theme_bw(base_size = 13)
  
  print(p_pca)
  ggsave("PCA_LC_variables_correlation_circle_185sp.jpg", p_pca,
         width = 6.5, height = 6.5, dpi = 300)
  ggsave("PCA_LC_variables_correlation_circle_185sp.pdf", p_pca,
         width = 6.5, height = 6.5)
}

# 备份全局自变量，用于各区域分析
global_dat <- dat
global_predictors <- predictors

# ---- 全局相关性可视化（可选）----
if (length(predictors) >= 2) {
  cor_mat_global <- cor(dat[, predictors], use = "pairwise.complete.obs")
  corrplot::corrplot(cor_mat_global, method = "color")
}

# ---- 1. 去除高度相关变量 ----
cat("▶ 检查并删除全局高度相关变量...\n")
if (length(predictors) >= 2) {
  cor_mat <- cor(as.matrix(dat[, predictors]), use = "pairwise.complete.obs")
  high_cor <- findCorrelation(cor_mat, cutoff = COR_CUTOFF, names = TRUE)
  if (length(high_cor) > 0) {
    cat("⚠️ 删除高度相关变量(|r|>", COR_CUTOFF, "): ",
        paste(high_cor, collapse = ", "), "\n", sep = "")
    predictors <- setdiff(predictors, high_cor)
  } else {
    cat("✅ 全局未检测到高度相关变量 (|r|>", COR_CUTOFF, ")。\n", sep = "")
  }
}

# ---- 2. 标准化自变量 ----
dat_scaled <- dat %>% mutate(across(all_of(predictors), scale))

# ---- 3. 初始全局模型 ----
form_global <- as.formula(paste("log_richness ~", paste(predictors, collapse = " + ")))
mod_global <- lm(form_global, data = dat_scaled)

# ---- 4. VIF 筛选 ----
vif_vals <- car::vif(mod_global)
keep_vars <- names(vif_vals[vif_vals <= VIF_CUTOFF])
cat("VIF 检查完成，移除高 VIF 变量(>", VIF_CUTOFF, "): ",
    paste(setdiff(names(vif_vals), keep_vars), collapse = ", "), "\n", sep = "")

# ---- 5. stepAIC ----
mod_step <- stepAIC(
  lm(as.formula(paste("log_richness ~", paste(keep_vars, collapse = " + "))),
     data = dat_scaled),
  direction = "both", trace = FALSE
)
cat("✅ 全局逐步回归完成。\n")

coef_global <- broom::tidy(mod_step) %>%
  filter(term != "(Intercept)") %>%
  arrange(desc(abs(estimate))) %>%
  rename(std_coef = estimate, p_value = p.value) %>%
  mutate(region = "Global")

write_csv(coef_global, OUT_GLOBAL)
cat("✅ 全局 OLS 系数表已保存：", OUT_GLOBAL, "\n")

# ============================================================
# 3️⃣ 生物地理分区定义（Asia 合并东方+东古北）
#     + 两大板块：Laurasia vs Gondwana
# ============================================================
dat <- dat %>%
  mutate(
    region = case_when(
      lat >= 0 & lon >= 60 ~ "Asia",
      lat >= 23 & lon >= -30 & lon < 60 ~ "West_Palearctic",
      lat >= 0 & lon < -30 ~ "Nearctic",
      lat < 0 & lon < -30 ~ "Neotropical",
      lat < 23 & lat > -35 & lon >= -30 & lon <= 60 ~ "Afrotropical",
      lat < 0 & lon > 110 & lon < 180 ~ "Australasian",
      TRUE ~ "Other"
    ),
    # -------- 生物地理两大板块：劳亚 vs 冈瓦纳 --------
    biogeo = case_when(
      region %in% c("Asia", "West_Palearctic", "Nearctic") ~ "Laurasia",
      region %in% c("Afrotropical", "Neotropical", "Australasian") ~ "Gondwana",
      TRUE ~ "Unclassified"
    )
  )

cat("✅ 区域划分完成：\n")
print(table(dat$region))
cat("✅ Biogeo 划分完成：\n")
print(table(dat$biogeo))

# ============================================================
# 4️⃣ 生物地理两大板块 OLS（Laurasia / Gondwana）
# ============================================================
biogeos <- unique(na.omit(dat$biogeo))
biogeo_results <- data.frame()

for (r in biogeos) {
  if (r == "Unclassified") next
  
  cat("\n\n==============================\n")
  cat("▶ 分析生物地理板块：", r, "\n")
  cat("==============================\n")
  
  dat_r <- dat %>%
    filter(biogeo == r, !is.na(species_richness))
  
  if (nrow(dat_r) < 5) {
    cat("⚠️ 板块", r, "有效样本数过少，跳过。\n")
    next
  }
  
  predictors_r <- global_predictors
  
  # 1) 去相关
  if (length(predictors_r) >= 2) {
    env_mat_r <- as.matrix(dat_r[, predictors_r, drop = FALSE])
    cor_mat_r <- stats::cor(env_mat_r, use = "pairwise.complete.obs")
    high_cor_r <- caret::findCorrelation(cor_mat_r, cutoff = COR_CUTOFF, names = TRUE)
    if (length(high_cor_r) > 0) {
      cat("⚠️ 板块", r, "删除高度相关变量(|r|>", COR_CUTOFF, "): ",
          paste(high_cor_r, collapse = ", "), "\n", sep = "")
      predictors_r <- setdiff(predictors_r, high_cor_r)
    } else {
      cat("✅ 板块", r, "未检测到高度相关变量 (|r|>", COR_CUTOFF, ")。\n", sep = "")
    }
  }
  
  if (length(predictors_r) == 0) {
    cat("⚠️ 板块", r, "去相关后无变量，跳过。\n")
    next
  }
  
  # 2) 标准化
  dat_r_scaled <- dat_r %>%
    mutate(across(all_of(predictors_r), scale))
  
  # 3) 构建初始模型 + 别名检查
  if (length(predictors_r) == 1) {
    form_r <- as.formula(paste("log_richness ~", predictors_r[1]))
    mod_r_init <- lm(form_r, data = dat_r_scaled)
  } else {
    form_r <- as.formula(paste("log_richness ~", paste(predictors_r, collapse = " + ")))
    mod_r_init <- lm(form_r, data = dat_r_scaled)
    
    al <- alias(mod_r_init)
    if (length(al$Complete) > 0) {
      aliased_vars_r <- unique(unlist(al$Complete))
      cat("⚠️ 板块", r, "检测到完全共线变量，已移除：",
          paste(aliased_vars_r, collapse = ", "), "\n")
      predictors_r <- setdiff(predictors_r, aliased_vars_r)
      if (length(predictors_r) == 0) {
        cat("⚠️ 板块", r, "移除别名变量后无可用变量，跳过。\n")
        next
      }
      if (length(predictors_r) == 1) {
        form_r <- as.formula(paste("log_richness ~", predictors_r[1]))
      } else {
        form_r <- as.formula(paste("log_richness ~", paste(predictors_r, collapse = " + ")))
      }
      mod_r_init <- lm(form_r, data = dat_r_scaled)
    }
  }
  
  # 4) VIF 筛选 + 确保每类至少保留一个变量
  if (length(predictors_r) >= 2) {
    vif_vals_r <- tryCatch(car::vif(mod_r_init), error = function(e) e)
    print(vif_vals_r)
    
    if (inherits(vif_vals_r, "error")) {
      cat("⚠️ 板块", r, "VIF 计算失败（可能仍有别名/奇异矩阵），暂不按 VIF 删变量。\n")
      keep_r <- predictors_r
    } else {
      keep_r <- names(vif_vals_r[vif_vals_r <= VIF_CUTOFF])
      dropped <- setdiff(names(vif_vals_r), keep_r)
      
      if (length(dropped) > 0) {
        cat("⚠️ 板块", r, "删除高VIF变量(>", VIF_CUTOFF, "): ",
            paste(dropped, collapse = ", "), "\n", sep = "")
      } else {
        cat("✅ 板块", r, "VIF 均 ≤ ", VIF_CUTOFF, "。\n", sep = "")
      }
      
      # 确保每个变量类别至少保留一个
      groups_all  <- term_to_group(names(vif_vals_r))
      groups_keep <- term_to_group(keep_r)
      missing_groups <- setdiff(unique(groups_all), unique(groups_keep))
      
      if (length(missing_groups) > 0) {
        cat("⚠️ 板块", r, "以下类别被完全删除，尝试补回代表变量：",
            paste(missing_groups, collapse = ", "), "\n")
        for (g in missing_groups) {
          vars_in_group <- names(vif_vals_r)[groups_all == g]
          if (length(vars_in_group) > 1) {
            rep_var <- names(sort(vif_vals_r[vars_in_group]))[1]
          } else {
            rep_var <- vars_in_group
          }
          keep_r <- union(keep_r, rep_var)
        }
      }
      
      if (length(keep_r) == 0) {
        cat("⚠️ 板块", r, "VIF 筛选后无变量，跳过。\n")
        next
      }
    }
  } else {
    keep_r <- predictors_r
  }
  
  # 5) stepAIC 或单变量模型
  if (length(keep_r) == 1) {
    form_final <- as.formula(paste("log_richness ~", keep_r[1]))
    mod_r_step <- lm(form_final, data = dat_r_scaled)
    cat("ℹ️ 板块", r, "仅 1 个自变量，使用单变量 OLS。\n")
  } else {
    form_reduced <- as.formula(paste("log_richness ~", paste(keep_r, collapse = " + ")))
    mod_r_step <- MASS::stepAIC(
      lm(form_reduced, data = dat_r_scaled),
      direction = "both", trace = FALSE
    )
    cat("✅ 板块", r, "逐步回归完成。\n")
  }
  
  # 6) 提取标准化系数
  coef_df_r <- broom::tidy(mod_r_step) %>%
    filter(term != "(Intercept)") %>%
    mutate(region = r) %>%
    rename(std_coef = estimate, p_value = p.value)
  biogeo_results <- bind_rows(biogeo_results, coef_df_r)
}

write_csv(biogeo_results, OUT_BIOGEO)
cat("✅ 生物地理板块 OLS 结果已保存：", OUT_BIOGEO, "\n")

# ============================================================
# 5️⃣ 区域 OLS（与全局一致流程）
# ============================================================
regions <- unique(na.omit(dat$region))
region_results <- data.frame()

for (r in regions) {
  if (r %in% c("Other", "Unclassified")) next
  
  cat("\n\n==============================\n")
  cat("▶ 分析区域：", r, "\n")
  cat("==============================\n")
  
  dat_r <- dat %>%
    filter(region == r, !is.na(species_richness))
  
  if (nrow(dat_r) < 5) {
    cat("⚠️ 区域", r, "有效样本数过少，跳过。\n")
    next
  }
  
  predictors_r <- global_predictors
  
  # 1) 去相关
  if (length(predictors_r) >= 2) {
    env_mat_r <- as.matrix(dat_r[, predictors_r, drop = FALSE])
    cor_mat_r <- stats::cor(env_mat_r, use = "pairwise.complete.obs")
    high_cor_r <- caret::findCorrelation(cor_mat_r, cutoff = COR_CUTOFF, names = TRUE)
    if (length(high_cor_r) > 0) {
      cat("⚠️ 区域", r, "删除高度相关变量(|r|>", COR_CUTOFF, "): ",
          paste(high_cor_r, collapse = ", "), "\n", sep = "")
      predictors_r <- setdiff(predictors_r, high_cor_r)
    } else {
      cat("✅ 区域", r, "未检测到高度相关变量 (|r|>", COR_CUTOFF, ")。\n", sep = "")
    }
  }
  
  if (length(predictors_r) == 0) {
    cat("⚠️ 区域", r, "去相关后无变量，跳过。\n")
    next
  }
  
  # 2) 标准化
  dat_r_scaled <- dat_r %>%
    mutate(across(all_of(predictors_r), scale))
  
  # 3) 初始模型 + 别名检查
  if (length(predictors_r) == 1) {
    form_r <- as.formula(paste("log_richness ~", predictors_r[1]))
    mod_r_init <- lm(form_r, data = dat_r_scaled)
  } else {
    form_r <- as.formula(paste("log_richness ~", paste(predictors_r, collapse = " + ")))
    mod_r_init <- lm(form_r, data = dat_r_scaled)
    
    al <- alias(mod_r_init)
    if (length(al$Complete) > 0) {
      aliased_vars_r <- unique(unlist(al$Complete))
      cat("⚠️ 区域", r, "检测到完全共线变量，已移除：",
          paste(aliased_vars_r, collapse = ", "), "\n")
      predictors_r <- setdiff(predictors_r, aliased_vars_r)
      if (length(predictors_r) == 0) {
        cat("⚠️ 区域", r, "移除别名变量后无可用变量，跳过。\n")
        next
      }
      if (length(predictors_r) == 1) {
        form_r <- as.formula(paste("log_richness ~", predictors_r[1]))
      } else {
        form_r <- as.formula(paste("log_richness ~", paste(predictors_r, collapse = " + ")))
      }
      mod_r_init <- lm(form_r, data = dat_r_scaled)
    }
  }
  
  # 4) VIF + 确保每类至少保留一个变量
  if (length(predictors_r) >= 2) {
    vif_vals_r <- tryCatch(car::vif(mod_r_init), error = function(e) e)
    print(vif_vals_r)
    
    if (inherits(vif_vals_r, "error")) {
      cat("⚠️ 区域", r, "VIF 计算失败（可能仍有别名/奇异矩阵），暂不按 VIF 删变量。\n")
      keep_r <- predictors_r
    } else {
      keep_r <- names(vif_vals_r[vif_vals_r <= VIF_CUTOFF])
      dropped <- setdiff(names(vif_vals_r), keep_r)
      
      if (length(dropped) > 0) {
        cat("⚠️ 区域", r, "删除高VIF变量(>", VIF_CUTOFF, "): ",
            paste(dropped, collapse = ", "), "\n", sep = "")
      } else {
        cat("✅ 区域", r, "VIF 均 ≤ ", VIF_CUTOFF, "。\n", sep = "")
      }
      
      # 确保每类至少保留一个变量
      groups_all  <- term_to_group(names(vif_vals_r))
      groups_keep <- term_to_group(keep_r)
      missing_groups <- setdiff(unique(groups_all), unique(groups_keep))
      
      if (length(missing_groups) > 0) {
        cat("⚠️ 区域", r, "以下类别被完全删除，尝试补回代表变量：",
            paste(missing_groups, collapse = ", "), "\n")
        for (g in missing_groups) {
          vars_in_group <- names(vif_vals_r)[groups_all == g]
          if (length(vars_in_group) > 1) {
            rep_var <- names(sort(vif_vals_r[vars_in_group]))[1]
          } else {
            rep_var <- vars_in_group
          }
          keep_r <- union(keep_r, rep_var)
        }
      }
      
      if (length(keep_r) == 0) {
        cat("⚠️ 区域", r, "VIF 筛选后无变量，跳过。\n")
        next
      }
    }
  } else {
    keep_r <- predictors_r
  }
  
  # 5) stepAIC 或单变量 OLS
  if (length(keep_r) == 1) {
    form_final <- as.formula(paste("log_richness ~", keep_r[1]))
    mod_r_step <- lm(form_final, data = dat_r_scaled)
    cat("ℹ️ 区域", r, "仅 1 个自变量，使用单变量 OLS。\n")
  } else {
    form_reduced <- as.formula(paste("log_richness ~", paste(keep_r, collapse = " + ")))
    mod_r_step <- MASS::stepAIC(
      lm(form_reduced, data = dat_r_scaled),
      direction = "both", trace = FALSE
    )
    cat("✅ 区域", r, "逐步回归完成。\n")
  }
  
  # 6) 提取系数
  coef_df_r <- broom::tidy(mod_r_step) %>%
    filter(term != "(Intercept)") %>%
    mutate(region = r) %>%
    rename(std_coef = estimate, p_value = p.value)
  
  region_results <- bind_rows(region_results, coef_df_r)
}

write_csv(region_results, OUT_REGION)
cat("✅ 区域 OLS 结果已保存：", OUT_REGION, "\n")

# ============================================================
# 6️⃣ 合并 Global + Biogeo + Region & 统一作图
# ============================================================
all_results <- bind_rows(coef_global, biogeo_results, region_results) %>%
  dplyr::select(region, term, std_coef, p_value)




# ---- 统一变量全集 ----
all_terms   <- sort(unique(all_results$term))
all_regions <- unique(all_results$region)

plot_df <- complete(all_results, region = all_regions, term = all_terms)

# ---- 分类 + 符号 ----
plot_df <- plot_df %>%
  mutate(
    group   = term_to_group(term),
    abs_coef = abs(std_coef),
    sign_lab = ifelse(std_coef > 0, "+",
                      ifelse(std_coef < 0, "–", ""))
  )

# 指定 facet 的顺序：Global → Laurasia / Gondwana → 各区域
region_order <- c(
  "Global",
  "Laurasia", "Gondwana",
  "Asia", "Nearctic", "West_Palearctic",
  "Australasian", "Neotropical", "Afrotropical"
)

plot_df$region <- factor(plot_df$region, levels = region_order)

# 颜色方案
group_colors <- c(
  "Temperature"   = "#E64B35FF",  # 温度相关：红橙
  "Precipitation" = "#3C5488FF",  # 降水相关：蓝
  "Topography"    = "#4DBBD5FF",  # 地形
  "Soil"          = "#00A087FF",  # 土壤
  "Landcover"     = "#F39B7FFF",  # 地表覆盖
  "Other"         = "grey70"      # 其他
)

# 为 term 设定一个在各区域中统一的顺序（按类别 + 全局平均 abs_coef 排）
plot_df <- plot_df %>%
  filter(!is.na(std_coef)) %>%
  group_by(term) %>%
  mutate(mean_abs_coef = mean(abs_coef, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(group, desc(mean_abs_coef), term) %>%
  mutate(term = fct_inorder(term))

# 绘图
p_all <- ggplot(plot_df %>% filter(!is.na(std_coef)),
                aes(x = term, y = abs_coef, fill = group)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = sign_lab),
            hjust = -0.25, size = 3.8, fontface = "bold") +
  coord_flip(clip = "off") +
  facet_wrap(~ region, ncol = 3) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(
    x = NULL,
    y = "Standardized coefficient",
    title = "OLS standardized coefficients across Global, Biogeographic Realms and Regions"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 9),
    strip.text  = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(5.5, 20, 5.5, 5.5)
  )

print(p_all)
ggsave(OUT_PLOT, p_all, width = 14, height = 8, dpi = 600)
ggsave(OUT_PDF,  p_all, width = 14, height = 8)

cat("✅ 图像已输出：", OUT_PLOT, "和", OUT_PDF, "\n")
