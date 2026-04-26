# =====================================================
# 过滤坐标精度低的点（例如没有小数或仅一位小数）
# 读取 CSV → 过滤 → 写出新 CSV
# =====================================================

# 安装一次即可
# install.packages(c("dplyr", "readr", "stringr"))

library(dplyr)
library(readr)
library(stringr)

# ------------ 参数设置 ------------
setwd("F:/实验室/project/Cupressaceae/taxa204/10.location/整理/")

INPUT_CSV  <- "Cupressaceae_171sp_locations.csv"     # 原始文件
OUTPUT_CSV <- "Cupressaceae_171sp_locations_filtered.csv"  # 输出文件

# 列名（可按实际情况修改）
COL_SPECIES <- "species"
COL_LON     <- "lon"
COL_LAT     <- "lat"

# ------------ 读入数据 ------------
dat_raw <- read_csv(INPUT_CSV, show_col_types = FALSE)

# ------------ 计算小数位数并过滤 ------------
dat_filtered <- dat_raw %>%
  rename(
    species = all_of(COL_SPECIES),
    lon     = all_of(COL_LON),
    lat     = all_of(COL_LAT)
  ) %>%
  mutate(
    lon = as.character(lon),
    lat = as.character(lat)
  ) %>%
  # 计算小数位数
  mutate(
    lon_decimals = ifelse(grepl("\\.", lon),
                          nchar(sub("^[^.]*\\.", "", lon)), 0),
    lat_decimals = ifelse(grepl("\\.", lat),
                          nchar(sub("^[^.]*\\.", "", lat)), 0)
  ) %>%
  # 筛选：经纬度都至少有两位小数
  filter(lon_decimals >= 2, lat_decimals >= 2) %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    species = str_squish(as.character(species))
  ) %>%
  select(species, lon, lat)

# ------------ 输出结果 ------------
write_csv(dat_filtered, OUTPUT_CSV)

cat("✅ 原始数据点数：", nrow(dat_raw),
    "\n✅ 过滤后剩余：", nrow(dat_filtered),
    "\n✅ 已输出文件：", OUTPUT_CSV, "\n")
