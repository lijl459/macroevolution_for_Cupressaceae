

# =====================================================
# 使用本地中国底图绘制两个来源的物种分布图
# 黑色点：dat1（主数据）
# 红色点：dat2（仅当该物种在dat2中存在时叠加）
# 每个物种输出一张 PNG 图片
# =====================================================

# install.packages(c("sf","ggplot2","dplyr","readr"))

library(sf)
library(ggplot2)
library(dplyr)
library(readr)

# ------------ 参数设置 ------------
setwd("F:/实验室/project/Cupressaceae/taxa204/10.location/整理/")

INPUT_CSV1 <- "Cupressaceae_171sp_locations_filtered.csv"            # 主数据
INPUT_CSV2 <- "world_conifer_dataset.csv"  # 可能部分物种有匹配
LOCAL_SHAPEFILE <- "F:/arcGIS/cn_countries/cn_countries.shp"     # 本地底图
OUT_DIR <- "species_maps_cn"                                   # 输出目录

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)

# ------------ 1) 读取两个文件 ------------
dat1 <- read_csv(INPUT_CSV1, show_col_types = FALSE) %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  filter(between(lat, -90, 90), between(lon, -180, 180))

dat2 <- read_csv(INPUT_CSV2, show_col_types = FALSE) %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  filter(between(lat, -90, 90), between(lon, -180, 180))

# ------------ 2) 统一列名 ------------
if (!"species" %in% names(dat1)) {
  if ("species" %in% names(dat1)) dat1 <- rename(dat1, normalized_species = species)
}
if (!"species" %in% names(dat2)) {
  if ("species" %in% names(dat2)) dat2 <- rename(dat2, normalized_species = species)
}

# ------------ 3) 读取底图 ------------
china <- st_read(LOCAL_SHAPEFILE, quiet = TRUE)
cat("✅ 成功读取底图，共", nrow(china), "个面\n")

# ------------ 4) 绘图循环 ------------
species_list <- sort(unique(dat1$species))

for (sp in species_list) {
  sub1 <- dat1 %>% filter(species == sp)
  sub2 <- dat2 %>% filter(species == sp)
  
  if (nrow(sub1) == 0) next  # 没有主数据就跳过
  
  # 动态调整显示范围（左右上下各扩10°）
  lon_range <- range(sub1$lon, na.rm = TRUE)
  lat_range <- range(sub1$lat, na.rm = TRUE)
  lon_range <- lon_range + c(-5, 5)
  lat_range <- lat_range + c(-5, 5)
  
  # ---------- 绘图 ----------
  p <- ggplot() +
    geom_sf(data = china, fill = "grey95", color = "grey70", linewidth = 0.3) +
    geom_point(data = sub1, aes(x = lon, y = lat),
               color = "black", alpha = 0.6, size = 0.7) +
    {
      if (nrow(sub2) > 0)
        geom_point(data = sub2, aes(x = lon, y = lat),
                   color = "red", alpha = 0.8, size = 1.0)
    } +
    coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
    labs(title = paste0(
      sp, " (dat1: ", nrow(sub1),
      if (nrow(sub2) > 0) paste0(", dat2: ", nrow(sub2)) else "",
      ")"
    ),
    x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 12)
  
  outfile <- file.path(OUT_DIR, paste0(sp, ".png"))
  ggsave(outfile, p, width = 6, height = 4, dpi = 300)
  message("✅ Saved map for: ", sp)
}

cat("✅ 所有物种地图已生成！输出路径：", normalizePath(OUT_DIR), "\n")


















