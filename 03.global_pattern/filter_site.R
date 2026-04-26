# ==============================================
# 用 Natural Earth 数据：去海上点 + 去城市点
# 输出三个文件：
#   1) *_onland.csv           仅陆地点
#   2) *_onland_nonurban.csv  陆地且非城市（最终保留）
#   3) *_urban.csv            城市内点（被剔除）
# ==============================================

# 一次性安装：
# install.packages(c("sf","lwgeom","dplyr","readr","rnaturalearth","rnaturalearthdata","ggplot2"))

library(sf)
library(lwgeom)
library(dplyr)
library(readr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)

# ------------ 参数 ------------
setwd("F:/实验室/project/Cupressaceae/taxa204/10.location/整理/")
INPUT_CSV <- "Cupressaceae.csv"

# 城市缓冲半径（米）。0 = 不缓冲；建议 500~2000 可减少贴城边误差
URBAN_BUFFER_M <- 2000

# ------------ 读数据 ------------
dat <- read_csv(INPUT_CSV, show_col_types = FALSE) %>%
  filter(between(lon, -180, 180), between(lat, -90, 90)) %>%
  distinct()

pts <- st_as_sf(dat, coords = c("lon","lat"), crs = 4326, remove = FALSE)

# ------------ 去海上点（ne_countries）------------
sf_use_s2(FALSE)  # 关闭 s2，避免 NE 多边形自交报错

land <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid() %>%
  st_union() %>%
  st_buffer(0)

keep_land <- st_intersects(pts, land, sparse = FALSE)[,1]
on_land   <- pts[ keep_land, ]
off_sea   <- pts[!keep_land, ]

# 备份仅陆地点（可选导出）
write_csv(st_drop_geometry(on_land), "Cupressaceae_onland.csv")

# ------------ 去城市点（NE Urban areas）------------
urban <- rnaturalearth::ne_download(
  scale = "medium",
  type = "urban_areas",
  category = "cultural",
  returnclass = "sf"
) %>%
  st_make_valid() %>%
  st_collection_extract("POLYGON") %>%
  st_union() %>%
  st_buffer(0)

# 只取覆盖点分布范围的城市面，加速
bbox <- st_as_sfc(st_bbox(on_land))
urban <- suppressWarnings(st_intersection(urban, bbox))

# 可选：城市缓冲（用米制坐标做缓冲更合理）
if (URBAN_BUFFER_M > 0) {
  urban <- urban %>%
    st_transform(3857) %>%
    st_buffer(URBAN_BUFFER_M) %>%
    st_transform(4326)
}

in_urban <- st_intersects(on_land, urban, sparse = FALSE)[,1]
on_land_nonurban <- on_land[!in_urban, ]  # ✅ 最终保留：陆地且非城市
urban_points     <- on_land[ in_urban, ]  # 被剔除的城市点

sf_use_s2(TRUE)  # 可恢复 s2

# ------------ 导出结果 ------------
write_csv(st_drop_geometry(on_land_nonurban), "Cupressaceae_onland_nonurban.csv")
write_csv(st_drop_geometry(urban_points),     "Cupressaceae_urban.csv")

cat(sprintf("总点数: %d\n陆地点: %d（去海）\n城市点: %d（剔除）\n最终保留: %d（陆地&非城市）\n",
            nrow(pts), nrow(on_land), nrow(urban_points), nrow(on_land_nonurban)))

# ------------ 可选：快速检查图 ------------
p <- ggplot() +
  geom_sf(data = land, fill = "grey95", color = "grey75", linewidth = 0.2) +
  geom_sf(data = pts, alpha = 0.3, size = 0.1, color = "blue") +       # 原始点
  geom_sf(data = on_land_nonurban, alpha = 0.3, size = 0.1, color = "red") +
  geom_sf(data = urban, alpha = 0.3, size = 0.1, color = "orange") +
  labs(title = "Filtered occurrences (Natural Earth): on-land & non-urban",
       subtitle = sprintf("Kept: %d  |  Urban removed: %d", nrow(on_land_nonurban), nrow(urban_points)),
       x = NULL, y = NULL) +
  theme_minimal()
print(p)


# =====================================================
# 7) 物种内去除 5 km 内重复点
# =====================================================

library(sf)
library(dplyr)

# 若未定义 on_land_nonurban，请替换为你要处理的 sf 对象
dat_sf <- on_land_nonurban  

# 1) 投影到米制坐标（Web Mercator, EPSG:3857 或 UTM 坐标更精准）
dat_sf_merc <- st_transform(dat_sf, 3857)

# 2) 定义去重函数（按物种分组）
dedup_species <- function(gdf, dist_km = 5) {
  if (nrow(gdf) <= 1) return(gdf)
  keep <- rep(TRUE, nrow(gdf))
  for (i in seq_len(nrow(gdf))) {
    if (!keep[i]) next
    # 计算该点与其他点的距离
    dists <- as.numeric(st_distance(gdf[i, ], gdf))
    # 将5km内的点（除自身）标记为重复
    dup_idx <- which(dists <= dist_km * 1000 & dists > 0)
    keep[dup_idx] <- FALSE
  }
  gdf[keep, ]
}

# 3) 分物种去重
deduped <- dat_sf_merc %>%
  group_by(species) %>%
  group_modify(~ dedup_species(.x, dist_km = 5)) %>%
  ungroup()

# 4) 转回经纬度坐标
deduped <- st_transform(deduped, 4326)

# 5) 导出结果
write_csv(st_drop_geometry(deduped),
          "Cupressaceae_onland_nonurban_dedup5km.csv")

cat(sprintf("去除前: %d  去除后: %d  删除比例: %.2f%%\n",
            nrow(on_land_nonurban), nrow(deduped),
            100 * (1 - nrow(deduped) / nrow(on_land_nonurban))))



# =====================================================
# 8) 统计每个物种在不同清洗阶段的样点数量
# =====================================================

library(dplyr)

# 检查数据是否存在
stopifnot(exists("pts"), exists("on_land"), exists("on_land_nonurban"), exists("deduped"))

# 各阶段计数
count_raw <- pts %>% 
  st_drop_geometry() %>% 
  count(species, name = "n_raw")

count_land <- on_land %>% 
  st_drop_geometry() %>% 
  count(species, name = "n_onland")

count_nonurban <- on_land_nonurban %>% 
  st_drop_geometry() %>% 
  count(species, name = "n_onland_nonurban")

count_dedup <- deduped %>% 
  st_drop_geometry() %>% 
  count(species, name = "n_final")

# 合并所有阶段
summary_table <- count_raw %>%
  full_join(count_land, by = "species") %>%
  full_join(count_nonurban, by = "species") %>%
  full_join(count_dedup, by = "species") %>%
  mutate(across(starts_with("n_"), ~ replace_na(.x, 0))) %>%
  mutate(
    loss_land = n_raw - n_onland,
    loss_urban = n_onland - n_onland_nonurban,
    loss_dedup = n_onland_nonurban - n_final,
    total_loss = n_raw - n_final,
    perc_retained = round(100 * n_final / n_raw, 1)
  ) %>%
  arrange(desc(n_raw))

# 导出统计结果
write_csv(summary_table, "Cupressaceae_species_cleaning_summary.csv")

# 打印示例
print(head(summary_table, 10))
cat("\n✅ 统计表已导出：Cupressaceae_species_cleaning_summary.csv\n")
