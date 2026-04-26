# =========================
# 100 km × 100 km 物种数热图（sf + ggplot2）
# 去重口径：同一栅格内按 species 去重统计
# 投影：Equal Earth (EPSG:8857)，面积保持，便于按距离建格
# =========================

# 安装一次即可：
# install.packages(c("sf", "dplyr", "readr", "ggplot2", "viridis", "stringr"))

library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(viridis)
library(stringr)

# ------------ 参数设置（按需修改）------------
# ⚠️ Windows 路径建议使用正斜杠 / 或双反斜杠 \\
setwd("F:/实验室/project/Cupressaceae/taxa204/11.gloable_pattern/diversity")

INPUT_CSV       <- "../Cupressaceae_171sp_locations.csv"        # 主数据
LOCAL_SHAPEFILE <- "F:/arcGIS/cn_countries/cn_countries.shp"     # 本地底图


# CSV 列名（按你的文件更改）
COL_SPECIES <- "species"
COL_LON     <- "lon"
COL_LAT     <- "lat"

# 栅格边长（米）
GRID_SIZE_M <- 100000  # 100 km

# ------------ 读入数据 ------------
dat_raw <- read_csv(INPUT_CSV, show_col_types = FALSE)

# 基本清洗：去掉经纬度缺失/非数值
dat <- dat_raw %>%
  rename(
    species = all_of(COL_SPECIES),
    lon     = all_of(COL_LON),
    lat     = all_of(COL_LAT)
  ) %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    species = str_squish(as.character(species))
  ) %>%
  filter(!is.na(lon), !is.na(lat))

# 构建点数据（WGS84，经纬度）
pts_wgs84 <- st_as_sf(dat, coords = c("lon", "lat"), crs = 4326, remove = FALSE)


# ------------ 读取全球底图 ------------
world_sf <- st_read(LOCAL_SHAPEFILE, quiet = TRUE)


# ------------ 投影转换（Equal Earth 等积投影）------------
crs_equal_area <- 8857
pts_proj   <- st_transform(pts_wgs84, crs_equal_area)
world_proj <- NULL
if (!is.null(world_sf)) {
  world_proj <- st_transform(world_sf, crs_equal_area)
}

# ------------ 构建 100km × 100km 栅格 ------------
# 自动使用点数据范围，并向外扩展 1 个格子，防止边界截断
bb <- st_bbox(pts_proj)
grid_extent <- c(
  xmin = bb["xmin"] - GRID_SIZE_M,
  ymin = bb["ymin"] - GRID_SIZE_M,
  xmax = bb["xmax"] + GRID_SIZE_M,
  ymax = bb["ymax"] + GRID_SIZE_M
)

# pts_proj 是你的投影后点（EPSG:8857）
grid <- st_make_grid(
  pts_proj,
  cellsize = c(GRID_SIZE_M, GRID_SIZE_M),
  what = "polygons",
  square = TRUE
)

grid <- st_sf(grid_id = seq_along(grid), geometry = grid)


# ------------ 点与格网空间匹配，计算物种数（去重）------------
idx <- st_intersects(pts_proj, grid)
pt2grid <- sapply(idx, function(x) if (length(x) == 0) NA_integer_ else x[1])
pts_proj$grid_id <- pt2grid

rich_tbl <- pts_proj %>%
  st_drop_geometry() %>%
  filter(!is.na(grid_id)) %>%
  distinct(grid_id, species) %>%
  count(grid_id, name = "species_richness")

grid_rich <- grid %>%
  left_join(rich_tbl, by = "grid_id")



world_proj <- st_read(LOCAL_SHAPEFILE, quiet = TRUE)

# ------------ 绘图（全球底图 + 热图）------------
p <- ggplot() +
  geom_sf(data = world_proj, fill = "grey95", color = "grey60", linewidth = 0.2) +  # 全球底图
  geom_sf(data = grid_rich, aes(fill = species_richness), color = NA) +             # 热图
  #scale_fill_viridis(
  #  option = "plasma",
  #  na.value = NA,
  #  name = "Species\nper 100×100 km"
  #) +
  scale_fill_gradientn(
    colours = c("#4575B4", "#74ADD1", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D73027"),
    #colours = c("#4575B4", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D73027"),
    trans = "sqrt",                # √变换让高值不压低
    na.value = NA,
    name = "Species richness"
  ) +
  guides(fill = guide_colorbar(barheight = unit(60, "pt"))) +
  coord_sf(crs = st_crs(crs_equal_area)) +
  theme_void(base_size = 12) +
  theme(
    legend.position = c(0.25, 0.25),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Global Species Richness (100 km × 100 km)")

print(p)





# ------------ 保存输出 ------------
OUT_JPG <- "richness_100km_global_heatmap_185sp.jpg"           # 导出的热图文件名
OUT_PDF <- "richness_100km_global_heatmap_185sp.pdf"           # 导出的热图文件名

ggsave(OUT_JPG, p, width = 10, height = 6.5, dpi = 300)
ggsave(OUT_PDF, p, width = 10, height = 6.5)

cat("✅ 已输出全球热图：", OUT_PNG, "\n")



# =========================
# 计算每格物种的两两平均分化时间（MPD）并绘制热图
# 直接接在上一步（grid, pts_proj, world_proj, phy 已存在）之后
# =========================

library(ape)
library(dplyr)
library(ggplot2)
library(viridis)
library(scales)

TREE_FILE <- "../Cup171.Ma.tre"
OUT_JPG_MPD <- "mean_pairwise_divergence_100km_heatmap_185sp.jpg"
OUT_PDF_MPD <- "mean_pairwise_divergence_100km_heatmap_185sp.pdf"

cat("▶ 开始计算格内物种两两平均分化时间 (MPD)...\n")

phy <- read.tree(TREE_FILE)
# ---- 1) 构建 cophenetic 距离矩阵（物种两两分化时间）----
cat("  计算物种间分化时间矩阵...\n")
dist_matrix <- cophenetic.phylo(phy)   # 对称矩阵；单位与树分支长度一致

# ---- 2) 点与格网匹配，获取格内物种列表 ----
if (!("grid_id" %in% names(pts_proj))) {
  idx <- st_intersects(pts_proj, grid)
  pts_proj$grid_id <- sapply(idx, function(x) if (length(x)==0) NA_integer_ else x[1])
}

sp_in_cell <- pts_proj |>
  st_drop_geometry() |>
  filter(!is.na(grid_id), !is.na(species)) |>
  distinct(grid_id, species)

# ---- 3) 建立物种名匹配表 ----
norm_space <- function(x) str_squish(gsub("_", " ", x))
norm_under <- function(x) str_squish(gsub(" ", "_", x))

tips_norm_space <- norm_space(phy$tip.label)
tips_norm_under <- norm_under(phy$tip.label)

map_df <- tibble(
  species = as.character(sp_in_cell$species),
  s_space = norm_space(as.character(sp_in_cell$species)),
  s_under = norm_under(as.character(sp_in_cell$species))
) |>
  distinct() |>
  mutate(
    tip_match = case_when(
      s_space %in% tips_norm_space ~ phy$tip.label[match(s_space, tips_norm_space)],
      s_under %in% tips_norm_under ~ phy$tip.label[match(s_under, tips_norm_under)],
      TRUE ~ NA_character_
    )
  )


# 1) 物种间“分化时间矩阵”（而不是路径长度）
dist_matrix <- cophenetic.phylo(phy) / 2

# 2) 逐格 MPD（mean pairwise divergence time）
mean_age_tbl <- sp_in_cell |>
  left_join(map_df, by = "species") |>
  group_by(grid_id) |>
  summarise(
    mean_div_time = {
      tips_here <- unique(na.omit(tip_match))
      if (length(tips_here) < 2) NA_real_ else {
        d <- dist_matrix[tips_here, tips_here]
        mean(d[upper.tri(d)], na.rm = TRUE)
      }
    },
    n_species = n_distinct(species),
    .groups = "drop"
  )

# 树高度（从根到任一叶的距离，单位与树分支长度一致）
tree_height <- max(node.depth.edgelength(phy))

# 下面两个应该大致相近（数值单位一致时）
max_tMRCA <- max(dist_matrix, na.rm = TRUE)      # 已经 /2 了
cat("Tree height ~", tree_height, 
    " | max pairwise tMRCA ~", max_tMRCA, "\n")




# ---- 5) 合并到格网并绘图 ----
grid_mpd <- grid |>
  left_join(mean_age_tbl, by = "grid_id")

grid_plot_mpd <- filter(grid_mpd, !is.na(mean_div_time))
min_val <- min(grid_plot_mpd$mean_div_time, na.rm = TRUE)
max_val <- max(grid_plot_mpd$mean_div_time, na.rm = TRUE)

p_mpd <- ggplot() +
  geom_sf(data = world_proj, fill = "grey95", color = "grey60", linewidth = 0.2) +
  geom_sf(data = grid_plot_mpd, aes(fill = mean_div_time), color = NA) +
  #scale_fill_viridis_c(
  #  option = "plasma",
  #  limits = c(min_val, max_val),
  #  oob = squish,
  #  name = "Mean pairwise divergence"
  #) +
  scale_fill_gradientn(
    colours = c("#4575B4", "#74ADD1", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D73027"),
    #colours = c("#4575B4", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D73027"),
    trans = "sqrt",                # √变换让高值不压低
    na.value = NA,
    name = "Mean Divergence \nTime (Ma)"
  ) +
  guides(fill = guide_colorbar(
    title.position = "top",
    barheight = unit(60, "pt")
  )) +
  coord_sf(crs = st_crs(crs_equal_area)) +
  theme_void(base_size = 12) +
  theme(
    legend.position = c(0.25, 0.25),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Global Mean Pairwise Divergence (100 km × 100 km)")

print(p_mpd)

ggsave(OUT_JPG_MPD, p_mpd, width = 10, height = 6.5, dpi = 300)
ggsave(OUT_PDF_MPD, p_mpd, width = 10, height = 6.5)
cat("✅ 已输出格内平均分化时间（MPD）热图：", OUT_JPG_MPD, "\n")




# ============================================================
# 高效计算全球环境变量的100km格网平均值
# ============================================================

library(terra)
library(sf)
library(dplyr)
library(readr)

# ========== 参数设置 ==========
# （根据你前面的定义自行调整）
setwd("F:/实验室/project/Cupressaceae/taxa204/11.gloable_pattern/climate/")

ENV <- rast("ENV_stack.tif")
GRID_SIZE_M <- 100000                       # 100 km
OUT_RASTER  <- "ENV_mean_100km.tif"
OUT_CSV     <- "grid_100km_env_values_185sp.csv"

# 假设你已有:
#   grid       —— sf 格网 (EPSG:8857)
#   pts_proj   —— 物种分布点 (sf)，带 grid_id 列
#   mean_age_tbl —— 含 grid_id 和 mean_div_time
# ============================================================


# ---- 1️⃣ 聚合环境栅格：按100 km格计算平均值 ----
cat("▶ 开始聚合全球环境栅格 (约100 km 分辨率)...\n")

# 如果数据是等面积投影 (如EPSG:8857)，直接按米聚合；
# 若仍是经纬度(WGS84)，可使用下方近似：100 km ≈ 1°。
res_m <- res(ENV)[1]

fact <- round(GRID_SIZE_M / res_m)  # 聚合倍数
cat("  原始分辨率:", res_m, " → 聚合倍数:", fact, "\n")


# 第一步：先聚合到约10km
ENV_10km <- aggregate(ENV, fact = 10, fun = mean, na.rm = TRUE,
                      filename = "ENV_10km.tif", overwrite = TRUE)

# 第二步：再从10km聚合到100km
ENV_mean <- aggregate(ENV_10km, fact = 10, fun = mean, na.rm = TRUE,
                      filename = "ENV_mean_100km.tif", overwrite = TRUE)


ENV_mean <- rast("F:\\实验室\\project\\Cupressaceae\\taxa204\\11.gloable_pattern\\climate\\ENV_mean_100km.tif")

# ---- 2️⃣ 提取仅有物种分布的格网 ----
cat("▶ 提取有物种分布的格网...\n")
grid_has_sp <- pts_proj %>%
  st_drop_geometry() %>%
  filter(!is.na(grid_id)) %>%
  distinct(grid_id) %>%
  pull(grid_id)

grid_sub <- grid %>% filter(grid_id %in% grid_has_sp)
cat("  格网数:", nrow(grid_sub), "\n")


# ---- 3️⃣ 计算每格中心点经纬度（WGS84）----
grid_cent <- grid_sub %>%
  st_centroid(of_largest_polygon = TRUE) %>%
  st_transform(4326)

cent_xy <- st_coordinates(grid_cent)
cent_df <- tibble(
  grid_id = grid_sub$grid_id,
  lon = cent_xy[, "X"],
  lat = cent_xy[, "Y"]
)


# ---- 4️⃣ 从聚合后的环境栅格提取格内均值 ----
cat("▶ 提取环境因子均值（来自已聚合栅格）...\n")

grid_envcrs <- st_transform(grid_sub, crs(ENV_mean))
vgrid <- vect(grid_envcrs)
env_mean_tbl <- terra::extract(ENV_mean, vgrid, fun = mean, na.rm = TRUE, touches = TRUE)

env_mean_tbl <- env_mean_tbl %>%
  dplyr::rename(row_id = ID) %>%
  dplyr::mutate(grid_id = grid_sub$grid_id[row_id]) %>%
  dplyr::select(-row_id)



# ---- 5️⃣ 整合物种丰富度、分化时间、环境均值 ----
cat("▶ 整合物种丰度与分化时间...\n")

rich_tbl <- pts_proj %>%
  st_drop_geometry() %>%
  filter(!is.na(grid_id)) %>%
  distinct(grid_id, species) %>%
  count(grid_id, name = "species_richness")

mpd_tbl <- mean_age_tbl %>%
  dplyr::select(grid_id, mean_div_time)


grid_env_final <- grid_sub %>%
  st_drop_geometry() %>%
  left_join(cent_df, by = "grid_id") %>%
  left_join(rich_tbl, by = "grid_id") %>%
  left_join(mpd_tbl, by = "grid_id") %>%
  left_join(env_mean_tbl, by = "grid_id") %>%
  arrange(desc(species_richness))


# ---- 6️⃣ 导出结果 ----
OUT_CSV     <- "grid_100km_env_values_185sp.csv"
write_csv(grid_env_final, OUT_CSV)
cat("✅ 已输出CSV文件：", OUT_CSV, "\n")

# 预览前几行
print(head(grid_env_final))










