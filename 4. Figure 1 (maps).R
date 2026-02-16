########################################################################
#################### Map: Mozambique + ONE zoom box (Maputo Province)
## Change (requested):
##  - NO faceting by year
##  - Plot ONE Maputo zoom panel with MEAN incidence across 2017–2019
##  - Color scale MAX fixed at 500 (values >500 are squished to top)
##  - Add incidence value tags on the points (+ keep US numeric IDs + key)
##
## Keeps:
##  - Rivers (Natural Earth; if download fails, continues without rivers)
##  - True metric scale bar (UTM) in zoom
##  - True lon/lat graticule + degree labels + custom 0.5° breaks
########################################################################

setwd("C:/HRP2 analysis/PhD_HRP2_qPCR")

pkgs <- c(
  "dplyr","sf","ggplot2","viridis","patchwork","ggspatial","ggrepel",
  "stringr","grid","scales","readxl",
  "rnaturalearth","rnaturalearthdata"
)
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install)

library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
library(patchwork)
library(ggspatial)
library(ggrepel)
library(stringr)
library(grid)
library(scales)
library(readxl)
library(rnaturalearth)
library(rnaturalearthdata)

# ---- 1) Data load ----
inc    <- read_excel("df_incidencia_anual_us_GPS_adjust_Maputo_2017_2019.xlsx")
moz_db <- st_read("kontur_boundaries_MZ_20230628.gpkg", quiet = TRUE)

prov_name <- "Maputo Province"
target_years <- c(2017L, 2018L, 2019L)

stopifnot(all(c("US","Distrito","Provincia","Ano",
                "Adjusted_Annual_Incidence_per_1000","Latitude","Longitude") %in% names(inc)))
stopifnot(all(c("name_en","geom") %in% names(moz_db)))

# ---- helper: make valid robustly ----
make_valid_safe <- function(x) {
  if ("st_make_valid" %in% getNamespaceExports("sf")) {
    suppressWarnings(sf::st_make_valid(x))
  } else {
    if (!requireNamespace("lwgeom", quietly = TRUE)) install.packages("lwgeom")
    suppressWarnings(lwgeom::st_make_valid(x))
  }
}

## ---- 2) Polygons as sf; active geometry = geom ----
moz_sf <- if (inherits(moz_db, "sf")) st_set_geometry(moz_db, "geom") else st_as_sf(moz_db, sf_column_name = "geom")
if (is.na(st_crs(moz_sf))) st_crs(moz_sf) <- 4326
moz_sf <- st_transform(moz_sf, 4326)
moz_sf <- make_valid_safe(moz_sf)

## ---- 3) Province outline ----
prov_sf <- moz_sf %>% filter(name_en == prov_name)
if (nrow(prov_sf) == 0) stop("No polygons found for name_en == 'Maputo Province'. Check: sort(unique(moz_sf$name_en)).")

prov_union   <- st_union(prov_sf) %>% make_valid_safe()
prov_outline <- st_as_sf(data.frame(id = 1), geometry = st_sfc(prov_union), crs = st_crs(moz_sf))

## ---- 4) District borders inside Maputo (spatial select + clip) ----
idx <- st_intersects(moz_sf, prov_outline, sparse = FALSE)[, 1]
dist_maputo <- moz_sf[idx, , drop = FALSE]
suppressWarnings({
  dist_maputo <- st_intersection(dist_maputo, prov_union)
  dist_maputo <- make_valid_safe(dist_maputo)
})

## ---- 5) Incidence points (Maputo Province; restrict to 2017–2019) ----
inc_maputo <- inc %>%
  mutate(
    Ano = as.integer(Ano),
    Adjusted_Annual_Incidence_per_1000 = as.numeric(Adjusted_Annual_Incidence_per_1000),
    Latitude  = as.numeric(Latitude),
    Longitude = as.numeric(Longitude),
    US_clean  = str_squish(US)
  ) %>%
  filter(
    Provincia == prov_name,
    !is.na(Ano),
    !is.na(Latitude), !is.na(Longitude),
    !is.na(Adjusted_Annual_Incidence_per_1000)
  )

years_found <- sort(unique(inc_maputo$Ano))
if (!all(target_years %in% years_found)) {
  stop(paste0(
    "Expected years 2017–2019 not all present after filtering. Years found: ",
    paste(years_found, collapse = ", ")
  ))
}

inc_maputo_3y <- inc_maputo %>% filter(Ano %in% target_years)

# ---- 6) Mean incidence across the 3 years (per health centre) ----
inc_mean <- inc_maputo_3y %>%
  group_by(US_clean) %>%
  summarise(
    US        = dplyr::first(US),
    Distrito  = dplyr::first(Distrito),
    Provincia = dplyr::first(Provincia),
    Longitude = mean(Longitude, na.rm = TRUE),
    Latitude  = mean(Latitude,  na.rm = TRUE),
    mean_incidence_per_1000 = mean(Adjusted_Annual_Incidence_per_1000, na.rm = TRUE),
    n_years_nonmissing = n_distinct(Ano[!is.na(Adjusted_Annual_Incidence_per_1000)]),
    .groups = "drop"
  ) %>%
  filter(
    !is.na(Longitude), !is.na(Latitude),
    !is.na(mean_incidence_per_1000)
  )

if (nrow(inc_mean) == 0) stop("No rows in inc_mean after summarising mean incidence across 2017–2019.")

missing_3y <- inc_mean %>% filter(n_years_nonmissing < 3) %>% pull(US_clean)
if (length(missing_3y) > 0) {
  message("Note: some health centres have <3 non-missing year values; mean computed over available years for: ",
          paste(utils::head(missing_3y, 15), collapse = ", "),
          ifelse(length(missing_3y) > 15, " ...", ""))
}

# ---- 7) Stable numeric IDs + key ----
us_key <- inc_mean %>%
  distinct(US_clean) %>%
  arrange(US_clean) %>%
  mutate(us_id = row_number(), US_label = stringr::str_to_title(US_clean))

inc_mean <- inc_mean %>% left_join(us_key, by = "US_clean")

write.csv(us_key, "US_ID_key_MaputoProvince_mean2017_2019.csv", row.names = FALSE)

# incidence tag text
inc_mean <- inc_mean %>%
  mutate(inc_tag = scales::number(mean_incidence_per_1000, accuracy = 1))

## ---- 8) Zoom bbox (padded) + crop layers ----
bb <- st_bbox(prov_outline)
xpad <- as.numeric(bb["xmax"] - bb["xmin"]) * 0.06
ypad <- as.numeric(bb["ymax"] - bb["ymin"]) * 0.06

bbox_zoom <- st_bbox(c(
  xmin = as.numeric(bb["xmin"]) - xpad,
  ymin = as.numeric(bb["ymin"]) - ypad,
  xmax = as.numeric(bb["xmax"]) + xpad,
  ymax = as.numeric(bb["ymax"]) + ypad
), crs = st_crs(prov_outline))

zoom_rect_sf <- st_sf(geometry = st_as_sfc(bbox_zoom))

dist_maputo_c  <- st_crop(dist_maputo,  bbox_zoom)
prov_outline_c <- st_crop(prov_outline, bbox_zoom)

inc_mean_c <- inc_mean %>%
  filter(
    Longitude >= bbox_zoom["xmin"], Longitude <= bbox_zoom["xmax"],
    Latitude  >= bbox_zoom["ymin"], Latitude  <= bbox_zoom["ymax"]
  )

## ---- 9) Project zoom polygons to UTM (for nicer geometry + true scale bar) ----
zoom_crs <- 32736  # UTM 36S
dist_maputo_utm  <- st_transform(dist_maputo_c,  zoom_crs)
prov_outline_utm <- st_transform(prov_outline_c, zoom_crs)

## ---- 9a) Force 0.5° breaks for BOTH lon/lat axes ----
step_deg <- 0.5

lon_breaks <- seq(
  floor(as.numeric(bbox_zoom["xmin"]) / step_deg) * step_deg,
  ceiling(as.numeric(bbox_zoom["xmax"]) / step_deg) * step_deg,
  by = step_deg
)

lat_breaks <- seq(
  floor(as.numeric(bbox_zoom["ymin"]) / step_deg) * step_deg,
  ceiling(as.numeric(bbox_zoom["ymax"]) / step_deg) * step_deg,
  by = step_deg
)

## ---- 9b) Rivers (Natural Earth) -> crop to zoom -> UTM ----
rivers_utm <- NULL
try({
  rivers_ne <- rnaturalearth::ne_download(
    scale = "medium",
    type = "rivers_lake_centerlines",
    category = "physical",
    returnclass = "sf"
  )
  
  rivers_ne <- st_transform(rivers_ne, 4326)
  rivers_zoom <- st_crop(rivers_ne, bbox_zoom)
  
  if (nrow(rivers_zoom) > 0) {
    rivers_utm <- st_transform(rivers_zoom, zoom_crs)
  } else {
    message("Rivers layer loaded but none fall inside the zoom extent; continuing without rivers.")
  }
}, silent = TRUE)

if (is.null(rivers_utm)) {
  message("Could not load rivers (Natural Earth download may have failed). Continuing without rivers.")
}

## ---- 9c) Color scale: fixed 0–500 ----
col_limits <- c(0, 500)
col_breaks <- seq(0, 500, by = 100)

## ---- 10) MAIN MAP ----
p_main <- ggplot() +
  geom_sf(data = moz_sf, fill = "grey97", color = "grey70", linewidth = 0.12) +
  geom_sf(data = prov_outline, fill = NA, color = "grey30", linewidth = 0.5) +
  geom_sf(data = zoom_rect_sf, fill = NA, color = "black", linewidth = 0.7) +
  coord_sf(datum = NA) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text  = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

## ---- 11) ZOOM MAP (single panel) ----
p_zoom_noleg <- ggplot() +
  geom_sf(data = dist_maputo_utm, fill = "grey95", color = "grey55", linewidth = 0.25) +
  geom_sf(data = prov_outline_utm, fill = NA, color = "grey30", linewidth = 0.5)

if (!is.null(rivers_utm)) {
  p_zoom_noleg <- p_zoom_noleg +
    geom_sf(data = rivers_utm, color = "steelblue4", linewidth = 0.35, alpha = 0.85)
}

p_zoom_noleg <- p_zoom_noleg +
  geom_point(
    data = inc_mean_c,
    aes(x = Longitude, y = Latitude, color = mean_incidence_per_1000),
    size = 2.9, alpha = 0.95
  ) +
  
  # numeric ID (links to key)
  ggrepel::geom_text_repel(
    data = inc_mean_c,
    aes(x = Longitude, y = Latitude, label = us_id),
    size = 3.0,
    seed = 1,
    min.segment.length = 0,
    box.padding = 0.25,
    point.padding = 0.15,
    max.overlaps = Inf
  ) +
  coord_sf(
    crs = zoom_crs,
    datum = sf::st_crs(4326),
    default_crs = sf::st_crs(4326),
    expand = FALSE
  ) +
  scale_x_continuous(breaks = lon_breaks) +
  scale_y_continuous(breaks = lat_breaks) +
  ggspatial::annotation_scale(
    location = "bl",
    unit_category = "metric",
    width_hint = 0.25,
    pad_x = unit(0.25, "cm"),
    pad_y = unit(0.20, "cm"),
    text_cex = 0.65,
    line_width = 0.35
  ) +
  scale_color_viridis_c(
    option = "C",
    limits = col_limits,
    breaks = col_breaks,
    oob = squish
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_text(size = 7),
    axis.ticks = element_line(linewidth = 0.2),
    legend.position = "none",
    panel.border = element_rect(fill = NA, linewidth = 0.5),
    panel.grid.major = element_line(linewidth = 0.25),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

## ---- 12) Legend column (incidence legend + compact key) ----
p_leg <- ggplot(inc_mean_c, aes(x = Longitude, y = Latitude, color = mean_incidence_per_1000)) +
  geom_point(alpha = 0) +
  scale_color_viridis_c(
    option = "C",
    limits = col_limits,
    breaks = col_breaks,
    oob = squish,
    name = "Mean annual\nincidence (2017–2019)\nper 1,000"
  ) +
  guides(color = guide_colorbar(
    barheight = unit(4.2, "cm"),
    barwidth  = unit(0.55, "cm"),
    ticks = TRUE
  )) +
  theme_void(base_size = 12) +
  theme(legend.position = "left", plot.margin = margin(0, 0, 0, 0))

key_text <- paste0(us_key$us_id, " = ", us_key$US_label, collapse = "\n")

p_key <- ggplot() +
  labs(title = "Health centre ID") +
  annotate("text", x = 0, y = 1, label = key_text, hjust = 0, vjust = 1, size = 5) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.margin = margin(0, 0, 0, 0)
  )

legend_col <- p_leg / p_key + plot_layout(heights = c(0.58, 0.42))

## ---- 13) Combine layout ----
right_side <- (p_zoom_noleg | legend_col) +
  plot_layout(widths = c(1, 0.33))

final_plot <- p_main + right_side +
  plot_layout(widths = c(0.65, 2.55)) +
  plot_annotation(
    title = "Mozambique + Maputo Province with study health centres",
    subtitle = "Zoom panel: Maputo Province. Points show mean annual incidence across 2017–2019 (cases per 1,000)."
  )

print(final_plot)

## ---- 14) Save ----
ggsave("Mozambique_main_plus_Maputo_zoom_MEAN2017_2019_inc0_500_tags.png",
       final_plot, width = 15, height = 7, dpi = 600)
ggsave("Mozambique_main_plus_Maputo_zoom_MEAN2017_2019_inc0_500_tags.pdf",
       final_plot, width = 15, height = 7)
