## ============================================================
## PCR+ RESTRICTIONS + COMBINATIONS + FLOWCHART + SCATTER
## + COMBINED FIGURE: Flowchart (A) next to Scatter (B)
##
## - First ANC only
## - Positivity: HRP2 detected, PfLDH detected, PCR positive
## - MAIN regression set: PCR+ AND (HRP2 detected OR PfLDH detected)
## - Summarise ALL combinations among PCR+ (includes "PCR only")
## - Summarise combinations within regression set (excludes "PCR only")
## - Flowchart reflects: First ANC -> PCR+ -> (exclude PCR-only) -> regression set -> combos
## - Scatter: saturated HRP2 has burgundy fill + black outline AND same glyph in legend
## ============================================================

setwd("C:/HRP2 analysis/PhD_HRP2_qPCR")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr, tidyr, tibble,
  ggplot2, scales, grid,
  DiagrammeR, DiagrammeRsvg, rsvg,
  patchwork, png
)

## --------- colours (KEEP burgundy) ----------
burgundy <- "#800020"
blue     <- "#1f78b4"

## --------- input ----------
data <- read.csv("ANC_recoded.csv", stringsAsFactors = FALSE)

## -------------------------------------------------------------------
## Clean key fields (prevents NA propagation + handles trailing spaces)
## -------------------------------------------------------------------
data <- data %>%
  mutate(
    visit_clean            = tolower(trimws(visit)),
    pfhrp2_detection_clean  = tolower(trimws(pfhrp2_detection)),
    pfldh_detection_clean   = tolower(trimws(pfldh_detection)),
    pfhrp2_saturation_clean = tolower(trimws(pfhrp2_saturation)),
    pcrpos_num              = suppressWarnings(as.integer(as.character(pcrpos)))
  )

## Keep only first ANC visits
data_first <- data %>% filter(visit_clean == "first anc")
n_first_anc <- nrow(data_first)

## -------------------------------------------------------------------
## Positivity flags (explicit NA handling)
## -------------------------------------------------------------------
data_first <- data_first %>%
  mutate(
    hrp2_pos = !is.na(pfhrp2_detection_clean) & pfhrp2_detection_clean == "detected",
    pldh_pos = !is.na(pfldh_detection_clean)  & pfldh_detection_clean  == "detected",
    pcr_pos  = !is.na(pcrpos_num)             & pcrpos_num == 1,
    hrp2_saturated = !is.na(pfhrp2_saturation_clean) & pfhrp2_saturation_clean == "yes"
  )

## -------------------------------------------------------------------
## 1) PCR-positive dataset (for "all combinations" table)
## -------------------------------------------------------------------
data_pcrpos <- data_first %>% filter(pcr_pos)
n_pcrpos <- nrow(data_pcrpos)

## -------------------------------------------------------------------
## 2) Regression dataset: PCR+ with ≥1 other detectable marker (HRP2 or PfLDH)
## -------------------------------------------------------------------
data_regression <- data_pcrpos %>% filter(hrp2_pos | pldh_pos)
n_regression <- nrow(data_regression)

cat("\n================ SAMPLE COUNTS ================\n")
cat("First ANC visits:", n_first_anc, "\n")
cat("PCR-positive:", n_pcrpos, "\n")
cat("PCR-positive + ≥1 antigen detected (HRP2 and/or PfLDH):", n_regression, "\n")

## -------------------------------------------------------------------
## 3) Detection combinations among PCR+ (includes PCR-only)
## -------------------------------------------------------------------
data_pcrpos_summary <- data_pcrpos %>%
  mutate(
    combo = case_when(
      hrp2_pos & pldh_pos      ~ "HRP2 + PfLDH + PCR",
      hrp2_pos & !pldh_pos     ~ "HRP2 + PCR",
      !hrp2_pos & pldh_pos     ~ "PfLDH + PCR",
      !hrp2_pos & !pldh_pos    ~ "PCR only",
      TRUE                     ~ "Check/NA"
    )
  )

combo_summary_pcrpos <- data_pcrpos_summary %>%
  count(combo, name = "n") %>%
  tidyr::complete(
    combo = c("PCR only", "HRP2 + PCR", "PfLDH + PCR", "HRP2 + PfLDH + PCR"),
    fill = list(n = 0)
  ) %>%
  mutate(
    total_n = sum(n),
    prop    = n / total_n
  )

cat("\n================ COMBINATIONS AMONG PCR+ ================\n")
print(combo_summary_pcrpos)

sum_boxes_pcrpos <- sum(combo_summary_pcrpos$n, na.rm = TRUE)
if (sum_boxes_pcrpos != n_pcrpos) {
  message("WARNING: sum of combo boxes (", sum_boxes_pcrpos, ") != n_pcrpos (", n_pcrpos, ").")
  message("Rows in 'Check/NA' = ", sum(data_pcrpos_summary$combo == "Check/NA", na.rm = TRUE))
} else {
  message("OK: sum of combo boxes matches n_pcrpos (", n_pcrpos, ").")
}

marker_totals_pcrpos <- tibble::tibble(
  marker     = c("HRP2 detected", "PfLDH detected"),
  n_positive = c(
    sum(data_pcrpos_summary$hrp2_pos, na.rm = TRUE),
    sum(data_pcrpos_summary$pldh_pos, na.rm = TRUE)
  )
) %>%
  mutate(
    total_n = n_pcrpos,
    prop    = n_positive / total_n
  )

cat("\n================ MARKER TOTALS AMONG PCR+ ================\n")
print(marker_totals_pcrpos)

## -------------------------------------------------------------------
## 4) Combinations within regression dataset (PCR+ & ≥1 antigen)
## -------------------------------------------------------------------
data_reg_summary <- data_regression %>%
  mutate(
    combo = case_when(
      hrp2_pos & pldh_pos  ~ "HRP2 + PfLDH + PCR",
      hrp2_pos & !pldh_pos ~ "HRP2 + PCR",
      !hrp2_pos & pldh_pos ~ "PfLDH + PCR",
      TRUE                 ~ "Check/NA"
    )
  )

combo_summary_reg <- data_reg_summary %>%
  count(combo, name = "n") %>%
  tidyr::complete(
    combo = c("HRP2 + PCR", "PfLDH + PCR", "HRP2 + PfLDH + PCR"),
    fill = list(n = 0)
  ) %>%
  mutate(
    total_n = sum(n),
    prop    = n / total_n
  )

cat("\n================ COMBINATIONS IN REGRESSION SET (PCR+ & ≥1 ANTIGEN) ================\n")
print(combo_summary_reg)

sum_boxes_reg <- sum(combo_summary_reg$n, na.rm = TRUE)
if (sum_boxes_reg != n_regression) {
  message("WARNING: sum of regression combo boxes (", sum_boxes_reg, ") != n_regression (", n_regression, ").")
  message("Rows in 'Check/NA' = ", sum(data_reg_summary$combo == "Check/NA", na.rm = TRUE))
} else {
  message("OK: sum of regression combo boxes matches n_regression (", n_regression, ").")
}

## -------------------------------------------------------------------
## 5) Flowchart
## -------------------------------------------------------------------
n_pcr_only   <- combo_summary_pcrpos %>% filter(combo == "PCR only") %>% pull(n)
n_hrp2_pcr   <- combo_summary_reg    %>% filter(combo == "HRP2 + PCR") %>% pull(n)
n_pfldh_pcr  <- combo_summary_reg    %>% filter(combo == "PfLDH + PCR") %>% pull(n)
n_triple_pos <- combo_summary_reg    %>% filter(combo == "HRP2 + PfLDH + PCR") %>% pull(n)

n_pcr_only   <- ifelse(length(n_pcr_only)   == 0, 0, n_pcr_only)
n_hrp2_pcr   <- ifelse(length(n_hrp2_pcr)   == 0, 0, n_hrp2_pcr)
n_pfldh_pcr  <- ifelse(length(n_pfldh_pcr)  == 0, 0, n_pfldh_pcr)
n_triple_pos <- ifelse(length(n_triple_pos) == 0, 0, n_triple_pos)

flow_plot <- DiagrammeR::grViz(sprintf("
digraph flow {
  graph [rankdir = TB]
  node  [shape = box, style = rounded, fontname = Helvetica, fontsize = 10]

  A [label = 'First ANC visits\\n(n = %d)']
  B [label = 'PCR-positive infections\\n(n = %d)']
  X [label = 'Excluded: PCR only\\n(no detectable antigen)\\n(n = %d)']
  C [label = 'Regression dataset\\nPCR+ with ≥1 detectable antigen\\n(n = %d)']

  E1 [label = 'HRP2 + PCR\\n(n = %d)']
  E2 [label = 'PfLDH + PCR\\n(n = %d)']
  E3 [label = 'HRP2 + PfLDH + PCR\\n(n = %d)']

  A -> B
  B -> X
  B -> C

  C -> E1
  C -> E2
  C -> E3
}
",
n_first_anc,
n_pcrpos,
n_pcr_only,
n_regression,
n_hrp2_pcr,
n_pfldh_pcr,
n_triple_pos
))

svg_txt <- DiagrammeRsvg::export_svg(flow_plot)

rsvg::rsvg_png(
  charToRaw(svg_txt),
  file   = "flowchart_PCRpos_regression_profiles.png",
  width  = 1600,
  height = 1200
)

rsvg::rsvg_pdf(
  charToRaw(svg_txt),
  file   = "flowchart_PCRpos_regression_profiles.pdf",
  width  = 11.69,
  height = 8.27
)

cat("\nSaved flowchart:\n",
    "- flowchart_PCRpos_regression_profiles.png\n",
    "- flowchart_PCRpos_regression_profiles.pdf\n", sep = "")

## -------------------------------------------------------------------
## 6) Scatter plot (restricted to regression dataset)
##    Saturated HRP2 shown in legend with SAME glyph:
##    burgundy fill + black outline (shape 21)
## -------------------------------------------------------------------
plot_df <- data_regression %>%
  dplyr::select(density, PfLDH_pg.ml, PfHRP2_pg.ml, hrp2_saturated) %>%
  tidyr::pivot_longer(
    cols = c(PfLDH_pg.ml, PfHRP2_pg.ml),
    names_to = "marker_raw",
    values_to = "pg_ml"
  ) %>%
  dplyr::filter(
    is.finite(density), !is.na(density),
    is.finite(pg_ml), !is.na(pg_ml),
    density > 0, pg_ml > 0
  ) %>%
  dplyr::mutate(
    marker = dplyr::case_when(
      marker_raw == "PfLDH_pg.ml" ~ "PfLDH",
      marker_raw == "PfHRP2_pg.ml" & hrp2_saturated ~ "HRP2 (saturated)",
      marker_raw == "PfHRP2_pg.ml" ~ "HRP2",
      TRUE ~ marker_raw
    ),
    marker = factor(marker, levels = c("PfLDH", "HRP2", "HRP2 (saturated)"))
  )

scatter_plot <- ggplot(plot_df, aes(x = density, y = pg_ml,
                                    color = marker, shape = marker, fill = marker)) +
  geom_point(alpha = 0.6, size = 2.7, stroke = 0.7) +
  
  scale_shape_manual(
    values = c("PfLDH" = 16, "HRP2" = 16, "HRP2 (saturated)" = 21),
    limits = c("PfLDH", "HRP2", "HRP2 (saturated)"),
    drop = FALSE
  ) +
  scale_color_manual(
    values = c("PfLDH" = blue, "HRP2" = burgundy, "HRP2 (saturated)" = "black"),
    limits = c("PfLDH", "HRP2", "HRP2 (saturated)"),
    drop = FALSE
  ) +
  scale_fill_manual(
    values = c("PfLDH" = blue, "HRP2" = burgundy, "HRP2 (saturated)" = burgundy),
    limits = c("PfLDH", "HRP2", "HRP2 (saturated)"),
    drop = FALSE
  ) +
  
  # single legend (shape) with correct saturated glyph
  guides(
    color = "none",
    fill  = "none",
    shape = guide_legend(
      title = "Marker",
      override.aes = list(
        color  = c(blue, burgundy, "black"),
        fill   = c(NA,   NA,       burgundy),
        alpha  = 1,
        size   = 3.6,
        stroke = c(0,    0,        0.9)
      )
    )
  ) +
  
  labs(
    x = "Parasite density by qPCR (parasites/μL)",
    y = "Antigen concentration (pg/mL)",
    title = "Parasite density vs antigen concentrations\n(PCR+ with ≥1 detectable antigen)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.margin    = margin(10, 10, 10, 10),
    legend.text      = element_text(size = 25),
    legend.title     = element_text(size = 20, face = "bold")
  )+
  theme_classic()

ggsave(
  filename = "scatter_plot_density_vs_antigens_saturated.png",
  plot = scatter_plot,
  width = 10,
  height = 8,
  dpi = 600,
  units = "in"
)

cat("\nSaved scatter plot:\n- scatter_plot_density_vs_antigens_saturated.png\n")

## -------------------------------------------------------------------
## 7) Combine Flowchart (A) + Scatter (B) side-by-side
##    Robust: reads the saved PNGs
## -------------------------------------------------------------------
png_to_gg <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  img  <- png::readPNG(path)
  grob <- grid::rasterGrob(img, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), interpolate = TRUE)
  ggplot() +
    annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

f_A <- "flowchart_PCRpos_regression_profiles.png"
f_B <- "scatter_plot_density_vs_antigens_saturated.png"

pA <- png_to_gg(f_A) + labs(tag = "A")
pB <- png_to_gg(f_B) + labs(tag = "B")

fig_AB <- (pA | pB) &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.98)
  )

print(fig_AB)

ggsave("FIG_AB_flowchart_plus_scatter.png", fig_AB, width = 18, height = 8, dpi = 1300)

cat("\nSaved combined figure:\n",
    "- FIG_AB_flowchart_plus_scatter.png\n",
    "- FIG_AB_flowchart_plus_scatter.pdf\n", sep = "")

# -------------------------------------------------------------------
# 8) HRP2 deletion candidates among PCR+ (PCR+ AND no HRP2 detected)
#    (includes both "PCR only" and "PfLDH + PCR")
# -------------------------------------------------------------------
deletion_candidates <- data_pcrpos_summary %>%
  filter(!hrp2_pos) %>%
  mutate(hrp2_deletion_candidate = TRUE) %>%
  select(
    nida, combo, hrp2_deletion_candidate,
    pfhrp2_detection, pfldh_detection, pcrpos, density, PfHRP2_pg.ml, PfLDH_pg.ml,
    everything()
  )

cat("\n================ HRP2-DELETION CANDIDATES (PCR+ & HRP2 not detected) ================\n")
print(deletion_candidates %>% select(nida, combo, density, PfHRP2_pg.ml, PfLDH_pg.ml) %>% head(20))

