############################################################
## HRP2 vs qPCR: identify samples with lowest HRP2 residuals
##  - First ANC only
##  - PCR-positive samples
##  - Uses PfHRP2_pg.ml (HRP2) and PfLDH_pg.ml (PfLDH)
############################################################

## ---------- Setup ----------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, ggplot2, broom, ggrepel)

## ---------- Paths ----------
# Change this if needed
setwd("C:/HRP2 analysis/PhD_HRP2_qPCR")

infile <- "ANC_recoded.csv"
outfile_candidates <- "low_HRP2_residual_candidates.csv"

## ---------- 1) Read data and keep first ANC ----------
anc <- read.csv(infile, stringsAsFactors = FALSE)%>%
  filter(!is.na(count_pfhrp2), count_pfhrp2 >= 20)

anc <- anc %>%
  filter(visit == "first ANC")

## Quick sanity check
# str(anc)
# names(anc)

## ---------- 2) Keep PCR-positive samples and relevant columns ----------

pcr_df <- anc %>%
  filter(pcrpos == 1) %>%
  select(
    nida,
    pcrpos,
    density,
    PfHRP2_pg.ml,
    PfLDH_pg.ml,
    pfhrp2_detection,
    pfldh_detection
  ) %>%
  rename(
    hrp2_value  = PfHRP2_pg.ml,
    pfldh_value = PfLDH_pg.ml
  )

## ---------- 3) Define pseudo-LODs and log-transform ----------

# Minimum detected HRP2 among PCR+ to define pseudo-LOD
min_detected_hrp2 <- pcr_df %>%
  filter(
    pfhrp2_detection == "detected",
    !is.na(hrp2_value),
    hrp2_value > 0
  ) %>%
  summarise(min_val = min(hrp2_value)) %>%
  pull(min_val)

if (length(min_detected_hrp2) == 0 || is.na(min_detected_hrp2)) {
  stop("No detected HRP2 values found among PCR+ samples. Cannot define pseudo LOD.")
}

pseudo_lod_hrp2 <- min_detected_hrp2 / 2

# Minimum positive density among PCR+
min_positive_density <- pcr_df %>%
  filter(!is.na(density), density > 0) %>%
  summarise(min_val = min(density)) %>%
  pull(min_val)

if (length(min_positive_density) == 0 || is.na(min_positive_density)) {
  stop("No positive density values found among PCR+ samples. Cannot define pseudo LOD for density.")
}

pseudo_lod_density <- min_positive_density / 2

# Add model-ready variables
pcr_df <- pcr_df %>%
  mutate(
    hrp2_for_model = dplyr::case_when(
      pfhrp2_detection == "detected" & !is.na(hrp2_value) & hrp2_value > 0 ~ hrp2_value,
      TRUE ~ pseudo_lod_hrp2   # non-detects / zeros -> LOD/2
    ),
    density_for_model = dplyr::case_when(
      !is.na(density) & density > 0 ~ density,
      TRUE ~ pseudo_lod_density
    ),
    log_hrp2    = log10(hrp2_for_model),
    log_density = log10(density_for_model)
  )

## ---------- 4) Fit log–log regression and get residuals ----------

model_hrp2_density <- lm(log_hrp2 ~ log_density, data = pcr_df)

pcr_aug <- augment(model_hrp2_density, pcr_df)
# pcr_aug includes:
#  - original columns
#  - .fitted (expected log_hrp2)
#  - .resid  (observed - expected)

## ---------- 5) Identify lowest-residual samples ----------

# How many extreme low-HRP2 samples you want to inspect
k <- 20   # change to 30, 50, etc., if you want

low_hrp2_dataset <- pcr_aug %>%
  arrange(.resid) %>%       # most negative residuals first
  slice(1:k) %>%            # take k worst
  select(
    nida,
    pcrpos,
    density,
    hrp2_value,
    pfhrp2_detection,
    pfldh_value,
    pfldh_detection,
    hrp2_for_model,
    density_for_model,
    log_hrp2,
    log_density,
    .fitted,
    .resid
  )

print(low_hrp2_dataset)

# Just the NIDAs if you want to copy/paste somewhere
low_hrp2_nidas <- low_hrp2_dataset %>% pull(nida)
print(low_hrp2_nidas)

## ---------- 6) Plot HRP2 vs density with lowest residuals highlighted ----------
lowest_residuals <- ggplot(pcr_aug, aes(x = density_for_model, y = hrp2_for_model)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(
    data = low_hrp2_dataset,
    aes(x = density_for_model, y = hrp2_for_model),
    color = "red",
    size = 2
  ) +
  geom_text_repel(
    data = low_hrp2_dataset,
    aes(x = density_for_model, y = hrp2_for_model, label = nida),
    size = 4
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "qPCR parasite density (log10, parasites/uL)",
    y = "PfHRP2 concentration (log10, pg/mL)",
    title = "PCR-positive samples: HRP2 vs qPCR density",
    subtitle = "Red points = lowest HRP2 residuals (possible hrp2 deletions / low expression)"
  ) +
  theme_bw() +
  theme(
    # keep title as-is (don’t set plot.title)
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 14),
    plot.subtitle = element_text(size = 14),
    legend.title  = element_text(size = 14),
    legend.text   = element_text(size = 13),
    strip.text    = element_text(size = 14)  # harmless if no facets
  )

lowest_residuals

ggsave(
  "C:/HRP2 analysis/PhD_HRP2_qPCR/hrp23_deletion candidates.png",
  lowest_residuals,
  width = 14.0, height = 8, dpi = 400
)


## ---------- 7) Save candidate dataset ----------

write.csv(
  low_hrp2_dataset,
  file = outfile_candidates,
  row.names = FALSE
)

cat("Saved lowest-residual candidates to:", outfile_candidates, "\n")
