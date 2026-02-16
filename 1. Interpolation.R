## ============================================================
## qSAT-style calibration with 5PL + LLOD/LLOQ (from pooled blanks, RAW MFI)
## + ULOQ (mean curve max - 10%)
## ROBUST (does NOT rename ANC columns; creates canonical numeric pvldh/pfhrp2/pfldh)
##
## KEY SETTINGS (per your clarification):
##   - Standards: top well IN ASSAY WELL is:
##       PfHRP2 = 100,000 pg/mL  (0.1 µg/mL)
##       PfLDH  = 1,000,000 pg/mL (1 µg/mL)
##       PvLDH  = 1,000,000 pg/mL (1 µg/mL)
##     and standard labels represent 1:100 then 3-fold serial:
##       100, 300, 900, ...  (15 points, Unknown74..88)
##     => conc_std_pgml_assay = top_pgml_assay * (100 / dilution_num)
##
##   - Samples are diluted 1:50 in the assay:
##     => Conc_pgml_sample = Conc_pgml_assay * 50
##
## OUTPUTS (same style as your other pipeline):
##   - GLOBAL_blank_stats_LLOD_LLOQ_rawMFI.csv
##   - MEAN_standard_points_RAWMFI_by_analyte.csv
##   - MEAN_curve_5PL_parameters_RAWMFI_by_analyte.csv
##   - GLOBAL_thresholds_rawMFI_and_pgml_assay_and_sample.csv
##   - ULOQ_meanCurve_by_analyte.csv
##   - ANC_interpolated_long_meanCurve.csv
##   - ANC_nonDetected_by_analyte.csv
##   - ANC_belowLLOQ_by_analyte.csv
##   - ANC_saturated_by_analyte.csv
##   - Final merged ANC with concentrations + flags:
##       C:/HRP2 analysis/ANC_interpolated_5PL_with_detection_and_saturation.csv
##       C:/HRP2 analysis/PhD_HRP2_qPCR/ANC_concentrations.csv
## ============================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr, tidyr, readr, stringr, purrr, ggplot2,
  minpack.lm, tibble
)

## ---------- Paths ----------
input_directory  <- "C:/HRP2 analysis/ag_results-only"  # plates
output_directory <- "C:/HRP2 analysis/"                 # outputs + ANC file
anc_file         <- file.path(output_directory, "ANC_data.csv")

final_dir <- "C:/HRP2 analysis/PhD_HRP2_qPCR"
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)

dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
setwd(input_directory)

qc_dir <- file.path(output_directory, "QC_5PL")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- Analytes ----------
analytes <- c("PvLDH","PfHRP2","PfLDH")

## ---------- Dilutions ----------
STD_TOP_DILUTION <- 100   # standards start at 1:100
SAMPLE_DILUTION  <- 50    # samples diluted 1:50 (report back-calculated concentrations)

## ---------- Top well concentrations in ASSAY WELL (pg/mL) ----------
top_ug_ml <- c(PfHRP2 = 0.1, PfLDH = 1, PvLDH = 1)
top_pg_ml_assay <- top_ug_ml * 1e6  # 1 µg/mL = 1e6 pg/mL

## ============================================================
## Helpers
## ============================================================

msg <- function(...) cat(sprintf(...), "\n")
clean_name <- function(x) tolower(gsub("[^a-z0-9]", "", x))

five_pl <- function(x, A, B, C, D, E) {
  A + (B - A) / (1 + (x / C)^D)^E
}

inv_5pl_vec <- function(MFI, A, B, C, D, E, min_c, max_c, n_grid = 3000) {
  if (anyNA(c(A, B, C, D, E, min_c, max_c))) return(rep(NA_real_, length(MFI)))
  if (!is.finite(min_c) || !is.finite(max_c) || min_c <= 0 || max_c <= 0 || min_c >= max_c) {
    return(rep(NA_real_, length(MFI)))
  }
  
  grid_c  <- 10^seq(log10(min_c), log10(max_c), length.out = n_grid)
  grid_mf <- suppressWarnings(five_pl(grid_c, A, B, C, D, E))
  
  keep <- is.finite(grid_mf) & is.finite(grid_c)
  grid_c  <- grid_c[keep]
  grid_mf <- grid_mf[keep]
  if (length(grid_mf) < 50) return(rep(NA_real_, length(MFI)))
  
  o <- order(grid_mf)
  x <- grid_mf[o]
  y <- grid_c[o]
  
  df <- tibble(x = x, y = y) %>%
    group_by(x) %>%
    summarise(y = mean(y), .groups = "drop") %>%
    arrange(x)
  
  if (nrow(df) < 50) return(rep(NA_real_, length(MFI)))
  
  approx(x = df$x, y = df$y, xout = MFI, rule = 2)$y
}

standardise_plate_names <- function(df) {
  nm_orig  <- names(df)
  nm_clean <- clean_name(nm_orig)
  
  nm_orig[str_detect(nm_clean, "pvldh")]  <- "PvLDH"
  nm_orig[str_detect(nm_clean, "pfhrp2")] <- "PfHRP2"
  nm_orig[str_detect(nm_clean, "pfldh")]  <- "PfLDH"
  
  names(df) <- make.unique(nm_orig)
  df
}

# identify standards by label; fallback to rows 74:88
get_std_idx <- function(df) {
  if (!"Sample" %in% names(df)) {
    if (nrow(df) >= 88) return(74:88)
    return(integer(0))
  }
  s <- tolower(trimws(as.character(df$Sample)))
  
  idx_lab <- which(
    stringr::str_detect(s, "^unknown(74|75|76|77|78|79|80|81|82|83|84|85|86|87|88)$") |
      stringr::str_detect(s, "recombinant_\\s*1\\s*[:/]\\s*[0-9]+") |
      stringr::str_detect(s, "sc_point") |
      stringr::str_detect(s, "^sc$")
  )
  if (length(idx_lab) > 0) return(idx_lab)
  if (nrow(df) >= 88) return(74:88)
  integer(0)
}

# parse dilution denominators (100, 300, 900, ...)
parse_dilution_num <- function(sample_vec) {
  s <- as.character(sample_vec)
  s_trim <- trimws(s)
  
  denom_unknown <- STD_TOP_DILUTION * 3^(0:14)  # 100, 300, ..., 100*3^14
  map_unknown <- setNames(denom_unknown, paste0("Unknown", 74:88))
  
  out <- suppressWarnings(as.numeric(map_unknown[s_trim]))
  
  need <- is.na(out)
  if (any(need)) {
    tmp <- s_trim[need]
    m <- stringr::str_match(tmp, "1\\s*[:/]\\s*([0-9]+)")
    out2 <- suppressWarnings(as.numeric(m[,2]))
    out[need] <- out2
  }
  
  need2 <- is.na(out)
  if (any(need2)) {
    out3 <- suppressWarnings(as.numeric(s_trim[need2]))
    out[need2] <- out3
  }
  
  out
}

extract_blanks <- function(df, plate_name, n_blank_default = 2) {
  analytes_present <- intersect(analytes, names(df))
  if (length(analytes_present) == 0) return(tibble())
  
  out <- NULL
  if ("Sample" %in% names(df)) {
    s <- tolower(as.character(df$Sample))
    is_blank <- str_detect(s, "blank|buffer|assaybuffer|assay buffer")
    if (any(is_blank, na.rm = TRUE)) out <- df %>% dplyr::filter(is_blank)
  }
  if (is.null(out)) {
    if (nrow(df) < n_blank_default) return(tibble())
    out <- df %>% dplyr::slice_tail(n = n_blank_default)
  }
  
  out %>%
    mutate(plate = plate_name) %>%
    select(plate, any_of(analytes_present))
}

fit_5pl_for_group <- function(df_group, xcol, ycol) {
  df_use <- df_group %>%
    filter(is.finite(.data[[xcol]]), is.finite(.data[[ycol]]),
           .data[[xcol]] > 0, .data[[ycol]] > 0) %>%
    arrange(.data[[xcol]])
  
  if (nrow(df_use) < 6) {
    return(tibble(A=NA_real_, B=NA_real_, C=NA_real_, D=NA_real_, E=NA_real_, converged=FALSE))
  }
  
  x <- df_use[[xcol]]
  y <- df_use[[ycol]]
  
  start_list <- list(
    A = min(y, na.rm = TRUE),
    B = max(y, na.rm = TRUE),
    C = stats::median(x, na.rm = TRUE),
    D = 1,
    E = 1
  )
  
  lower <- c(A = 0,   B = 0,   C = min(x, na.rm = TRUE) / 100, D = 0.01, E = 0.01)
  upper <- c(A = Inf, B = Inf, C = max(x, na.rm = TRUE) * 100, D = 20,   E = 20)
  
  fit <- tryCatch(
    minpack.lm::nlsLM(
      y ~ five_pl(x, A, B, C, D, E),
      start   = start_list,
      lower   = lower,
      upper   = upper,
      control = minpack.lm::nls.lm.control(maxiter = 1024, ftol = 1e-10, ptol = 1e-10)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble(A=NA_real_, B=NA_real_, C=NA_real_, D=NA_real_, E=NA_real_, converged=FALSE))
  }
  
  co <- coef(fit)
  tibble(A=unname(co[["A"]]), B=unname(co[["B"]]), C=unname(co[["C"]]),
         D=unname(co[["D"]]), E=unname(co[["E"]]), converged=TRUE)
}

compute_uloq <- function(A, B, C, D, E, min_c, max_c, n_grid = 4000) {
  grid_c  <- 10^seq(log10(min_c), log10(max_c), length.out = n_grid)
  grid_mf <- suppressWarnings(five_pl(grid_c, A, B, C, D, E))
  grid_mf <- grid_mf[is.finite(grid_mf)]
  if (length(grid_mf) < 50) return(tibble(ULOQ_MFI = NA_real_, ULOQ_pgml_assay = NA_real_))
  
  max_mfi  <- max(grid_mf, na.rm = TRUE)
  uloq_mfi <- 0.9 * max_mfi
  uloq_pg  <- inv_5pl_vec(uloq_mfi, A, B, C, D, E, min_c, max_c)
  
  tibble(ULOQ_MFI = uloq_mfi, ULOQ_pgml_assay = uloq_pg)
}

## --- ANC: pick MFI columns safely (NO RENAMING) ---
pick_best_mfi_col <- function(df, target_clean) {
  nm <- names(df)
  nm_clean <- clean_name(nm)
  
  exclude <- str_detect(nm_clean, "pgml|pg|conc|concentration|detect|satur|llod|lloq|uloq|range|interpolat")
  
  score <- rep(0, length(nm))
  score[!exclude & nm_clean == target_clean] <- 3
  score[!exclude & str_detect(nm_clean, paste0("^", target_clean, "mfi$"))] <- 2
  score[!exclude & str_detect(nm_clean, target_clean)] <- pmax(score[!exclude & str_detect(nm_clean, target_clean)], 1)
  
  if (max(score) == 0) return(NA_character_)
  nm[which.max(score)]
}

ensure_anc_mfi_columns <- function(df) {
  if (!"nida" %in% names(df)) stop("Column 'nida' not found in ANC_data.csv")
  if (any(duplicated(names(df)))) {
    dups <- unique(names(df)[duplicated(names(df))])
    stop(paste0("ANC_data.csv has duplicated column names: ", paste(dups, collapse=", "),
                "\nFix CSV headers to be unique, then rerun."))
  }
  
  col_pvldh  <- pick_best_mfi_col(df, "pvldh")
  col_pfhrp2 <- pick_best_mfi_col(df, "pfhrp2")
  col_pfldh  <- pick_best_mfi_col(df, "pfldh")
  
  msg("ANC column picked for pvldh  MFI:  %s", ifelse(is.na(col_pvldh),  "NOT FOUND", col_pvldh))
  msg("ANC column picked for pfhrp2 MFI:  %s", ifelse(is.na(col_pfhrp2), "NOT FOUND", col_pfhrp2))
  msg("ANC column picked for pfldh  MFI:  %s", ifelse(is.na(col_pfldh),  "NOT FOUND", col_pfldh))
  
  if (!is.na(col_pvldh))  df$pvldh  <- suppressWarnings(as.numeric(as.character(df[[col_pvldh]])))
  if (!is.na(col_pfhrp2)) df$pfhrp2 <- suppressWarnings(as.numeric(as.character(df[[col_pfhrp2]])))
  if (!is.na(col_pfldh))  df$pfldh  <- suppressWarnings(as.numeric(as.character(df[[col_pfldh]])))
  
  df
}

## ============================================================
## 1) READ PLATES
## ============================================================

csv_files <- list.files(path = input_directory, pattern = "\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("No CSV files found in input_directory.")

plate_list_raw <- lapply(csv_files, function(f) read.csv(f, check.names = TRUE))
names(plate_list_raw) <- tools::file_path_sans_ext(basename(csv_files))
plate_list_raw <- lapply(plate_list_raw, standardise_plate_names)

msg("Loaded %d plate CSV(s).", length(plate_list_raw))

## ============================================================
## 2) GLOBAL BLANKS (RAW MFI): LLOD/LLOQ per analyte
## ============================================================

blank_raw_long <- imap_dfr(plate_list_raw, function(df, plate) {
  extract_blanks(df, plate) %>%
    pivot_longer(cols = any_of(analytes), names_to = "Analyte", values_to = "MFI_blank_raw") %>%
    mutate(MFI_blank_raw = suppressWarnings(as.numeric(MFI_blank_raw)))
})

global_blank_stats <- blank_raw_long %>%
  group_by(Analyte) %>%
  summarise(
    n_blank_total  = sum(is.finite(MFI_blank_raw)),
    mean_blank_raw = mean(MFI_blank_raw, na.rm = TRUE),
    sd_blank_raw   = sd(MFI_blank_raw,   na.rm = TRUE),
    LLOD_MFI_raw   = mean_blank_raw + 3*sd_blank_raw,
    LLOQ_MFI_raw   = mean_blank_raw + 6*sd_blank_raw,
    .groups = "drop"
  )

write.csv(global_blank_stats,
          file.path(qc_dir, "GLOBAL_blank_stats_LLOD_LLOQ_rawMFI.csv"),
          row.names = FALSE)

## ============================================================
## 3) STANDARDS: build mean standard points (RAW MFI) and mean 5PL per analyte
## ============================================================

std_long <- imap_dfr(plate_list_raw, function(df, plate) {
  std_idx <- get_std_idx(df)
  if (length(std_idx) == 0) return(tibble())
  
  df_std <- df %>%
    slice(std_idx) %>%
    mutate(
      plate = plate,
      dilution_num = if ("Sample" %in% names(.)) parse_dilution_num(Sample) else NA_real_
    ) %>%
    select(plate, dilution_num, any_of(analytes)) %>%
    pivot_longer(cols = any_of(analytes), names_to = "Analyte", values_to = "MFI_std_raw") %>%
    mutate(MFI_std_raw = suppressWarnings(as.numeric(MFI_std_raw)))
  
  df_std
})

if (nrow(std_long) == 0) stop("No standards extracted. Check standards labels/positions.")

std_long <- std_long %>%
  filter(is.finite(dilution_num), dilution_num > 0, is.finite(MFI_std_raw), MFI_std_raw > 0) %>%
  mutate(
    conc_std_pgml_assay = case_when(
      Analyte == "PfHRP2" ~ top_pg_ml_assay[["PfHRP2"]] / dilution_num,
      Analyte == "PfLDH"  ~ top_pg_ml_assay[["PfLDH"]]  / dilution_num,
      Analyte == "PvLDH"  ~ top_pg_ml_assay[["PvLDH"]]  / dilution_num,
      TRUE ~ NA_real_
    )
    
  ) %>%
  filter(is.finite(conc_std_pgml_assay), conc_std_pgml_assay > 0)

write.csv(std_long, file.path(qc_dir, "QC_standards_long_RAWMFI_with_conc.csv"), row.names = FALSE)

msg("STANDARD max concentrations (ASSAY WELL units) — should be ~100,000 HRP2 and ~1,000,000 LDHs:")
print(std_long %>% group_by(Analyte) %>% summarise(max_conc = max(conc_std_pgml_assay), min_conc = min(conc_std_pgml_assay), .groups="drop"))

mean_std_points <- std_long %>%
  group_by(Analyte, conc_std_pgml_assay) %>%
  summarise(
    mean_MFI_raw = mean(MFI_std_raw, na.rm = TRUE),
    n_plates     = n_distinct(plate),
    .groups = "drop"
  ) %>%
  arrange(Analyte, conc_std_pgml_assay)

write.csv(mean_std_points,
          file.path(qc_dir, "MEAN_standard_points_RAWMFI_by_analyte.csv"),
          row.names = FALSE)

mean_curve_params <- mean_std_points %>%
  group_by(Analyte) %>%
  group_modify(~ fit_5pl_for_group(.x, xcol = "conc_std_pgml_assay", ycol = "mean_MFI_raw")) %>%
  ungroup()

write.csv(mean_curve_params,
          file.path(qc_dir, "MEAN_curve_5PL_parameters_RAWMFI_by_analyte.csv"),
          row.names = FALSE)

conc_ranges <- mean_std_points %>%
  group_by(Analyte) %>%
  summarise(
    min_c = min(conc_std_pgml_assay, na.rm = TRUE) / 3,
    max_c = max(conc_std_pgml_assay, na.rm = TRUE) * 3,
    .groups = "drop"
  )

## ============================================================
## 4) Convert GLOBAL LLOD/LLOQ (RAW MFI) -> pg/mL using MEAN curve
## ============================================================

global_thresholds <- global_blank_stats %>%
  left_join(mean_curve_params, by = "Analyte") %>%
  left_join(conc_ranges,       by = "Analyte") %>%
  rowwise() %>%
  mutate(
    LLOD_pgml_assay = if (isTRUE(converged)) inv_5pl_vec(LLOD_MFI_raw, A,B,C,D,E, min_c, max_c) else NA_real_,
    LLOQ_pgml_assay = if (isTRUE(converged)) inv_5pl_vec(LLOQ_MFI_raw, A,B,C,D,E, min_c, max_c) else NA_real_,
    LLOD_pgml_sample = LLOD_pgml_assay * SAMPLE_DILUTION,
    LLOQ_pgml_sample = LLOQ_pgml_assay * SAMPLE_DILUTION
  ) %>%
  ungroup() %>%
  select(Analyte, n_blank_total, mean_blank_raw, sd_blank_raw,
         LLOD_MFI_raw, LLOQ_MFI_raw,
         LLOD_pgml_assay, LLOQ_pgml_assay,
         LLOD_pgml_sample, LLOQ_pgml_sample)

write.csv(global_thresholds,
          file.path(qc_dir, "GLOBAL_thresholds_rawMFI_and_pgml_assay_and_sample.csv"),
          row.names = FALSE)

## ============================================================
## 5) ULOQ from mean curve: 0.9 * max predicted MFI
## ============================================================

uloq_table <- mean_curve_params %>%
  left_join(conc_ranges, by = "Analyte") %>%
  rowwise() %>%
  do(bind_cols(tibble(Analyte = .$Analyte),
               compute_uloq(.$A, .$B, .$C, .$D, .$E, .$min_c, .$max_c))) %>%
  ungroup() %>%
  mutate(
    ULOQ_pgml_sample = ULOQ_pgml_assay * SAMPLE_DILUTION
  )

write.csv(uloq_table,
          file.path(qc_dir, "ULOQ_meanCurve_by_analyte.csv"),
          row.names = FALSE)

## ============================================================
## 6) APPLY MEAN 5PL TO ANC DATA (RAW MFI) + censor/cap + dilution correction
## ============================================================

if (!file.exists(anc_file)) stop(sprintf("ANC file not found: %s", anc_file))
anc_data <- read.csv(anc_file, check.names = TRUE)
anc_data <- ensure_anc_mfi_columns(anc_data)

mfi_cols <- intersect(c("pvldh", "pfhrp2", "pfldh"), names(anc_data))
if (length(mfi_cols) == 0) stop("No usable MFI columns found after canonicalization (pvldh/pfhrp2/pfldh).")

anc_long <- anc_data %>%
  mutate(nida = as.character(nida)) %>%
  select(nida, all_of(mfi_cols)) %>%
  pivot_longer(cols = all_of(mfi_cols), names_to = "Analyte_raw", values_to = "MFI_raw") %>%
  mutate(
    Analyte = recode(
      Analyte_raw,
      "pvldh"  = "PvLDH",
      "pfhrp2" = "PfHRP2",
      "pfldh"  = "PfLDH"
    ),
    MFI_raw = suppressWarnings(as.numeric(MFI_raw))
  ) %>%
  left_join(mean_curve_params, by = "Analyte") %>%
  left_join(conc_ranges,       by = "Analyte") %>%
  left_join(global_thresholds %>% select(Analyte, LLOD_MFI_raw, LLOQ_MFI_raw, LLOD_pgml_assay, LLOQ_pgml_assay), by="Analyte") %>%
  left_join(uloq_table       %>% select(Analyte, ULOQ_MFI, ULOQ_pgml_assay), by="Analyte") %>%
  group_by(Analyte) %>%
  mutate(
    Conc_pgml_assay_uncensored = if (isTRUE(first(converged))) inv_5pl_vec(MFI_raw, first(A), first(B), first(C), first(D), first(E), first(min_c), first(max_c)) else NA_real_
  ) %>%
  ungroup() %>%
  mutate(
    Detection = case_when(
      !is.finite(MFI_raw) | !is.finite(LLOD_MFI_raw) ~ NA_character_,
      MFI_raw >= LLOD_MFI_raw ~ "detected",
      TRUE ~ "non detected"
    ),
    below_LLOQ = case_when(
      !is.finite(MFI_raw) | !is.finite(LLOQ_MFI_raw) ~ NA_character_,
      MFI_raw < LLOQ_MFI_raw ~ "Yes",
      TRUE ~ "No"
    ),
    analyte_saturation = case_when(
      !is.finite(MFI_raw) | !is.finite(ULOQ_MFI) ~ NA_character_,
      MFI_raw > ULOQ_MFI ~ "Yes",
      TRUE ~ "No"
    ),
    within_analytical_range = case_when(
      !is.finite(MFI_raw) | !is.finite(LLOQ_MFI_raw) | !is.finite(ULOQ_MFI) ~ NA_character_,
      (MFI_raw >= LLOQ_MFI_raw) & (MFI_raw <= ULOQ_MFI) ~ "Yes",
      TRUE ~ "No"
    ),
    # censor/cap on MFI thresholds, in ASSAY WELL concentration units
    Conc_pgml_assay = case_when(
      !is.finite(MFI_raw) ~ NA_real_,
      is.finite(LLOQ_MFI_raw) & MFI_raw < LLOQ_MFI_raw ~ NA_real_,
      is.finite(ULOQ_MFI)     & MFI_raw > ULOQ_MFI     ~ ULOQ_pgml_assay,
      TRUE ~ Conc_pgml_assay_uncensored
    ),
    # report concentration corrected back to original (pre 1:50 dilution)
    Conc_pgml_sample = Conc_pgml_assay * SAMPLE_DILUTION
  )

write.csv(anc_long,
          file.path(qc_dir, "ANC_interpolated_long_meanCurve.csv"),
          row.names = FALSE)

## --- Save side lists (like your other pipeline) ---
anc_non_detected <- anc_long %>%
  filter(Detection == "non detected") %>%
  select(nida, Analyte, MFI_raw, LLOD_MFI_raw, LLOQ_MFI_raw)

write.csv(anc_non_detected,
          file.path(qc_dir, "ANC_nonDetected_by_analyte.csv"),
          row.names = FALSE)

anc_below_lloq <- anc_long %>%
  filter(below_LLOQ == "Yes") %>%
  select(nida, Analyte, MFI_raw, LLOQ_MFI_raw)

write.csv(anc_below_lloq,
          file.path(qc_dir, "ANC_belowLLOQ_by_analyte.csv"),
          row.names = FALSE)

anc_saturated <- anc_long %>%
  filter(analyte_saturation == "Yes") %>%
  select(nida, Analyte, MFI_raw, ULOQ_MFI, ULOQ_pgml_assay, Conc_pgml_sample)

write.csv(anc_saturated,
          file.path(qc_dir, "ANC_saturated_by_analyte.csv"),
          row.names = FALSE)

## ============================================================
## 7) WIDE MERGE BACK TO ANC + SAVE FINAL DATABASES
## ============================================================

conc_wide <- anc_long %>%
  mutate(conc_col = case_when(
    Analyte=="PfHRP2" ~ "PfHRP2_pg.ml",
    Analyte=="PfLDH"  ~ "PfLDH_pg.ml",
    Analyte=="PvLDH"  ~ "PvLDH_pg.ml",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(conc_col)) %>%
  select(nida, conc_col, Conc_pgml_sample) %>%
  pivot_wider(names_from = conc_col, values_from = Conc_pgml_sample)

det_wide <- anc_long %>%
  filter(Analyte %in% c("PfHRP2","PfLDH")) %>%
  mutate(det_col = ifelse(Analyte=="PfHRP2","pfhrp2_detection","pfldh_detection")) %>%
  select(nida, det_col, Detection) %>%
  pivot_wider(names_from = det_col, values_from = Detection)

sat_wide <- anc_long %>%
  mutate(sat_col = case_when(
    Analyte=="PfHRP2" ~ "pfhrp2_saturation",
    Analyte=="PfLDH"  ~ "pfldh_saturation",
    Analyte=="PvLDH"  ~ "pvldh_saturation",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(sat_col)) %>%
  select(nida, sat_col, analyte_saturation) %>%
  pivot_wider(names_from = sat_col, values_from = analyte_saturation)

anc_with_conc <- anc_data %>%
  mutate(nida = as.character(nida)) %>%
  left_join(conc_wide, by = "nida") %>%
  left_join(det_wide,  by = "nida") %>%
  left_join(sat_wide,  by = "nida")

# Save in output_directory (as before)
out_main <- file.path(output_directory, "ANC_interpolated_5PL_with_detection_and_saturation.csv")
write.csv(anc_with_conc, out_main, row.names = FALSE)

# Save final requested file in PhD_HRP2_qPCR
out_final <- file.path(final_dir, "ANC_concentrations.csv")
write.csv(anc_with_conc, out_final, row.names = FALSE)

## ============================================================
## 8) QC plots (optional but useful)
## ============================================================

curve_grid <- mean_curve_params %>%
  left_join(conc_ranges, by="Analyte") %>%
  filter(converged, is.finite(min_c), is.finite(max_c), min_c > 0, max_c > 0) %>%
  rowwise() %>%
  do({
    conc <- 10^seq(log10(.$min_c), log10(.$max_c), length.out = 400)
    tibble(
      Analyte = .$Analyte,
      conc_pgml_assay = conc,
      MFI_pred = five_pl(conc, .$A, .$B, .$C, .$D, .$E)
    )
  }) %>%
  ungroup()

p_std_fit <- ggplot() +
  geom_point(data = mean_std_points, aes(x = conc_std_pgml_assay, y = mean_MFI_raw),
             alpha = 0.6, size = 1) +
  geom_line(data = curve_grid, aes(x = conc_pgml_assay, y = MFI_pred), linewidth = 0.9) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~Analyte, scales="free") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Mean standard points (RAW MFI) + mean 5PL fit",
    x = "Concentration in assay well (pg/mL, log10)",
    y = "MFI (log10)"
  )

ggsave(file.path(qc_dir, "QC_meanStandards_plus_mean5PL.png"),
       p_std_fit, width = 10, height = 6, dpi = 300)

p_mfi_conc <- anc_long %>%
  filter(is.finite(MFI_raw), is.finite(Conc_pgml_sample), Conc_pgml_sample > 0) %>%
  ggplot(aes(x = MFI_raw, y = Conc_pgml_sample, color = within_analytical_range)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~Analyte, scales="free") +
  theme_minimal(base_size = 12) +
  labs(
    title = sprintf("ANC: RAW MFI vs concentration (sample-corrected; ×%d)", SAMPLE_DILUTION),
    x = "RAW MFI (log10)",
    y = "Concentration in original sample (pg/mL, log10)",
    color = "Within range"
  )

ggsave(file.path(qc_dir, "QC_ANC_MFI_vs_Conc_sampleCorrected.png"),
       p_mfi_conc, width = 10, height = 6, dpi = 300)

msg("DONE.")
msg("Main output: %s", out_main)
msg("Final output: %s", out_final)
msg("QC folder: %s", qc_dir)
