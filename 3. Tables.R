# ============================================================
# Tables by centro (FIRST ANC) — aligned to regression cleaning
#  - Table A: ALL first ANC samples (descriptive cohort)
#  - Table B: REGRESSION-MATCHED POS SET:
#       PCR+ AND (HRP2 detected OR PfLDH detected)
#       AND complete-case for regression covariates (as in lm drop)
#
# Continuous handling:
#   - For PfHRP2 and PfLDH: if *_detection != "detected",
#     concentration set to NA (excluded) to avoid tiny values.
#   - Zeros in continuous vars excluded from ALL calculations + tests (0 -> NA).
#   - PfLDH (pg/mL), PfHRP2 (pg/mL), density:
#       GEOMETRIC MEAN (min–max) via log10.
#   - VAR2CSA + cumulative antibodies (already log-transformed):
#       GEOMETRIC MEAN (min–max) on ORIGINAL scale (back-transformed).
#
# PfHRP2 saturation handling:
#   - Raw > pfhrp2_sat_thr = saturated
#   - For summaries, saturated values are capped to pfhrp2_cap_val
#   - If saturation occurred in a centro, displayed max becomes >threshold
# ============================================================

setwd("C:/HRP2 analysis/PhD_HRP2_qPCR")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr, tidyr, tibble, stringr,
  gtsummary, broom, gt,
  officer, flextable
)

data_raw <- read.csv("ANC_recoded.csv", stringsAsFactors = FALSE)

# ---- PfHRP2 saturation threshold (pg/mL) ----
pfhrp2_sat_thr <- 10000
pfhrp2_cap_val <- 10000
pfhrp2_sat_thr_fmt <- formatC(pfhrp2_sat_thr, format = "f", digits = 2, big.mark = ",")

# ---- Antibody back-transform assumption ----
# TRUE  = antibodies stored as log10(x + 1)  -> original = 10^log - 1
# FALSE = antibodies stored as log10(x)      -> original = 10^log
antibody_log_plus1 <- TRUE

# ---- Helper: Fisher for 2x2 if any expected < 5, else Pearson chi-square ----
fisher_if_needed <- function(data, variable, by, ...) {
  tab <- table(data[[variable]], data[[by]], useNA = "no")
  if (length(tab) == 0L || any(dim(tab) < 2L)) return(tibble::tibble(p.value = NA_real_))
  if (all(dim(tab) == 2L)) {
    exp_counts <- suppressWarnings(chisq.test(tab, correct = FALSE)$expected)
    if (any(exp_counts < 5)) return(broom::tidy(fisher.test(tab)))
    return(broom::tidy(chisq.test(tab, correct = FALSE)))
  }
  broom::tidy(chisq.test(tab, correct = FALSE))
}

# ---- Helpers: detection parsing + numeric hygiene ----
is_detected <- function(x) {
  x2 <- tolower(trimws(as.character(x)))
  x2 %in% c("detected", "positive", "pos")
}
zero_to_na <- function(x) dplyr::if_else(!is.na(x) & x == 0, NA_real_, x)
nonfinite_to_na <- function(x) dplyr::if_else(!is.finite(x), NA_real_, x)

# ---- Standardize key fields like regression scripts ----
standardize_core_fields <- function(df) {
  df %>%
    mutate(
      # centre recode seen in your regression scripts
      centro = as.character(centro),
      centro = dplyr::recode(centro, "ManhiÃ§a" = "Manhiça", .default = centro),
      
      # numeric coercions consistent with your models
      pcrpos_num = suppressWarnings(as.numeric(as.character(pcrpos))),
      density    = suppressWarnings(as.numeric(as.character(density))),
      gestweeks  = suppressWarnings(as.numeric(as.character(gestweeks))),
      
      PfLDH_pg.ml  = suppressWarnings(as.numeric(as.character(PfLDH_pg.ml))),
      PfHRP2_pg.ml = suppressWarnings(as.numeric(as.character(PfHRP2_pg.ml))),
      
      antibody_pregnancy  = suppressWarnings(as.numeric(antibody_pregnancy)),
      antibody_cumulative = suppressWarnings(as.numeric(antibody_cumulative)),
      
      # net/fumig missing code handling aligned to regression scripts (999 -> NA)
      net   = dplyr::na_if(as.character(net),   "999"),
      fumig = dplyr::na_if(as.character(fumig), "999"),
      
      # detection flags
      pfhrp2_is_detected = is_detected(pfhrp2_detection),
      pfldh_is_detected  = is_detected(pfldh_detection)
    ) %>%
    mutate(
      # optional recode if net/fumig stored as 0/1 strings
      net = dplyr::case_when(
        net %in% c("0") ~ "Yes",
        net %in% c("1") ~ "No",
        TRUE ~ net
      ),
      fumig = dplyr::case_when(
        fumig %in% c("0") ~ "Yes",
        fumig %in% c("1") ~ "No",
        TRUE ~ fumig
      )
    ) %>%
    mutate(
      across(any_of(c("age_group","Season","Period","net","fumig","sede","hiv","centro","gravidity")),
             as.factor),
      # keep net/fumig levels but DO NOT create a fake "NA" category
      net   = if ("net" %in% names(.))   forcats::fct_drop(net)   else net,
      fumig = if ("fumig" %in% names(.)) forcats::fct_drop(fumig) else fumig
    )
}

# ---- Build Table B population matched to regression eligibility ----
# This mirrors the *analysis set* that enters models:
#  - first ANC
#  - PCR+
#  - AND ≥1 antigen detected (HRP2 or PfLDH)
#  - AND complete-case for regression covariates used
build_regression_matched_pos_set <- function(df) {
  
  df0 <- df %>%
    filter(visit == "first ANC") %>%
    standardize_core_fields() %>%
    filter(!is.na(pcrpos_num) & pcrpos_num == 1) %>%
    filter(pfhrp2_is_detected | pfldh_is_detected)
  
  # covariates used in regression scripts (keep only those present)
  covars_all <- c("gestweeks","density","gravidity","centro",
                  "age_group","Season","Period","net","fumig","sede","hiv")
  covars_used <- covars_all[covars_all %in% names(df0)]
  
  # make Table B reflect the same complete-case removal that lm will do
  df0 %>%
    filter(if_all(all_of(covars_used), ~ !is.na(.x)))
}

# ---- Cleaning function for TABLE values (your rules kept) ----
clean_for_table_values <- function(df_std) {
  df_std %>%
    mutate(
      centro = factor(as.character(centro),
                      levels = intersect(c("Ilha Josina","Manhiça","Magude"), unique(as.character(centro)))),
      
      # EXCLUDE concentrations when corresponding marker is NON detected
      PfLDH_pg.ml_tbl  = if_else(pfldh_is_detected,  PfLDH_pg.ml,  NA_real_),
      PfHRP2_pg.ml_tbl = if_else(pfhrp2_is_detected, PfHRP2_pg.ml, NA_real_)
    ) %>%
    # PfHRP2 saturation capping (applies only to detected values because others are NA)
    mutate(
      PfHRP2_saturated = !is.na(PfHRP2_pg.ml_tbl) & PfHRP2_pg.ml_tbl > pfhrp2_sat_thr,
      PfHRP2_pg.ml_tbl = if_else(PfHRP2_saturated, pfhrp2_cap_val, PfHRP2_pg.ml_tbl)
    ) %>%
    # Remove zeros from ALL continuous vars AFTER detection filtering + capping
    mutate(
      across(
        c(PfLDH_pg.ml_tbl, PfHRP2_pg.ml_tbl, density, gestweeks,
          antibody_pregnancy, antibody_cumulative),
        ~ zero_to_na(nonfinite_to_na(.x))
      )
    )
}

# ---- Force displayed max to >threshold when saturation exists ----
fix_pfhrp2_max_display <- function(tbl, df_clean) {
  lvl <- levels(df_clean$centro)
  stat_cols <- paste0("stat_", seq_along(lvl))
  
  any_sat <- df_clean %>%
    group_by(centro) %>%
    summarise(any_sat = any(PfHRP2_saturated, na.rm = TRUE), .groups = "drop")
  
  tbl %>%
    gtsummary::modify_table_body(~{
      tb <- .
      for (i in seq_along(lvl)) {
        centro_i <- lvl[i]
        col_i <- stat_cols[i]
        has_sat <- isTRUE(any_sat$any_sat[match(centro_i, any_sat$centro)])
        if (!has_sat || !col_i %in% names(tb)) next
        
        tb <- tb %>%
          mutate(
            !!col_i := ifelse(
              variable == "PfHRP2_pg.ml_tbl" & row_type == "label" & !is.na(.data[[col_i]]),
              sub("(\\([^\\-]+-)([^\\)]+)(\\))",
                  paste0("\\1>", pfhrp2_sat_thr_fmt, "\\3"),
                  .data[[col_i]]),
              .data[[col_i]]
            )
          )
      }
      tb
    })
}

# ---- Geometric mean helpers (zeros already NA; use >0 safety) ----
gmean_raw_log10 <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (!length(x)) return(NA_real_)
  10^(mean(log10(x)))
}
min_pos <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (!length(x)) return(NA_real_)
  min(x)
}
max_pos <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (!length(x)) return(NA_real_)
  max(x)
}

ab_back <- function(x_log10) if (antibody_log_plus1) 10^x_log10 - 1 else 10^x_log10
ab_gmean <- function(x_log10) { x_log10 <- x_log10[!is.na(x_log10)]; if (!length(x_log10)) return(NA_real_); ab_back(mean(x_log10)) }
ab_min   <- function(x_log10) { x_log10 <- x_log10[!is.na(x_log10)]; if (!length(x_log10)) return(NA_real_); ab_back(min(x_log10)) }
ab_max   <- function(x_log10) { x_log10 <- x_log10[!is.na(x_log10)]; if (!length(x_log10)) return(NA_real_); ab_back(max(x_log10)) }

fmt_num <- function(x, digits = 2) {
  if (is.na(x)) return(NA_character_)
  thr <- 10^-digits
  if (x > 0 && x < thr) return(paste0("<", formatC(thr, format = "f", digits = digits)))
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

replace_stats_with_geomean <- function(tbl, df_clean) {
  lvl <- levels(df_clean$centro)
  stat_cols <- paste0("stat_", seq_along(lvl))
  
  summ <- df_clean %>%
    group_by(centro) %>%
    summarise(
      PfLDH_gm  = gmean_raw_log10(PfLDH_pg.ml_tbl),
      PfLDH_min = min_pos(PfLDH_pg.ml_tbl),
      PfLDH_max = max_pos(PfLDH_pg.ml_tbl),
      
      PfHRP2_gm  = gmean_raw_log10(PfHRP2_pg.ml_tbl),
      PfHRP2_min = min_pos(PfHRP2_pg.ml_tbl),
      PfHRP2_max = max_pos(PfHRP2_pg.ml_tbl),
      
      dens_gm  = gmean_raw_log10(density),
      dens_min = min_pos(density),
      dens_max = max_pos(density),
      
      ab_preg_gm  = ab_gmean(antibody_pregnancy),
      ab_preg_min = ab_min(antibody_pregnancy),
      ab_preg_max = ab_max(antibody_pregnancy),
      
      ab_cum_gm  = ab_gmean(antibody_cumulative),
      ab_cum_min = ab_min(antibody_cumulative),
      ab_cum_max = ab_max(antibody_cumulative),
      .groups = "drop"
    )
  
  tbl %>%
    gtsummary::modify_table_body(~{
      tb <- .
      for (i in seq_along(lvl)) {
        centro_i <- lvl[i]
        col_i <- stat_cols[i]
        if (!col_i %in% names(tb)) next
        s <- summ[summ$centro == centro_i, , drop = FALSE]
        if (nrow(s) == 0) next
        
        tb <- tb %>%
          mutate(
            !!col_i := dplyr::case_when(
              variable == "PfLDH_pg.ml_tbl" & row_type == "label" ~
                paste0(fmt_num(s$PfLDH_gm), " (", fmt_num(s$PfLDH_min), "-", fmt_num(s$PfLDH_max), ")"),
              
              variable == "PfHRP2_pg.ml_tbl" & row_type == "label" ~
                paste0(fmt_num(s$PfHRP2_gm), " (", fmt_num(s$PfHRP2_min), "-", fmt_num(s$PfHRP2_max), ")"),
              
              variable == "density" & row_type == "label" ~
                paste0(fmt_num(s$dens_gm), " (", fmt_num(s$dens_min), "-", fmt_num(s$dens_max), ")"),
              
              variable == "antibody_pregnancy" & row_type == "label" ~
                paste0(fmt_num(s$ab_preg_gm), " (", fmt_num(s$ab_preg_min), "-", fmt_num(s$ab_preg_max), ")"),
              
              variable == "antibody_cumulative" & row_type == "label" ~
                paste0(fmt_num(s$ab_cum_gm), " (", fmt_num(s$ab_cum_min), "-", fmt_num(s$ab_cum_max), ")"),
              
              TRUE ~ .data[[col_i]]
            )
          )
      }
      tb
    })
}

# ---- Build datasets ----
data_first_all <- data_raw %>%
  filter(visit == "first ANC") %>%
  standardize_core_fields()

data_pos_regression_matched <- build_regression_matched_pos_set(data_raw)

# Apply table-value cleaning rules (detection -> NA, zeros -> NA, HRP2 cap)
data_all_clean <- clean_for_table_values(data_first_all)
data_pos_clean <- clean_for_table_values(data_pos_regression_matched)

cat("\n================ TABLE POPULATIONS ================\n")
cat("Table A (First ANC, all): n =", nrow(data_all_clean), "\n")
cat("Table B (Regression-matched: PCR+ & ≥1 antigen detected & complete-case): n =", nrow(data_pos_clean), "\n")

# ---- Variables ----
vars_cont <- c("PfLDH_pg.ml_tbl", "PfHRP2_pg.ml_tbl", "density", "gestweeks",
               "antibody_pregnancy", "antibody_cumulative")
vars_cat  <- c("gravidity", "hiv", "age_group", "net", "fumig",
               "Season", "Period", "sede")

# ---- Labels ----
var_labels <- list(
  PfLDH_pg.ml_tbl  ~ "Pf-LDH (pg/mL; only if detected)",
  PfHRP2_pg.ml_tbl ~ paste0("PfHRP2 (pg/mL; only if detected; capped at ", pfhrp2_cap_val,
                            " if > ", pfhrp2_sat_thr_fmt, ")"),
  density ~ "Parasite density (parasites/µL)",
  gestweeks ~ "Gestational week",
  antibody_pregnancy  ~ "VAR2CSA IgG (pregnancy-specific)",
  antibody_cumulative ~ "Cumulative IgG",
  gravidity ~ "Gravidity",
  hiv   ~ "HIV status",
  age_group ~ "Age group",
  net       ~ "Bed net use",
  fumig     ~ "IRS fumigation",
  Season    ~ "Season",
  Period    ~ "Period",
  sede      ~ "Type of setting"
)

gtsummary::theme_gtsummary_journal(journal = "jama")
gtsummary::theme_gtsummary_compact()

build_tbl <- function(df_clean) {
  tbl <- df_clean %>%
    select(centro, all_of(vars_cont), all_of(vars_cat)) %>%
    gtsummary::tbl_summary(
      by = centro,
      type = list(all_continuous() ~ "continuous", all_categorical() ~ "categorical"),
      statistic = list(
        all_continuous()  ~ "{median} ({min}-{max})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      percent = "column",
      missing = "no",
      label = var_labels,
      digits = list(all_continuous() ~ 2, all_categorical() ~ 1)
    ) %>%
    gtsummary::add_p(
      test = list(all_continuous() ~ "kruskal.test", all_categorical() ~ fisher_if_needed),
      pvalue_fun = ~ gtsummary::style_pvalue(.x, digits = 3)
    ) %>%
    gtsummary::bold_labels() %>%
    gtsummary::bold_p(t = 0.05)
  
  tbl %>%
    replace_stats_with_geomean(df_clean) %>%
    fix_pfhrp2_max_display(df_clean)
}

tbl_all <- build_tbl(data_all_clean)
tbl_pos <- build_tbl(data_pos_clean)

# ---- Export to editable Word with BOTH tables ----
ft_all <- gtsummary::as_flex_table(tbl_all) %>%
  flextable::autofit() %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::align(align = "center", part = "header")

ft_pos <- gtsummary::as_flex_table(tbl_pos) %>%
  flextable::autofit() %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::align(align = "center", part = "header")

doc <- officer::read_docx()

ab_note <- if (antibody_log_plus1) {
  "Antibodies back-transformed assuming log10(x+1): original = 10^log - 1."
} else {
  "Antibodies back-transformed assuming log10(x): original = 10^log."
}

note_txt <- paste0(
  "For PfHRP2 and PfLDH, concentrations were excluded (set to missing) when the corresponding *_detection was not detected. ",
  "Zeros in continuous variables were excluded from summaries and tests (0 -> missing). ",
  "PfLDH, PfHRP2, and parasite density are shown as geometric mean (min–max) computed via log10. ",
  "VAR2CSA and cumulative antibodies are shown as geometric mean (min–max) after back-transformation. ",
  ab_note, " ",
  "PfHRP2 saturation: values > ", pfhrp2_sat_thr_fmt,
  " were capped to ", pfhrp2_cap_val,
  " for summary stats; if saturation occurred in a centro, the maximum is shown as >",
  pfhrp2_sat_thr_fmt, "."
)

doc <- officer::body_add_par(doc, "Supplementary Table A. First ANC — all samples", style = "heading 1")
doc <- officer::body_add_par(doc, note_txt, style = "Normal")
doc <- flextable::body_add_flextable(doc, value = ft_all)
doc <- officer::body_add_par(doc, "", style = "Normal")

doc <- officer::body_add_par(doc, "Supplementary Table B. First ANC — regression-matched positive set", style = "heading 1")
doc <- officer::body_add_par(doc, "Subset: PCR+ and ≥1 antigen detected (HRP2 and/or PfLDH); complete-case for regression covariates.", style = "Normal")
doc <- officer::body_add_par(doc, note_txt, style = "Normal")
doc <- flextable::body_add_flextable(doc, value = ft_pos)

print(doc, target = "tables_by_centro_firstANC_ALL_and_REGRESSION_POS.docx")

# ---- Optional: HTML + PNG exports ----
tbl_all_gt <- tbl_all %>%
  gtsummary::as_gt() %>%
  gt::tab_header(
    title = gt::md("**Supplementary Table A. First ANC — all samples**"),
    subtitle = gt::md("Table-value rules: antigens excluded if not detected; zeros excluded; HRP2 capped if saturated.")
  )

tbl_pos_gt <- tbl_pos %>%
  gtsummary::as_gt() %>%
  gt::tab_header(
    title = gt::md("**Supplementary Table B. First ANC — regression-matched positive set**"),
    subtitle = gt::md("PCR+ and ≥1 antigen detected; complete-case for regression covariates; table-value rules applied.")
  )

gt::gtsave(tbl_all_gt, "table_A_firstANC_all.html")
gt::gtsave(tbl_pos_gt, "table_B_firstANC_regression_pos.html")

if (!requireNamespace("webshot2", quietly = TRUE)) install.packages("webshot2")
gt::gtsave(tbl_all_gt, "table_A_firstANC_all.png")
gt::gtsave(tbl_pos_gt, "table_B_firstANC_regression_pos.png")

message("Saved Word: tables_by_centro_firstANC_ALL_and_REGRESSION_POS.docx (editable). Also saved HTML/PNG versions.")
