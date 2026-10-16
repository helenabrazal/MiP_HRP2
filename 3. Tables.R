# ============================================================
# Tables by centro (FIRST ANC)
#  - Table A: ALL first ANC samples (descriptive cohort)
#  - Table B: ALL qPCR-POSITIVE infections, which is the set
#             entering the regression models
#
# TWO CHANGES FROM THE PREVIOUS VERSION
#
# 1. Bed net and fumigation are no longer recoded here. The
#    previous version had 0 as "Yes" and 1 as "No", the reverse
#    of the data preparation script, and converted 999 to NA so
#    that 158 infections went missing from the table while the
#    models retained them. vector_control now comes from the
#    data preparation script, so there is one definition.
#
# 2. Table B excluded infections with no detectable antigen,
#    leaving 343 of 394. The models include all 394, with
#    undetectable antigen handled by left censoring, so that
#    filter is removed and the table now describes the set the
#    models actually use.
#
# Continuous handling:
#  - For PfHRP2 and PfLDH: if *_detection != "detected",
#    concentration set to NA (excluded) to avoid tiny values.
#  - Zeros in continuous vars excluded from ALL calculations
#    and tests (0 -> NA).
#  - PfLDH, PfHRP2, density: GEOMETRIC MEAN (min-max) via log10.
#  - VAR2CSA and cumulative antibodies (already log-transformed):
#    GEOMETRIC MEAN (min-max) back-transformed.
#
# PfHRP2 saturation handling:
#  - Raw > pfhrp2_sat_thr = saturated
#  - For summaries, saturated values capped to pfhrp2_cap_val
#  - If saturation occurred in a centro, displayed max becomes
#    >threshold
# ============================================================

## ------------------------------------------------------------
## Paths
##
## DATA_DIR is where ANC_recoded.csv is read from; OUT_DIR is
## where the tables are written. Every output below uses an
## explicit path, so nothing depends on the working directory
## and files cannot end up wherever the session happened to be.
## ------------------------------------------------------------

DATA_DIR <- "C:/HRP2 analysis/MiP_HRP2_manuscriptcodes/MiP_HRP2"
OUT_DIR  <- "C:/HRP2 analysis/PhD_HRP2_qPCR"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

out <- function(f) file.path(OUT_DIR, f)


# ============================================================
# Packages
# ============================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  dplyr, tidyr, tibble, stringr, gtsummary,
  broom, gt, officer, flextable, forcats
)


# ============================================================
# Read data
# ============================================================

data_raw <- read.csv(
  file.path(DATA_DIR, "ANC_recoded.csv"),
  stringsAsFactors = FALSE
)

if (!"vector_control" %in% names(data_raw)) {
  stop(
    "vector_control is missing from ANC_recoded.csv. Run the data ",
    "preparation script first; this script does not build it."
  )
}


# ============================================================
# PfHRP2 saturation threshold
# ============================================================

pfhrp2_sat_thr <- 10000
pfhrp2_cap_val <- 10000

pfhrp2_sat_thr_fmt <- formatC(
  pfhrp2_sat_thr, format = "f", digits = 2, big.mark = ","
)


# ============================================================
# Antibody back-transform assumption
# ============================================================
# TRUE  = antibodies stored as log10(x + 1); original = 10^log - 1
# FALSE = antibodies stored as log10(x);     original = 10^log
# ============================================================

antibody_log_plus1 <- TRUE


# ============================================================
# Helper: Fisher for 2x2 if any expected < 5, else chi-square
# ============================================================

fisher_if_needed <- function(data, variable, by, ...) {
  
  tab <- table(data[[variable]], data[[by]], useNA = "no")
  
  if (length(tab) == 0L || any(dim(tab) < 2L)) {
    return(tibble::tibble(p.value = NA_real_))
  }
  
  if (all(dim(tab) == 2L)) {
    
    exp_counts <- suppressWarnings(
      chisq.test(tab, correct = FALSE)$expected
    )
    
    if (any(exp_counts < 5)) {
      return(broom::tidy(fisher.test(tab)))
    }
    
    return(broom::tidy(chisq.test(tab, correct = FALSE)))
  }
  
  broom::tidy(chisq.test(tab, correct = FALSE))
}


# ============================================================
# Helpers: detection parsing + numeric hygiene
# ============================================================

is_detected <- function(x) {
  x2 <- tolower(trimws(as.character(x)))
  x2 %in% c("detected", "positive", "pos")
}


zero_to_na <- function(x) {
  dplyr::if_else(!is.na(x) & x == 0, NA_real_, x)
}


nonfinite_to_na <- function(x) {
  dplyr::if_else(!is.finite(x), NA_real_, x)
}


# ============================================================
# Standardize key fields
# ============================================================

standardize_core_fields <- function(df) {
  
  df %>%
    
    mutate(
      
      # ------------------------------------------------------
      # Centre recode
      # ------------------------------------------------------
      
      centro = as.character(centro),
      
      centro = dplyr::recode(
        centro, "ManhiÃ§a" = "Manhiça", .default = centro
      ),
      
      
      # ------------------------------------------------------
      # Numeric coercions consistent with models
      # ------------------------------------------------------
      
      pcrpos_num = suppressWarnings(as.numeric(as.character(pcrpos))),
      density    = suppressWarnings(as.numeric(as.character(density))),
      gestweeks  = suppressWarnings(as.numeric(as.character(gestweeks))),
      
      PfLDH_pg.ml  = suppressWarnings(as.numeric(as.character(PfLDH_pg.ml))),
      PfHRP2_pg.ml = suppressWarnings(as.numeric(as.character(PfHRP2_pg.ml))),
      
      antibody_pregnancy  = suppressWarnings(as.numeric(antibody_pregnancy)),
      antibody_cumulative = suppressWarnings(as.numeric(antibody_cumulative)),
      
      
      # ------------------------------------------------------
      # Detection flags
      # ------------------------------------------------------
      
      pfhrp2_is_detected = is_detected(pfhrp2_detection),
      pfldh_is_detected  = is_detected(pfldh_detection)
    ) %>%
    
    # ======================================================
  # VECTOR CONTROL
  #
  # Built in the data preparation script, not here. That is
  # where net and fumig are recoded (1 = Yes, 0 = No) and
  # combined into a count of measures reported, with 999 kept
  # as "Undetermined". Rebuilding it here would be a second
  # place for the two to drift apart.
  # ======================================================
  
  mutate(
    
    # ------------------------------------------------------
    # Convert categorical variables to factors
    # ------------------------------------------------------
    
    across(
      any_of(c(
        "age_group", "Season", "Period",
        "sede", "hiv", "centro", "gravidity", "vector_control"
      )),
      as.factor
    ),
    
    
    # ------------------------------------------------------
    # Explicit order for vector control
    # ------------------------------------------------------
    
    vector_control = factor(
      vector_control,
      levels = c("0 measures", "1 measure", "2 measures", "Undetermined")
    )
  )
}


# ============================================================
# Build the qPCR-positive set
#
# This is the population entering the models: first ANC and
# qPCR-positive. Infections with no detectable antigen are NOT
# excluded; they are informative observations handled by left
# censoring in the models, and excluding them here would make
# the table describe a different sample from the analysis.
#
# No complete-case step is applied either: vector_control now
# carries Undetermined as a level rather than NA, so no record
# is dropped for it, and the models fit on the same rows.
# ============================================================

build_pcrpos_set <- function(df) {
  
  df %>%
    filter(visit == "first ANC") %>%
    standardize_core_fields() %>%
    filter(!is.na(pcrpos_num) & pcrpos_num == 1)
}


# ============================================================
# Cleaning function for TABLE values
# ============================================================

clean_for_table_values <- function(df_std) {
  
  df_std %>%
    
    mutate(
      
      centro = factor(
        as.character(centro),
        levels = intersect(
          c("Ilha Josina", "Manhiça", "Magude"),
          unique(as.character(centro))
        )
      ),
      
      
      # ------------------------------------------------------
      # Exclude concentrations when the marker is NOT detected
      # ------------------------------------------------------
      
      PfLDH_pg.ml_tbl = if_else(
        pfldh_is_detected, PfLDH_pg.ml, NA_real_
      ),
      
      PfHRP2_pg.ml_tbl = if_else(
        pfhrp2_is_detected, PfHRP2_pg.ml, NA_real_
      )
    ) %>%
    
    mutate(
      
      PfHRP2_saturated =
        !is.na(PfHRP2_pg.ml_tbl) & PfHRP2_pg.ml_tbl > pfhrp2_sat_thr,
      
      PfHRP2_pg.ml_tbl =
        if_else(PfHRP2_saturated, pfhrp2_cap_val, PfHRP2_pg.ml_tbl)
    ) %>%
    
    mutate(
      
      across(
        c(
          PfLDH_pg.ml_tbl, PfHRP2_pg.ml_tbl, density,
          gestweeks, antibody_pregnancy, antibody_cumulative
        ),
        ~ zero_to_na(nonfinite_to_na(.x))
      )
    )
}


# ============================================================
# Force displayed maximum to > threshold when saturation exists
# ============================================================

fix_pfhrp2_max_display <- function(tbl, df_clean) {
  
  lvl <- levels(df_clean$centro)
  stat_cols <- paste0("stat_", seq_along(lvl))
  
  any_sat <- df_clean %>%
    group_by(centro) %>%
    summarise(any_sat = any(PfHRP2_saturated, na.rm = TRUE), .groups = "drop")
  
  tbl %>%
    
    gtsummary::modify_table_body(~ {
      
      tb <- .
      
      for (i in seq_along(lvl)) {
        
        centro_i <- lvl[i]
        col_i <- stat_cols[i]
        
        has_sat <- isTRUE(
          any_sat$any_sat[match(centro_i, any_sat$centro)]
        )
        
        if (!has_sat || !col_i %in% names(tb)) next
        
        tb <- tb %>%
          mutate(
            !!col_i := ifelse(
              variable == "PfHRP2_pg.ml_tbl" &
                row_type == "label" &
                !is.na(.data[[col_i]]),
              sub(
                "(\\([^\\-]+-)([^\\)]+)(\\))",
                paste0("\\1>", pfhrp2_sat_thr_fmt, "\\3"),
                .data[[col_i]]
              ),
              .data[[col_i]]
            )
          )
      }
      
      tb
    })
}


# ============================================================
# Geometric mean helpers
# ============================================================

gmean_raw_log10 <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (!length(x)) return(NA_real_)
  10 ^ mean(log10(x))
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


# ============================================================
# Antibody back-transformation helpers
# ============================================================

ab_back <- function(x_log10) {
  if (antibody_log_plus1) 10 ^ x_log10 - 1 else 10 ^ x_log10
}

ab_gmean <- function(x_log10) {
  x_log10 <- x_log10[!is.na(x_log10)]
  if (!length(x_log10)) return(NA_real_)
  ab_back(mean(x_log10))
}

ab_min <- function(x_log10) {
  x_log10 <- x_log10[!is.na(x_log10)]
  if (!length(x_log10)) return(NA_real_)
  ab_back(min(x_log10))
}

ab_max <- function(x_log10) {
  x_log10 <- x_log10[!is.na(x_log10)]
  if (!length(x_log10)) return(NA_real_)
  ab_back(max(x_log10))
}


# ============================================================
# Number formatting
# ============================================================

fmt_num <- function(x, digits = 2) {
  
  if (is.na(x)) return(NA_character_)
  
  thr <- 10 ^ -digits
  
  if (x > 0 && x < thr) {
    return(paste0("<", formatC(thr, format = "f", digits = digits)))
  }
  
  formatC(x, format = "f", digits = digits, big.mark = ",")
}


# ============================================================
# Replace table statistics with geometric means
# ============================================================

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
    
    gtsummary::modify_table_body(~ {
      
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
                paste0(
                  fmt_num(s$PfLDH_gm), " (",
                  fmt_num(s$PfLDH_min), "-",
                  fmt_num(s$PfLDH_max), ")"
                ),
              
              variable == "PfHRP2_pg.ml_tbl" & row_type == "label" ~
                paste0(
                  fmt_num(s$PfHRP2_gm), " (",
                  fmt_num(s$PfHRP2_min), "-",
                  fmt_num(s$PfHRP2_max), ")"
                ),
              
              variable == "density" & row_type == "label" ~
                paste0(
                  fmt_num(s$dens_gm), " (",
                  fmt_num(s$dens_min), "-",
                  fmt_num(s$dens_max), ")"
                ),
              
              variable == "antibody_pregnancy" & row_type == "label" ~
                paste0(
                  fmt_num(s$ab_preg_gm), " (",
                  fmt_num(s$ab_preg_min), "-",
                  fmt_num(s$ab_preg_max), ")"
                ),
              
              variable == "antibody_cumulative" & row_type == "label" ~
                paste0(
                  fmt_num(s$ab_cum_gm), " (",
                  fmt_num(s$ab_cum_min), "-",
                  fmt_num(s$ab_cum_max), ")"
                ),
              
              TRUE ~ .data[[col_i]]
            )
          )
      }
      
      tb
    })
}


# ============================================================
# Build datasets
# ============================================================

data_first_all <- data_raw %>%
  filter(visit == "first ANC") %>%
  standardize_core_fields()

data_pcrpos <- build_pcrpos_set(data_raw)

data_all_clean <- clean_for_table_values(data_first_all)
data_pos_clean <- clean_for_table_values(data_pcrpos)


# ============================================================
# Display table populations
#
# Table B must come to 394. If it does not, the filters here
# have drifted from the analysis scripts.
# ============================================================

cat("\n================ TABLE POPULATIONS ================\n")
cat("Table A (First ANC, all)        : n =", nrow(data_all_clean), "\n")
cat("Table B (qPCR-positive)         : n =", nrow(data_pos_clean), "\n")

cat("\nBy centre:\n")
print(table(data_pos_clean$centro, useNA = "ifany"))

cat("\nVector control:\n")
print(table(data_pos_clean$vector_control, useNA = "ifany"))

if (nrow(data_pos_clean) != 394) {
  warning(
    "Table B has ", nrow(data_pos_clean), " rows, not 394. ",
    "Check the filters against the analysis scripts."
  )
}

n_undet <- sum(data_pos_clean$vector_control == "Undetermined", na.rm = TRUE)
cat("\nUndetermined vector control:", n_undet, "of", nrow(data_pos_clean), "\n")
cat("This level is retained rather than treated as missing; it is\n")
cat("commoner in multigravidae at every site and is noted in the\n")
cat("limitations.\n")


# ============================================================
# Variables
# ============================================================

vars_cont <- c(
  "PfLDH_pg.ml_tbl",
  "PfHRP2_pg.ml_tbl",
  "density",
  "gestweeks",
  "antibody_pregnancy",
  "antibody_cumulative"
)

vars_cat <- c(
  "gravidity",
  "hiv",
  "age_group",
  "vector_control",
  "Season",
  "Period",
  "sede"
)


# ============================================================
# Labels
# ============================================================

var_labels <- list(
  
  PfLDH_pg.ml_tbl ~ "Pf-LDH (pg/mL; only if detected)",
  
  PfHRP2_pg.ml_tbl ~ paste0(
    "PfHRP2 (pg/mL; only if detected; capped at ",
    pfhrp2_cap_val, " if > ", pfhrp2_sat_thr_fmt, ")"
  ),
  
  density ~ "Parasite density (parasites/µL)",
  
  gestweeks ~ "Gestational week",
  
  antibody_pregnancy ~ "VAR2CSA IgG (pregnancy-specific)",
  
  antibody_cumulative ~ "Cumulative IgG",
  
  gravidity ~ "Gravidity",
  
  hiv ~ "HIV status",
  
  age_group ~ "Age group",
  
  vector_control ~ "Vector control measures reported",
  
  Season ~ "Season",
  
  Period ~ "Period",
  
  sede ~ "Type of setting"
)


# ============================================================
# gtsummary themes
# ============================================================

gtsummary::theme_gtsummary_journal(journal = "jama")
gtsummary::theme_gtsummary_compact()


# ============================================================
# Build table
# ============================================================

build_tbl <- function(df_clean) {
  
  tbl <- df_clean %>%
    
    select(centro, all_of(vars_cont), all_of(vars_cat)) %>%
    
    gtsummary::tbl_summary(
      
      by = centro,
      
      type = list(
        all_continuous() ~ "continuous",
        all_categorical() ~ "categorical"
      ),
      
      statistic = list(
        all_continuous() ~ "{median} ({min}-{max})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      
      percent = "column",
      
      missing = "no",
      
      label = var_labels,
      
      digits = list(
        all_continuous() ~ 2,
        all_categorical() ~ c(0, 1)
      )
    ) %>%
    
    gtsummary::add_p(
      
      test = list(
        all_continuous() ~ "kruskal.test",
        all_categorical() ~ fisher_if_needed
      ),
      
      pvalue_fun = ~ gtsummary::style_pvalue(.x, digits = 3)
    ) %>%
    
    gtsummary::bold_labels() %>%
    
    gtsummary::bold_p(t = 0.05)
  
  tbl %>%
    replace_stats_with_geomean(df_clean) %>%
    fix_pfhrp2_max_display(df_clean)
}


# ============================================================
# Generate both tables
# ============================================================

tbl_all <- build_tbl(data_all_clean)
tbl_pos <- build_tbl(data_pos_clean)


# ============================================================
# Export to editable Word with BOTH tables
# ============================================================

ft_all <- gtsummary::as_flex_table(tbl_all) %>%
  flextable::autofit() %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::align(align = "center", part = "header")

ft_pos <- gtsummary::as_flex_table(tbl_pos) %>%
  flextable::autofit() %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::align(align = "center", part = "header")


# ============================================================
# Create Word document
# ============================================================

doc <- officer::read_docx()

ab_note <- if (antibody_log_plus1) {
  "Antibodies back-transformed assuming log10(x+1): original = 10^log - 1."
} else {
  "Antibodies back-transformed assuming log10(x): original = 10^log."
}

note_txt <- paste0(
  "For PfHRP2 and PfLDH, concentrations were excluded ",
  "(set to missing) when the corresponding detection flag was ",
  "not detected; such infections are retained in the table and ",
  "in the models, where undetectable antigen is handled by left ",
  "censoring. ",
  "Zeros in continuous variables were excluded from summaries ",
  "and tests. ",
  "PfLDH, PfHRP2 and parasite density are shown as geometric ",
  "mean (min-max) computed via log10. ",
  "VAR2CSA and cumulative antibodies are shown as geometric ",
  "mean (min-max) after back-transformation. ",
  ab_note, " ",
  "Vector control counts how many measures were reported ",
  "(insecticide-treated net, indoor residual spraying), with ",
  "Undetermined where neither was established; this is the ",
  "variable used in the regression models. ",
  "PfHRP2 saturation: values > ", pfhrp2_sat_thr_fmt,
  " were capped to ", pfhrp2_cap_val,
  " for summary statistics; where saturation occurred in a ",
  "centre the maximum is shown as >", pfhrp2_sat_thr_fmt, "."
)


# ---- Table A ----

doc <- officer::body_add_par(
  doc, "Table A. First ANC - all samples", style = "heading 1"
)

doc <- officer::body_add_par(doc, note_txt, style = "Normal")

doc <- flextable::body_add_flextable(doc, value = ft_all)

doc <- officer::body_add_par(doc, "", style = "Normal")


# ---- Table B ----

doc <- officer::body_add_par(
  doc, "Table B. First ANC - qPCR-positive infections", style = "heading 1"
)

doc <- officer::body_add_par(
  doc,
  paste0(
    "All qPCR-positive infections at first antenatal care (n = ",
    nrow(data_pos_clean), "), which is the set entering the ",
    "regression models."
  ),
  style = "Normal"
)

doc <- officer::body_add_par(doc, note_txt, style = "Normal")

doc <- flextable::body_add_flextable(doc, value = ft_pos)


# ============================================================
# Save Word
# ============================================================

print(doc, target = out("tables_by_centro_firstANC_ALL_and_PCRPOS.docx"))


# ============================================================
# Optional: HTML + PNG exports
# ============================================================

tbl_all_gt <- tbl_all %>%
  gtsummary::as_gt() %>%
  gt::tab_header(
    title = gt::md("**Table A. First ANC - all samples**"),
    subtitle = gt::md(
      "Antigens excluded if not detected; zeros excluded; HRP2 capped if saturated."
    )
  )

tbl_pos_gt <- tbl_pos %>%
  gtsummary::as_gt() %>%
  gt::tab_header(
    title = gt::md("**Table B. First ANC - qPCR-positive infections**"),
    subtitle = gt::md(
      "All qPCR-positive infections, the set entering the regression models."
    )
  )

gt::gtsave(tbl_all_gt, out("table_A_firstANC_all.html"))
gt::gtsave(tbl_pos_gt, out("table_B_firstANC_pcrpos.html"))

if (!requireNamespace("webshot2", quietly = TRUE)) {
  install.packages("webshot2")
}

gt::gtsave(tbl_all_gt, out("table_A_firstANC_all.png"))
gt::gtsave(tbl_pos_gt, out("table_B_firstANC_pcrpos.png"))


# ============================================================
# Finished
# ============================================================

cat("\n================ FILES WRITTEN ================\n")
cat("Folder:", OUT_DIR, "\n")
for (f in c("tables_by_centro_firstANC_ALL_and_PCRPOS.docx",
            "table_A_firstANC_all.html", "table_B_firstANC_pcrpos.html",
            "table_A_firstANC_all.png",  "table_B_firstANC_pcrpos.png")) {
  cat(sprintf("  %-50s %s\n", f,
              if (file.exists(out(f))) "ok" else "NOT WRITTEN"))
}
cat("\nThe file:/// line printed by gtsave is webshot2's temporary\n")
cat("HTML, not the saved image. The PNG is in the folder above.\n")