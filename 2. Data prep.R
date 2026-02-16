## The aim of this script will be to rename variables of epi data
## for data interpretation before merging with interpolated variables

# Set the working directory
setwd("C:/HRP2 analysis/PhD_HRP2_qPCR")

# Install pacman if not already installed
if (!require("pacman")) install.packages("pacman")

# Load the necessary libraries using pacman
pacman::p_load(
  dplyr, ggplot2, gridExtra, gtsummary, purrr, corrplot, tidyr, haven
)

# Read the data
data <- read.csv("ANC_concentrations.csv", stringsAsFactors = FALSE)

# Read Stata file and recode 'sede'
sede <- read_dta("sede_Helena.dta") %>%
  mutate(
    sede = as.numeric(sede),
    sede = case_when(
      sede == 0 ~ "rural",
      sede == 1 ~ "urban",
      TRUE      ~ NA_character_
    )
  )

# Keep only nida + sede and harmonise type of nida
sede_join <- sede %>%
  select(nida, sede) %>%
  mutate(nida = as.character(nida))

# Make sure nida in 'data' is the same type and join
data <- data %>%
  mutate(nida = as.character(nida)) %>%
  left_join(sede_join, by = "nida")

min_count <- 20

## Coerce counts to numeric (robust to character import)
data <- data %>%
  mutate(
    count_pfhrp2 = suppressWarnings(as.numeric(count_pfhrp2)),
    count_pfldh  = suppressWarnings(as.numeric(count_pfldh))
  )

## Detect value columns robustly (kept for compatibility)
pick_value_col <- function(df, preferred, pattern, count_col, exclude_pattern = NULL) {
  if (preferred %in% names(df)) return(preferred)
  
  hits <- names(df)[grepl(pattern, names(df), ignore.case = TRUE)]
  hits <- setdiff(hits, count_col)
  hits <- hits[!grepl("count", hits, ignore.case = TRUE)]
  
  if (!is.null(exclude_pattern)) {
    hits <- hits[!grepl(exclude_pattern, hits, ignore.case = TRUE)]
  }
  
  if (length(hits) == 1) return(hits)
  
  if (length(hits) > 1) {
    message(
      "NOTE: Multiple candidate value columns for ", preferred, ": ",
      paste(hits, collapse = ", "),
      ". Using: ", hits[1]
    )
    return(hits[1])
  }
  
  stop(
    "Could not find a value column for ", preferred,
    ". Expected '", preferred, "' or something matching pattern: ", pattern
  )
}

pfhrp2_val <- pick_value_col(
  data,
  preferred = "PfHRP2_pg.ml",
  pattern   = "pfhrp2|hrp2",
  count_col = "count_pfhrp2"
)

pfldh_val <- pick_value_col(
  data,
  preferred = "PfLDH_pg.ml",
  pattern   = "pfldh|\\bpfldh\\b|ldh",
  count_col = "count_pfldh",
  exclude_pattern = "pvldh"
)

## ---------- QC: KEEP ROWS, set antigen value(s) to NA when count < min_count ----------
before_n <- nrow(data)

will_low_pfhrp2 <- !is.na(data$count_pfhrp2) & data$count_pfhrp2 < min_count
will_low_pfldh  <- !is.na(data$count_pfldh)  & data$count_pfldh  < min_count

n_low_pfhrp2 <- sum(will_low_pfhrp2)
n_low_pfldh  <- sum(will_low_pfldh)
n_low_both   <- sum(will_low_pfhrp2 & will_low_pfldh)

## Optional: keep affected rows for audit (NOT removed)
affected_rows <- data[will_low_pfhrp2 | will_low_pfldh, ]

## Set corresponding antigen concentration values to NA (rows are kept)
data <- data %>%
  mutate(
    # optional flags to track what happened
    qc_pfhrp2_count_lt20 = if_else(will_low_pfhrp2, "Yes", "No"),
    qc_pfldh_count_lt20  = if_else(will_low_pfldh,  "Yes", "No"),
    
    # ensure value cols are numeric and set NA only where count < min_count
    "{pfhrp2_val}" := if_else(
      will_low_pfhrp2,
      NA_real_,
      suppressWarnings(as.numeric(.data[[pfhrp2_val]]))
    ),
    "{pfldh_val}" := if_else(
      will_low_pfldh,
      NA_real_,
      suppressWarnings(as.numeric(.data[[pfldh_val]]))
    )
  )

message(
  "QC applied (min_count=", min_count, "): antigen values set to NA (rows kept). ",
  "PfHRP2 low-count rows: ", n_low_pfhrp2 - n_low_both,
  "; PfLDH low-count rows: ", n_low_pfldh - n_low_both,
  "; both: ", n_low_both,
  "; total rows in data: ", before_n, "."
)

# ----------- Create antibody composite score -----------

safe_log10p1 <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  suppressWarnings(ifelse(is.na(x_num) | x_num < 0, NA_real_, log10(x_num + 1)))
}

# df: data frame
# antigen_vec: vector of antigen column names
# out_name: name of composite output column
# min_count: minimum per-antigen count to be included in the composite for that row
# center: if TRUE, mean-center the composite (subtract global mean); no scaling
# keep_log_cols: if TRUE, keep the temporary log_ columns in output
make_composite <- function(df, antigen_vec, out_name, min_count = 20,
                           center = FALSE, keep_log_cols = FALSE) {
  
  v <- intersect(antigen_vec, names(df))
  if (!length(v)) {
    df[[out_name]] <- NA_real_
    return(df)
  }
  
  count_cols <- paste0("count_", v)
  
  missing_counts <- setdiff(count_cols, names(df))
  if (length(missing_counts) > 0) {
    df[missing_counts] <- NA_real_
  }
  
  log_cols <- paste0("log_", v)
  
  out <- df %>%
    mutate(across(all_of(v), safe_log10p1, .names = "log_{.col}")) %>%
    rowwise() %>%
    mutate(
      "{out_name}" := {
        vals   <- c_across(all_of(log_cols))
        counts <- c_across(all_of(count_cols))
        
        use <- !is.na(vals) & !is.na(counts) & counts >= min_count
        if (any(use)) mean(vals[use]) else NA_real_
      }
    ) %>%
    ungroup()
  
  if (isTRUE(center)) {
    m <- mean(out[[out_name]], na.rm = TRUE)
    if (is.finite(m)) out[[out_name]] <- out[[out_name]] - m
  }
  
  if (!isTRUE(keep_log_cols)) {
    out <- out %>% select(-all_of(log_cols))
  }
  
  out
}

# ---- Specify your antigen panels ----

pregnancy_antigens <- c("dbl34", "dbl6e")
if (!"antibody_pregnancy" %in% names(data)) {
  data <- make_composite(data, pregnancy_antigens, "antibody_pregnancy", min_count = 20)
}

protection_antigens <- c("eba175", "pfrh2", "pfrh5")
if (!"antibody_protection" %in% names(data)) {
  data <- make_composite(data, protection_antigens, "antibody_protection", min_count = 20)
}

cumulative_antigens <- c("pama1", "msp1")
if (!"antibody_cumulative" %in% names(data)) {
  data <- make_composite(data, cumulative_antigens, "antibody_cumulative", min_count = 20)
}

###### recode variables to show real names

# Gravidity + gestational trimester (First trimester as reference)
data_recoded <- data %>%
  mutate(
    gravidity = ifelse(
      gestnum == 1, "primigravidae",
      ifelse(gestnum > 1, "multigravidae", NA_character_)
    ),
    
    gestweeks_num = suppressWarnings(as.numeric(as.character(gestweeks))),
    
    gest_trimester = case_when(
      is.na(gestweeks_num)                      ~ NA_character_,
      gestweeks_num >= 1  & gestweeks_num <= 12 ~ "First trimester",
      gestweeks_num >= 13 & gestweeks_num <= 27 ~ "Second trimester",
      gestweeks_num >= 28                       ~ "Third trimester",
      TRUE                                      ~ NA_character_
    ),
    
    gest_trimester = factor(
      gest_trimester,
      levels = c("First trimester", "Second trimester", "Third trimester")
    )
  ) %>%
  select(-gestweeks_num)

# Recode variable "centro"
data_recoded <- data_recoded %>%
  mutate(
    centro = dplyr::recode_factor(
      as.character(centro),
      `1` = "Manhiça",
      `2` = "Ilha Josina",
      `3` = "Magude",
      .default = NA_character_
    )
  )

# Recode HIV (0/1)
data_recoded <- data_recoded %>%
  mutate(
    hiv = dplyr::recode_factor(
      as.character(hiv),
      `0` = "Negative",
      `1` = "Positive",
      .default = NA_character_
    )
  )

# Recode IPTp use
data_recoded <- data_recoded %>%
  mutate(
    tiptime_numeric = suppressWarnings(as.numeric(as.character(tiptime))),
    IPTp_use = case_when(
      is.na(tiptime_numeric)          ~ NA_character_,
      tiptime_numeric %in% c(0, 1, 2) ~ "IPTp3-",
      tiptime_numeric > 2             ~ "IPTp3+",
      TRUE                            ~ NA_character_
    )
  ) %>%
  select(-tiptime_numeric)

# Replace NA values in IPTp_use with blank strings
data_recoded <- data_recoded %>%
  mutate(IPTp_use = replace_na(IPTp_use, ""))

# Recode seasonality
data_recoded <- data_recoded %>%
  mutate(
    Season = dplyr::recode_factor(
      as.character(rainy),
      `0` = "Dry season",
      `1` = "Rainy season",
      .default = NA_character_
    )
  )

# Recode period
data_recoded <- data_recoded %>%
  mutate(
    Period = dplyr::recode_factor(
      as.character(period_proj),
      `0` = "Nov. 2016–Oct. 2017",
      `1` = "Nov. 2017-Oct. 2018 ",
      `2` = "Nov. 2018-Oct. 2019",
      .default = NA_character_
    )
  )

# Recode hemoglobin (anemia)
data_recoded <- data_recoded %>%
  mutate(
    anemia = ifelse(!is.na(hemog) & hemog < 11, "1",
                    ifelse(!is.na(hemog), "0", NA_character_))
  )

# Recode age into a categorical variable
data_recoded <- data_recoded %>%
  mutate(
    age_group = ifelse(!is.na(age) & age < 18, "<18",
                       ifelse(!is.na(age), "18 and older", NA_character_))
  )

## add saturation vars (kept commented as in your original script)
# data_recoded <- data_recoded %>%
#   dplyr::mutate(
#     pfhrp2_saturation = dplyr::if_else(is.na(pfhrp2), NA_character_,
#                                       dplyr::if_else(pfhrp2 > 18235.8, "Yes", "No")),
#     pfldh_saturation  = dplyr::if_else(is.na(pfldh), NA_character_,
#                                       dplyr::if_else(pfldh > 20637.9, "Yes", "No")),
#     pvldh_saturation  = dplyr::if_else(is.na(pvldh), NA_character_,
#                                       dplyr::if_else(pvldh > 21472.2, "Yes", "No"))
#   )

write.csv(
  data_recoded,
  file = "C:/HRP2 analysis/PhD_HRP2_qPCR/ANC_recoded.csv",
  row.names = FALSE
)
