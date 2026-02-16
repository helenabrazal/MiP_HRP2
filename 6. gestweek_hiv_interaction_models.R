## ============================================================
## HRP2 / PfLDH models + forest plots
##  - Reference group in models: multigravidae, Ilha Josina
##  - Contrasts plotted: ALL primigravidae vs multigravidae
##    combinations across centres
##    (only opposite gravidity pairs; no primi–primi or multi–multi)
##  - Adjusted p-values: Sidak for these contrast families
## ============================================================

## ---------- Setup ----------
setwd("C:/HRP2 analysis/PhD_HRP2_qPCR")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr, readr, tidyr,
  ggplot2, ggtext, patchwork,
  broom, emmeans, openxlsx, car
)


## ---------- Read data ----------
data <- read.csv("ANC_recoded.csv") 

## ============================================================
## 1. Data preparation
##    - restrict to first ANC
##    - log-transform density and antigens
##    - center log density
##    - recode 999 in net / fumig to NA
##    - set factor reference levels: centro = Ilha Josina, gravidity = multigravidae
##    - keep observations with any evidence of infection (antigen detected or PCR+)
## ============================================================
data <- data %>%
  filter(visit == "first ANC") %>%
  mutate(
    ## outcomes: include zeros
    log_pfhrp2 = log1p(PfHRP2_pg.ml),
    log_pfldh  = log1p(PfLDH_pg.ml),
    
    ## log density + centering
    log_density   = log1p(density),
    c_log_density = as.numeric(scale(log_density, center = TRUE, scale = FALSE)),
    
    ## gestational age numeric + centered
    gestweeks   = suppressWarnings(as.numeric(as.character(gestweeks))),
    gestweeks_c = as.numeric(scale(gestweeks, center = TRUE, scale = FALSE)),
    
    ## age numeric + centered (if needed later)
    age   = suppressWarnings(as.numeric(as.character(age))),
    age_c = as.numeric(scale(age, center = TRUE, scale = FALSE)),
    
    ## recode 999 in net / fumig to NA
    net   = as.character(net),
    net   = dplyr::na_if(net, "999"),
    fumig = as.character(fumig),
    fumig = dplyr::na_if(fumig, "999"),
    
    ## categorical covariates
    hiv   = as.factor(hiv),
    age_group = as.factor(age_group),
    net       = as.factor(net),
    fumig     = as.factor(fumig),
    centro    = as.factor(centro),
    gravidity = as.factor(gravidity),
    Season    = as.factor(Season),
    Period    = as.factor(Period),
    sede      = as.factor(sede)
  ) %>%
  ## EDIT 1) PCR+ only
  dplyr::filter(pcrpos == 1) %>%
  ## set reference: centro = Ilha Josina, gravidity = primigravidae
  mutate(
    centro    = stats::relevel(centro, ref = "Ilha Josina"),
    gravidity = stats::relevel(gravidity, ref = "primigravidae")
  ) %>%
  ## infection subset
  filter(
    pfhrp2_detection == "detected" |
      pfldh_detection == "detected" |
      pcrpos == 1
  )


## quick sanity check if needed:
# dplyr::count(data, gravidity, centro)

## ============================================================
## 2. Models (ANCOVA-style)
##    log antigen ~ c_log_density * gravidity * centro + covariates
## ============================================================
data_hrp2 <- data %>% filter(pfhrp2_saturation != "Yes")

model_hiv <- lm(
  log_pfhrp2 ~  hiv * c_log_density  +
    Season + Period + gestweeks  + net + fumig + sede + centro + gravidity + age_group,
  data = data_hrp2
)

model_getweeks <- lm(
  log_pfhrp2 ~  gestweeks * c_log_density  +
    Season + Period + hiv  + net + fumig + sede + centro + gravidity+ age_group,
  data = data_hrp2
)


model_hiv_ldh <- lm(
  log_pfldh ~  hiv * c_log_density  +
    Season + Period + gestweeks  + net + fumig + sede + centro + gravidity+ age_group,
  data = data
)

model_getweeks_ldh <- lm(
  log_pfldh ~  gestweeks * c_log_density  +
    Season + Period + hiv  + net + fumig + sede + centro + gravidity+ age_group,
  data = data
)
# Save model summaries and AIC to a tidy dataframe for further inspection
model_hiv_summary <- broom::tidy(model_hiv) %>%
  mutate(model = "hiv_pfhrp2")


model_getweeks_total <- broom::tidy(model_getweeks) %>%
  mutate(model = "gestweeks_hrp2")

model_hiv_ldh <- broom::tidy(model_hiv_ldh) %>%
  mutate(model = "hiv_ldh")


model_getweeks_ldh <- broom::tidy(model_getweeks_ldh) %>%
  mutate(model = "gestweeks_ldh")

# Combine the results into a single dataframe
combined_model_summaries_hiv_gestweeks <- bind_rows(model_hiv_summary, model_getweeks_total, model_hiv_ldh, model_getweeks_ldh)

# Save the results to a CSV file
write.csv(combined_model_summaries_hiv_gestweeks, "model_summaries_hiv_gestweeks.csv", row.names = FALSE)

