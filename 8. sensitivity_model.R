## ============================================================
## SCRIPT 1 — Antigens as outcomes (SENSITIVITY) — EDITED
## HRP2 / PfLDH models + forest plots + “hasta/desde” curves
##
## CHANGES REQUESTED:
##   1) KEEP ONLY PCR+  (pcrpos == 1)
##   2) USE ONLY DENSITY VALUES COMMON ACROSS GROUPS
##      (common support of c_log_density across gravidity × centro)
##   3) HRP2 SATURATION: DO NOT DROP SATURATED VALUES
##      -> Treat pfhrp2_saturation == "Yes" as RIGHT-CENSORED at a threshold (pg/mL)
##
## Everything else: identical to your script.
## ============================================================

## ---------- Setup ----------
setwd("C:/HRP2 analysis/PhD_HRP2_qPCR")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr, readr, tidyr,
  ggplot2, ggtext, patchwork,
  broom, emmeans, car, scales,
  survival
)

burgundy <- "#800020"
blue     <- "#1f78b4"

## ---------- Read data ----------
data <- read.csv("ANC_recoded.csv")

## ============================================================
## IMPORTANT: HRP2 saturation threshold in pg/mL
## - This MUST be the concentration value that corresponds to pfhrp2_saturation == "Yes"
## - Saturated samples are treated as: true HRP2 >= HRP2_SAT_PGML
## ============================================================
HRP2_SAT_PGML <- 5829.7016  # <-- REPLACE with your exact saturation pg/mL threshold

## ============================================================
## 0) Helper: common support filter (c_log_density) across groups
## ============================================================
apply_common_support <- function(df,
                                 x = "c_log_density",
                                 group_vars = c("gravidity", "centro"),
                                 probs = c(0.05, 0.95)) {
  stopifnot(x %in% names(df))
  stopifnot(all(group_vars %in% names(df)))
  stopifnot(length(probs) == 2)
  
  rng <- df %>%
    dplyr::filter(is.finite(.data[[x]]), !is.na(.data[[x]])) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      q_low  = as.numeric(stats::quantile(.data[[x]], probs[1], na.rm = TRUE)),
      q_high = as.numeric(stats::quantile(.data[[x]], probs[2], na.rm = TRUE)),
      n_grp  = dplyr::n(),
      .groups = "drop"
    )
  
  if (nrow(rng) < 2) {
    warning("Common support: <2 groups after filtering; skipping common-support restriction.")
    return(list(df = df, c_low = NA_real_, c_high = NA_real_, rng = rng))
  }
  
  c_low  <- max(rng$q_low,  na.rm = TRUE)
  c_high <- min(rng$q_high, na.rm = TRUE)
  
  if (!is.finite(c_low) || !is.finite(c_high) || c_low >= c_high) {
    warning("Common support empty (c_low >= c_high); skipping common-support restriction.")
    return(list(df = df, c_low = c_low, c_high = c_high, rng = rng))
  }
  
  df2 <- df %>% dplyr::filter(.data[[x]] >= c_low, .data[[x]] <= c_high)
  list(df = df2, c_low = c_low, c_high = c_high, rng = rng)
}

## ============================================================
## 1) Data preparation
## ============================================================
data <- data %>%
  dplyr::filter(visit == "first ANC") %>%
  dplyr::mutate(
    ## fix encoding issues
    centro = as.character(centro),
    centro = dplyr::recode(centro, "ManhiÃ§a" = "Manhiça", .default = centro),
    
    ## ensure PCR flag numeric then filter later
    pcrpos = suppressWarnings(as.numeric(as.character(pcrpos))),
    
    ## ensure numeric outcomes/predictors
    PfHRP2_pg.ml = suppressWarnings(as.numeric(PfHRP2_pg.ml)),
    PfLDH_pg.ml  = suppressWarnings(as.numeric(PfLDH_pg.ml)),
    density      = suppressWarnings(as.numeric(density)),
    
    ## saturation flag (robust to case)
    hrp2_sat = tolower(as.character(pfhrp2_saturation)) %in% c("yes","y","true","1"),
    
    ## cap saturated HRP2 at threshold for modeling on log-scale
    PfHRP2_pg.ml_cens = ifelse(hrp2_sat, HRP2_SAT_PGML, PfHRP2_pg.ml),
    
    ## outcomes
    log_pfhrp2 = log1p(PfHRP2_pg.ml_cens),
    log_pfldh  = log1p(PfLDH_pg.ml),
    
    ## censoring indicator for Surv():
    ## 1 = observed (reliable); 0 = right-censored (saturated)
    hrp2_event = as.integer(ifelse(hrp2_sat, 0, 1)),
    
    ## density + centering
    log_density   = log1p(density),
    c_log_density = log_density - mean(log_density, na.rm = TRUE),
    
    ## covariates
    gestweeks = suppressWarnings(as.numeric(as.character(gestweeks))),
    age       = suppressWarnings(as.numeric(as.character(age))),
    
    net   = dplyr::na_if(as.character(net),   "999"),
    fumig = dplyr::na_if(as.character(fumig), "999")
  ) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::any_of(c("age_group","Season","Period","net","fumig","centro","gravidity","sede","hiv")),
      as.factor
    )
  ) %>%
  ## references (kept as you had; comparisons do not depend on reference)
  dplyr::mutate(
    centro    = if ("Ilha Josina" %in% levels(centro)) stats::relevel(centro, ref = "Ilha Josina") else centro,
    gravidity = if ("multigravidae" %in% levels(gravidity)) stats::relevel(gravidity, ref = "multigravidae") else gravidity
  ) %>%
  ## CHANGE 1) PCR+ only
  dplyr::filter(pcrpos == 1)

if (!("gestweeks" %in% names(data))) stop("gestweeks column is missing; cannot run sensitivity model.")
if (all(is.na(data$gestweeks)))      stop("gestweeks is all NA after recoding; cannot run sensitivity model.")
if (nrow(data) == 0)                 stop("No rows after filtering to PCR+ (pcrpos == 1).")

## CHANGE 2) restrict to COMMON SUPPORT of c_log_density across gravidity × centro
cs <- apply_common_support(data, x = "c_log_density", group_vars = c("gravidity","centro"), probs = c(0.05, 0.95))
data <- cs$df
COMMON_C_LOW  <- cs$c_low
COMMON_C_HIGH <- cs$c_high

cat("\n================ DATA PREP =================\n")
cat("Rows (PCR+ & common support):", nrow(data), "\n")
cat("Centres:", paste(levels(droplevels(data$centro)), collapse = ", "), "\n")
cat("Gravidity:", paste(levels(droplevels(data$gravidity)), collapse = ", "), "\n")
if (is.finite(COMMON_C_LOW) && is.finite(COMMON_C_HIGH) && COMMON_C_LOW < COMMON_C_HIGH) {
  cat(sprintf("Common support (c_log_density): [%.4f, %.4f] (5–95%% per group)\n", COMMON_C_LOW, COMMON_C_HIGH))
} else {
  cat("Common support (c_log_density): not applied (empty/invalid)\n")
}

## ============================================================
## 2) Models (ANCOVA-style)
##    log antigen ~ c_log_density * gravidity * centro + (gestweeks:c_log_density) + covariates
## ============================================================

covars <- c("age_group", "Season", "Period", "gestweeks", "net", "fumig", "sede", "hiv")
covars <- covars[covars %in% names(data)]
covar_str <- if (length(covars)) paste(covars, collapse = " + ") else "1"

## --- SENSITIVITY TERM (only change in original script) ---
sens_int <- "gestweeks:c_log_density"

## Build RHS once (used by both models)
rhs <- paste0("c_log_density * gravidity * centro + ", sens_int, " + ", covar_str)

## PfLDH: standard LM
form_ldh <- stats::as.formula(paste0("log_pfldh ~ ", rhs))

## HRP2: censored regression (right-censored at HRP2_SAT_PGML for saturated samples)
form_hrp2_cens <- stats::as.formula(
  paste0("survival::Surv(log_pfhrp2, hrp2_event, type = 'right') ~ ", rhs)
)

## Fit models
model_pfhrp2_total <- survival::survreg(form_hrp2_cens, data = data, dist = "gaussian")
model_pfldh_total  <- lm(form_ldh, data = data)

## ============================================================
## 3) Print model info + ANOVA
## ============================================================
n_hrp2 <- nrow(model.frame(model_pfhrp2_total))
cat("\n================ MODELS =================\n")
cat("N HRP2 (censored model):", n_hrp2, "\n")
cat("N PfLDH (lm):", nobs(model_pfldh_total), "\n")
cat("AIC HRP2 :", AIC(model_pfhrp2_total), "\n")
cat("AIC PfLDH:", AIC(model_pfldh_total), "\n\n")

cat("\nType II ANOVA HRP2 (survreg):\n");  print(car::Anova(model_pfhrp2_total, type = 2))
cat("\nType II ANOVA PfLDH (lm):\n");     print(car::Anova(model_pfldh_total,  type = 2))

## ============================================================
## 4) Tidy + save CSV
## ============================================================

tidy_survreg <- function(fit) {
  s <- summary(fit)$table
  out <- as.data.frame(s)
  out$term <- rownames(out)
  rownames(out) <- NULL
  out %>%
    dplyr::rename(
      estimate  = Value,
      std.error = `Std. Error`,
      statistic = z,
      p.value   = p
    ) %>%
    dplyr::select(term, estimate, std.error, statistic, p.value)
}

out <- dplyr::bind_rows(
  tidy_survreg(model_pfhrp2_total) %>%
    mutate(model = "log_pfhrp2_censored", AIC = AIC(model_pfhrp2_total), n = n_hrp2),
  broom::tidy(model_pfldh_total) %>%
    mutate(model = "log_pfldh", AIC = AIC(model_pfldh_total), n = nobs(model_pfldh_total))
)

write.csv(
  out,
  "SENSITIVITY_MODEL_RESULTS_HRP2_PfLDH_PCRpos_CommonSupport.csv",
  row.names = FALSE
)
cat("Saved: SENSITIVITY_MODEL_RESULTS_HRP2_PfLDH_PCRpos_CommonSupport.csv\n")

## ============================================================
## 3) CONTRASTS: Primi-centre anchored vs ALL groups (no duplicates)
## ============================================================

grav_short <- function(g) {
  g <- as.character(g)
  ifelse(g == "primigravidae", "Primi",
         ifelse(g == "multigravidae", "Multi", g))
}

make_methods_primi_vs_all_nodup <- function(emm_obj,
                                            grav_var   = "gravidity",
                                            centro_var = "centro",
                                            primi_lvl  = "primigravidae") {
  emm_df <- as.data.frame(emm_obj)
  K <- nrow(emm_df)
  if (K < 2) stop("emm_obj has <2 rows; cannot build contrasts.")
  
  key <- paste(as.character(emm_df[[grav_var]]), as.character(emm_df[[centro_var]]), sep = "||")
  idx_map <- seq_len(K); names(idx_map) <- key
  
  primi_centres <- sort(unique(as.character(emm_df[[centro_var]][as.character(emm_df[[grav_var]]) == primi_lvl])))
  if (!length(primi_centres)) stop("No Primi centres found in emm_obj.")
  
  methods <- list()
  labels  <- character()
  idx <- 1L
  
  for (c1 in primi_centres) {
    a_key <- paste(primi_lvl, c1, sep = "||")
    if (!a_key %in% names(idx_map)) next
    a_idx <- idx_map[[a_key]]
    
    for (j in seq_len(K)) {
      if (j == a_idx) next
      
      g2 <- as.character(emm_df[[grav_var]][j])
      c2 <- as.character(emm_df[[centro_var]][j])
      
      ## avoid repetition for Primi vs Primi
      if (g2 == primi_lvl) {
        if (!(c1 < c2)) next
      }
      
      v <- numeric(K)
      v[a_idx] <-  1
      v[j]     <- -1
      
      labels[idx]    <- paste0("Primi ", c1, " - ", grav_short(g2), " ", c2)
      methods[[idx]] <- v
      idx <- idx + 1L
    }
  }
  
  if (!length(methods)) stop("No contrasts built (check levels).")
  names(methods) <- labels
  methods
}

tidy_contrasts <- function(pw_sidak, pw_raw, label) {
  as.data.frame(summary(pw_sidak, infer = c(TRUE, TRUE))) %>%
    dplyr::rename(p_sidak = p.value) %>%
    dplyr::left_join(
      as.data.frame(summary(pw_raw)) %>%
        dplyr::select(contrast, p.value) %>%
        dplyr::rename(p_raw = p.value),
      by = "contrast"
    ) %>%
    dplyr::mutate(model = label)
}

## ============================================================
## 4) LEVELS at a COMMON density point
##    (use c_log_density = 0 only if it lies within common support)
## ============================================================

slope_var <- "c_log_density"

c_at <- 0
if (is.finite(COMMON_C_LOW) && is.finite(COMMON_C_HIGH) && COMMON_C_LOW < COMMON_C_HIGH) {
  if (!(c_at >= COMMON_C_LOW && c_at <= COMMON_C_HIGH)) {
    c_at <- mean(c(COMMON_C_LOW, COMMON_C_HIGH))
  }
}
cat(sprintf("\nEMMEANS evaluated at c_log_density = %.4f\n", c_at))

emm_pfhrp2 <- emmeans(model_pfhrp2_total, ~ gravidity * centro, at = setNames(list(c_at), slope_var))
emm_pfldh  <- emmeans(model_pfldh_total,  ~ gravidity * centro, at = setNames(list(c_at), slope_var))

meth_levels_pfhrp2 <- make_methods_primi_vs_all_nodup(emm_pfhrp2)
meth_levels_pfldh  <- make_methods_primi_vs_all_nodup(emm_pfldh)

pw_pfhrp2_raw   <- contrast(emm_pfhrp2, method = meth_levels_pfhrp2, adjust = "none")
pw_pfhrp2_sidak <- contrast(emm_pfhrp2, method = meth_levels_pfhrp2, adjust = "sidak")

pw_pfldh_raw    <- contrast(emm_pfldh,  method = meth_levels_pfldh,  adjust = "none")
pw_pfldh_sidak  <- contrast(emm_pfldh,  method = meth_levels_pfldh,  adjust = "sidak")

levels_hrp2_df <- tidy_contrasts(pw_pfhrp2_sidak, pw_pfhrp2_raw, "HRP2 level")
levels_ldh_df  <- tidy_contrasts(pw_pfldh_sidak,  pw_pfldh_raw,  "PfLDH level")

write.csv(levels_hrp2_df, "LEVELS_contrasts_HRP2_PrimiAnchorVsAll_Sidak_PCRpos_CommonSupport.csv", row.names = FALSE)
write.csv(levels_ldh_df,  "LEVELS_contrasts_PfLDH_PrimiAnchorVsAll_Sidak_PCRpos_CommonSupport.csv", row.names = FALSE)
cat("Saved: LEVELS_contrasts_*_PCRpos_CommonSupport.csv\n")

## ============================================================
## 5) SLOPES (emtrends) + same comparisons
## ============================================================

tr_hrp2  <- emtrends(model_pfhrp2_total, specs = ~ gravidity * centro, var = slope_var)
tr_pfldh <- emtrends(model_pfldh_total,  specs = ~ gravidity * centro, var = slope_var)

meth_slope_hrp2  <- make_methods_primi_vs_all_nodup(tr_hrp2)
meth_slope_pfldh <- make_methods_primi_vs_all_nodup(tr_pfldh)

pw_hrp2_raw   <- contrast(tr_hrp2,  method = meth_slope_hrp2,  adjust = "none")
pw_hrp2_sidak <- contrast(tr_hrp2,  method = meth_slope_hrp2,  adjust = "sidak")

pw_pfldh_raw   <- contrast(tr_pfldh, method = meth_slope_pfldh, adjust = "none")
pw_pfldh_sidak <- contrast(tr_pfldh, method = meth_slope_pfldh, adjust = "sidak")

slopes_hrp2_df <- tidy_contrasts(pw_hrp2_sidak,  pw_hrp2_raw,  "HRP2 slope")
slopes_ldh_df  <- tidy_contrasts(pw_pfldh_sidak, pw_pfldh_raw, "PfLDH slope")

slopes_all <- bind_rows(slopes_hrp2_df, slopes_ldh_df)
write.csv(slopes_all, "SLOPES_contrasts_PrimiAnchorVsAll_Sidak_PCRpos_CommonSupport.csv", row.names = FALSE)
cat("Saved: SLOPES_contrasts_PrimiAnchorVsAll_Sidak_PCRpos_CommonSupport.csv\n")

## ============================================================
## 6) Forest plot helper (unchanged)
## ============================================================

fmt_p <- function(p, lbl) {
  base <- ifelse(is.na(p), paste0(lbl, " p = NA"), sprintf("%s p = %.3f", lbl, p))
  ifelse(!is.na(p) & p < 0.05, paste0("<b>", base, "</b>"), base)
}

make_panel_forest <- function(df, title, colour_hex, show_y_labels = TRUE, ylim_max = NULL) {
  g <- ggplot(df, aes(x = contrast, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.22, linewidth = 0.85, color = colour_hex) +
    geom_point(shape = 22, size = 4.1, stroke = 0.95, fill = colour_hex, color = "black") +
    ggtext::geom_richtext(aes(y = label_y, label = p_label),
                          hjust = 0, size = 3, fill = NA, label.color = NA, color = "black") +
    coord_flip(clip = "off") +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black"),
      plot.title  = element_text(face = "bold"),
      plot.margin = margin(8, 46, 8, 8)
    )
  
  if (!is.null(ylim_max)) g <- g + scale_y_continuous(limits = c(-ylim_max, ylim_max))
  if (!show_y_labels) g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  g
}

prep_forest_df <- function(df_in) {
  df <- df_in %>%
    mutate(
      contrast = as.character(contrast),
      p_label  = paste0(fmt_p(p_raw, "raw"), "<br/>", fmt_p(p_sidak, "Sidak"))
    )
  ylim_max <- max(abs(c(df$lower.CL, df$upper.CL)), na.rm = TRUE) * 1.20
  span <- 2 * ylim_max
  df <- df %>%
    mutate(label_y = pmin(upper.CL + 0.04 * span, ylim_max - 0.02 * span))
  list(df = df, ylim_max = ylim_max)
}

## ---- Forest: LEVELS ----
lev_hrp2 <- prep_forest_df(levels_hrp2_df)$df
lev_ldh  <- prep_forest_df(levels_ldh_df)$df

contrast_order <- lev_hrp2 %>% arrange(estimate) %>% pull(contrast)
lev_hrp2 <- lev_hrp2 %>% mutate(contrast = factor(contrast, levels = contrast_order))
lev_ldh  <- lev_ldh  %>% mutate(contrast = factor(contrast, levels = contrast_order))

ylim_max_levels <- max(
  abs(c(lev_hrp2$lower.CL, lev_hrp2$upper.CL, lev_ldh$lower.CL, lev_ldh$upper.CL)),
  na.rm = TRUE
) * 1.20

p_levels <- make_panel_forest(lev_hrp2, "HRP2 levels (qPCR based densities as predictor)", burgundy, TRUE, ylim_max_levels) +
  make_panel_forest(lev_ldh,  "PfLDH levels (qPCR based densities as predictor)", blue, FALSE, ylim_max_levels) +
  patchwork::plot_layout(widths = c(1, 1))

cat("\n--- Printing FOREST (LEVELS) ---\n")
print(p_levels)

ggsave("SENSITIVITY_FOREST_LEVELS_HRP2_PfLDH_PrimiAnchorVsAll_PCRpos_CommonSupport.png",
       p_levels, width = 12, height = 6, dpi = 300)

## ---- Forest: SLOPES ----
sl_hrp2 <- slopes_hrp2_df %>% mutate(p_label = paste0(fmt_p(p_raw,"raw"),"<br/>",fmt_p(p_sidak,"Sidak")))
sl_ldh  <- slopes_ldh_df  %>% mutate(p_label = paste0(fmt_p(p_raw,"raw"),"<br/>",fmt_p(p_sidak,"Sidak")))

contrast_order_s <- sl_hrp2 %>% arrange(estimate) %>% pull(contrast)
sl_hrp2 <- sl_hrp2 %>% mutate(contrast = factor(as.character(contrast), levels = contrast_order_s))
sl_ldh  <- sl_ldh  %>% mutate(contrast = factor(as.character(contrast), levels = contrast_order_s))

ylim_max_slopes <- max(
  abs(c(sl_hrp2$lower.CL, sl_hrp2$upper.CL, sl_ldh$lower.CL, sl_ldh$upper.CL)),
  na.rm = TRUE
) * 1.20
span_s <- 2 * ylim_max_slopes

sl_hrp2 <- sl_hrp2 %>% mutate(label_y = pmin(upper.CL + 0.04 * span_s, ylim_max_slopes - 0.02 * span_s))
sl_ldh  <- sl_ldh  %>% mutate(label_y = pmin(upper.CL + 0.04 * span_s, ylim_max_slopes - 0.02 * span_s))

p_slopes <- make_panel_forest(sl_hrp2, "HRP2 slopes (qPCR based densities as predictor)", burgundy, TRUE, ylim_max_slopes) +
  make_panel_forest(sl_ldh,  "PfLDH slopes (qPCR based densities as predictor)", blue, FALSE, ylim_max_slopes) +
  patchwork::plot_layout(widths = c(1, 1))

cat("\n--- Printing FOREST (SLOPES) ---\n")
print(p_slopes)

ggsave("SENSITIVITY_FOREST_SLOPES_HRP2_PfLDH_PrimiAnchorVsAll_PCRpos_CommonSupport.png",
       p_slopes, width = 12, height = 6, dpi = 300)

## ============================================================
## 7) Lines: adjusted antigen–density curves (COMMON SUPPORT)
## ============================================================

mu_logdens <- mean(data$log_density, na.rm = TRUE)

# use a grid strictly within common support
if (is.finite(COMMON_C_LOW) && is.finite(COMMON_C_HIGH) && COMMON_C_LOW < COMMON_C_HIGH) {
  c_seq <- seq(COMMON_C_LOW, COMMON_C_HIGH, length.out = 25)
} else {
  c_vals <- na.omit(model.frame(model_pfhrp2_total)[[slope_var]])
  c_seq  <- unique(quantile(c_vals, probs = seq(0.05, 0.95, by = 0.05), na.rm = TRUE))
}

predict_lines <- function(model, marker_label) {
  out <- lapply(c_seq, function(cc) {
    em <- emmeans(model, specs = ~ gravidity * centro, at = list(c_log_density = cc))
    as.data.frame(em) %>% mutate(c_val = cc)
  }) %>% bind_rows()
  
  out %>%
    mutate(
      marker  = marker_label,
      stratum = interaction(gravidity, centro, sep = " × "),
      density = pmax(exp(c_val + mu_logdens) - 1, 1e-6),
      antigen       = exp(emmean) - 1,
      antigen_lower = exp(lower.CL) - 1,
      antigen_upper = exp(upper.CL) - 1
    )
}

lines_hrp2  <- predict_lines(model_pfhrp2_total, "HRP2 (pg/mL)")
lines_pfldh <- predict_lines(model_pfldh_total,  "PfLDH (pg/mL)")
lines_all   <- bind_rows(lines_hrp2, lines_pfldh)

p_lines <- ggplot(lines_all, aes(x = density, y = antigen, group = stratum)) +
  geom_ribbon(aes(ymin = antigen_lower, ymax = antigen_upper, fill = stratum), alpha = 0.10, color = NA) +
  geom_line(aes(color = stratum), linewidth = 0.9) +
  scale_x_log10(labels = label_number(accuracy = 0.01, big.mark = ",")) +
  scale_y_log10(labels = label_number(accuracy = 0.01, big.mark = ",")) +
  facet_wrap(~ marker, ncol = 2, scales = "free_y") +
  labs(x = "qPCR density (parasites/µL)", y = "Antigen (pg/mL)", color = NULL) +
  guides(fill = "none") +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.key.width = grid::unit(1.4, "lines")
  )

cat("\n--- Printing LINES plot ---\n")
print(p_lines)

ggsave("SENSITIVITY_LINES_adjusted_antigen_density_PCRpos_CommonSupport.png",
       p_lines, width = 11, height = 6.5, dpi = 300)

## ============================================================
## 8) HASTA/DESDE curves — edited to use COMMON SUPPORT grid
## ============================================================

## ---------- Excel-safe text + CSV writers ----------
clean_text <- function(x) {
  if (is.null(x)) return(x)
  x <- as.character(x)
  x <- enc2utf8(x)
  
  x <- gsub("ManhiÃ§a", "Manhiça", x, fixed = TRUE)
  x <- gsub("â€“", "-", x, fixed = TRUE)
  x <- gsub("\u2013", "-", x, fixed = TRUE)
  x <- gsub("\u2014", "-", x, fixed = TRUE)
  x <- gsub("\u2212", "-", x, fixed = TRUE)
  
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- gsub("[[:space:]]+", " ", x)
  x <- gsub("[[:cntrl:]]", "", x)
  
  trimws(x)
}

write_csv_utf8_bom <- function(df, file) {
  con <- file(file, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw("\ufeff"), con)
  utils::write.csv(df, con, row.names = FALSE)
}

find_segments <- function(x, flag) {
  flag <- as.logical(flag)
  flag[is.na(flag)] <- FALSE
  if (!any(flag)) return(data.frame(start = numeric(0), end = numeric(0)))
  r <- rle(flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  idx_true <- which(r$values)
  data.frame(start = x[starts[idx_true]], end = x[ends[idx_true]])
}

segments_to_string <- function(seg) {
  if (nrow(seg) == 0) return("none")
  paste(sprintf("%g-%g", seg$start, seg$end), collapse = "; ")
}

build_primi_anchor_pairs <- function(groups_df, primi_lvl = "primigravidae") {
  groups_df <- groups_df %>% distinct(gravidity, centro)
  groups_df$gravidity <- as.character(groups_df$gravidity)
  groups_df$centro    <- as.character(groups_df$centro)
  
  primi_centres <- sort(unique(groups_df$centro[groups_df$gravidity == primi_lvl]))
  if (!length(primi_centres)) return(data.frame())
  
  out <- list(); k <- 1L
  for (c1 in primi_centres) {
    for (i in seq_len(nrow(groups_df))) {
      g2 <- groups_df$gravidity[i]
      c2 <- groups_df$centro[i]
      if (g2 == primi_lvl && c2 == c1) next
      if (g2 == primi_lvl && !(c1 < c2)) next
      
      out[[k]] <- data.frame(
        g1 = primi_lvl, c1 = c1,
        g2 = g2,        c2 = c2,
        contrast = paste0("Primi ", c1, " - ", grav_short(g2), " ", c2),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  bind_rows(out) %>% distinct(contrast, .keep_all = TRUE)
}

make_newdata_fixed <- function(model, grav, cent, c_var_name, c_value) {
  mf <- model.frame(model)
  nd <- as.list(rep(NA, ncol(mf))); names(nd) <- names(mf)
  
  if (!c_var_name %in% names(nd)) stop("Model missing predictor: ", c_var_name)
  
  nd[[c_var_name]] <- c_value
  if ("gravidity" %in% names(nd)) nd$gravidity <- grav
  if ("centro"    %in% names(nd)) nd$centro    <- cent
  
  for (nm in names(nd)) {
    if (nm == "(Intercept)") next
    if (!is.na(nd[[nm]])) next
    if (!nm %in% names(mf)) next
    
    if (is.factor(mf[[nm]])) {
      tab <- table(mf[[nm]])
      nd[[nm]] <- names(tab)[which.max(tab)]
    } else if (is.numeric(mf[[nm]]) || is.integer(mf[[nm]])) {
      nd[[nm]] <- mean(mf[[nm]], na.rm = TRUE)
    } else {
      nd[[nm]] <- mf[[nm]][which(!is.na(mf[[nm]]))[1]]
    }
  }
  
  nd <- as.data.frame(nd, stringsAsFactors = FALSE)
  
  for (nm in names(nd)) {
    if (nm %in% names(mf) && is.factor(mf[[nm]])) {
      nd[[nm]] <- factor(nd[[nm]], levels = levels(mf[[nm]]))
    }
  }
  nd
}

contrast_curve_exact <- function(model, raw_grid, mu_log_raw, c_var_name, g1, c1, g2, c2, alpha = 0.05) {
  
  beta <- coef(model)
  V0   <- vcov(model)
  
  # ---- FIX CLAVE para survreg(): vcov incluye log(scale); lo eliminamos ----
  # Nos quedamos solo con la parte de V correspondiente a coef(model)
  if (!all(names(beta) %in% colnames(V0))) {
    stop("vcov(model) no contiene todos los nombres de coef(model). Revisa nombres/contrastes.")
  }
  V <- V0[names(beta), names(beta), drop = FALSE]
  # ------------------------------------------------------------------------
  
  tt_noy <- delete.response(terms(model))
  
  diff_hat <- rep(NA_real_, length(raw_grid))
  se_diff  <- rep(NA_real_, length(raw_grid))
  
  for (i in seq_along(raw_grid)) {
    x0 <- raw_grid[i]
    xc <- log1p(x0) - mu_log_raw
    
    nd1 <- make_newdata_fixed(model, grav = g1, cent = c1, c_var_name = c_var_name, c_value = xc)
    nd2 <- make_newdata_fixed(model, grav = g2, cent = c2, c_var_name = c_var_name, c_value = xc)
    
    X1 <- tryCatch(model.matrix(tt_noy, nd1), error = function(e) NULL)
    X2 <- tryCatch(model.matrix(tt_noy, nd2), error = function(e) NULL)
    if (is.null(X1) || is.null(X2)) next
    
    # Alinear exactamente con beta (y con V ya recortada)
    if (!all(names(beta) %in% colnames(X1)) || !all(names(beta) %in% colnames(X2))) next
    X1 <- X1[, names(beta), drop = FALSE]
    X2 <- X2[, names(beta), drop = FALSE]
    
    xd <- X1 - X2
    
    diff_hat[i] <- as.numeric(xd %*% beta)
    se_diff[i]  <- as.numeric(sqrt(xd %*% V %*% t(xd)))
  }
  
  z <- qnorm(1 - alpha/2)
  lower <- diff_hat - z * se_diff
  upper <- diff_hat + z * se_diff
  
  out <- data.frame(x_raw = raw_grid, diff = diff_hat, lower = lower, upper = upper) %>%
    mutate(
      sig_pos = !is.na(lower) & lower > 0,
      sig_neg = !is.na(upper) & upper < 0
    )
  
  seg_pos <- find_segments(out$x_raw, out$sig_pos)
  seg_neg <- find_segments(out$x_raw, out$sig_neg)
  
  list(
    curve = out,
    summary = data.frame(
      last_pos_x  = if (any(out$sig_pos, na.rm = TRUE)) max(out$x_raw[out$sig_pos], na.rm = TRUE) else NA_real_,
      first_neg_x = if (any(out$sig_neg, na.rm = TRUE)) min(out$x_raw[out$sig_neg], na.rm = TRUE) else NA_real_,
      seg_pos     = segments_to_string(seg_pos),
      seg_neg     = segments_to_string(seg_neg),
      x_min_used  = min(raw_grid),
      x_max_used  = max(raw_grid),
      n_grid      = length(raw_grid),
      note        = "sig_pos: lower>0; sig_neg: upper<0; SE via vcov(model) on linear predictor."
    )
  )
}


## ---------- boundary labels: safe when there are NO sig segments ----------
make_shade_boundary_labels <- function(df_one_contrast, unit_suffix) {
  df <- df_one_contrast
  if (!nrow(df)) return(NULL)
  
  sp <- df$sig_pos; sp[is.na(sp)] <- FALSE
  sn <- df$sig_neg; sn[is.na(sn)] <- FALSE
  
  seg_pos <- find_segments(df$x_raw, sp)
  seg_neg <- find_segments(df$x_raw, sn)
  
  ymin <- min(c(df$lower, df$diff, 0), na.rm = TRUE)
  ymax <- max(c(df$upper, df$diff, 0), na.rm = TRUE)
  span <- if (is.finite(ymax - ymin) && (ymax - ymin) > 0) (ymax - ymin) else 1
  
  y_top_base <- ymax - 0.08 * span
  y_bot_base <- ymin + 0.08 * span
  step <- 0.14 * span  # bigger separation
  
  labfun <- scales::label_number(accuracy = 0.01, big.mark = ",")
  
  build_labels_from_segments <- function(seg, y_base, direction = c("down", "up")) {
    direction <- match.arg(direction)
    if (!nrow(seg)) return(NULL)
    
    labs <- dplyr::bind_rows(
      data.frame(x_raw = seg$start, boundary = "start", stringsAsFactors = FALSE),
      data.frame(x_raw = seg$end,   boundary = "end",   stringsAsFactors = FALSE)
    ) %>%
      dplyr::distinct(x_raw, boundary) %>%
      dplyr::arrange(x_raw, boundary)
    
    if (!nrow(labs)) return(NULL)
    
    labs %>%
      dplyr::mutate(
        idx = dplyr::row_number(),
        y   = if (direction == "down") (y_base - (idx - 1) * step) else (y_base + (idx - 1) * step),
        y   = pmin(pmax(y, ymin + 0.03 * span), ymax - 0.03 * span),
        label = paste0(labfun(x_raw), "\n", unit_suffix),
        hjust = 0.5
      ) %>%
      dplyr::select(x_raw, y, label, hjust)
  }
  
  out_pos <- build_labels_from_segments(seg_pos, y_top_base, direction = "down")
  out_neg <- build_labels_from_segments(seg_neg, y_bot_base, direction = "up")
  
  out <- dplyr::bind_rows(out_pos, out_neg)
  
  # IMPORTANT: if no significant segments -> no labels
  if (is.null(out) || !nrow(out)) return(NULL)
  
  out %>%
    dplyr::mutate(contrast = df$contrast[1]) %>%
    dplyr::select(contrast, x_raw, y, label, hjust)
}

## ---------- plot: safe when labels_df is NULL/empty ----------
plot_contrast_curves <- function(curves_df, labels_df, title_txt, colour_hex, x_label, y_expr) {
  
  if (is.null(labels_df) || !nrow(labels_df)) {
    vlines_df <- data.frame(contrast = character(0), x_raw = numeric(0))
    labels_df <- data.frame(contrast = character(0), x_raw = numeric(0), y = numeric(0),
                            label = character(0), hjust = numeric(0))
  } else {
    vlines_df <- labels_df %>% dplyr::distinct(contrast, x_raw)
  }
  
  ggplot(curves_df, aes(x = x_raw)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = colour_hex, alpha = 0.18, colour = NA) +
    geom_ribbon(data = curves_df %>% dplyr::filter(sig_pos),
                aes(ymin = lower, ymax = upper), fill = colour_hex, alpha = 0.30, colour = NA) +
    geom_ribbon(data = curves_df %>% dplyr::filter(sig_neg),
                aes(ymin = lower, ymax = upper), fill = colour_hex, alpha = 0.30, colour = NA) +
    geom_line(aes(y = diff), colour = colour_hex, linewidth = 0.9) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    
    geom_vline(
      data = vlines_df,
      aes(xintercept = x_raw),
      inherit.aes = FALSE,
      linetype = "dotted",
      linewidth = 0.45,
      colour = "black",
      alpha = 0.65
    ) +
    
    geom_text(
      data = labels_df,
      aes(x = x_raw, y = y, label = label, hjust = hjust),
      inherit.aes = FALSE,
      size  = 3.5,
      vjust = 0.5,
      color = "black",
      lineheight = 0.92,
      check_overlap = FALSE
    ) +
    
    scale_x_log10(labels = scales::label_number(accuracy = 0.01, big.mark = ",")) +
    labs(title = title_txt, x = x_label, y = y_expr) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
      strip.text = element_text(face = "bold", size = 11, color = "black"),
      plot.title = element_text(face = "bold", size = 14, color = "black"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

run_hasta_desde_primi_anchor_vs_all <- function(model_obj, marker_tag, colour_hex, y_expr, unit_suffix,
                                                common_c_low = COMMON_C_LOW, common_c_high = COMMON_C_HIGH) {
  cat("\n================ Parasite-density ranges with significant between-group differences:", marker_tag, "(qPCR based densities as predictor) ================\n")
  
  mf_used <- model.frame(model_obj)
  idx <- suppressWarnings(as.integer(rownames(mf_used)))
  d_used <- data[idx, , drop = FALSE]
  
  # exact mu for c_log_density transform (consistent with your previous logic)
  mu_logdens_center <- mean(d_used$log_density - d_used$c_log_density, na.rm = TRUE)
  
  # COMMON SUPPORT grid in raw density space
  if (is.finite(common_c_low) && is.finite(common_c_high) && common_c_low < common_c_high) {
    x_min <- pmax(exp(common_c_low  + mu_logdens_center) - 1, 1e-8)
    x_max <- pmax(exp(common_c_high + mu_logdens_center) - 1, 1e-8)
  } else {
    # fallback if common support not available
    raw_vec <- d_used$density
    raw_vec <- raw_vec[is.finite(raw_vec) & !is.na(raw_vec) & raw_vec > 0]
    x_min <- as.numeric(quantile(raw_vec, 0.01, na.rm = TRUE))
    x_max <- as.numeric(quantile(raw_vec, 0.99, na.rm = TRUE))
  }
  
  if (!is.finite(x_min) || !is.finite(x_max) || x_min >= x_max) {
    cat("WARNING: invalid density grid bounds; skipping.\n")
    return(invisible(NULL))
  }
  
  x_grid <- exp(seq(log(x_min), log(x_max), length.out = 200))
  
  groups <- mf_used %>% distinct(gravidity, centro)
  pairs_df <- build_primi_anchor_pairs(groups)
  pairs_df$contrast <- clean_text(pairs_df$contrast)
  
  cat("Grid:", signif(min(x_grid), 4), "to", signif(max(x_grid), 4), "(n=200)\n")
  cat("Total contrasts:", nrow(pairs_df), "\n")
  
  res_list <- lapply(seq_len(nrow(pairs_df)), function(ii) {
    rr <- contrast_curve_exact(
      model      = model_obj,
      raw_grid   = x_grid,
      mu_log_raw = mu_logdens_center,
      c_var_name = "c_log_density",
      g1 = pairs_df$g1[ii], c1 = pairs_df$c1[ii],
      g2 = pairs_df$g2[ii], c2 = pairs_df$c2[ii]
    )
    curve <- rr$curve %>% mutate(contrast = pairs_df$contrast[ii])
    summ  <- rr$summary %>% mutate(contrast = pairs_df$contrast[ii])
    list(curve = curve, summary = summ)
  })
  
  curves_all <- bind_rows(lapply(res_list, `[[`, "curve"))
  summ_all   <- bind_rows(lapply(res_list, `[[`, "summary"))
  
  curves_all <- curves_all %>% mutate(contrast = clean_text(contrast))
  summ_all   <- summ_all %>%
    mutate(
      contrast = clean_text(contrast),
      seg_pos  = clean_text(seg_pos),
      seg_neg  = clean_text(seg_neg),
      note     = clean_text(note)
    )
  
  curves_file <- paste0("CURVES_", marker_tag, "_hasta_desde_PrimiAnchorVsAll_PCRpos_CommonSupport.csv")
  summ_file   <- paste0("SUMMARY_", marker_tag, "_hasta_desde_PrimiAnchorVsAll_PCRpos_CommonSupport.csv")
  write_csv_utf8_bom(curves_all, curves_file)
  write_csv_utf8_bom(summ_all,   summ_file)
  
  cat("Saved (UTF-8 BOM):\n  -", curves_file, "\n  -", summ_file, "\n")
  
  labels_all <- bind_rows(lapply(split(curves_all, curves_all$contrast),
                                 make_shade_boundary_labels,
                                 unit_suffix = unit_suffix))
  
  cat("\n================ KEY (Ilha Josina) —", marker_tag, "================\n")
  print(summ_all %>% filter(grepl("Ilha Josina", contrast)) %>% arrange(contrast))
  
  curves_IJ <- curves_all %>%
    filter(grepl("Ilha Josina", contrast)) %>%
    mutate(contrast = factor(contrast, levels = unique(contrast)))
  
  labels_IJ <- labels_all %>%
    filter(grepl("Ilha Josina", contrast)) %>%
    mutate(contrast = factor(contrast, levels = levels(curves_IJ$contrast)))
  
  p_IJ <- plot_contrast_curves(
    curves_IJ, labels_IJ,
    title_txt  = paste0("Parasite-density ranges with significant between-group differences — ", marker_tag, " outcome (qPCR based densities as predictor)"),
    colour_hex = colour_hex,
    x_label    = "qPCR parasite density (parasites/µL; log scale)",
    y_expr     = y_expr
  ) + facet_wrap(~ contrast, ncol = 2, scales = "free_y")
  
  cat("\n--- Printing HASTA/DESDE plot:", marker_tag, "---\n")
  print(p_IJ)
  
  out_png <- paste0("HASTA_DESDE_", marker_tag, "_qPCRdensity_PrimiAnchorVsAll_IJ_facets_PCRpos_CommonSupport.png")
  out_pdf <- paste0("HASTA_DESDE_", marker_tag, "_qPCRdensity_PrimiAnchorVsAll_IJ_facets_PCRpos_CommonSupport.pdf")
  ggsave(out_png, p_IJ, width = 11, height = 10, dpi = 300)
  ggsave(out_pdf, p_IJ, width = 11, height = 10)
  
  cat("Saved:\n  -", out_png, "\n  -", out_pdf, "\n")
  
  invisible(list(curves = curves_all, summary = summ_all, labels = labels_all, plot_IJ = p_IJ))
}

## ---- RUN hasta/desde for BOTH outcomes (HRP2 and PfLDH) ----
res_hasta_hrp2  <- run_hasta_desde_primi_anchor_vs_all(
  model_obj   = model_pfhrp2_total,
  marker_tag  = "HRP2",
  colour_hex  = burgundy,
  y_expr      = expression("HRP2 pg/mL; log scale"),
  unit_suffix = "p/µL"
)

res_hasta_pfldh <- run_hasta_desde_primi_anchor_vs_all(
  model_obj   = model_pfldh_total,
  marker_tag  = "PfLDH",
  colour_hex  = blue,
  y_expr      = expression("PfLDH pg/mL; log scale"),
  unit_suffix = "p/µL"
)

## ============================================================
## 9) MAGNITUDE of significant differences for HASTA/DESDE curves
##    (intervals where 95% CI does NOT cross 0)
##
## FIXED:
##   - avoids group_modify() dropping 'contrast'
##   - robust to BOM/encoding header issues (e.g. "ï..contrast")
##
## Saves:
##   MAGNITUDE_HASTA_DESDE_HRP2_PCRpos_CommonSupport.csv
##   MAGNITUDE_HASTA_DESDE_PfLDH_PCRpos_CommonSupport.csv
##   MAGNITUDE_HASTA_DESDE_ALL_PCRpos_CommonSupport.csv
## ============================================================

## ---- helpers (only if not defined earlier) ----
if (!exists("clean_text", inherits = TRUE)) {
  clean_text <- function(x) {
    if (is.null(x)) return(x)
    x <- as.character(x)
    x <- enc2utf8(x)
    x <- gsub("ManhiÃ§a", "Manhiça", x, fixed = TRUE)
    x <- gsub("â€“", "-", x, fixed = TRUE)
    x <- gsub("\u2013", "-", x, fixed = TRUE)
    x <- gsub("\u2014", "-", x, fixed = TRUE)
    x <- gsub("\u2212", "-", x, fixed = TRUE)
    x <- gsub("\u00A0", " ", x, fixed = TRUE)
    x <- gsub("[[:space:]]+", " ", x)
    x <- gsub("[[:cntrl:]]", "", x)
    trimws(x)
  }
}

if (!exists("write_csv_utf8_bom", inherits = TRUE)) {
  write_csv_utf8_bom <- function(df, file) {
    con <- file(file, open = "wb")
    on.exit(close(con), add = TRUE)
    writeBin(charToRaw("\ufeff"), con)  # UTF-8 BOM (Excel friendly)
    utils::write.csv(df, con, row.names = FALSE)
  }
}

trapz <- function(x, y) {
  ok <- is.finite(x) & is.finite(y) & !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  o <- order(x)
  x <- x[o]; y <- y[o]
  sum((y[-1] + y[-length(y)]) / 2 * diff(x))
}

add_segment_id <- function(flag) {
  flag <- as.logical(flag)
  flag[is.na(flag)] <- FALSE
  r <- rle(flag)
  seg_id <- rep(NA_integer_, length(flag))
  idx_end <- cumsum(r$lengths)
  idx_start <- idx_end - r$lengths + 1
  k <- 0L
  for (i in seq_along(r$values)) {
    if (isTRUE(r$values[i])) {
      k <- k + 1L
      seg_id[idx_start[i]:idx_end[i]] <- k
    }
  }
  seg_id
}

## ---- robust header normaliser (BOM + name variants) ----
fix_hasta_curves_cols <- function(df) {
  stopifnot(is.data.frame(df))
  nm <- names(df)
  
  # strip BOM if attached to header
  nm <- gsub("^\ufeff", "", nm)
  nm <- gsub("^ï\\.\\.", "", nm)
  
  nm <- trimws(nm)
  nm_low <- tolower(nm)
  
  nm_low <- gsub("\\.+", "_", nm_low)
  nm_low <- gsub("[^a-z0-9_]", "_", nm_low)
  nm_low <- gsub("_+", "_", nm_low)
  
  names(df) <- nm_low
  
  rename_first <- function(target, candidates) {
    if (target %in% names(df)) return(invisible(NULL))
    cand <- intersect(candidates, names(df))
    if (length(cand) > 0) names(df)[match(cand[1], names(df))] <<- target
    invisible(NULL)
  }
  
  rename_first("contrast", c("contrast","comparison","pair","label","facet"))
  rename_first("x_raw",    c("x_raw","x","xraw","density","dens","qpcr_density","parasite_density"))
  rename_first("diff",     c("diff","estimate","delta","difference"))
  rename_first("lower",    c("lower","lower_cl","lowercl","lwr","lower_ci","lower_cl_","lower_cl__"))
  rename_first("upper",    c("upper","upper_cl","uppercl","upr","upper_ci","upper_cl_","upper_cl__"))
  rename_first("sig_pos",  c("sig_pos","sigpos","pos_sig","sig_positive"))
  rename_first("sig_neg",  c("sig_neg","signeg","neg_sig","sig_negative"))
  
  df
}

## ---- magnitude summariser (split-based, no group_modify) ----
summarise_magnitude_from_curves <- function(curves_df, marker_tag) {
  
  df <- fix_hasta_curves_cols(curves_df)
  
  req <- c("contrast","x_raw","diff","lower","upper","sig_pos","sig_neg")
  miss <- setdiff(req, names(df))
  if (length(miss)) {
    stop(
      "Curves are missing required columns: ", paste(miss, collapse = ", "),
      "\nColumns I see are:\n  ", paste(names(df), collapse = ", ")
    )
  }
  
  df <- df %>%
    dplyr::mutate(
      marker   = marker_tag,
      contrast = clean_text(.data$contrast),
      x_raw    = as.numeric(.data$x_raw),
      diff     = as.numeric(.data$diff),
      lower    = as.numeric(.data$lower),
      upper    = as.numeric(.data$upper),
      sig_pos  = as.logical(.data$sig_pos),
      sig_neg  = as.logical(.data$sig_neg)
    ) %>%
    dplyr::filter(is.finite(x_raw), !is.na(x_raw), x_raw > 0) %>%
    dplyr::arrange(contrast, x_raw)
  
  if (!nrow(df)) return(tibble::tibble())
  
  summarise_one_direction <- function(dsub, contrast_key, flag_col, direction_label) {
    dsub <- dsub %>%
      dplyr::arrange(x_raw) %>%
      dplyr::mutate(
        flag    = as.logical(.data[[flag_col]]),
        log10x  = log10(x_raw),
        ratio_1p = exp(diff),                 # fold-change in (1+marker)
        pct_1p   = (exp(diff) - 1) * 100      # % change in (1+marker)
      )
    
    dsub$flag[is.na(dsub$flag)] <- FALSE
    
    if (!any(dsub$flag, na.rm = TRUE)) {
      return(tibble::tibble())  # no sig region in this direction
    }
    
    # segment IDs for contiguous TRUE runs
    dsub <- dsub %>%
      dplyr::mutate(
        direction = direction_label,
        seg_id = add_segment_id(.data$flag)
      )
    
    # OVERALL (pooled over all sig points)
    overall <- dsub %>%
      dplyr::filter(flag) %>%
      dplyr::summarise(
        marker = marker_tag,
        contrast = contrast_key,
        direction = direction_label,
        level = "overall",
        seg_id = NA_integer_,
        n_points = dplyr::n(),
        density_min = min(x_raw, na.rm = TRUE),
        density_max = max(x_raw, na.rm = TRUE),
        log10_range = max(log10x, na.rm = TRUE) - min(log10x, na.rm = TRUE),
        
        mean_diff_log1p   = mean(diff, na.rm = TRUE),
        median_diff_log1p = stats::median(diff, na.rm = TRUE),
        
        mean_ratio_1p   = mean(ratio_1p, na.rm = TRUE),
        median_ratio_1p = stats::median(ratio_1p, na.rm = TRUE),
        
        mean_pct_1p   = mean(pct_1p, na.rm = TRUE),
        median_pct_1p = stats::median(pct_1p, na.rm = TRUE),
        
        avg_diff_log1p_per_log10 = {
          w <- max(log10x, na.rm = TRUE) - min(log10x, na.rm = TRUE)
          if (is.finite(w) && w > 0) trapz(log10x, diff) / w else NA_real_
        },
        avg_pct_1p_per_log10 = {
          w <- max(log10x, na.rm = TRUE) - min(log10x, na.rm = TRUE)
          if (is.finite(w) && w > 0) trapz(log10x, pct_1p) / w else NA_real_
        }
      )
    
    # PER SEGMENT
    segs <- dsub %>%
      dplyr::filter(!is.na(seg_id)) %>%
      dplyr::group_by(seg_id) %>%
      dplyr::summarise(
        marker = marker_tag,
        contrast = contrast_key,
        direction = direction_label,
        level = "segment",
        n_points = dplyr::n(),
        density_min = min(x_raw, na.rm = TRUE),
        density_max = max(x_raw, na.rm = TRUE),
        log10_range = max(log10x, na.rm = TRUE) - min(log10x, na.rm = TRUE),
        
        mean_diff_log1p   = mean(diff, na.rm = TRUE),
        median_diff_log1p = stats::median(diff, na.rm = TRUE),
        
        mean_ratio_1p   = mean(ratio_1p, na.rm = TRUE),
        median_ratio_1p = stats::median(ratio_1p, na.rm = TRUE),
        
        mean_pct_1p   = mean(pct_1p, na.rm = TRUE),
        median_pct_1p = stats::median(pct_1p, na.rm = TRUE),
        
        avg_diff_log1p_per_log10 = {
          w <- max(log10x, na.rm = TRUE) - min(log10x, na.rm = TRUE)
          if (is.finite(w) && w > 0) trapz(log10x, diff) / w else NA_real_
        },
        avg_pct_1p_per_log10 = {
          w <- max(log10x, na.rm = TRUE) - min(log10x, na.rm = TRUE)
          if (is.finite(w) && w > 0) trapz(log10x, pct_1p) / w else NA_real_
        },
        .groups = "drop"
      )
    
    dplyr::bind_rows(overall, segs)
  }
  
  # split by contrast (keeps contrast key safely)
  spl <- split(df, df$contrast, drop = TRUE)
  
  out <- lapply(names(spl), function(cc) {
    dsub <- spl[[cc]]
    
    dplyr::bind_rows(
      summarise_one_direction(dsub, cc, "sig_pos", "sig_pos"),
      summarise_one_direction(dsub, cc, "sig_neg", "sig_neg")
    )
  })
  
  dplyr::bind_rows(out) %>%
    dplyr::mutate(
      direction = factor(direction, levels = c("sig_pos","sig_neg")),
      contrast  = clean_text(contrast)
    )
}

## ---- load curves (prefer res_* objects; else read CSV) ----
get_curves_safe <- function(res_obj, fallback_file) {
  if (!is.null(res_obj) && is.list(res_obj) &&
      "curves" %in% names(res_obj) && is.data.frame(res_obj$curves)) {
    return(res_obj$curves)
  }
  if (file.exists(fallback_file)) {
    return(read.csv(fallback_file, stringsAsFactors = FALSE, check.names = FALSE))
  }
  stop("Could not find curves: res object missing and file not found: ", fallback_file)
}

curves_hrp2  <- get_curves_safe(res_hasta_hrp2,
                                "CURVES_HRP2_hasta_desde_PrimiAnchorVsAll_PCRpos_CommonSupport.csv")
curves_pfldh <- get_curves_safe(res_hasta_pfldh,
                                "CURVES_PfLDH_hasta_desde_PrimiAnchorVsAll_PCRpos_CommonSupport.csv")

## ---- run + save ----
mag_hrp2  <- summarise_magnitude_from_curves(curves_hrp2,  "HRP2")
mag_pfldh <- summarise_magnitude_from_curves(curves_pfldh, "PfLDH")

write_csv_utf8_bom(mag_hrp2,  "MAGNITUDE_HASTA_DESDE_HRP2_PCRpos_CommonSupport.csv")
write_csv_utf8_bom(mag_pfldh, "MAGNITUDE_HASTA_DESDE_PfLDH_PCRpos_CommonSupport.csv")

mag_all <- dplyr::bind_rows(mag_hrp2, mag_pfldh)
write_csv_utf8_bom(mag_all, "MAGNITUDE_HASTA_DESDE_ALL_PCRpos_CommonSupport.csv")

cat("\nSaved magnitude summaries:\n")
cat(" - MAGNITUDE_HASTA_DESDE_HRP2_PCRpos_CommonSupport.csv\n")
cat(" - MAGNITUDE_HASTA_DESDE_PfLDH_PCRpos_CommonSupport.csv\n")
cat(" - MAGNITUDE_HASTA_DESDE_ALL_PCRpos_CommonSupport.csv\n")

cat("\n--- OVERALL significant intervals (preview) ---\n")
print(
  mag_all %>%
    dplyr::filter(level == "overall") %>%
    dplyr::arrange(marker, contrast, direction) %>%
    utils::head(20)
)

## ============================================================
## Combine SLOPES forest (A) + LINES plot (B) stacked
## Robust: uses the saved PNGs
## ============================================================

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("png", quietly = TRUE)) install.packages("png")

library(ggplot2)
library(patchwork)
library(png)
library(grid)

png_to_gg <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  img  <- png::readPNG(path)
  grob <- grid::rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
  
  ggplot() +
    annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

f_A <- "SENSITIVITY_FOREST_SLOPES_HRP2_PfLDH_PrimiAnchorVsAll_PCRpos_CommonSupport.png"
f_B <- "SENSITIVITY_LINES_adjusted_antigen_density_PCRpos_CommonSupport.png"

pA <- png_to_gg(f_A) + labs(tag = "A")
pB <- png_to_gg(f_B) + labs(tag = "B")

p_AB <- (pA / pB) &
  theme(
    plot.tag = element_text(face = "bold", size = 18),
    plot.tag.position = c(0.02, 0.98)
  )

print(p_AB)

ggsave("SENSITIVITY_FIG_AB_SLOPES_FOREST_plus_LINES.png", p_AB, width = 14, height = 14, dpi = 300)
ggsave("SENSITIVITY_FIG_AB_SLOPES_FOREST_plus_LINES.pdf", p_AB, width = 14, height = 14)

cat("Saved:\n - SENSITIVITY_FIG_AB_SLOPES_FOREST_plus_LINES.png\n - SENSITIVITY_FIG_AB_SLOPES_FOREST_plus_LINES.pdf\n")

