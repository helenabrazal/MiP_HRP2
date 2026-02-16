## ============================================================
## SCRIPT 2 — VAR2CSA / Cumulative antibodies as outcomes (ROBUST)
##
## 1) PCR+ only (pcrpos == 1)
## 2) Common support PER predictor across gravidity × centro:
##      - c_log_density  (qPCR models)
##      - c_log_pfhrp2   (HRP2 models)
##      - c_log_pfldh    (PfLDH models)
##    using per-group 5–95% quantiles → intersection
##
## Comparisons (NO repetitions; Primi ALWAYS group 1):
##   - For each Primi centre c1:
##       Primi(c1) - [every other group] (Primi other centres + Multi any centre)
##   - Primi vs Primi only once (alphabetical centre order)
##
## Forest:
##   - ONLY VAR2CSA outcome + HRP2 predictor slopes
##   - uses SAME comparison set as hasta/desde
##
## Hasta/Desde (VAR2CSA only):
##   - for BOTH: HRP2 and qPCR density predictors
##   - MATCHES SCRIPT 1 formatting:
##       * stacked boundary labels (2-line: value \n unit)
##       * dotted vertical lines at boundaries
##       * labels only if significant segments exist
##       * theme_classic, no grid, square borders
##       * Excel-friendly CSV (UTF-8 BOM)
## ============================================================

## ---------- Setup ----------
setwd("C:/HRP2 analysis/PhD_HRP2_qPCR")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  dplyr, readr, tidyr,
  ggplot2, ggtext, patchwork,
  broom, emmeans, car, scales,
  grid
)

ambar <- "#e2943a"
blue  <- "#1f78b4"

## ---------- Read data ----------
data <- read.csv("ANC_recoded.csv", stringsAsFactors = FALSE, check.names = FALSE)

## ============================================================
## 0) Helpers: safe logs + common support + clean text + BOM CSV
## ============================================================
safe_log1p <- function(x) {
  suppressWarnings(ifelse(is.na(x) | x < 0, NA_real_, log1p(x)))
}

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
    warning("Common support: <2 groups with finite values; skipping restriction for ", x)
    return(list(df = df, c_low = NA_real_, c_high = NA_real_, rng = rng))
  }
  
  c_low  <- suppressWarnings(max(rng$q_low,  na.rm = TRUE))
  c_high <- suppressWarnings(min(rng$q_high, na.rm = TRUE))
  
  if (!is.finite(c_low) || !is.finite(c_high) || c_low >= c_high) {
    warning("Common support empty/invalid for ", x, " (c_low >= c_high). Skipping restriction.")
    return(list(df = df, c_low = c_low, c_high = c_high, rng = rng))
  }
  
  df2 <- df %>% dplyr::filter(.data[[x]] >= c_low, .data[[x]] <= c_high)
  list(df = df2, c_low = c_low, c_high = c_high, rng = rng)
}

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
  writeBin(charToRaw("\ufeff"), con)  # UTF-8 BOM
  utils::write.csv(df, con, row.names = FALSE)
}

## ============================================================
## 1) Data preparation (PCR+ only; center predictors on PCR+ base)
## ============================================================
req <- c(
  "visit","centro","gravidity","density","gestweeks","pcrpos",
  "PfHRP2_pg.ml","PfLDH_pg.ml","pfhrp2_saturation",
  "antibody_pregnancy","antibody_cumulative"
)
miss <- setdiff(req, names(data))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

data_base <- data %>%
  dplyr::filter(visit == "first ANC") %>%
  dplyr::mutate(
    row_id = dplyr::row_number(),
    
    centro = as.character(centro),
    centro = dplyr::recode(centro, "ManhiÃ§a" = "Manhiça", .default = centro),
    
    pcrpos    = suppressWarnings(as.numeric(as.character(pcrpos))),
    density   = suppressWarnings(as.numeric(as.character(density))),
    gestweeks = suppressWarnings(as.numeric(as.character(gestweeks))),
    
    log_density = safe_log1p(density),
    log_pfhrp2  = safe_log1p(PfHRP2_pg.ml),
    log_pfldh   = safe_log1p(PfLDH_pg.ml),
    
    net   = dplyr::na_if(as.character(net),   "999"),
    fumig = dplyr::na_if(as.character(fumig), "999"),
    
    hrp2_sat = factor(
      ifelse(tolower(as.character(pfhrp2_saturation)) %in% c("yes","y","true","1"), "Yes", "No"),
      levels = c("No","Yes")
    ),
    
    antibody_pregnancy  = suppressWarnings(as.numeric(antibody_pregnancy)),
    antibody_cumulative = suppressWarnings(as.numeric(antibody_cumulative))
  ) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::any_of(c("age_group","Season","Period","net","fumig","centro","gravidity","sede","hiv")),
      as.factor
    )
  ) %>%
  dplyr::mutate(
    centro    = if ("Ilha Josina" %in% levels(centro)) stats::relevel(centro, ref = "Ilha Josina") else centro,
    gravidity = if ("primigravidae" %in% levels(gravidity)) stats::relevel(gravidity, ref = "primigravidae") else gravidity
  ) %>%
  dplyr::filter(pcrpos == 1) %>%
  dplyr::mutate(
    c_log_density = log_density - mean(log_density, na.rm = TRUE),
    c_log_pfhrp2  = log_pfhrp2  - mean(log_pfhrp2,  na.rm = TRUE),
    c_log_pfldh   = log_pfldh   - mean(log_pfldh,   na.rm = TRUE)
  )

if (nrow(data_base) == 0) stop("No rows after filtering to PCR+ (pcrpos == 1).")

rownames(data_base) <- as.character(data_base$row_id)

cat("\n================ DATA PREP (BASE) =================\n")
cat("Rows (PCR+):", nrow(data_base), "\n")
cat("Centres:", paste(levels(droplevels(data_base$centro)), collapse = ", "), "\n")
cat("Gravidity:", paste(levels(droplevels(data_base$gravidity)), collapse = ", "), "\n")

## ============================================================
## 1b) Apply common support PER predictor (gravidity × centro)
## ============================================================
cs_den  <- apply_common_support(data_base, x = "c_log_density", group_vars = c("gravidity","centro"), probs = c(0.05, 0.95))
cs_hrp2 <- apply_common_support(data_base, x = "c_log_pfhrp2",  group_vars = c("gravidity","centro"), probs = c(0.05, 0.95))
cs_ldh  <- apply_common_support(data_base, x = "c_log_pfldh",   group_vars = c("gravidity","centro"), probs = c(0.05, 0.95))

data_den  <- cs_den$df
data_hrp2 <- cs_hrp2$df %>% dplyr::filter(hrp2_sat == "No")  # keep as you had in SCRIPT 2
data_ldh  <- cs_ldh$df

cat("\n================ COMMON SUPPORT =================\n")
if (is.finite(cs_den$c_low) && is.finite(cs_den$c_high) && cs_den$c_low < cs_den$c_high) {
  cat(sprintf("qPCR common support (c_log_density): [%.4f, %.4f] | rows kept: %d\n",
              cs_den$c_low, cs_den$c_high, nrow(data_den)))
} else {
  cat("qPCR common support not applied (empty/invalid) | rows kept:", nrow(data_den), "\n")
}
if (is.finite(cs_hrp2$c_low) && is.finite(cs_hrp2$c_high) && cs_hrp2$c_low < cs_hrp2$c_high) {
  cat(sprintf("HRP2 common support (c_log_pfhrp2):  [%.4f, %.4f] | rows kept: %d (after dropping sat: %d)\n",
              cs_hrp2$c_low, cs_hrp2$c_high, nrow(cs_hrp2$df), nrow(data_hrp2)))
} else {
  cat("HRP2 common support not applied (empty/invalid) | rows kept:", nrow(data_hrp2), "\n")
}
if (is.finite(cs_ldh$c_low) && is.finite(cs_ldh$c_high) && cs_ldh$c_low < cs_ldh$c_high) {
  cat(sprintf("PfLDH common support (c_log_pfldh):  [%.4f, %.4f] | rows kept: %d\n",
              cs_ldh$c_low, cs_ldh$c_high, nrow(data_ldh)))
} else {
  cat("PfLDH common support not applied (empty/invalid) | rows kept:", nrow(data_ldh), "\n")
}

## ============================================================
## 2) Comparison builder (Primi anchored vs all; no duplicate Primi-Primi)
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
      
      ## avoid duplicate Primi–Primi (keep only c1 < c2)
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

fmt_p <- function(p, lbl) {
  base <- ifelse(is.na(p), paste0(lbl, " p = NA"), sprintf("%s p = %.3f", lbl, p))
  ifelse(!is.na(p) & p < 0.05, paste0("<b>", base, "</b>"), base)
}

make_panel_slopes <- function(df, title, colour_hex, ylim_max = NULL) {
  g <- ggplot(df, aes(x = contrast, y = estimate)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.22, linewidth = 0.85, color = colour_hex) +
    geom_point(shape = 22, size = 4.3, stroke = 0.95, fill = colour_hex, color = "black") +
    ggtext::geom_richtext(aes(y = label_y, label = p_label), hjust = 0, size = 3,
                          fill = NA, label.color = NA, color = "black") +
    coord_flip(clip = "off") +
    labs(
      title = title,
      x = "Contrast (HRP2 pg/mL as predictor)",
      y = "Difference in slopes (group 1 − group 2)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black"),
      plot.title  = element_text(face = "bold"),
      plot.margin = margin(8, 46, 8, 8)
    )
  if (!is.null(ylim_max)) g <- g + scale_y_continuous(limits = c(-ylim_max, ylim_max))
  g
}

## ============================================================
## 3) Hasta/Desde helpers (MATCH SCRIPT 1 formatting)
## ============================================================
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

build_primi_anchor_pairs_from_mf <- function(mf_used, primi_lvl = "primigravidae") {
  groups_df <- mf_used %>% dplyr::distinct(gravidity, centro)
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
  dplyr::bind_rows(out) %>% dplyr::distinct(contrast, .keep_all = TRUE)
}

make_newdata_fixed <- function(model, grav, cent, c_var_name, c_value) {
  mf <- model.frame(model)
  nd <- as.list(rep(NA, ncol(mf))); names(nd) <- names(mf)
  
  if (!c_var_name %in% names(nd)) stop("Model missing predictor: ", c_var_name)
  
  nd[[c_var_name]] <- c_value
  if ("gravidity" %in% names(nd)) nd$gravidity <- grav
  if ("centro"    %in% names(nd)) nd$centro    <- cent
  
  resp_nm <- names(mf)[1]  # response column in model.frame (lm)
  
  for (nm in names(nd)) {
    if (nm == "(Intercept)") next
    if (nm == resp_nm) next
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
  beta   <- coef(model)
  V      <- vcov(model)
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
    
    if (!all(names(beta) %in% colnames(X1))) next
    if (!all(names(beta) %in% colnames(X2))) next
    
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
    dplyr::mutate(
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
  step <- 0.14 * span
  
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
  if (is.null(out) || !nrow(out)) return(NULL)
  
  out %>%
    dplyr::mutate(contrast = df$contrast[1]) %>%
    dplyr::select(contrast, x_raw, y, label, hjust)
}

plot_contrast_curves <- function(curves_df, labels_df, title_txt, colour_hex, x_label, y_label_txt) {
  
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
    labs(title = title_txt, x = x_label, y = y_label_txt) +
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

run_hasta_desde_primi_anchor_vs_all <- function(model_obj,
                                                predictor_tag,
                                                c_var_name,
                                                mu_log_raw,
                                                c_low,
                                                c_high,
                                                unit_suffix,
                                                colour_hex,
                                                file_stub,
                                                key_centre = "Ilha Josina") {
  cat("\n================ Parasite-density ranges with significant between-group differences:", predictor_tag, "(Primi anchored vs all) ================\n")
  
  mf_used <- model.frame(model_obj)
  
  ## grid restricted to common support if available, else fallback 1–99% of model-frame approx
  if (is.finite(c_low) && is.finite(c_high) && c_low < c_high) {
    x_min <- pmax(exp(c_low  + mu_log_raw) - 1, 1e-8)
    x_max <- pmax(exp(c_high + mu_log_raw) - 1, 1e-8)
  } else {
    x_centered <- suppressWarnings(as.numeric(mf_used[[c_var_name]]))
    x_centered <- x_centered[is.finite(x_centered) & !is.na(x_centered)]
    if (length(x_centered) < 20) {
      cat("WARNING: Too few finite values in model.frame for ", c_var_name, ". Skipping.\n", sep = "")
      return(invisible(NULL))
    }
    x_raw_approx <- pmax(exp(x_centered + mu_log_raw) - 1, 1e-8)
    x_min <- as.numeric(stats::quantile(x_raw_approx, 0.01, na.rm = TRUE))
    x_max <- as.numeric(stats::quantile(x_raw_approx, 0.99, na.rm = TRUE))
  }
  
  if (!is.finite(x_min) || !is.finite(x_max) || x_min >= x_max) {
    cat("WARNING: invalid grid bounds; skipping.\n")
    return(invisible(NULL))
  }
  
  x_grid <- exp(seq(log(x_min), log(x_max), length.out = 200))
  
  pairs_df <- build_primi_anchor_pairs_from_mf(mf_used)
  if (nrow(pairs_df) == 0) {
    cat("WARNING: No contrasts available (check groups present in model.frame). Skipping.\n")
    return(invisible(NULL))
  }
  pairs_df$contrast <- clean_text(pairs_df$contrast)
  
  cat("Grid:", signif(min(x_grid), 4), "to", signif(max(x_grid), 4), "(n=200)\n")
  cat("Total contrasts:", nrow(pairs_df), "\n")
  
  res_list <- lapply(seq_len(nrow(pairs_df)), function(ii) {
    rr <- contrast_curve_exact(
      model      = model_obj,
      raw_grid   = x_grid,
      mu_log_raw = mu_log_raw,
      c_var_name = c_var_name,
      g1 = pairs_df$g1[ii], c1 = pairs_df$c1[ii],
      g2 = pairs_df$g2[ii], c2 = pairs_df$c2[ii]
    )
    curve <- rr$curve %>% dplyr::mutate(contrast = pairs_df$contrast[ii])
    summ  <- rr$summary %>% dplyr::mutate(contrast = pairs_df$contrast[ii])
    list(curve = curve, summary = summ)
  })
  
  curves_all <- dplyr::bind_rows(lapply(res_list, `[[`, "curve"))
  summ_all   <- dplyr::bind_rows(lapply(res_list, `[[`, "summary"))
  
  curves_all <- curves_all %>% dplyr::mutate(contrast = clean_text(contrast))
  summ_all <- summ_all %>%
    dplyr::mutate(
      contrast = clean_text(contrast),
      seg_pos  = clean_text(seg_pos),
      seg_neg  = clean_text(seg_neg),
      note     = clean_text(note)
    )
  
  curves_file <- paste0("CURVES_VAR2CSA_", predictor_tag, "_hasta_desde_PrimiAnchorVsAll_", file_stub, ".csv")
  summ_file   <- paste0("SUMMARY_VAR2CSA_", predictor_tag, "_hasta_desde_PrimiAnchorVsAll_", file_stub, ".csv")
  write_csv_utf8_bom(curves_all, curves_file)
  write_csv_utf8_bom(summ_all,   summ_file)
  
  cat("Saved (UTF-8 BOM):\n  -", curves_file, "\n  -", summ_file, "\n")
  
  labels_all <- dplyr::bind_rows(
    lapply(split(curves_all, curves_all$contrast),
           make_shade_boundary_labels,
           unit_suffix = unit_suffix)
  )
  
  cat("\n================ KEY (", key_centre, ") — ", predictor_tag, "================\n", sep = "")
  print(summ_all %>% dplyr::filter(grepl(key_centre, contrast)) %>% dplyr::arrange(contrast))
  
  curves_IJ <- curves_all %>%
    dplyr::filter(grepl(key_centre, contrast)) %>%
    dplyr::mutate(contrast = factor(contrast, levels = unique(contrast)))
  
  labels_IJ <- labels_all %>%
    dplyr::filter(grepl(key_centre, contrast)) %>%
    dplyr::mutate(contrast = factor(contrast, levels = levels(curves_IJ$contrast)))
  
  xlab <- if (predictor_tag == "HRP2") "HRP2 (pg/mL; log scale)" else "qPCR parasite density (parasites/µL; log scale)"
  ylab <- "VAR2CSA IgG MFIs, log scale"
  
  p_IJ <- plot_contrast_curves(
    curves_IJ, labels_IJ,
    title_txt = paste0("Parasite-density ranges with significant between-group differences — VAR2CSA — ", predictor_tag,
                       " predictor"),
    colour_hex = colour_hex,
    x_label = xlab,
    y_label_txt = ylab
  ) + facet_wrap(~ contrast, ncol = 2, scales = "free_y")
  
  print(p_IJ)
  
  out_png <- paste0("HASTA_DESDE_VAR2CSA_", predictor_tag, "_PrimiAnchorVsAll_IJ_facets_", file_stub, ".png")
  out_pdf <- paste0("HASTA_DESDE_VAR2CSA_", predictor_tag, "_PrimiAnchorVsAll_IJ_facets_", file_stub, ".pdf")
  ggsave(out_png, p_IJ, width = 11, height = 10, dpi = 300)
  ggsave(out_pdf, p_IJ, width = 11, height = 10)
  
  cat("Saved:\n  -", out_png, "\n  -", out_pdf, "\n")
  
  invisible(list(curves = curves_all, summary = summ_all, labels = labels_all, plot_IJ = p_IJ))
}

## ============================================================
## 4) Runner per outcome (predictor-specific common-support datasets)
## ============================================================
run_antibody_suite <- function(df_den, df_hrp2, df_ldh,
                               outcome_var, outcome_label, file_stub,
                               cs_den, cs_hrp2, cs_ldh,
                               mu_logdens_center, mu_loghrp2_center) {
  
  cat("\n\n============================================================\n")
  cat("RUNNING OUTCOME:", outcome_label, "(", outcome_var, ")\n")
  cat("============================================================\n")
  
  covars_all <- c("Season", "Period", "gestweeks", "age_group", "net", "fumig", "sede", "hiv")
  covars_den  <- covars_all[covars_all %in% names(df_den)]
  covars_hrp2 <- covars_all[covars_all %in% names(df_hrp2)]
  covars_ldh  <- covars_all[covars_all %in% names(df_ldh)]
  
  covar_str_den  <- if (length(covars_den)  > 0) paste(covars_den,  collapse = " + ") else "1"
  covar_str_hrp2 <- if (length(covars_hrp2) > 0) paste(covars_hrp2, collapse = " + ") else "1"
  covar_str_ldh  <- if (length(covars_ldh)  > 0) paste(covars_ldh,  collapse = " + ") else "1"
  
  f_qpcr <- as.formula(paste0(outcome_var, " ~ c_log_density * gravidity * centro + ", covar_str_den))
  f_hrp2 <- as.formula(paste0(outcome_var, " ~ c_log_pfhrp2 * gravidity * centro + ", covar_str_hrp2))
  f_ldh  <- as.formula(paste0(outcome_var, " ~ c_log_pfldh  * gravidity * centro + ", covar_str_ldh))
  
  m_qpcr <- lm(f_qpcr, data = df_den)
  m_hrp2 <- lm(f_hrp2, data = df_hrp2)
  m_ldh  <- lm(f_ldh,  data = df_ldh)
  
  cat("\nModel N (after common support):\n")
  cat("  qPCR :", nobs(m_qpcr), "\n")
  cat("  HRP2 :", nobs(m_hrp2), "\n")
  cat("  PfLDH:", nobs(m_ldh),  "\n")
  
  ## ---- Save model results ----
  write.csv(
    dplyr::bind_rows(
      broom::tidy(m_qpcr) %>% dplyr::mutate(model = "qpcr", AIC = AIC(m_qpcr), n = nobs(m_qpcr), outcome = outcome_label),
      broom::tidy(m_hrp2) %>% dplyr::mutate(model = "hrp2", AIC = AIC(m_hrp2), n = nobs(m_hrp2), outcome = outcome_label),
      broom::tidy(m_ldh)  %>% dplyr::mutate(model = "ldh",  AIC = AIC(m_ldh),  n = nobs(m_ldh),  outcome = outcome_label)
    ),
    paste0("MODEL_RESULTS_", file_stub, "_PCRpos_CommonSupport.csv"),
    row.names = FALSE
  )
  
  ## ---- Slopes and contrasts (SAME comparison set) ----
  tr_qpcr <- emtrends(m_qpcr, specs = ~ gravidity * centro, var = "c_log_density")
  tr_hrp2 <- emtrends(m_hrp2, specs = ~ gravidity * centro, var = "c_log_pfhrp2")
  tr_ldh  <- emtrends(m_ldh,  specs = ~ gravidity * centro, var = "c_log_pfldh")
  
  meth_qpcr <- make_methods_primi_vs_all_nodup(tr_qpcr)
  meth_hrp2 <- make_methods_primi_vs_all_nodup(tr_hrp2)
  meth_ldh  <- make_methods_primi_vs_all_nodup(tr_ldh)
  
  slopes_qpcr <- tidy_contrasts(
    contrast(tr_qpcr, method = meth_qpcr, adjust = "sidak"),
    contrast(tr_qpcr, method = meth_qpcr, adjust = "none"),
    "Antibodies ~ qPCR slope"
  )
  slopes_hrp2 <- tidy_contrasts(
    contrast(tr_hrp2, method = meth_hrp2, adjust = "sidak"),
    contrast(tr_hrp2, method = meth_hrp2, adjust = "none"),
    "Antibodies ~ HRP2 slope"
  )
  slopes_ldh <- tidy_contrasts(
    contrast(tr_ldh, method = meth_ldh, adjust = "sidak"),
    contrast(tr_ldh, method = meth_ldh, adjust = "none"),
    "Antibodies ~ PfLDH slope"
  )
  
  slopes_all <- dplyr::bind_rows(slopes_qpcr, slopes_hrp2, slopes_ldh) %>%
    dplyr::mutate(outcome = outcome_label)
  
  write.csv(slopes_all,
            paste0("SLOPES_contrasts_PrimiAnchorVsAll_Sidak_", file_stub, "_PCRpos_CommonSupport.csv"),
            row.names = FALSE)
  
  ## ---- Forest plot: ONLY VAR2CSA + HRP2 slopes ----
  p_forest <- NULL
  if (outcome_var == "antibody_pregnancy") {
    forest_df <- slopes_hrp2 %>%
      dplyr::mutate(
        contrast = as.character(contrast),
        p_label  = paste0(fmt_p(p_raw, "raw"), "<br/>", fmt_p(p_sidak, "Sidak"))
      )
    
    contrast_order <- forest_df %>% dplyr::arrange(estimate) %>% dplyr::pull(contrast)
    forest_df <- forest_df %>% dplyr::mutate(contrast = factor(contrast, levels = contrast_order))
    
    ylim_max <- max(abs(c(forest_df$lower.CL, forest_df$upper.CL)), na.rm = TRUE) * 1.20
    span <- 2 * ylim_max
    forest_df <- forest_df %>% dplyr::mutate(label_y = pmin(upper.CL + 0.04 * span, ylim_max - 0.02 * span))
    
    p_forest <- make_panel_slopes(
      forest_df,
      title = "Forest — VAR2CSA | HRP2 slopes (Primi anchored vs all) — PCR+ + common support",
      colour_hex = ambar,
      ylim_max = ylim_max
    )
    
    print(p_forest)
    
    h_forest <- max(6, min(30, 0.20 * nrow(forest_df) + 2))
    ggsave(paste0("FOREST_VAR2CSA_HRP2_PrimiAnchorVsAll_", file_stub, "_PCRpos_CommonSupport.png"),
           p_forest, width = 15, height = h_forest, dpi = 300)
    ggsave(paste0("FOREST_VAR2CSA_HRP2_PrimiAnchorVsAll_", file_stub, "_PCRpos_CommonSupport.pdf"),
           p_forest, width = 15, height = h_forest)
  }
  
  ## ---- Line plot (adjusted antibody–HRP2) on HRP2 common support ----
  if (is.finite(cs_hrp2$c_low) && is.finite(cs_hrp2$c_high) && cs_hrp2$c_low < cs_hrp2$c_high) {
    c_seq <- seq(cs_hrp2$c_low, cs_hrp2$c_high, length.out = 25)
  } else {
    mf_h <- model.frame(m_hrp2)
    c_vals <- na.omit(mf_h$c_log_pfhrp2)
    c_seq <- as.numeric(unique(stats::quantile(c_vals, probs = seq(0.05, 0.95, by = 0.05), na.rm = TRUE)))
  }
  
  lines_df <- lapply(c_seq, function(cc) {
    em <- emmeans(m_hrp2, specs = ~ gravidity * centro, at = list(c_log_pfhrp2 = cc))
    as.data.frame(em) %>% dplyr::mutate(c_val = cc)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(
      stratum = interaction(gravidity, centro, sep = " × "),
      hrp2_pgml = pmax(exp(c_val + mu_loghrp2_center) - 1, 1e-6),
      pred  = emmean,
      lower = lower.CL,
      upper = upper.CL
    )
  
  p_lines <- ggplot(lines_df, aes(x = hrp2_pgml, y = pred, group = stratum)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = stratum), alpha = 0.10, color = NA) +
    geom_line(aes(color = stratum), linewidth = 0.9) +
    scale_x_log10(labels = scales::label_number(accuracy = 0.01, big.mark = ",")) +
    labs(
      title = paste0(outcome_label, ": adjusted antibody–HRP2 by gravidity × centro"),
      x = "HRP2 (pg/mL; log scale)",
      y = "Predicted antibody level",
      color = NULL
    ) +
    guides(fill = "none") +
    theme_classic(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.4, "lines")
    )
  
  print(p_lines)
  
  ggsave(paste0("LINES_adjusted_", file_stub, "_HRP2_PCRpos_CommonSupport.png"),
         p_lines, width = 11, height = 6.5, dpi = 300)
  ggsave(paste0("LINES_adjusted_", file_stub, "_HRP2_PCRpos_CommonSupport.pdf"),
         p_lines, width = 11, height = 6.5)
  
  ## ---- Hasta/Desde: VAR2CSA only (HRP2 + qPCR) with common-support grids ----
  hasta_out <- NULL
  if (outcome_var == "antibody_pregnancy") {
    
    hrp2_out <- run_hasta_desde_primi_anchor_vs_all(
      model_obj     = m_hrp2,
      predictor_tag = "HRP2",
      c_var_name    = "c_log_pfhrp2",
      mu_log_raw    = mu_loghrp2_center,
      c_low         = cs_hrp2$c_low,
      c_high        = cs_hrp2$c_high,
      unit_suffix   = "pg/mL",
      colour_hex    = ambar,
      file_stub     = paste0(file_stub, "_PCRpos_CommonSupport")
    )
    
    qpcr_out <- run_hasta_desde_primi_anchor_vs_all(
      model_obj     = m_qpcr,
      predictor_tag = "qPCRdensity",
      c_var_name    = "c_log_density",
      mu_log_raw    = mu_logdens_center,
      c_low         = cs_den$c_low,
      c_high        = cs_den$c_high,
      unit_suffix   = "parasites/µL",
      colour_hex    = blue,
      file_stub     = paste0(file_stub, "_PCRpos_CommonSupport")
    )
    
    hasta_out <- list(HRP2 = hrp2_out, qPCRdensity = qpcr_out)
  }
  
  invisible(list(
    models = list(qpcr = m_qpcr, hrp2 = m_hrp2, ldh = m_ldh),
    slopes = slopes_all,
    forest_plot = p_forest,
    lines_plot = p_lines,
    hasta_cuando = hasta_out
  ))
}

## ============================================================
## 5) Run both outcomes
## ============================================================

## exact centering means implied by your c_ variables (PCR+ base)
MU_LOG_DENS <- mean(data_base$log_density - data_base$c_log_density, na.rm = TRUE)
MU_LOG_HRP2 <- mean(data_base$log_pfhrp2  - data_base$c_log_pfhrp2,  na.rm = TRUE)

res_preg <- run_antibody_suite(
  df_den  = data_den,
  df_hrp2 = data_hrp2,
  df_ldh  = data_ldh,
  outcome_var   = "antibody_pregnancy",
  outcome_label = "VAR2CSA IgG",
  file_stub     = "antibody_pregnancy",
  cs_den = cs_den, cs_hrp2 = cs_hrp2, cs_ldh = cs_ldh,
  mu_logdens_center = MU_LOG_DENS,
  mu_loghrp2_center = MU_LOG_HRP2
)

res_cum <- run_antibody_suite(
  df_den  = data_den,
  df_hrp2 = data_hrp2,
  df_ldh  = data_ldh,
  outcome_var   = "antibody_cumulative",
  outcome_label = "Cumulative blood-stage IgG",
  file_stub     = "antibody_cumulative",
  cs_den = cs_den, cs_hrp2 = cs_hrp2, cs_ldh = cs_ldh,
  mu_logdens_center = MU_LOG_DENS,
  mu_loghrp2_center = MU_LOG_HRP2
)

cat("\nDONE SCRIPT 2 (PCR+ + common support per predictor) — with SCRIPT 1 hasta/desde formatting.\n")

######================patwork with all hasta cuando plots (with the condition to have ran sensiticity script)
## ============================================================
## Combine 3 HASTA/DESDE PNGs into one figure (a–c)
## Layout: (a | b) / c
## ============================================================

## ============================================================
## Combine 3 HASTA/DESDE PNGs into one figure
## Layout: (a | b) on top, and c centered below (same width as a/b)
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

## --- filenames (edit if yours differ) ---
f_a <- "HASTA_DESDE_HRP2_qPCRdensity_PrimiAnchorVsAll_IJ_facets_PCRpos_CommonSupport.png"
f_b <- "HASTA_DESDE_PfLDH_qPCRdensity_PrimiAnchorVsAll_IJ_facets_PCRpos_CommonSupport.png"
f_c <- "HASTA_DESDE_VAR2CSA_HRP2_PrimiAnchorVsAll_IJ_facets_antibody_pregnancy_PCRpos_CommonSupport.png"

p_a <- png_to_gg(f_a) + labs(tag = "a")
p_b <- png_to_gg(f_b) + labs(tag = "b")
p_c <- png_to_gg(f_c) + labs(tag = "c")

top_row <- wrap_plots(p_a, p_b, ncol = 2)

## Bottom row: spacer | C | spacer, with widths 1:2:1
## => C takes 2/(1+2+1) = 0.5 of the total width (same as A or B)
bottom_row <- wrap_plots(plot_spacer(), p_c, plot_spacer(),
                         ncol = 3, widths = c(1, 2, 1))

p_3panel <- top_row / bottom_row &
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    plot.tag.position = c(0.02, 0.98)
  )

print(p_3panel)

ggsave("HASTA_DESDE_3PANEL_a_b_c_centered.png", p_3panel, width = 16, height = 14, dpi = 300)
ggsave("HASTA_DESDE_3PANEL_a_b_c_centered.pdf", p_3panel, width = 16, height = 14)

cat("Saved:\n - HASTA_DESDE_3PANEL_a_b_c_centered.png\n - HASTA_DESDE_3PANEL_a_b_c_centered.pdf\n")

## ============================================================
## 6) MAGNITUDE (VAR2CSA) within significant hasta/desde regions
##     - Computes effect size restricted to the DARK regions:
##       sig_pos (CI fully > 0) and sig_neg (CI fully < 0)
##     - Outputs both model-scale differences AND fold-change
##       assuming ln OR log10 (since base is not explicit in script)
##
## Saves (Excel-friendly UTF-8 BOM):
##   MAGNITUDE_VAR2CSA_HRP2predictor_PCRpos_CommonSupport.csv
##   MAGNITUDE_VAR2CSA_qPCRpredictor_PCRpos_CommonSupport.csv
##   MAGNITUDE_VAR2CSA_ALLpredictors_PCRpos_CommonSupport.csv
## ============================================================

clean_text <- function(x) {
  if (is.null(x)) return(x)
  x <- as.character(x)
  x <- enc2utf8(x)
  x <- gsub("ManhiÃ§a", "Manhiça", x, fixed = TRUE)
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

fix_hasta_curves_cols <- function(df) {
  stopifnot(is.data.frame(df))
  nm <- names(df)
  nm <- gsub("^\ufeff", "", nm)   # strip BOM if present
  nm <- gsub("^ï\\.\\.", "", nm)  # strip weird BOM artefact
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
  rename_first("lower",    c("lower","lower_cl","lowercl","lwr","lower_ci"))
  rename_first("upper",    c("upper","upper_cl","uppercl","upr","upper_ci"))
  rename_first("sig_pos",  c("sig_pos","sigpos","pos_sig","sig_positive"))
  rename_first("sig_neg",  c("sig_neg","signeg","neg_sig","sig_negative"))
  
  df
}

summarise_magnitude_from_curves_var2csa <- function(curves_df, predictor_tag) {
  
  df <- fix_hasta_curves_cols(curves_df)
  
  req <- c("contrast","x_raw","diff","lower","upper","sig_pos","sig_neg")
  miss <- setdiff(req, names(df))
  if (length(miss)) {
    stop(
      "Curves missing required columns: ", paste(miss, collapse = ", "),
      "\nColumns seen:\n  ", paste(names(df), collapse = ", ")
    )
  }
  
  df <- df %>%
    dplyr::mutate(
      predictor = predictor_tag,
      contrast  = clean_text(.data$contrast),
      x_raw     = as.numeric(.data$x_raw),
      diff      = as.numeric(.data$diff),
      lower     = as.numeric(.data$lower),
      upper     = as.numeric(.data$upper),
      sig_pos   = as.logical(.data$sig_pos),
      sig_neg   = as.logical(.data$sig_neg)
    ) %>%
    dplyr::filter(is.finite(x_raw), !is.na(x_raw), x_raw > 0) %>%
    dplyr::arrange(contrast, x_raw)
  
  if (!nrow(df)) return(tibble::tibble())
  
  summarise_one_direction <- function(dsub, contrast_key, flag_col, direction_label) {
    
    dsub <- dsub %>%
      dplyr::arrange(x_raw) %>%
      dplyr::mutate(
        flag     = as.logical(.data[[flag_col]]),
        log10x   = log10(x_raw),
        
        ## Fold-change options (choose the correct one later)
        ratio_if_ln    = exp(diff),     # valid if outcome is ln-transformed
        ratio_if_log10 = 10^diff        # valid if outcome is log10-transformed
      )
    
    dsub$flag[is.na(dsub$flag)] <- FALSE
    if (!any(dsub$flag, na.rm = TRUE)) return(tibble::tibble())
    
    dsub <- dsub %>%
      dplyr::mutate(
        direction = direction_label,
        seg_id    = add_segment_id(.data$flag)
      )
    
    overall <- dsub %>%
      dplyr::filter(flag) %>%
      dplyr::summarise(
        predictor = predictor_tag,
        contrast  = contrast_key,
        direction = direction_label,
        level     = "overall",
        seg_id    = NA_integer_,
        n_points  = dplyr::n(),
        x_min     = min(x_raw, na.rm = TRUE),
        x_max     = max(x_raw, na.rm = TRUE),
        log10_range = max(log10x, na.rm = TRUE) - min(log10x, na.rm = TRUE),
        
        ## Model-scale magnitude (always valid)
        mean_diff   = mean(diff, na.rm = TRUE),
        median_diff = stats::median(diff, na.rm = TRUE),
        
        ## Fold-change summaries (interpret depending on log base)
        mean_ratio_if_ln    = mean(ratio_if_ln, na.rm = TRUE),
        median_ratio_if_ln  = stats::median(ratio_if_ln, na.rm = TRUE),
        geom_mean_ratio_if_ln = exp(mean(diff, na.rm = TRUE)),
        
        mean_ratio_if_log10   = mean(ratio_if_log10, na.rm = TRUE),
        median_ratio_if_log10 = stats::median(ratio_if_log10, na.rm = TRUE),
        geom_mean_ratio_if_log10 = 10^(mean(diff, na.rm = TRUE)),
        
        ## Average effect per 1-log10 increase in predictor (on model scale)
        avg_diff_per_log10 = {
          w <- max(log10x, na.rm = TRUE) - min(log10x, na.rm = TRUE)
          if (is.finite(w) && w > 0) trapz(log10x, diff) / w else NA_real_
        }
      )
    
    segs <- dsub %>%
      dplyr::filter(!is.na(seg_id)) %>%
      dplyr::group_by(seg_id) %>%
      dplyr::summarise(
        predictor = predictor_tag,
        contrast  = contrast_key,
        direction = direction_label,
        level     = "segment",
        n_points  = dplyr::n(),
        x_min     = min(x_raw, na.rm = TRUE),
        x_max     = max(x_raw, na.rm = TRUE),
        log10_range = max(log10x, na.rm = TRUE) - min(log10x, na.rm = TRUE),
        
        mean_diff   = mean(diff, na.rm = TRUE),
        median_diff = stats::median(diff, na.rm = TRUE),
        
        mean_ratio_if_ln    = mean(ratio_if_ln, na.rm = TRUE),
        median_ratio_if_ln  = stats::median(ratio_if_ln, na.rm = TRUE),
        geom_mean_ratio_if_ln = exp(mean(diff, na.rm = TRUE)),
        
        mean_ratio_if_log10   = mean(ratio_if_log10, na.rm = TRUE),
        median_ratio_if_log10 = stats::median(ratio_if_log10, na.rm = TRUE),
        geom_mean_ratio_if_log10 = 10^(mean(diff, na.rm = TRUE)),
        
        avg_diff_per_log10 = {
          w <- max(log10x, na.rm = TRUE) - min(log10x, na.rm = TRUE)
          if (is.finite(w) && w > 0) trapz(log10x, diff) / w else NA_real_
        },
        .groups = "drop"
      )
    
    dplyr::bind_rows(overall, segs)
  }
  
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

get_curves_safe <- function(obj, fallback_file) {
  if (!is.null(obj) && is.list(obj) && "curves" %in% names(obj) && is.data.frame(obj$curves)) {
    return(obj$curves)
  }
  if (file.exists(fallback_file)) {
    return(read.csv(fallback_file, stringsAsFactors = FALSE, check.names = FALSE))
  }
  stop("Could not find curves: object missing and file not found: ", fallback_file)
}

## ---- grab VAR2CSA curves from res_preg (preferred), else read CSV ----
curves_var2_hrp2 <- get_curves_safe(
  res_preg$hasta_cuando$HRP2,
  "CURVES_VAR2CSA_HRP2_hasta_desde_PrimiAnchorVsAll_antibody_pregnancy_PCRpos_CommonSupport.csv"
)

curves_var2_qpcr <- get_curves_safe(
  res_preg$hasta_cuando$qPCRdensity,
  "CURVES_VAR2CSA_qPCRdensity_hasta_desde_PrimiAnchorVsAll_antibody_pregnancy_PCRpos_CommonSupport.csv"
)

## ---- compute magnitude summaries ----
mag_var2_hrp2 <- summarise_magnitude_from_curves_var2csa(curves_var2_hrp2, "VAR2CSA | HRP2 predictor")
mag_var2_qpcr <- summarise_magnitude_from_curves_var2csa(curves_var2_qpcr, "VAR2CSA | qPCR density predictor")

write_csv_utf8_bom(mag_var2_hrp2, "MAGNITUDE_VAR2CSA_HRP2predictor_PCRpos_CommonSupport.csv")
write_csv_utf8_bom(mag_var2_qpcr, "MAGNITUDE_VAR2CSA_qPCRpredictor_PCRpos_CommonSupport.csv")

mag_var2_all <- dplyr::bind_rows(mag_var2_hrp2, mag_var2_qpcr)
write_csv_utf8_bom(mag_var2_all, "MAGNITUDE_VAR2CSA_ALLpredictors_PCRpos_CommonSupport.csv")

cat("\nSaved VAR2CSA magnitude summaries:\n",
    "- MAGNITUDE_VAR2CSA_HRP2predictor_PCRpos_CommonSupport.csv\n",
    "- MAGNITUDE_VAR2CSA_qPCRpredictor_PCRpos_CommonSupport.csv\n",
    "- MAGNITUDE_VAR2CSA_ALLpredictors_PCRpos_CommonSupport.csv\n", sep = "")

cat("\nPreview (overall rows):\n")
print(
  mag_var2_all %>%
    dplyr::filter(level == "overall") %>%
    dplyr::arrange(predictor, contrast, direction) %>%
    utils::head(20)
)

