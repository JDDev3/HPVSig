evaluate_threshold <- function(x, group, thr, method_name = NA_character_) {
  if (!is.factor(group) || nlevels(group) != 2L) {
    stop("`group` must be a factor with exactly 2 levels.")
  }
  g_levels <- levels(group)
  pos <- g_levels[2]  # treat second level as "positive"
  
  # Decide which side of the threshold should be "positive"
  mean_pos <- mean(x[group == pos], na.rm = TRUE)
  mean_neg <- mean(x[group != pos], na.rm = TRUE)
  
  if (mean_pos >= mean_neg) {
    pred_pos <- x >= thr
  } else {
    pred_pos <- x <= thr
  }
  
  pred <- ifelse(pred_pos, pos, g_levels[1])
  pred <- factor(pred, levels = g_levels)
  
  tab <- table(Observed = group, Predicted = pred)
  
  TP <- tab[pos, pos]
  TN <- tab[g_levels[1], g_levels[1]]
  FP <- tab[g_levels[1], pos]
  FN <- tab[pos, g_levels[1]]
  
  acc  <- (TP + TN) / sum(tab)
  sens <- TP / (TP + FN)
  spec <- TN / (TN + FP)
  
  list(
    method      = method_name,
    threshold   = thr,
    accuracy    = acc,
    sensitivity = sens,
    specificity = spec,
    confusion   = tab
  )
}

# -------------------------------------------------------------------------
# Method 1: ROC + Youden's index (supervised)
# -------------------------------------------------------------------------
#threshold_youden <- function(x, group) {
#  if (!is.factor(group) || nlevels(group) != 2L) {
#    stop("`group` must be a factor with exactly 2 levels.")
#  }
#  roc_obj <- roc(response = group, predictor = x, direction = "<")
#  best    <- coords(roc_obj, x = "best", best.method = "youden", transpose = FALSE)
#  print(best)
#  thr     <- as.numeric(best["threshold"])
#  
#  eval <- evaluate_threshold(x, group, thr, method_name = "ROC+Youden")
#  
#  list(
#    method     = "ROC+Youden",
#    roc        = roc_obj,
#    coords     = best,
#    evaluation = eval
#  )
#}

threshold_youden <- function(x,
                             group,
                             ci_level = 0.95,
                             roc_test_method = "delong") {
  # x: numeric predictor (e.g., PC1)
  # group: factor with exactly 2 levels: c(neg, pos)
  # Returns: list(method, threshold, roc_plot, auc, auc_ci, p_value, evaluation, ...)
  
  if (is.null(group)) stop("`group` is required for Youden ROC thresholding.")
  if (!is.factor(group)) group <- factor(group)
  
  # Pairwise complete cases
  x <- as.numeric(x)
  ok <- is.finite(x) & !is.na(group)
  x <- x[ok]
  group <- droplevels(group[ok])
  
  if (nlevels(group) != 2L) {
    stop("`group` must have exactly 2 levels after dropping NA/unknown.")
  }
  if (length(x) < 4) stop("Need at least 4 finite values for ROC/Youden.")
  
  levs <- levels(group)   # levs[1]=negative/control, levs[2]=positive/case
  
  # ROC: pROC treats the 2nd level as "case" by default when levels are provided
  roc_obj <- pROC::roc(
    response  = group,
    predictor = x,
    levels    = levs,
    direction = "auto",
    quiet     = TRUE
  )
  
  # AUC + CI (CI may fail for tiny n; guard with tryCatch)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  auc_ci <- tryCatch({
    ci <- pROC::ci.auc(roc_obj, conf.level = ci_level)
    as.numeric(ci)[c(1, 2, 3)]  # lower, median, upper
  }, error = function(e) {
    c(NA_real_, NA_real_, NA_real_)
  })
  
  # DeLong p-value for H0: AUC = 0.5 (guard for failures)
  p_val <- tryCatch({
    pROC::roc.test(roc_obj, auc = 0.5, method = roc_test_method)$p.value
  }, error = function(e) {
    NA_real_
  })
  
  # Youden threshold (best.method="youden")
  best <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  # coords() can return a named vector or 1-row matrix/data.frame depending on pROC version
  thr <- if (is.matrix(best) || is.data.frame(best)) as.numeric(best[1, "threshold"]) else as.numeric(best["threshold"])
  sens <- if (is.matrix(best) || is.data.frame(best)) as.numeric(best[1, "sensitivity"]) else as.numeric(best["sensitivity"])
  spec <- if (is.matrix(best) || is.data.frame(best)) as.numeric(best[1, "specificity"]) else as.numeric(best["specificity"])
  
  # Youden J
  youden_J <- sens + spec - 1
  
  # Build ROC plot data
  roc_coords <- pROC::coords(
    roc_obj,
    x = "all",
    ret = c("specificity", "sensitivity"),
    transpose = FALSE
  )
  roc_df <- as.data.frame(roc_coords)
  roc_df$fpr <- 1 - roc_df$specificity
  roc_df$tpr <- roc_df$sensitivity
  
  auc_ci_txt <- if (all(is.finite(auc_ci))) {
    sprintf("AUC = %.3f (%.1f%% CI %.3f–%.3f)", auc_val, 100 * ci_level, auc_ci[1], auc_ci[3])
  } else {
    sprintf("AUC = %.3f (CI unavailable)", auc_val)
  }
  p_txt <- if (is.finite(p_val)) sprintf("DeLong p = %.3g", p_val) else "DeLong p = NA"
  
  roc_plot <- ggplot2::ggplot(roc_df, ggplot2::aes(x = fpr, y = tpr)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "ROC curve (Youden)",
      subtitle = paste(auc_ci_txt, p_txt, sep = " • "),
      x = "False positive rate (1 - specificity)",
      y = "True positive rate (sensitivity)"
    ) +
    ggplot2::theme_classic()
  
  # Keep your existing evaluation path
  eval <- evaluate_threshold(
    x,
    group,
    thr,
    method_name = "Youden (ROC)"
  )
  
  # Add diagnostics into evaluation too (handy for tables/export if you want later)
  eval$auc <- auc_val
  eval$auc_ci <- auc_ci
  eval$p_value_auc_gt_0.5 <- p_val
  eval$roc_sensitivity <- sens
  eval$roc_specificity <- spec
  eval$youden_J <- youden_J
  
  list(
    method     = "Youden (ROC)",
    threshold  = thr,
    roc_obj    = roc_obj,
    roc_plot   = roc_plot,
    auc        = auc_val,
    auc_ci     = auc_ci,    # c(lower, median, upper)
    p_value    = p_val,
    sensitivity = sens,
    specificity = spec,
    youden_J    = youden_J,
    evaluation  = eval
  )
}



# -------------------------------------------------------------------------
# Method 2: Gaussian (auto) – decide equal vs unequal variance via var.test
# -------------------------------------------------------------------------
threshold_gaussian_auto <- function(x,
                                    group,
                                    alpha = 0.05,
                                    unequal_boundary = c("density", "bayes")) {
  unequal_boundary <- match.arg(unequal_boundary)
  
  # Ensure consistent, reproducible filtering (Shiny vs manual)
  x <- as.numeric(x)
  group <- factor(group)
  
  ok <- is.finite(x) & !is.na(group)
  x <- x[ok]
  group <- droplevels(group[ok])
  
  if (nlevels(group) != 2L) {
    stop("`group` must be a factor with exactly 2 levels after dropping NA/invalid.")
  }
  
  g_levels <- levels(group)
  g1 <- g_levels[1]
  g2 <- g_levels[2]
  
  x1 <- x[group == g1]
  x2 <- x[group == g2]
  
  mu1 <- mean(x1)
  mu2 <- mean(x2)
  s1  <- stats::sd(x1)
  s2  <- stats::sd(x2)
  
  n1 <- length(x1)
  n2 <- length(x2)
  n_tot <- n1 + n2
  
  if (n1 < 1 || n2 < 1) stop("Both groups must have at least one observation.")
  
  # var.test can error if a group has <2 obs or zero variance; guard it
  vt <- tryCatch(stats::var.test(x1, x2), error = function(e) NULL)
  p_var <- if (!is.null(vt)) vt$p.value else NA_real_
  
  use_equal <- is.finite(p_var) && p_var >= alpha
  chosen <- if (!is.finite(p_var)) "not_tested" else if (use_equal) "equal" else "unequal"
  
  mid <- (mu1 + mu2) / 2
  
  # ------------------------------------------------------------------
  # Rule requested:
  #   - equal variance -> midpoint
  #   - unequal variance -> density intersection by default
  #       (option: Bayes boundary with empirical priors)
  # ------------------------------------------------------------------
  thr <- mid
  boundary_used <- "midpoint"
  
  if (!use_equal) {
    # if SDs are unusable, fall back to midpoint
    if (is.finite(s1) && s1 > 0 && is.finite(s2) && s2 > 0) {
      
      # weights for boundary equation w1 f1(x) = w2 f2(x)
      if (unequal_boundary == "bayes") {
        w1 <- n1 / n_tot
        w2 <- n2 / n_tot
      } else {
        w1 <- 0.5
        w2 <- 0.5
      }
      boundary_used <- unequal_boundary
      
      eps <- sqrt(.Machine$double.eps)
      
      A <- 1 / s1^2 - 1 / s2^2
      B <- -2 * mu1 / s1^2 + 2 * mu2 / s2^2
      C <- (mu1^2 / s1^2) - (mu2^2 / s2^2) -
        2 * log((w2 * s1) / (w1 * s2))
      
      if (abs(A) < eps) {
        # equal-variance-ish linear solution
        sigma2 <- (s1^2 + s2^2) / 2
        if (mu1 == mu2) {
          thr <- mid
        } else {
          thr <- mid + (sigma2 / (mu2 - mu1)) * log(w1 / w2)
        }
      } else {
        disc <- B^2 - 4 * A * C
        if (is.finite(disc) && disc >= 0) {
          sqrt_disc <- sqrt(disc)
          r1 <- (-B + sqrt_disc) / (2 * A)
          r2 <- (-B - sqrt_disc) / (2 * A)
          candidates <- c(r1, r2)
          
          # Prefer root between means; otherwise closest to midpoint
          between <- candidates >= min(mu1, mu2) & candidates <= max(mu1, mu2)
          if (any(between)) {
            thr <- candidates[between][which.min(abs(candidates[between] - mid))]
          } else {
            thr <- candidates[which.min(abs(candidates - mid))]
          }
        } else {
          thr <- mid
          boundary_used <- "midpoint_fallback"
        }
      }
    } else {
      boundary_used <- "midpoint_fallback"
      thr <- mid
    }
  }
  
  method_label <- if (use_equal) {
    "Gaussian(midpoint; equal var)"
  } else {
    paste0("Gaussian(", boundary_used, "; unequal var)")
  }
  
  method_name <- if (use_equal) {
    "Gaussian (equal var)"
  } else {
    paste0("Gaussian (unequal var)")
  }
  
  eval <- evaluate_threshold(x, group, thr, method_name = method_name)
  
  sd_pooled <- if (n1 >= 2 && n2 >= 2 && is.finite(s1) && is.finite(s2)) {
    sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  } else {
    NA_real_
  }
  
  list(
    method     = method_label,
    chosen     = chosen,
    p_var_test = p_var,
    threshold  = thr,
    params     = list(
      mu1 = mu1, mu2 = mu2,
      sd1 = s1,  sd2 = s2,
      sd_pooled = sd_pooled,
      n1 = n1, n2 = n2,
      pi1 = n1 / n_tot, pi2 = n2 / n_tot,
      unequal_boundary = unequal_boundary
    ),
    evaluation = eval
  )
}




# -------------------------------------------------------------------------
# Method 3: Gaussian Mixture Model (unsupervised fit, optional supervised eval)
# -------------------------------------------------------------------------
threshold_gmm <- function(x,
                          group = NULL,
                          method = c("equal_component", "posterior_0.5", "valley", "discrete"),
                          fallback = c("discrete", "posterior_0.5", "equal_component", "valley", "none"),
                          grid_n = 512) {
  method   <- match.arg(method)
  fallback <- match.arg(fallback)
  
  x <- as.numeric(x)
  
  if (!is.null(group)) {
    group <- factor(group)
    ok <- is.finite(x) & !is.na(group)
    x <- x[ok]
    group <- droplevels(group[ok])
    if (nlevels(group) != 2L) {
      stop("If provided, `group` must be a factor with exactly 2 levels after dropping NA/invalid.")
    }
  } else {
    x <- x[is.finite(x)]
  }
  
  if (length(x) < 4) stop("`x` must have at least 4 finite values.")
  
  fit <- mclust::Mclust(x, G = 2, verbose = FALSE)
  means <- fit$parameters$mean
  vars <- as.numeric(fit$parameters$variance$sigmasq)
  if (length(vars) == 1L) vars <- rep(vars, 2L)
  if (length(vars) != 2L) stop("Unexpected variance structure from Mclust.")
  props <- fit$parameters$pro
  
  if (length(means) != 2L) {
    stop("GMM did not fit exactly 2 components; got ", length(means), ".")
  }
  
  # Order components by mean (lo, hi) for consistent bracketing
  ordc <- order(means)
  m_lo <- means[ordc[1]]
  m_hi <- means[ordc[2]]
  s_lo <- sqrt(vars[ordc[1]])
  s_hi <- sqrt(vars[ordc[2]])
  pi_lo <- props[ordc[1]]
  pi_hi <- props[ordc[2]]
  
  # Map back to "high mean component" index for posterior-based discrete fallback
  pos_comp <- which.max(means)
  
  # Helper: find a root of g(x)=0 near the midpoint between means.
  # If endpoint sign change fails, scan a grid to find a bracketing interval.
  solve_root <- function(g, lower, upper, target = (m_lo + m_hi) / 2) {
    if (!is.finite(lower) || !is.finite(upper) || lower >= upper) return(NA_real_)
    
    gl <- g(lower)
    gu <- g(upper)
    
    if (is.finite(gl) && is.finite(gu) && gl == 0) return(lower)
    if (is.finite(gl) && is.finite(gu) && gu == 0) return(upper)
    
    if (is.finite(gl) && is.finite(gu) && (gl * gu < 0)) {
      return(stats::uniroot(g, lower = lower, upper = upper)$root)
    }
    
    # Grid scan for a sign change bracket
    xs <- seq(lower, upper, length.out = grid_n)
    gs <- vapply(xs, g, numeric(1))
    ok <- is.finite(gs)
    xs <- xs[ok]; gs <- gs[ok]
    if (length(xs) < 4) return(NA_real_)
    
    sgn <- sign(gs)
    flip <- which(sgn[-1] != sgn[-length(sgn)] & sgn[-1] != 0 & sgn[-length(sgn)] != 0)
    if (!length(flip)) return(NA_real_)
    
    # Choose bracket whose midpoint is closest to target (typically mean midpoint)
    mids <- (xs[flip] + xs[flip + 1]) / 2
    i <- flip[which.min(abs(mids - target))]
    stats::uniroot(g, lower = xs[i], upper = xs[i + 1])$root
  }
  
  # Continuous methods ---------------------------------------------------------
  thr_continuous <- function(which_method) {
    # bracket region: between means, with a modest extension
    mid  <- (m_lo + m_hi) / 2
    span <- max(abs(m_hi - m_lo), 1e-6)
    ext  <- span * 0.5 + 3 * max(s_lo, s_hi, na.rm = TRUE)

    
    lower <- m_lo - ext

    upper <- m_hi + ext

    
    if (which_method == "equal_component") {
      g <- function(z) stats::dnorm(z, mean = m_lo, sd = s_lo) - stats::dnorm(z, mean = m_hi, sd = s_hi)
      # Prefer the root between means if possible
      thr <- solve_root(g, lower = m_lo, upper = m_hi, target = mid)
      if (!is.finite(thr)) thr <- solve_root(g, lower = lower, upper = upper, target = mid)
      return(thr)
    }
    
    if (which_method == "posterior_0.5") {
      # MAP boundary: pi_hi*f_hi - pi_lo*f_lo = 0
      g <- function(z) (pi_hi * stats::dnorm(z, mean = m_hi, sd = s_hi)) -
        (pi_lo * stats::dnorm(z, mean = m_lo, sd = s_lo))
      thr <- solve_root(g, lower = m_lo, upper = m_hi, target = mid)
      if (!is.finite(thr)) thr <- solve_root(g, lower = lower, upper = upper, target = mid)
      return(thr)
    }
    
    if (which_method == "valley") {
      print("working")
      # Minimize mixture density between the means (even if not strongly bimodal)
      p <- function(z) (pi_lo * stats::dnorm(z, mean = m_lo, sd = s_lo)) +
        (pi_hi * stats::dnorm(z, mean = m_hi, sd = s_hi))
      if (!is.finite(m_lo) || !is.finite(m_hi) || m_lo >= m_hi) return(NA_real_)
      out <- stats::optimize(p, interval = c(m_lo, m_hi))
      return(out$minimum)
    }
    
    NA_real_
  }
  
  # Discrete fallback ----------------------------------------------------------
  thr_discrete <- function() {
    post_pos <- fit$z[, pos_comp]
    ord <- order(x)
    x_ord <- x[ord]
    post_ord <- post_pos[ord]
    
    is_hi <- post_ord >= 0.5
    flip_idx <- which(is_hi[-1] != is_hi[-length(is_hi)])
    
    if (!length(flip_idx)) {
      # No crossing: use midpoint of component means (stable, interpretable)
      return((m_lo + m_hi) / 2)
    }
    
    cand <- vapply(flip_idx, function(i) mean(c(x_ord[i], x_ord[i + 1])), numeric(1))
    mid <- (m_lo + m_hi) / 2
    
    cand_between <- cand[cand >= m_lo & cand <= m_hi]
    if (length(cand_between)) {
      cand_between[which.min(abs(cand_between - mid))]
    } else {
      cand[which.min(abs(cand - mid))]
    }
  }
  
  # Choose threshold with fallback logic --------------------------------------
  thr <- switch(
    method,
    equal_component = thr_continuous("equal_component"),
    posterior_0.5   = thr_continuous("posterior_0.5"),
    valley          = thr_continuous("valley"),
    discrete        = thr_discrete()
  )
  
  used <- method
  if (!is.finite(thr) && fallback != "none") {
    thr <- switch(
      fallback,
      discrete        = thr_discrete(),
      posterior_0.5   = thr_continuous("posterior_0.5"),
      equal_component = thr_continuous("equal_component"),
      valley          = thr_continuous("valley"),
      none            = NA_real_
    )
    used <- paste0(method, " (fallback→", fallback, ")")
  }
  
  eval <- NULL
  if (!is.null(group) && is.finite(thr)) {
    if (!is.factor(group) || nlevels(group) != 2L) {
      stop("If provided, `group` must be a factor with exactly 2 levels.")
    }
    eval <- evaluate_threshold(x, group, thr, method_name = "GMM")
  }
  
  list(
    method     = paste0("GMM: ", used),
    mclust_fit = fit,
    threshold  = thr,
    params     = list(
      means    = means,
      vars     = vars,
      props    = props,
      pos_comp = which.max(means),
      neg_comp = setdiff(1:2, which.max(means)),
      ordered  = list(m_lo = m_lo, m_hi = m_hi, s_lo = s_lo, s_hi = s_hi, pi_lo = pi_lo, pi_hi = pi_hi)
    ),
    evaluation = eval
  )
}




# -------------------------------------------------------------------------
# Build raw summary table from a list of method result objects
# -------------------------------------------------------------------------
threshold_summary_table <- function(res_list) {
  rows <- lapply(res_list, function(res) {
    ev <- res$evaluation
    if (is.null(ev)) return(NULL)
    
    tab <- ev$confusion
    if (is.null(tab)) return(NULL)
    
    g_levels <- rownames(tab)
    neg <- g_levels[1]
    pos <- g_levels[2]
    
    TP <- tab[pos, pos]
    TN <- tab[neg, neg]
    FP <- tab[neg, pos]
    FN <- tab[pos, neg]
    
    tibble(
      method      = if (!is.null(ev$method) && !is.na(ev$method)) ev$method else res$method,
      threshold   = ev$threshold,
      accuracy    = ev$accuracy,
      sensitivity = ev$sensitivity,
      specificity = ev$specificity,
      TP = TP, TN = TN, FP = FP, FN = FN
    )
  })
  
  bind_rows(rows)
}

# -------------------------------------------------------------------------
# Format summary into a publication-ready table (variable-agnostic)
# -------------------------------------------------------------------------
format_threshold_summary <- function(summary_df,
                                     var_label = "Value",
                                     digits = 3,
                                     include_confusion = FALSE) {
  if (nrow(summary_df) == 0) return(summary_df)
  
  df <- summary_df %>%
    dplyr::mutate(
      threshold   = round(.data$threshold,   digits = digits),
      accuracy    = round(.data$accuracy,    digits = digits),
      sensitivity = round(.data$sensitivity, digits = digits),
      specificity = round(.data$specificity, digits = digits)
    )
  
  if (!include_confusion) {
    df <- df %>%
      dplyr::select(method, threshold, accuracy, sensitivity, specificity)
  }
  
  thr_col_name <- paste0("Threshold (", var_label, ")")
  
  df %>%
    dplyr::rename(
      Method          = method,
      !!thr_col_name := threshold,
      Accuracy        = accuracy,
      Sensitivity     = sensitivity,
      Specificity     = specificity
    )
}


# -------------------------------------------------------------------------
# UPDATED: Plot distributions by group with threshold lines (generic variable)
#          - Manual colors for groups via group_colors
#          - Palette-based colors for methods via method_palette (Brewer)
# -------------------------------------------------------------------------
plot_thresholds <- function(x,
                            group,
                            summary_table,
                            var_label      = "Value",
                            group_colors   = NULL,
                            method_palette = NULL) {
  dat <- tibble(value = x, group = group)
  
  # Identify the threshold column in the summary table
  thr_col <- grep("^Threshold \\(", names(summary_table), value = TRUE)
  if (length(thr_col) != 1L) {
    stop("Could not uniquely identify the threshold column in summary_table.")
  }
  
  p <- ggplot(dat, aes(x = value, fill = group)) +
    geom_density(alpha = 0.5, color = NA) +
    geom_vline(
      data = summary_table,
      aes(xintercept = .data[[thr_col]], color = Method),
      linetype = "dashed",
      linewidth = 2
    ) +
    theme_classic() +
    labs(
      x = var_label,
      y = "Density",
      fill = "Group",
      color = "Method"
    )+
    theme(legend.position = "right")+
    theme(legend.text = element_text(size = 24))+
    theme(legend.title = element_text(size = 24))+
    theme(axis.text.x = element_text(size = 36))+
    theme(axis.text.y = element_text(size = 36))+
    theme(axis.title = element_text(size = 36))
  
  # Optional: manual colors for groups (fill)
  if (!is.null(group_colors)) {
    lev_g <- levels(group)
    cols_g <- group_colors
    
    if (is.null(names(cols_g))) {
      # Unnamed vector: assume order corresponds to factor levels
      if (length(cols_g) < length(lev_g)) {
        stop("`group_colors` has fewer colors than there are group levels.")
      }
      cols_g <- cols_g[seq_along(lev_g)]
      names(cols_g) <- lev_g
    } else {
      # Named vector: align by level name
      if (!all(lev_g %in% names(cols_g))) {
        stop("Names of `group_colors` must include all group levels: ",
             paste(lev_g, collapse = ", "))
      }
      cols_g <- cols_g[lev_g]
    }
    
    p <- p + scale_fill_manual(values = cols_g)
  }
  
  # Optional: palette-based colors for methods (lines)
  # method_palette should be a Brewer palette name, e.g. "Set1", "Dark2", "Paired"
  if (!is.null(method_palette)) {
    p <- p + scale_color_brewer(palette = method_palette)
  }
  
  p
}

# -------------------------------------------------------------------------
# NEW: Build a DT datatable from the publication-ready summary table
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# NEW: Build a gt table from the publication-ready summary table
#      Suitable for PDF/Word/HTML via gt::gtsave()
# -------------------------------------------------------------------------
threshold_gt_table <- function(summary_table,
                               digits  = 3,
                               title   = NULL,
                               subtitle = NULL) {
  if (!requireNamespace("gt", quietly = TRUE)) {
    warning("Package 'gt' is not installed; returning NULL for summary_gt.")
    return(NULL)
  }
  
  # Identify the threshold column name (e.g. "Threshold (PC1 score)")
  thr_col <- grep("^Threshold \\(", names(summary_table), value = TRUE)
  if (length(thr_col) != 1L) {
    thr_col <- NULL
  }
  
  # Columns we want to format as numbers
  numeric_cols <- intersect(
    names(summary_table),
    c(thr_col, "Accuracy", "Sensitivity", "Specificity")
  )
  
  tab <- gt::gt(summary_table)
  
  # Optional title/subtitle
  if (!is.null(title) || !is.null(subtitle)) {
    tab <- gt::tab_header(
      tab,
      title    = if (!is.null(title)) title else gt::md(""),
      subtitle = if (!is.null(subtitle)) subtitle else gt::md("")
    )
  }
  
  # Number formatting for threshold + performance metrics
  if (length(numeric_cols) > 0) {
    tab <- gt::fmt_number(
      tab,
      columns  = dplyr::all_of(numeric_cols),
      decimals = digits
    )
  }
  
  # Some light styling (tweak to taste)
  tab <- tab |>
    gt::tab_options(
      table.font.size = gt::px(20),
      data_row.padding = gt::px(4)
    ) |>
    gt::cols_align(
      align   = "center",
      columns = dplyr::everything()
    )
  
  tab
}



determine_threshold_two_groups <- function(x,
                                           group,
                                           methods           = c("youden", "gaussian_auto", "gmm"),
                                           var_label         = NULL,
                                           digits            = 3,
                                           include_confusion = FALSE,
                                           group_colors      = NULL,
                                           method_palette    = NULL) {
  if (is.null(var_label)) {
    var_label <- deparse(substitute(x))
  }
  
  group <- factor(group)
  if (nlevels(group) != 2L) {
    stop("`group` must have exactly 2 levels.")
  }
  
  valid_methods <- c("youden", "gaussian_auto", "gmm")
  if (!all(methods %in% valid_methods)) {
    stop("`methods` must be a subset of: ",
         paste(valid_methods, collapse = ", "))
  }
  
  results <- list()
  
  if ("youden" %in% methods) {
    results$youden <- threshold_youden(x, group)
  }
  if ("gaussian_auto" %in% methods) {
    results$gaussian_auto <- threshold_gaussian_auto(x, group)
  }
  if ("gmm" %in% methods) {
    results$gmm <- threshold_gmm(x, group)
  }
  
  summary_raw   <- threshold_summary_table(results)
  summary_table <- format_threshold_summary(
    summary_raw,
    var_label         = var_label,
    digits            = digits,
    include_confusion = include_confusion
  )
  
  plt <- plot_thresholds(
    x             = x,
    group         = group,
    summary_table = summary_table,
    var_label     = var_label,
    group_colors  = group_colors,
    method_palette = method_palette
  )
  
  # New: gt table for static output (PDF/Word/HTML via gtsave)
  summary_gt <- threshold_gt_table(
    summary_table,
    digits  = digits,
    subtitle = NULL
  )
  
  # Keep kable option for Rmd PDF/Word if you want it
  summary_kable <- NULL
  if (requireNamespace("knitr", quietly = TRUE)) {
    summary_kable <- knitr::kable(
      summary_table,
      digits = digits,
      align = paste0(c("l", rep("c", ncol(summary_table) - 1)), collapse = "")
    )
  }
  
  list(
    results       = results,
    summary_raw   = summary_raw,
    summary_table = summary_table,
    summary_kable = summary_kable,
    summary_gt    = summary_gt,  # <--- gt object to save with gt::gtsave
    plot          = plt
  )
}

#align new pca to old pca
align_pcatools_pc_sign <- function(ref_pca,
                                   new_pca,
                                   pc = "PC1",
                                   min_common_genes = 50) {
  if (is.null(ref_pca$loadings) || is.null(new_pca$loadings)) {
    stop("Both ref_pca and new_pca must have $loadings.")
  }
  if (!pc %in% colnames(ref_pca$loadings) || !pc %in% colnames(new_pca$loadings)) {
    stop("PC '", pc, "' not found in both $loadings.")
  }
  if (is.null(new_pca$rotated) || !pc %in% colnames(new_pca$rotated)) {
    stop("new_pca must have $rotated with column ", pc, " to flip.")
  }
  
  common <- intersect(rownames(ref_pca$loadings), rownames(new_pca$loadings))
  if (length(common) < 2L) stop("No overlapping genes between reference and new PCA loadings.")
  if (length(common) < min_common_genes) warning("Only ", length(common), " common genes; alignment may be unstable.")
  
  r <- suppressWarnings(stats::cor(
    ref_pca$loadings[common, pc],
    new_pca$loadings[common, pc],
    use = "pairwise.complete.obs"
  ))
  
  flipped <- FALSE
  if (is.finite(r) && !is.na(r) && r < 0) {
    new_pca$loadings[, pc] <- -new_pca$loadings[, pc]
    new_pca$rotated[, pc]  <- -new_pca$rotated[, pc]
    flipped <- TRUE
    r <- -r
  }
  if(flipped){
    message("New Signature Flipped")
  }
  
  list(
    new_pca_aligned = new_pca,
    cor_loadings    = r,
    flipped         = flipped,
    n_common_genes  = length(common),
    pc              = pc
  )
}

plot_pcatools_loadings_correlation <- function(ref_pca,
                                               new_pca,
                                               pc = "PC1",
                                               min_common_genes = 50,
                                               point_alpha = 0.5) {
  if (is.null(ref_pca$loadings) || is.null(new_pca$loadings)) {
    stop("Both ref_pca and new_pca must have $loadings.")
  }
  if (!pc %in% colnames(ref_pca$loadings) || !pc %in% colnames(new_pca$loadings)) {
    stop("PC '", pc, "' not found in both $loadings.")
  }
  
  common <- intersect(rownames(ref_pca$loadings), rownames(new_pca$loadings))
  if (length(common) < 2L) stop("No overlapping genes between reference and new PCA loadings.")
  if (length(common) < min_common_genes) warning("Only ", length(common), " common genes; plot may be unstable.")
  
  df <- tibble::tibble(
    ref = as.numeric(ref_pca$loadings[common, pc]),
    new = as.numeric(new_pca$loadings[common, pc])
  )
  
  r <- suppressWarnings(stats::cor(df$ref, df$new, use = "pairwise.complete.obs"))
  
  #ggpubr::ggscatter(
  #  df,
  #  x = "ref",
  #  y = "new"
  #) +
  #  ggpubr::stat_cor()
  
  ggplot2::ggplot(df, ggplot2::aes(x = ref, y = new)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey85") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
    ggplot2::geom_point(alpha = point_alpha) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    ggplot2::coord_equal() +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title    = paste0("Loadings correlation (", pc, ")"),
      subtitle = paste0("Pearson r = ", signif(r, 3), " | n = ", nrow(df), " common genes"),
      x = "Reference loadings",
      y = "New loadings"
    )
}
