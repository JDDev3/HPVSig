normalize_label <- function(x) tolower(trimws(as.character(x)))

is_unknown_label <- function(x) {
  is_na <- is.na(x)
  x_chr <- normalize_label(x)
  unknown_patterns <- c(
    "", "na", "n/a", "missing", "unknown", "unk",
    "not done", "notdone", "nd", "n.d.", "."
  )
  is_na | x_chr %in% unknown_patterns
}

infer_binary_group <- function(x,
                               flip = FALSE,
                               drop_unknown = TRUE,
                               allow_unknown_in_full = TRUE) {
  x0 <- x
  
  # Keep both raw and trimmed character versions
  x_chr_raw  <- as.character(x0)
  x_chr_trim <- trimws(x_chr_raw)
  
  # Single source of truth for "unknown"
  unknown_mask <- is_unknown_label(x0)
  
  # Define "known" for training/inference
  known <- rep(TRUE, length(x0))
  if (isTRUE(drop_unknown)) {
    known <- !unknown_mask
  }
  
  # Training labels used to infer the two levels
  x_train_chr <- x_chr_trim[known]
  x_train_chr <- x_train_chr[!is.na(x_train_chr)]
  x_train_chr <- x_train_chr[nzchar(x_train_chr)]
  
  levs_obs <- unique(x_train_chr)
  
  if (length(levs_obs) != 2L) {
    stop(sprintf(
      "infer_binary_group() expects exactly 2 distinct non-unknown levels; found %d: %s",
      length(levs_obs),
      paste(levs_obs, collapse = ", ")
    ))
  }
  
  reason <- NULL
  neg <- pos <- NULL
  
  # logical-like
  if (is.logical(x0) || all(levs_obs %in% c("TRUE", "FALSE"))) {
    neg <- "FALSE"
    pos <- "TRUE"
    reason <- "Logical detected: FALSE=negative, TRUE=positive."
    
  } else {
    # numeric-like (handles "0"/"1" too)
    suppressWarnings(num_vals <- as.numeric(levs_obs))
    if (all(is.finite(num_vals))) {
      ord <- order(num_vals)
      neg <- levs_obs[ord[1]]
      pos <- levs_obs[ord[2]]
      reason <- "Numeric-like detected: smaller value=negative, larger value=positive."
      
    } else {
      lev1 <- normalize_label(levs_obs[1])
      lev2 <- normalize_label(levs_obs[2])
      
      neg_keywords <- c("neg","negative","hpv-","hpv -","no","none","absent",
                        "wt","wildtype","control","nonhpv","non-hpv","hpv negative","0")
      pos_keywords <- c("pos","positive","hpv+","hpv +","yes","present",
                        "mut","mutant","case","hpv positive","1")
      
      has_any <- function(lbl, patterns) {
        any(vapply(patterns, function(p) grepl(p, lbl, fixed = TRUE), logical(1)))
      }
      
      lev1_is_neg <- has_any(lev1, neg_keywords)
      lev2_is_neg <- has_any(lev2, neg_keywords)
      lev1_is_pos <- has_any(lev1, pos_keywords)
      lev2_is_pos <- has_any(lev2, pos_keywords)
      
      if (lev1_is_neg && lev2_is_pos) {
        neg <- levs_obs[1]; pos <- levs_obs[2]
        reason <- "Label keywords matched: first=negative-like, second=positive-like."
      } else if (lev2_is_neg && lev1_is_pos) {
        neg <- levs_obs[2]; pos <- levs_obs[1]
        reason <- "Label keywords matched: second=negative-like, first=positive-like."
      } else if (lev1_is_neg && !lev2_is_neg) {
        neg <- levs_obs[1]; pos <- levs_obs[2]
        reason <- "First level looks negative-like; treating other as positive."
      } else if (lev2_is_neg && !lev1_is_neg) {
        neg <- levs_obs[2]; pos <- levs_obs[1]
        reason <- "Second level looks negative-like; treating other as positive."
      } else if (lev1_is_pos && !lev2_is_pos) {
        pos <- levs_obs[1]; neg <- levs_obs[2]
        reason <- "First level looks positive-like; treating other as negative."
      } else if (lev2_is_pos && !lev1_is_pos) {
        pos <- levs_obs[2]; neg <- levs_obs[1]
        reason <- "Second level looks positive-like; treating other as negative."
      } else {
        ord <- sort(levs_obs)
        neg <- ord[1]; pos <- ord[2]
        reason <- "Ambiguous labels; using alphabetical order (first=negative, second=positive)."
      }
    }
  }
  
  if (isTRUE(flip)) {
    tmp <- neg; neg <- pos; pos <- tmp
    reason <- paste0(reason, " (User flip applied.)")
  }
  
  # Use trimmed values for factor creation (avoids whitespace mismatch)
  train_vals <- x_chr_trim[known]
  train_group <- factor(train_vals, levels = c(as.character(neg), as.character(pos)))
  train_group <- droplevels(train_group)
  
  if (any(is.na(train_group))) {
    bad <- unique(train_vals[is.na(train_group)])
    stop(sprintf(
      "infer_binary_group(): training values did not match inferred levels. Unmatched: %s",
      paste(bad, collapse = ", ")
    ))
  }
  
  # Full factor: unknowns may become NA
  full_group <- factor(x_chr_trim, levels = c(as.character(neg), as.character(pos)))
  
  if (!isTRUE(allow_unknown_in_full) && any(is.na(full_group) & !unknown_mask)) {
    bad <- unique(x_chr_trim[is.na(full_group) & !unknown_mask])
    stop(sprintf(
      "infer_binary_group(): some non-unknown values did not match inferred levels. Unmatched: %s",
      paste(bad, collapse = ", ")
    ))
  }
  
  list(
    train_group     = train_group,
    full_group      = full_group,
    neg             = levels(train_group)[1],
    pos             = levels(train_group)[2],
    reason          = reason,
    observed_levels = levs_obs,
    known_idx       = known
  )
}

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
                             ci_level = 0.95) {
  # x: numeric predictor (e.g., PC1)
  # group: factor with exactly 2 levels: c(neg, pos)

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
    direction = if (
      mean(x[group == levs[2]]) >= mean(x[group == levs[1]])
    ) "<" else ">",
    quiet     = TRUE
  )
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  auc_ci <- tryCatch({
    ci <- pROC::ci.auc(roc_obj, conf.level = ci_level)
    as.numeric(ci)[c(1, 2, 3)]  # lower, median, upper
  }, error = function(e) {
    c(NA_real_, NA_real_, NA_real_)
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

  roc_plot <- ggplot2::ggplot(roc_df, ggplot2::aes(x = fpr, y = tpr)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "ROC curve (Youden)",
      subtitle = auc_ci_txt,
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
  
  # Add diagnostics into evaluation 
  eval$auc <- auc_val
  eval$auc_ci <- auc_ci
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
      C <- (mu1^2 / s1^2) - (mu2^2 / s2^2) +
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
                          fallback = c("none", "discrete", "posterior_0.5", "equal_component", "valley"),
                          grid_n = 512) {
  method   <- match.arg(method)
  fallback <- match.arg(fallback)
  
  x <- as.numeric(x)
  
  if (!is.null(group) && length(group) != length(x)) {
    stop("`group` and `x` must have the same length.")
  }
  
  # Keep every finite score for GMM fitting.
  ok <- is.finite(x)
  x <- x[ok]
  
  # Keep reference labels aligned with those scores for later evaluation.
  if (!is.null(group)) {
    group <- group[ok]
    group[is_unknown_label(group)] <- NA
    group <- droplevels(factor(group))
    
    if (nlevels(group) > 2L) {
      stop("`group` must contain no more than 2 known levels.")
    }
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
  if (!is.null(group) && nlevels(group) == 2L && is.finite(thr)) {
    known <- !is.na(group)
    
    eval <- evaluate_threshold(
      x = x[known],
      group = group[known],
      thr = thr,
      method_name = "GMM"
    )
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
  
  # Exact Clopper-Pearson 95% confidence interval
  format_ci <- function(successes, total) {
    if (is.na(successes) || is.na(total) || total <= 0) {
      return(NA_character_)
    }
    
    ci <- stats::binom.test(
      x = successes,
      n = total,
      conf.level = 0.95
    )$conf.int
    
    sprintf(
      paste0("%.", digits, "f–%.", digits, "f"),
      ci[1], ci[2]
    )
  }
  
  df <- summary_df
  
  # Unsupervised mode has no reference labels or confusion counts
  for (column in c("TP", "TN", "FP", "FN")) {
    if (!column %in% names(df)) {
      df[[column]] <- NA_real_
    }
  }
  
  df[["Accuracy 95% CI"]] <- vapply(
    seq_len(nrow(df)),
    function(i) format_ci(
      df$TP[i] + df$TN[i],
      df$TP[i] + df$TN[i] + df$FP[i] + df$FN[i]
    ),
    character(1)
  )
  
  df[["Sensitivity 95% CI"]] <- vapply(
    seq_len(nrow(df)),
    function(i) format_ci(df$TP[i], df$TP[i] + df$FN[i]),
    character(1)
  )
  
  df[["Specificity 95% CI"]] <- vapply(
    seq_len(nrow(df)),
    function(i) format_ci(df$TN[i], df$TN[i] + df$FP[i]),
    character(1)
  )
  
  for (column in c("threshold", "accuracy", "sensitivity", "specificity")) {
    df[[column]] <- round(df[[column]], digits = digits)
  }
  
  # Combine each estimate with its confidence interval
  for (metric in c("accuracy", "sensitivity", "specificity")) {
    label <- switch(
      metric,
      accuracy = "Accuracy",
      sensitivity = "Sensitivity",
      specificity = "Specificity"
    )
    
    ci <- df[[paste0(label, " 95% CI")]]
    
    df[[metric]] <- ifelse(
      is.na(df[[metric]]),
      NA_character_,
      ifelse(
        is.na(ci),
        sprintf("%.*f", digits, df[[metric]]),
        sprintf("%.*f (%s)", digits, df[[metric]], ci)
      )
    )
  }
  
  columns <- c(
    "method", "threshold", "accuracy", "sensitivity", "specificity"
  )
  
  if (include_confusion) {
    columns <- c(columns, "TP", "TN", "FP", "FN")
  }
  
  df <- dplyr::select(df, dplyr::all_of(columns))
  
  names(df)[match(
    c("method", "threshold", "accuracy", "sensitivity", "specificity"),
    names(df)
  )] <- c(
    "Method",
    paste0("Threshold (", var_label, ")"),
    "Accuracy (95% CI)",
    "Sensitivity (95% CI)",
    "Specificity (95% CI)"
  )
  
  df
}

plot_thresholds <- function(x,
                            group,
                            summary_table,
                            var_label = "Value",
                            method_palette = "Dark2") {
  
  group <- factor(group)
  
  # Density colors: negative = blue, positive = red, other = gray.
  group_colors <- setNames(
    rep("#A6A4A4", nlevels(group)),
    levels(group)
  )
  
  group_labels <- tolower(trimws(levels(group)))
  group_colors[group_labels == "negative"] <- "#3B4992"
  group_colors[group_labels == "positive"] <- "#EE0000"
  
  dat <- tibble::tibble(value = x, group = group)
  
  thr_col <- grep("^Threshold \\(", names(summary_table), value = TRUE)
  if (length(thr_col) != 1L) {
    stop("Could not uniquely identify the threshold column in summary_table.")
  }
  
  ggplot2::ggplot(dat, ggplot2::aes(x = value, fill = group)) +
    ggplot2::geom_density(alpha = 0.5, color = NA) +
    ggplot2::geom_vline(
      data = summary_table,
      ggplot2::aes(xintercept = .data[[thr_col]], color = Method),
      linetype = "dashed",
      linewidth = 2
    ) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::scale_color_brewer(palette = method_palette) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = var_label,
      y = "Density",
      fill = "Group",
      color = "Method"
    ) +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 24),
      legend.title = ggplot2::element_text(size = 24),
      axis.text.x = ggplot2::element_text(size = 36),
      axis.text.y = ggplot2::element_text(size = 36),
      axis.title = ggplot2::element_text(size = 36)
    )
}

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
  
  tab <- gt::text_transform(
    tab,
    locations = gt::cells_body(
      columns = dplyr::all_of(c(
        "Accuracy (95% CI)",
        "Sensitivity (95% CI)",
        "Specificity (95% CI)"
      ))
    ),
    fn = function(x) {
      sub(
        "^([0-9]+\\.[0-9]+)",
        "<strong>\\1</strong>",
        x
      )
    }
  )
  
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

# -------------------------------------------------------------------------
# Build the fixed HNSCC reference from the saved PCA and original expression
# -------------------------------------------------------------------------
build_fixed_hnscc <- function(
    pca_path = "App_files/hpv_signature_tcga_hnscc_pca_model.rds",
    expression_path = "App_files/tcga_primary_hnc.rds",
    hpv_column = "hpv_status",
    pca = NULL,
    tcga = NULL) {
  
  if (is.null(pca)) pca <- readRDS(pca_path)
  if (is.null(tcga)) tcga <- readRDS(expression_path)
  
  genes <- pca$xvars
  samples <- pca$yvars
  expression <- SummarizedExperiment::assay(tcga)
  
  if (!length(genes) || !length(samples) ||
      anyDuplicated(genes) || anyDuplicated(samples) ||
      anyDuplicated(rownames(expression)) ||
      anyDuplicated(colnames(expression)) ||
      !all(genes %in% rownames(expression)) ||
      !all(samples %in% colnames(expression))) {
    stop("Cannot uniquely match the saved PCA genes and samples to HNSCC expression.")
  }
  
  # The stored assay is already processed; do not normalize or log again.
  expression <- as.matrix(expression[genes, samples, drop = FALSE])
  
  if (!is.numeric(expression) || any(!is.finite(expression))) {
    stop("HNSCC reference expression must be numeric and finite.")
  }
  
  if (!all(genes %in% rownames(pca$loadings)) ||
      !all(samples %in% rownames(pca$rotated))) {
    stop("Saved PCA loadings or scores do not match the reference identifiers.")
  }
  
  centers <- rowMeans(expression)
  weights <- setNames(as.numeric(pca$loadings[genes, "PC1"]), genes)
  saved_scores <- as.numeric(pca$rotated[samples, "PC1"])
  
  if (any(!is.finite(weights)) || any(!is.finite(saved_scores))) {
    stop("Saved PCA weights and scores must be finite.")
  }
  
  # Verify that these means and weights reproduce the saved training scores.
  reconstructed <- as.numeric(crossprod(
    weights,
    sweep(expression, 1, centers, "-")
  ))
  
  if (!isTRUE(all.equal(
    reconstructed,
    saved_scores,
    tolerance = 1e-6
  ))) {
    stop(
      "The HNSCC expression does not reproduce the saved PC1 scores. ",
      "Check the reference data and PCA preprocessing before proceeding."
    )
  }
  
  metadata <- as.data.frame(SummarizedExperiment::colData(tcga))
  
  if (!hpv_column %in% names(metadata) ||
      !all(samples %in% rownames(metadata))) {
    stop("Cannot match HNSCC HPV labels to the saved PCA samples.")
  }
  
  labels <- tolower(trimws(as.character(metadata[samples, hpv_column])))
  positive <- which(labels == "positive")
  negative <- which(labels == "negative")
  
  if (!length(positive) || !length(negative)) {
    stop("HNSCC labels must include 'positive' and 'negative'.")
  }
  
  # Training labels establish which direction represents HPV positivity.
  difference <- mean(saved_scores[positive]) -
    mean(saved_scores[negative])
  
  if (!is.finite(difference) || difference == 0) {
    stop("Cannot establish the HPV-positive score direction.")
  }
  
  # Fit the threshold using HNSCC scores only, without status labels.
  gmm <- threshold_gmm(
    saved_scores,
    group = NULL,
    method = "equal_component",
    fallback = "none"
  )
  
  cutoff <- gmm$threshold
  component_means <- range(gmm$params$means)
  
  if (length(cutoff) != 1L || !is.finite(cutoff) ||
      cutoff <= component_means[1] ||
      cutoff >= component_means[2]) {
    stop("Could not identify a valid HNSCC GMM cutoff between component means.")
  }
  
  # Learn the two supervised cutoffs from HNSCC reference labels.
  known <- labels %in% c("negative", "positive")
  
  training_group <- factor(
    labels[known],
    levels = c("negative", "positive")
  )
  
  # Use the same HNSCC-derived direction for all three methods.
  roc_reference <- pROC::roc(
    response = training_group,
    predictor = saved_scores[known],
    levels = c("negative", "positive"),
    direction = if (difference > 0) "<" else ">",
    quiet = TRUE
  )
  
  youden_cutoffs <- as.numeric(unlist(
    pROC::coords(
      roc_reference,
      x = "best",
      best.method = "youden",
      ret = "threshold",
      transpose = FALSE
    ),
    use.names = FALSE
  ))
  
  youden_cutoffs <- youden_cutoffs[is.finite(youden_cutoffs)]
  
  if (!length(youden_cutoffs)) {
    stop("Could not identify a finite HNSCC ROC/Youden cutoff.")
  }
  
  gaussian_reference <- threshold_gaussian_auto(
    saved_scores[known],
    training_group
  )
  
  thresholds <- c(
    youden = youden_cutoffs[1],
    gaussian_auto = gaussian_reference$threshold,
    gmm = cutoff
  )
  
  if (any(!is.finite(thresholds))) {
    stop("One or more HNSCC reference thresholds are not finite.")
  }
  
  list(
    genes = genes,
    center = centers,
    loadings = weights,
    pca_loadings = as.matrix(
      pca$loadings[genes, c("PC1", "PC2"), drop = FALSE]
    ),
    thresholds = thresholds,
    method_labels = c(
      youden = "Fixed ROC/Youden",
      gaussian_auto = "Fixed Gaussian",
      gmm = "Fixed GMM"
    ),
    threshold = cutoff,
    positive_high = difference > 0,
    threshold_method = gmm$method,
    reconstruction_max_abs_error = max(abs(reconstructed - saved_scores))
  )
}


# -------------------------------------------------------------------------
# Apply the fixed HNSCC reference to already processed expression
# -------------------------------------------------------------------------
predict_fixed_hnscc <- function(expression, model) {
  
  expression <- as.matrix(expression)
  
  if (!is.numeric(expression) || ncol(expression) < 1L ||
      is.null(rownames(expression)) ||
      is.null(colnames(expression)) ||
      anyDuplicated(rownames(expression)) ||
      anyDuplicated(colnames(expression))) {
    stop("Expression must be numeric with unique gene and sample names.")
  }
  
  present <- intersect(model$genes, rownames(expression))
  absent <- setdiff(model$genes, rownames(expression))
  
  if (!length(present)) {
    stop("No HNSCC signature genes were found.")
  }
  
  observed <- expression[present, , drop = FALSE]
  
  if (any(!is.finite(observed))) {
    stop("Signature expression must contain finite values.")
  }
  
  # Missing genes contribute zero after HNSCC-mean centering.
  centered <- matrix(
    0,
    nrow = length(model$genes),
    ncol = ncol(expression),
    dimnames = list(model$genes, colnames(expression))
  )
  
  centered[present, ] <- sweep(
    observed,
    1,
    model$center[present],
    "-"
  )
  
  # Project onto the original HNSCC PC1 and PC2.
  scores <- as.data.frame(
    crossprod(
      centered,
      model$pca_loadings[model$genes, , drop = FALSE]
    )
  )
  
  calls <- data.frame(
    sample = rownames(scores),
    scores,
    signature_genes_present = length(present),
    signature_genes_total = length(model$genes),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  for (method in names(model$thresholds)) {
    
    positive <- if (model$positive_high) {
      scores$PC1 >= model$thresholds[[method]]
    } else {
      scores$PC1 <= model$thresholds[[method]]
    }
    
    calls[[paste0("call_fixed_", method)]] <- ifelse(
      positive,
      "Positive",
      "Negative"
    )
  }
  
  list(
    scores = scores,
    calls = calls,
    missing_genes = absent
  )
}

make_fixed_hnscc_report <- function(
    pca,
    clinical = NULL,
    sample_col = NULL,
    hpv_col = NULL,
    flip = FALSE) {
  
  model <- pca$fixed_reference
  calls <- pca$fixed_calls
  methods <- names(model$thresholds)
  has_clinical <- !is.null(clinical)
  
  observed <- rep(NA_character_, nrow(calls))
  
  color_options <- stats::setNames(
    paste0("call_fixed_", methods),
    unname(model$method_labels[methods])
  )
  
  if (has_clinical) {
    
    if (is.null(sample_col) || is.null(hpv_col) ||
        !all(c(sample_col, hpv_col) %in% names(clinical))) {
      stop("Select valid clinical sample-ID and HPV-status columns.")
    }
    
    ids <- as.character(clinical[[sample_col]])
    
    if (anyDuplicated(ids[!is.na(ids)])) {
      stop("Clinical sample identifiers must be unique.")
    }
    
    # Retain every expression sample in the output metadata.
    meta <- as.data.frame(
      clinical[match(calls$sample, ids), , drop = FALSE]
    )
    rownames(meta) <- NULL
    meta[[sample_col]] <- calls$sample
    
    raw_status <- meta[[hpv_col]]
    known <- !is_unknown_label(raw_status)
    
    if (any(known)) {
      
      # Establish label mapping from the supplied clinical table.
      # This does not change any fixed-model predictions.
      mapping <- infer_binary_group(
        clinical[[hpv_col]],
        flip = flip,
        drop_unknown = TRUE
      )
      
      observed[known] <- ifelse(
        trimws(as.character(raw_status[known])) ==
          as.character(mapping$pos),
        "Positive",
        "Negative"
      )
    }
    
    # Append all three fixed calls and the projected scores.
    for (nm in setdiff(names(calls), "sample")) {
      meta[[nm]] <- calls[[nm]]
    }
    
    plot_data <- meta
    plot_data$sample <- calls$sample
    plot_data$hpv_plot <- as.character(raw_status)
    plot_data$hpv_plot[!known] <- "Unknown"
    
    color_options <- c(
      "HPV status" = "hpv_plot",
      color_options
    )
    
  } else {
    meta <- calls
    plot_data <- calls
  }
  
  levels_hpv <- c("Negative", "Positive")
  observed <- factor(observed, levels = levels_hpv)
  evaluable <- !is.na(observed)
  
  safe_ratio <- function(numerator, denominator) {
    if (denominator > 0) numerator / denominator else NA_real_
  }
  
  # Evaluate the already assigned calls.
  # No target-cohort threshold fitting or direction selection occurs here.
  results <- stats::setNames(
    lapply(methods, function(method) {
      
      predicted <- factor(
        calls[[paste0("call_fixed_", method)]],
        levels = levels_hpv
      )
      
      confusion <- table(
        Observed = observed[evaluable],
        Predicted = predicted[evaluable]
      )
      
      tp <- as.numeric(confusion["Positive", "Positive"])
      tn <- as.numeric(confusion["Negative", "Negative"])
      fp <- as.numeric(confusion["Negative", "Positive"])
      fn <- as.numeric(confusion["Positive", "Negative"])
      
      method_label <- unname(model$method_labels[[method]])
      threshold <- unname(model$thresholds[[method]])
      
      list(
        method = method_label,
        threshold = threshold,
        evaluation = list(
          method = method_label,
          threshold = threshold,
          accuracy = safe_ratio(tp + tn, tp + tn + fp + fn),
          sensitivity = safe_ratio(tp, tp + fn),
          specificity = safe_ratio(tn, tn + fp),
          confusion = confusion
        )
      )
    }),
    methods
  )
  
  # Reuse the existing exact-CI and table-formatting functions.
  summary_raw <- threshold_summary_table(results)
  
  summary_table <- format_threshold_summary(
    summary_raw,
    var_label = "PC1 score",
    digits = 3,
    include_confusion = FALSE
  )
  
  summary_table[["Evaluable samples"]] <- sum(evaluable)
  
  summary_gt <- threshold_gt_table(
    summary_table,
    subtitle = sprintf(
      "%d/%d signature genes available; cutoffs were estimated in HNSCC.",
      calls$signature_genes_present[1],
      calls$signature_genes_total[1]
    )
  )
  
  # Threshold distributions:
  # With clinical labels, use only known reference classifications.
  # Without clinical labels, show one pooled sample distribution.
  if (has_clinical) {
    plot_scores <- calls$PC1[evaluable]
    plot_group <- droplevels(observed[evaluable])
  } else {
    plot_scores <- calls$PC1
    plot_group <- factor(rep("All samples", nrow(calls)))
  }
  
  if (!length(plot_scores)) {
    
    threshold_plot <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = "No known reference labels available for the threshold plot."
      ) +
      ggplot2::theme_void()
    
  } else if (all(table(droplevels(plot_group)) >= 2L)) {
    
    threshold_plot <- plot_thresholds(
      x = plot_scores,
      group = plot_group,
      summary_table = summary_table,
      var_label = "PC1 score",
      method_palette = "Dark2"
    )
    
  } else {
    
    # Density estimation requires more than one observation per group.
    threshold_plot <- ggplot2::ggplot(
      data.frame(
        PC1 = plot_scores,
        Group = plot_group
      ),
      ggplot2::aes(x = PC1, fill = Group)
    ) +
      ggplot2::geom_histogram(
        bins = 30,
        position = "identity",
        alpha = 0.5
      ) +
      ggplot2::geom_vline(
        data = summary_table,
        ggplot2::aes(
          xintercept = .data[["Threshold (PC1 score)"]],
          color = Method
        ),
        linetype = "dashed",
        linewidth = 1
      ) +
      ggplot2::scale_color_brewer(palette = "Dark2") +
      ggplot2::theme_classic() +
      ggplot2::labs(
        x = "PC1 score",
        y = "Samples",
        color = "Method"
      )
  }
  
  # A single ROC curve evaluates the shared projected PC1 scores.
  roc_result <- NULL
  
  if (length(unique(as.character(observed[evaluable]))) == 2L) {
    
    roc_object <- pROC::roc(
      response = observed[evaluable],
      predictor = calls$PC1[evaluable],
      levels = levels_hpv,
      direction = if (model$positive_high) "<" else ">",
      quiet = TRUE
    )
    
    auc_ci <- tryCatch(
      as.numeric(pROC::ci.auc(roc_object, method = "delong")),
      error = function(e) rep(NA_real_, 3)
    )
    
    roc_result <- list(
      auc = as.numeric(pROC::auc(roc_object)),
      auc_ci = auc_ci,
      roc_plot = pROC::ggroc(roc_object) +
        ggplot2::theme_classic() +
        ggplot2::labs(title = "Fixed HNSCC score ROC")
    )
  }
  
  threshold_res <- list(
    results = results,
    summary_raw = summary_raw,
    summary_table = summary_table,
    summary_gt = summary_gt,
    plot = threshold_plot
  )
  
  list(
    mode = if (has_clinical) "supervised" else "unsupervised",
    scoring_method = "fixed",
    pca_object = pca,
    pca_scores = calls[, c("sample", "PC1", "PC2")],
    threshold_res = threshold_res,
    threshold_plot = threshold_plot,
    pca_data = plot_data,
    pca_color_options = color_options,
    clin_augmented = if (has_clinical) meta else NULL,
    meta_unsupervised = if (!has_clinical) meta else NULL,
    gmm_threshold = unname(model$thresholds[["gmm"]]),
    gmm_summary_table = summary_table,
    gmm_summary_gt = summary_gt,
    hpv_col_name = if (has_clinical) hpv_col else NULL,
    fixed_roc = roc_result
  )
}

make_cohort_hpv_report <- function(pca, clinical = NULL,
                                   sample_col = NULL, hpv_col = NULL,
                                   flip = FALSE) {
  pc <- as.data.frame(pca$rotated)
  pc$sample <- rownames(pc)
  pc <- pc[order(pc$sample), , drop = FALSE]
  
  if (anyDuplicated(pc$sample) || any(!is.finite(pc$PC1))) {
    stop("PCA sample IDs must be unique and PC1 scores must be finite.")
  }
  
  if (!"PC2" %in% names(pc)) pc$PC2 <- NA_real_
  
  x <- setNames(pc$PC1, pc$sample)
  hpv_levels <- c("Negative", "Positive")
  
  # Reuse the GMM fitted with PCA, independently of clinical information.
  gmm <- pca$cohort_gmm
  if (is.null(gmm)) {
    stop("Cohort GMM fit is missing. Regenerate the PCA.")
  }
  
  cutoff <- gmm$threshold
  if (length(cutoff) != 1L || !is.finite(cutoff)) {
    stop("Could not determine an equal-component GMM cutoff.")
  }
  
  direction <- pca$hpv_positive_high
  oriented <- is.logical(direction) &&
    length(direction) == 1L &&
    !is.na(direction)
  
  gmm_calls <- if (oriented) {
    positive <- if (direction) x >= cutoff else x <= cutoff
    factor(
      ifelse(positive, "Positive", "Negative"),
      levels = hpv_levels
    )
  } else {
    factor(
      ifelse(x >= cutoff, "GMM_high", "GMM_low"),
      levels = c("GMM_low", "GMM_high")
    )
  }
  
  has_clinical <- !is.null(clinical)
  observed <- factor(
    rep(NA_character_, length(x)),
    levels = hpv_levels
  )
  meta <- pc[, c("sample", "PC1", "PC2"), drop = FALSE]
  mapping_available <- FALSE
  
  if (has_clinical) {
    clinical <- as.data.frame(clinical)
    
    if (length(sample_col) != 1L || length(hpv_col) != 1L ||
        !all(c(sample_col, hpv_col) %in% names(clinical))) {
      stop("Select valid clinical sample-ID and HPV-status columns.")
    }
    
    ids <- as.character(clinical[[sample_col]])
    if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
      stop("Clinical sample IDs must be nonmissing and unique.")
    }
    
    # Match metadata TO all expression samples.
    matched <- match(pc$sample, ids)
    meta <- clinical[matched, , drop = FALSE]
    meta[[sample_col]] <- pc$sample
    meta$sample <- pc$sample
    meta$PC1 <- pc$PC1
    meta$PC2 <- pc$PC2
    rownames(meta) <- NULL
    
    # Replace any calls from a previously exported report.
    meta <- meta[, !names(meta) %in% c(
      "call_youden", "call_gaussian_auto", "call_gmm"
    ), drop = FALSE]
    
    labels <- clinical[[hpv_col]]
    known_levels <- unique(
      trimws(as.character(labels[!is_unknown_label(labels)]))
    )
    
    if (length(known_levels) == 2L) {
      inf <- infer_binary_group(
        labels,
        flip = flip,
        drop_unknown = TRUE
      )
      
      mapped <- as.character(inf$full_group)[matched]
      observed <- factor(
        ifelse(
          is.na(mapped),
          NA_character_,
          ifelse(mapped == inf$pos, "Positive", "Negative")
        ),
        levels = hpv_levels
      )
      mapping_available <- TRUE
    }
    
    display_labels <- as.character(labels[matched])
    display_labels[is_unknown_label(display_labels)] <- "Unknown"
    meta$hpv_plot <- factor(display_labels)
  }
  
  known <- !is.na(observed)
  results <- list()
  group_n <- table(observed[known])
  
  # Only supervised thresholds learn from the known-label subset.
  if (sum(known) >= 4L && all(group_n > 0L)) {
    results$youden <- threshold_youden(
      x[known], observed[known]
    )
  }
  
  if (all(group_n >= 2L)) {
    results$gaussian_auto <- threshold_gaussian_auto(
      x[known], observed[known]
    )
  }
  
  # Apply supervised thresholds to every expression sample.
  if (length(results)) {
    positive_high <- mean(x[known & observed == "Positive"]) >=
      mean(x[known & observed == "Negative"])
    
    for (method in names(results)) {
      threshold <- results[[method]]$threshold
      positive <- if (positive_high) {
        x >= threshold
      } else {
        x <= threshold
      }
      
      meta[[paste0("call_", method)]] <- factor(
        ifelse(positive, "Positive", "Negative"),
        levels = hpv_levels
      )
    }
  }
  
  meta$call_gmm <- gmm_calls
  
  # Evaluate existing GMM calls without changing their direction.
  evaluable <- known & oriented
  confusion <- table(
    Observed = observed[evaluable],
    Predicted = factor(
      as.character(gmm_calls[evaluable]),
      levels = hpv_levels
    )
  )
  
  tp <- confusion["Positive", "Positive"]
  tn <- confusion["Negative", "Negative"]
  fp <- confusion["Negative", "Positive"]
  fn <- confusion["Positive", "Negative"]
  
  ratio <- function(a, b) {
    if (b > 0) as.numeric(a / b) else NA_real_
  }
  
  gmm$evaluation <- list(
    method = "GMM",
    threshold = cutoff,
    accuracy = ratio(tp + tn, tp + tn + fp + fn),
    sensitivity = ratio(tp, tp + fn),
    specificity = ratio(tn, tn + fp),
    confusion = confusion
  )
  results$gmm <- gmm
  
  summary_raw <- threshold_summary_table(results)
  summary_table <- format_threshold_summary(
    summary_raw,
    var_label = "PC1 score"
  )
  
  summary_table[["Fitting samples"]] <- c(
    rep(sum(known), length(results) - 1L),
    length(x)
  )
  summary_table[["Evaluable samples"]] <- c(
    rep(sum(known), length(results) - 1L),
    sum(evaluable)
  )
  
  note <- if (!oriented) {
    paste(
      "HNSCC direction unavailable: GMM calls are high/low groups;",
      "HPV performance is not calculated."
    )
  } else if (has_clinical && !mapping_available) {
    paste(
      "GMM uses all expression samples.",
      "Two clinical label levels are required for comparison."
    )
  } else {
    paste(
      "GMM uses all expression samples;",
      "performance uses known clinical labels only."
    )
  }
  
  summary_gt <- threshold_gt_table(
    summary_table,
    subtitle = note
  )
  
  # Keep unknown clinical labels out of the comparison density plot.
  plot_idx <- if (any(known)) known else rep(TRUE, length(x))
  plot_group <- if (any(known)) {
    droplevels(observed[known])
  } else {
    droplevels(gmm_calls)
  }
  
  threshold_plot <- plot_thresholds(
    x[plot_idx],
    plot_group,
    summary_table,
    var_label = "PC1 score",
    method_palette = "Dark2"
  )
  
  if (any(table(plot_group) < 2L)) {
    threshold_plot$layers[[1]] <- ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 20,
      position = "identity",
      alpha = 0.5,
      color = NA
    )
  }
  
  color_options <- if (has_clinical) {
    c("HPV status" = "hpv_plot")
  } else {
    character()
  }
  
  call_columns <- grep("^call_", names(meta), value = TRUE)
  color_options <- c(
    color_options,
    setNames(
      call_columns,
      paste0("Call: ", sub("^call_", "", call_columns))
    )
  )
  
  threshold_res <- list(
    results = results,
    summary_raw = summary_raw,
    summary_table = summary_table,
    summary_gt = summary_gt,
    plot = threshold_plot
  )
  
  list(
    mode = if (has_clinical) "supervised" else "unsupervised",
    scoring_method = "cohort",
    pca_scores = pc,
    threshold_res = threshold_res,
    threshold_plot = threshold_plot,
    pca_data = meta,
    pca_color_options = color_options,
    clin_augmented = if (has_clinical) meta else NULL,
    meta_unsupervised = if (!has_clinical) meta else NULL,
    gmm_threshold = cutoff,
    gmm_summary_table = summary_table,
    gmm_summary_gt = summary_gt,
    hpv_col_name = if (has_clinical) hpv_col else NULL
  )
}

make_signature_heatmap <- function(reference_expression, new_expression,
                                   reference_hpv, new_hpv = NULL,
                                   annotation_label = "HPV call") {
  
  # Match shared genes in the same order.
  genes <- intersect(
    rownames(reference_expression),
    rownames(new_expression)
  )
  
  if (length(genes) < 2L || ncol(new_expression) < 2L) {
    stop("The comparison heatmap requires at least two shared genes and two uploaded samples.")
  }
  
  ref <- as.matrix(reference_expression[genes, , drop = FALSE])
  new <- as.matrix(new_expression[genes, , drop = FALSE])
  
  # Gene-wise z-scores within each cohort, for display only.
  ref_z <- t(scale(t(ref)))
  new_z <- t(scale(t(new)))
  
  # Exclude genes whose z-scores are undefined in either cohort.
  keep <- rowSums(!is.finite(ref_z)) == 0 &
    rowSums(!is.finite(new_z)) == 0
  
  omitted <- genes[!keep]
  ref_z <- ref_z[keep, , drop = FALSE]
  new_z <- new_z[keep, , drop = FALSE]
  
  if (nrow(ref_z) < 2L) {
    stop("Fewer than two shared genes have finite, variable expression in both cohorts.")
  }
  
  # Prepare annotation labels and colors.
  make_annotation <- function(labels, title) {
    if (is.null(labels)) return(NULL)
    
    labels <- trimws(as.character(labels))
    labels[is_unknown_label(labels)] <- "Unknown"
    normalized <- tolower(labels)
    
    labels[normalized %in% c(
      "positive", "pos", "true", "1", "hpv+"
    )] <- "Positive"
    
    labels[normalized %in% c(
      "negative", "neg", "false", "0", "hpv-"
    )] <- "Negative"
    
    categories <- unique(labels)
    
    colors <- stats::setNames(
      grDevices::hcl.colors(length(categories), "Dark 3"),
      categories
    )
    
    hpv_colors <- c(
      Positive = "#EE0000",
      Negative = "#3B4992",
      Unknown = "#A6A4A4"
    )
    
    for (label in intersect(categories, names(hpv_colors))) {
      colors[label] <- hpv_colors[[label]]
    }
    
    ComplexHeatmap::HeatmapAnnotation(
      Status = labels,
      col = list(Status = colors),
      show_annotation_name = FALSE,
      annotation_legend_param = list(
        Status = list(title = title)
      )
    )
  }
  
  # Use the same expression color scale for both cohorts.
  color_scale <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("#3B4992", "white", "#EE0000")
  )
  
  reference_heatmap <- ComplexHeatmap::Heatmap(
    ref_z,
    name = "HNSCC",
    col = color_scale,
    column_title = paste0("TCGA HNSCC (n = ", ncol(ref_z), ")"),
    show_column_names = FALSE,
    show_row_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "complete",
    top_annotation = make_annotation(
      reference_hpv,
      "HNSCC HPV status"
    ),
    heatmap_legend_param = list(title = "Gene z-score"),
    use_raster = FALSE
  )
  
  # Shares the HNSCC row order when the two heatmaps are drawn together.
  uploaded_heatmap <- ComplexHeatmap::Heatmap(
    new_z,
    name = "Uploaded",
    col = color_scale,
    column_title = paste0("Uploaded cohort (n = ", ncol(new_z), ")"),
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 9),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "complete",
    top_annotation = make_annotation(
      new_hpv,
      annotation_label
    ),
    show_heatmap_legend = FALSE,
    use_raster = FALSE
  )
  
  list(
    heatmap = reference_heatmap + uploaded_heatmap,
    genes = rownames(ref_z),
    omitted_genes = omitted,
    reference_z = ref_z,
    uploaded_z = new_z
  )
}
