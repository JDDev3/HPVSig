
setwd("Figure_Generation")
penile_data_abbreviated <- readr::read_delim("data/PSCC_Multiomic_Clincal_Data.txt")
export_folder_name <- "Exports"

#### Generation of Figure3B
library(ComplexHeatmap)

anno_col <- list(
  HPV_ISH = c(
    "Negative" = "#3b4992ff",
    "Positive" = "#ee0000ff",
    "Unknown" = "#A6A4A4"
  ),
  p16_IHC = c(
    "Negative" = "#3b4992ff",
    "Positive" = "#ee0000ff",
    "Unknown" = "#A6A4A4"
  ),
  HPV_EM = c(
    "Negative" = "#3b4992ff",
    "Positive" = "#ee0000ff"
  ),
  HPV_Signature = c(
    "Negative" = "#3b4992ff",
    "Positive" = "#ee0000ff"
  )
)

generate_counted_column <- function(count_col_name, count_data, count_data_columns, count_value = "Positive", sample_id = 1){
  total_count_name <- paste0(count_col_name, "_total_count")
  value_count_name <- paste0(count_col_name, "_", count_value, "_count")
  result_count_name <- paste0(count_col_name, "_", "count_percent")
  
  count_data_results <- as.data.frame(rowSums(count_data[count_data_columns] != "Unknown"))
  colnames(count_data_results) <- total_count_name
  
  count_data_results[value_count_name] <- rowSums(count_data[count_data_columns] == count_value)
  count_data_results[result_count_name] <- (count_data_results[2]/count_data_results[1]) * 100
  count_data_results <- cbind(count_data, count_data_results)
  return(count_data_results)
}


penile_data_abbreviated_counted <- generate_counted_column("combined_HPV_EM", penile_data_abbreviated, c("HPV_ISH", "p16_IHC", "HPV_EM"))
penile_data_abbreviated_counted <- generate_counted_column("combined_HPV_Signature", penile_data_abbreviated_counted, c("HPV_ISH", "p16_IHC", "HPV_Signature"))
penile_data_abbreviated_counted <- generate_counted_column("p16IHC_HPV_Sig", penile_data_abbreviated_counted, c("p16_IHC", "HPV_Signature"))
penile_data_abbreviated_counted <- generate_counted_column("p16IHC_HPV_EM", penile_data_abbreviated_counted, c("p16_IHC", "HPV_EM"))
penile_data_abbreviated_counted <- generate_counted_column("HPVISH_HPV_Sig", penile_data_abbreviated_counted, c("HPV_ISH", "HPV_Signature"))
penile_data_abbreviated_counted <- generate_counted_column("HPVISH_HPV_EM", penile_data_abbreviated_counted, c("HPV_ISH", "HPV_EM"))
penile_data_abbreviated_counted <- generate_counted_column("Cumulative_HPV_Call", penile_data_abbreviated_counted, c("HPV_ISH", "p16_IHC", "HPV_Signature", "HPV_EM"))

penile_data_abbreviated_counted$disagreement_calls <- sapply(penile_data_abbreviated_counted$Cumulative_HPV_Call_count_percent, function(x) if(x == 0){"Negative"} else if(x == 100){"Positive"} else if (x > 0 & x < 100){"Disagreement"})
penile_data_abbreviated_counted$combined_HPV_sig <- sapply(penile_data_abbreviated_counted$combined_HPV_Signature_count_percent, function(x) if(x < 50){"Negative"} else if(x > 50){"Positive"}else if (x == 50){"Unknown"})
penile_data_abbreviated_counted$combined_HPV_EM <- sapply(penile_data_abbreviated_counted$combined_HPV_EM_count_percent, function(x) if(x < 50){"Negative"} else if(x > 50){"Positive"}else if (x == 50){"Unknown"})
penile_data_abbreviated_counted$p16_HPV_sig <- sapply(penile_data_abbreviated_counted$p16IHC_HPV_Sig_count_percent, function(x) if(x == 0){"Negative"} else if(x == 100){"Positive"}else{"Unknown"})
penile_data_abbreviated_counted$p16_HPV_EM <- sapply(penile_data_abbreviated_counted$p16IHC_HPV_EM_count_percent, function(x) if(x == 0){"Negative"} else if(x == 100){"Positive"}else{"Unknown"})
penile_data_abbreviated_counted$ISH_HPV_sig <- sapply(penile_data_abbreviated_counted$HPVISH_HPV_Sig_count_percent, function(x) if(x == 0){"Negative"} else if(x == 100){"Positive"}else{"Unknown"})
penile_data_abbreviated_counted$ISH_HPV_EM <- sapply(penile_data_abbreviated_counted$HPVISH_HPV_EM_count_percent, function(x) if(x == 0){"Negative"} else if(x == 100){"Positive"}else{"Unknown"})


annotation_coloring <- sapply(penile_data_abbreviated_counted$disagreement_calls, function(x) if(x == "Positive"){"#ee0000ff"}else if(x == "Negative"){"#3b4992ff"}else if(x == "Disagreement"){"#8E159E"})

top_annotation <- ComplexHeatmap::HeatmapAnnotation(
  Sample = ComplexHeatmap::anno_text(penile_data_abbreviated$PT_Num_Multiomic_Clin, rot = 90, gp = gpar(col = annotation_coloring, fontface = "bold")),
  HPV_ISH = penile_data_abbreviated$HPV_ISH,
  p16_IHC = penile_data_abbreviated$p16_IHC,
  HPV_EM = penile_data_abbreviated$HPV_EM,
  HPV_Signature = penile_data_abbreviated$HPV_Signature,
  col = anno_col, 
  show_legend = TRUE,
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)


dummy_mat <- matrix(0, nrow=1, ncol = nrow(penile_data_abbreviated))

pdf(here::here(export_folder_name, "HPV_Call_annotation.pdf"))

heatmap <- ComplexHeatmap::Heatmap(
  dummy_mat,
  top_annotation = top_annotation,
  show_heatmap_legend = FALSE,
  cluster_columns = FALSE
) |> draw()

dev.off()

#generation of agreement tables (Figure 3C)

library(irr)
library(tibble)
library(gt)
library(tidyverse)

style_agreement_contingency_gt <- function(gt_tbl,
                                           nameA,
                                           positive_levels = c("Positive", "Pos", "HPV+", "Yes", "1"),
                                           negative_levels = c("Negative", "Neg", "HPV-", "No", "0"),
                                           unknown_levels  = c("Unknown", "NA"),
                                           col_pos = "#ee0000ff",
                                           col_neg = "#3b4992ff",
                                           col_dis = "#8E159E",
                                           col_unk = "grey70") {
  pos_fill <- scales::alpha(col_pos, 0.3)
  neg_fill <- scales::alpha(col_neg, 0.3)
  dis_fill <- scales::alpha(col_dis, 0.3)
  unk_fill <- scales::alpha(col_unk, 0.3)
  tot_fill <- scales::alpha("grey92", 1)
  
  symA <- rlang::sym(nameA)
  
  gt_tbl <- gt_tbl |>
    gt::tab_style(
      style = gt::cell_fill(color = pos_fill),
      locations = gt::cells_stub(rows = !!symA %in% positive_levels)
    ) |>
    gt::tab_style(
      style = gt::cell_fill(color = neg_fill),
      locations = gt::cells_stub(rows = !!symA %in% negative_levels)
    ) |>
    gt::tab_style(
      style = gt::cell_fill(color = unk_fill),
      locations = gt::cells_stub(rows = !!symA %in% unknown_levels)
    ) |>
    gt::tab_style(
      style = list(gt::cell_fill(color = tot_fill), gt::cell_text(weight = "bold")),
      locations = gt::cells_stub(rows = !!symA == "Total")
    )
  
  # ---- Style column labels
  gt_tbl <- gt_tbl |>
    gt::tab_style(
      style = gt::cell_fill(color = pos_fill),
      locations = gt::cells_column_labels(columns = any_of(positive_levels))
    ) |>
    gt::tab_style(
      style = gt::cell_fill(color = neg_fill),
      locations = gt::cells_column_labels(columns = any_of(negative_levels))
    ) |>
    gt::tab_style(
      style = gt::cell_fill(color = unk_fill),
      locations = gt::cells_column_labels(columns = any_of(unknown_levels))
    ) |>
    gt::tab_style(
      style = list(gt::cell_fill(color = tot_fill), gt::cell_text(weight = "bold")),
      locations = gt::cells_column_labels(columns = any_of("Total"))
    )
  
  # ---- Disagreement blocks: Pos (rows) x Neg (cols) and Neg x Pos
  gt_tbl <- gt_tbl |>
    gt::tab_style(
      style = gt::cell_fill(color = dis_fill),
      locations = gt::cells_body(
        rows    = !!symA %in% positive_levels,
        columns = any_of(negative_levels)
      )
    ) |>
    gt::tab_style(
      style = gt::cell_fill(color = dis_fill),
      locations = gt::cells_body(
        rows    = !!symA %in% negative_levels,
        columns = any_of(positive_levels)
      )
    )
  
  # ---- Unknown/NA blocks (if present)
  gt_tbl <- gt_tbl |>
    gt::tab_style(
      style = gt::cell_fill(color = unk_fill),
      locations = gt::cells_body(
        rows    = !!symA %in% unknown_levels,
        columns = everything()
      )
    ) |>
    gt::tab_style(
      style = gt::cell_fill(color = unk_fill),
      locations = gt::cells_body(
        rows    = everything(),
        columns = any_of(unknown_levels)
      )
    )
  
  # ---- Diagonal agreements: tint cells where row level == column level
  diag_levels <- unique(c(positive_levels, negative_levels, unknown_levels))
  for (lvl in diag_levels) {
    if (lvl %in% c(positive_levels, negative_levels, unknown_levels)) {
      fill <- if (lvl %in% positive_levels) pos_fill else if (lvl %in% negative_levels) neg_fill else unk_fill
      gt_tbl <- gt_tbl |>
        gt::tab_style(
          style = list(gt::cell_fill(color = fill), gt::cell_text(weight = "bold")),
          locations = gt::cells_body(rows = !!symA == lvl, columns = any_of(lvl))
        )
    }
  }
  
  # ---- Total row/col styling in body
  gt_tbl <- gt_tbl |>
    gt::tab_style(
      style = list(gt::cell_fill(color = tot_fill), gt::cell_text(weight = "bold")),
      locations = gt::cells_body(rows = !!symA == "Total", columns = everything())
    ) |>
    gt::tab_style(
      style = list(gt::cell_fill(color = tot_fill), gt::cell_text(weight = "bold")),
      locations = gt::cells_body(rows = everything(), columns = any_of("Total"))
    )
  
  gt_tbl
}

agreement_summary_gt <- function(HPVCallMethodA,
                                 HPVCallMethodB,
                                 unknown_levels = "Unknown",
                                 correct = TRUE,
                                 mcnemar_method = c("exact", "chisq")) {
  
  mcnemar_method <- match.arg(mcnemar_method)
  
  # Accept a vector or a one-column data frame
  extract_col <- function(x, default_name) {
    if (is.data.frame(x)) {
      if (ncol(x) != 1L) {
        stop("If you pass a data.frame, it must have exactly one column.")
      }
      
      vec <- x[[1]]
      name <- names(x)[1]
    } else {
      vec <- x
      name <- default_name
    }
    
    list(vec = vec, name = name)
  }
  
  format_p <- function(x) {
    ifelse(
      !is.finite(x),
      NA_character_,
      ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
    )
  }
  
  A <- extract_col(HPVCallMethodA, "HPVCallMethodA")
  B <- extract_col(HPVCallMethodB, "HPVCallMethodB")
  
  a <- as.character(A$vec)
  b <- as.character(B$vec)
  nameA <- A$name
  nameB <- B$name
  
  if (length(a) != length(b)) {
    stop(
      "HPVCallMethodA and HPVCallMethodB must contain the same number of observations."
    )
  }
  
  # ------------------------------------------------------------------
  # 1) Raw contingency table: include Unknown and NA
  # ------------------------------------------------------------------
  a_raw <- a
  b_raw <- b
  
  a_raw[is.na(a_raw)] <- "NA"
  b_raw[is.na(b_raw)] <- "NA"
  
  levels_all <- sort(unique(c(a_raw, b_raw)))
  
  raw_tab <- table(
    factor(a_raw, levels = levels_all),
    factor(b_raw, levels = levels_all),
    dnn = c(nameA, nameB)
  )
  
  raw_tab_tot <- addmargins(raw_tab)
  
  dn <- dimnames(raw_tab_tot)
  dn[[1]][length(dn[[1]])] <- "Total"
  dn[[2]][length(dn[[2]])] <- "Total"
  dimnames(raw_tab_tot) <- dn
  
  # ------------------------------------------------------------------
  # 2) Exclude pairs with Unknown or NA from statistical calculations
  # ------------------------------------------------------------------
  is_bad <- is.na(a) |
    is.na(b) |
    a %in% unknown_levels |
    b %in% unknown_levels
  
  keep <- !is_bad
  n_total <- length(a)
  n_used <- sum(keep)
  
  if (n_used == 0L) {
    stop("No rows left after removing Unknown/NA.")
  }
  
  a_use <- a[keep]
  b_use <- b[keep]
  
  levels_binary <- sort(unique(c(a_use, b_use)))
  
  analysis_tab <- table(
    factor(a_use, levels = levels_binary),
    factor(b_use, levels = levels_binary),
    dnn = c(nameA, nameB)
  )
  
  if (nrow(analysis_tab) != 2L || ncol(analysis_tab) != 2L) {
    stop("McNemar's test requires a 2×2 table after removing Unknown/NA.")
  }
  
  # ------------------------------------------------------------------
  # 3) Cohen's kappa and approximate 95% confidence interval
  # ------------------------------------------------------------------
  kappa_res <- irr::kappa2(
    data.frame(a_use, b_use),
    weight = "unweighted"
  )
  
  n_pairs <- sum(analysis_tab)
  cell_prob <- analysis_tab / n_pairs
  row_prob <- rowSums(cell_prob)
  col_prob <- colSums(cell_prob)
  
  observed_agreement <- sum(diag(cell_prob))
  expected_agreement <- sum(row_prob * col_prob)
  
  kappa_ci_low <- NA_real_
  kappa_ci_high <- NA_real_
  
  if (n_pairs > 1 && expected_agreement < 1) {
    
    # Multinomial delta-method variance for the estimated kappa
    gradient <- diag(nrow(cell_prob)) / (1 - expected_agreement) -
      ((1 - observed_agreement) / (1 - expected_agreement)^2) *
      outer(col_prob, row_prob, "+")
    
    mean_gradient <- sum(cell_prob * gradient)
    
    kappa_variance <- sum(
      cell_prob * (gradient - mean_gradient)^2
    ) / n_pairs
    
    kappa_se <- sqrt(kappa_variance)
    margin <- stats::qnorm(0.975) * kappa_se
    
    # Restrict reported limits to the possible range of kappa
    kappa_ci_low <- max(-1, kappa_res$value - margin)
    kappa_ci_high <- min(1, kappa_res$value + margin)
  }
  
  # ------------------------------------------------------------------
  # 4) McNemar: exact binomial or chi-square approximation
  # ------------------------------------------------------------------
  discordant_ab <- unname(analysis_tab[1, 2])
  discordant_ba <- unname(analysis_tab[2, 1])
  n_discordant <- discordant_ab + discordant_ba
  
  mcnemar_chi2 <- NA_real_
  mcnemar_df <- NA_real_
  
  if (mcnemar_method == "exact") {
    
    mcnemar_p <- if (n_discordant == 0) {
      1
    } else {
      stats::binom.test(
        x = discordant_ab,
        n = n_discordant,
        p = 0.5,
        alternative = "two.sided"
      )$p.value
    }
    
    mcnemar_label <- "Exact two-sided McNemar"
    
  } else {
    
    mc_res <- stats::mcnemar.test(
      analysis_tab,
      correct = correct
    )
    
    mcnemar_p <- mc_res$p.value
    mcnemar_chi2 <- unname(mc_res$statistic)
    mcnemar_df <- unname(mc_res$parameter)
    
    mcnemar_label <- if (correct) {
      "McNemar chi-square (continuity corrected)"
    } else {
      "McNemar chi-square (uncorrected)"
    }
  }
  
  stats_tbl <- tibble::tibble(
    HPVCallMethodA = nameA,
    HPVCallMethodB = nameB,
    n_total = n_total,
    n_used = n_used,
    kappa = kappa_res$value,
    kappa_ci_low = kappa_ci_low,
    kappa_ci_high = kappa_ci_high,
    kappa_p = kappa_res$p.value,
    discordant_ab = discordant_ab,
    discordant_ba = discordant_ba,
    n_discordant = n_discordant,
    mcnemar_method = mcnemar_label,
    mcnemar_chi2 = mcnemar_chi2,
    mcnemar_df = mcnemar_df,
    mcnemar_p = mcnemar_p
  )
  
  # ------------------------------------------------------------------
  # 5) Raw contingency table with statistical annotation
  # ------------------------------------------------------------------
  raw_df_tot <- as.data.frame.matrix(
    raw_tab_tot,
    stringsAsFactors = FALSE
  )
  
  raw_df_tot <- data.frame(
    Category = rownames(raw_df_tot),
    raw_df_tot,
    row.names = NULL,
    check.names = FALSE
  )
  
  names(raw_df_tot)[1] <- nameA
  
  ci_part <- if (all(is.finite(c(kappa_ci_low, kappa_ci_high)))) {
    sprintf("95%% CI %.3f–%.3f", kappa_ci_low, kappa_ci_high)
  } else {
    "95% CI unavailable"
  }
  
  mcnemar_p_text <- if (!is.finite(mcnemar_p)) {
    "= NA"
  } else if (mcnemar_p < 0.001) {
    "< 0.001"
  } else {
    sprintf("= %.3f", mcnemar_p)
  }
  
  stats_text <- gt::html(sprintf(
    paste0(
      "Cohen's kappa = <strong>%.3f</strong> (%s)",
      "<br>%s: p %s",
      "<br>n = %d (used = %d)."
    ),
    stats_tbl$kappa,
    ci_part,
    mcnemar_label,
    mcnemar_p_text,
    n_total,
    n_used
  ))
  
  num_cols_counts <- names(raw_df_tot)[
    vapply(raw_df_tot, is.numeric, logical(1))
  ]
  
  table_gt <- gt::gt(
    raw_df_tot,
    rowname_col = names(raw_df_tot)[1]
  ) |>
    gt::tab_stubhead(label = names(raw_df_tot)[1]) |>
    gt::tab_spanner(
      label = nameB,
      columns = dplyr::all_of(num_cols_counts)
    ) |>
    gt::fmt_number(
      columns = dplyr::all_of(num_cols_counts),
      decimals = 0
    ) |>
    gt::tab_options(
      table.font.size = 60,
      heading.title.font.size = 48,
      column_labels.font.size = 48,
      stub.font.size = 48,
      data_row.padding = gt::px(2),
      table_body.hlines.width = gt::px(0.5),
      table.border.top.width = gt::px(1),
      table.border.bottom.width = gt::px(1)
    ) |>
    gt::tab_source_note(source_note = stats_text)
  
  table_gt <- table_gt |>
    style_agreement_contingency_gt(
      nameA = nameA,
      positive_levels = "Positive",
      negative_levels = "Negative",
      unknown_levels = c(unknown_levels, "NA"),
      col_pos = "#ee0000ff",
      col_neg = "#3b4992ff",
      col_dis = "#8E159E",
      col_unk = "grey70"
    )
  
  # ------------------------------------------------------------------
  # 6) Analysis contingency table: Unknown and NA excluded
  # ------------------------------------------------------------------
  analysis_df <- as.data.frame.matrix(
    analysis_tab,
    stringsAsFactors = FALSE
  )
  
  analysis_df <- data.frame(
    Category = rownames(analysis_df),
    analysis_df,
    row.names = NULL,
    check.names = FALSE
  )
  
  names(analysis_df)[1] <- nameA
  
  analysis_gt <- gt::gt(analysis_df) |>
    gt::tab_header(
      title = paste0(
        "Contingency used for stats (",
        nameA, " vs ", nameB, ")"
      )
    )
  
  # ------------------------------------------------------------------
  # 7) Statistics-only table
  # ------------------------------------------------------------------
  stats_gt <- gt::gt(stats_tbl) |>
    gt::fmt_number(
      columns = dplyr::all_of(c(
        "kappa",
        "kappa_ci_low",
        "kappa_ci_high",
        "mcnemar_chi2"
      )),
      decimals = 5
    ) |>
    gt::fmt(
      columns = dplyr::all_of(c("kappa_p", "mcnemar_p")),
      fns = format_p
    ) |>
    gt::tab_header(
      title = "Agreement statistics (Cohen's kappa & McNemar's test)"
    )
  
  list(
    raw_contingency = raw_tab,
    raw_contingency_with_tot = raw_tab_tot,
    analysis_contingency = analysis_tab,
    stats_tbl = stats_tbl,
    table_gt = table_gt,
    analysis_gt = analysis_gt,
    stats_gt = stats_gt
  )
}

HPVEM_vs_HPVISH <- agreement_summary_gt(penile_data_abbreviated["HPV_ISH"], penile_data_abbreviated["HPV_EM"], mcnemar_method = "exact")
HPVEM_vs_p16IHC <- agreement_summary_gt(penile_data_abbreviated["p16_IHC"], penile_data_abbreviated["HPV_EM"], mcnemar_method = "exact")

HPV_Biological_vs_HPVISH <- agreement_summary_gt(penile_data_abbreviated["HPV_ISH"], penile_data_abbreviated["HPV_Signature"], mcnemar_method = "exact")
HPV_Biological_vs_p16_IHC <- agreement_summary_gt(penile_data_abbreviated["p16_IHC"], penile_data_abbreviated["HPV_Signature"], mcnemar_method = "exact")

HPVEM_vs_HPV_Biological <- agreement_summary_gt(penile_data_abbreviated["HPV_EM"], penile_data_abbreviated["HPV_Signature"], mcnemar_method = "exact")


#export
gt::gtsave(HPVEM_vs_HPVISH$table_gt, here::here(export_folder_name, "HPV_EM_vs_HPV_ISH.png"), expand = 100)
gt::gtsave(HPVEM_vs_p16IHC$table_gt, here::here(export_folder_name, "HPV_EM_vs_p16_IHC.png"), expand = 100)

gt::gtsave(HPV_Biological_vs_HPVISH$table_gt, here::here(export_folder_name, "HPV_Biological_vs_HPV_ISH.png"), expand = 100)
gt::gtsave(HPV_Biological_vs_p16_IHC$table_gt, here::here(export_folder_name, "HPV_Biological_vs_p16_IHC.png"), expand = 100)

gt::gtsave(HPVEM_vs_HPV_Biological$table_gt, here::here(export_folder_name, "HPV_EM_vs_HPV_Signature.png"), expand = 100)

#Generation of Figure 4 Survival data

#Helper to ensure survival data is in appropriate format

prepare_surv_data <- function(data,
                               time_col,
                               status_col,
                               group_col,
                               unknown_levels = "Unknown",
                               event_status = 1,
                               include_unknown_group = FALSE,
                               time_unit_in = c("auto", "days", "months", "years"),
                               time_unit_out = c("years", "months", "days")) {
  df <- as.data.frame(data)
  needed <- c(time_col, status_col, group_col)
  missing_cols <- setdiff(needed, names(df))
  if (length(missing_cols) > 0) {
    stop("These columns are missing from 'data': ",
         paste(missing_cols, collapse = ", "))
  }
  
  time_raw   <- as.numeric(df[[time_col]])
  status_raw <- df[[status_col]]
  group_raw  <- as.character(df[[group_col]])
  
  # Decide input time unit
  time_unit_in  <- match.arg(time_unit_in)
  time_unit_out <- match.arg(time_unit_out)
  
  if (time_unit_in == "auto") {
    max_t <- max(time_raw, na.rm = TRUE)
    if (max_t > 1000) {
      time_unit_in <- "days"
    } else if (max_t > 40) {
      time_unit_in <- "months"
    } else {
      time_unit_in <- "years"
    }
  }
  
  # Conversion factors using days as reference
  to_days <- c(
    days   = 1,
    months = 30.4375,
    years  = 365.25
  )
  
  factor <- to_days[time_unit_in] / to_days[time_unit_out]
  time   <- time_raw * factor
  
  # Filter Unknown / NA
  if (isTRUE(include_unknown_group)) {
    drop_rows <- is.na(time) | is.na(status_raw) | is.na(group_raw)
  } else {
    drop_rows <- is.na(time) |
      is.na(status_raw) |
      is.na(group_raw) |
      group_raw %in% unknown_levels
  }
  
  keep <- !drop_rows
  if (!any(keep)) stop("No rows left after filtering Unknown/NA.")
  
  time_use   <- time[keep]
  status_use <- status_raw[keep]
  group_use  <- factor(group_raw[keep])
  
  # status == event_status -> event; everything else censored
  event_ind <- as.numeric(status_use == event_status)
  
  df_used <- data.frame(
    time   = time_use,
    status = event_ind,
    group  = group_use
  )
  
  list(
    data_used       = df_used,
    time_unit_in    = time_unit_in,
    time_unit_out   = time_unit_out,
    conversion_fact = as.numeric(factor)
  )
}

#KM Plot (Figure 4B, 4C)
km_km_plot <- function(data,
                       time_col,
                       status_col,
                       group_col,
                       unknown_levels = "Unknown",
                       event_status = 1,
                       include_unknown_group = FALSE,
                       conf_int = FALSE,
                       palette = NULL,
                       time_unit_in = c("auto", "days", "months", "years"),
                       time_unit_out = c("years", "months", "days"),
                       title = NULL) {
  
  prep <- prepare_surv_data(
    data                = data,
    time_col            = time_col,
    status_col          = status_col,
    group_col           = group_col,
    unknown_levels      = unknown_levels,
    event_status        = event_status,
    include_unknown_group = include_unknown_group,
    time_unit_in        = time_unit_in,
    time_unit_out       = time_unit_out
  )
  
  df_used       <- prep$data_used
  time_unit_in  <- prep$time_unit_in
  time_unit_out <- prep$time_unit_out
  
  fit <- survival::survfit(
    survival::Surv(time, status) ~ group,
    data = df_used
  )
  
  if (is.null(title)) {
    title <- paste("Overall survival by", group_col)
  }
  xlab <- paste0("Time (", time_unit_out, ")")
  
  color_factor <- c(
    "Negative" = "#3b4992ff",
    "Positive" = "#ee0000ff"
  )
  
  km_plot <- survminer::ggsurvplot(
    fit,
    data                  = df_used,
    risk.table            = TRUE,
    risk.table.height     = 0.3,
    color = "strata",
    palette = c("#3b4992ff", "#ee0000ff"),
    conf.int              = conf_int,
    pval                  = TRUE,
    pval.method           = TRUE,
    pval.size = 20,
    title                 = title,
    xlab                  = xlab,
    ylab                  = "Survival probability",
    legend.title          = group_col,
    legend.labs           = levels(df_used$group),
    font.main             = c(40, "bold"),
    font.x                = c(40, "plain"),
    font.y                = c(40, "plain"),
    risk.table.fontsize = 20,
    tables.col = "strata"
  ) 
  
  km_plot$plot <- km_plot$plot +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 40))+
    ggplot2::theme(legend.text = ggplot2::element_text(size = 40))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 50))+
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 50))+
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 50, face = "bold"))+
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = 50, face = "bold"))+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 50, face = "bold"))
  
  km_plot$table <- km_plot$table +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 40, face = "bold"))+
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = 40, face = "bold"))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 50))+
    ggplot2::theme(axis.text.y = ggtext::element_markdown(size = 40, face = "bold"))+
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 50, face = "bold"))
  
  km_plot$table@theme$axis.text.y$size <- 40
  #km_plot$table@theme$axis.text.y$face <- "bold"
  
  list(
    fit            = fit,
    plot           = km_plot,
    data_used      = df_used,
    time_unit_in   = time_unit_in,
    time_unit_out  = time_unit_out
  )
}


#Removal of duplicated clinical data (Pt_82)
penile_data_abbreviated_deduplicated <- penile_data_abbreviated_counted[!duplicated(penile_data_abbreviated_counted$PT_Num_Multiomic_Clin),]

survival_analysis_list_deduplicated <- list()
survival_analysis_list_deduplicated[["HPV_ISH_survival"]] <- km_km_plot(penile_data_abbreviated_deduplicated, time_col = "OS.time", status_col = "OS", group_col = "HPV_ISH", include_unknown_group = FALSE)
survival_analysis_list_deduplicated[["p16_IHC_survival"]] <- km_km_plot(penile_data_abbreviated_deduplicated, time_col = "OS.time", status_col = "OS", group_col = "p16_IHC", include_unknown_group = FALSE)
survival_analysis_list_deduplicated[["HPV_EM_survival"]] <- km_km_plot(penile_data_abbreviated_deduplicated, time_col = "OS.time", status_col = "OS", group_col = "HPV_EM")
survival_analysis_list_deduplicated[["HPV_signature_survival"]] <- km_km_plot(penile_data_abbreviated_deduplicated, time_col = "OS.time", status_col = "OS", group_col = "HPV_Signature")

survival_analysis_list_deduplicated[["combined_HPV_signature"]] <- km_km_plot(penile_data_abbreviated_deduplicated, time_col = "OS.time", status_col = "OS", group_col = "combined_HPV_sig")
survival_analysis_list_deduplicated[["combined_HPV_EM"]] <- km_km_plot(penile_data_abbreviated_deduplicated, time_col = "OS.time", status_col = "OS", group_col = "combined_HPV_EM")

survival_analysis_list_deduplicated[["p16_HPV_sig"]] <- km_km_plot(penile_data_abbreviated_deduplicated, time_col = "OS.time", status_col = "OS", group_col = "p16_HPV_sig", include_unknown_group = FALSE)
survival_analysis_list_deduplicated[["p16_HPV_EM"]] <- km_km_plot(penile_data_abbreviated_deduplicated, time_col = "OS.time", status_col = "OS", group_col = "p16_HPV_EM", include_unknown_group = FALSE)
survival_analysis_list_deduplicated[["ISH_HPV_sig"]] <- km_km_plot(penile_data_abbreviated_deduplicated, time_col = "OS.time", status_col = "OS", group_col = "ISH_HPV_sig", include_unknown_group = FALSE)
survival_analysis_list_deduplicated[["ISH_HPV_EM"]] <- km_km_plot(penile_data_abbreviated_deduplicated, time_col = "OS.time", status_col = "OS", group_col = "ISH_HPV_EM", include_unknown_group = FALSE)


for(item in names(survival_analysis_list_deduplicated)){
  svg(here::here(export_folder_name, paste0("surival_analysis_deduplicated_", Sys.Date(), "_", item, ".svg")), width = 16, height = 16)
  print(survival_analysis_list_deduplicated[[item]]$plot)
  dev.off()
}

#Forest Plot (Figure 4A)

library(survival)
library(survminer)

multivar_forest_plot <- function(data,
                                 time_col,
                                 status_col,
                                 var_cols,             # character vector of variable names
                                 event_status = 1,
                                 unknown_levels = "Unknown",
                                 time_unit_in  = c("auto", "days", "months", "years"),
                                 time_unit_out = c("years", "months", "days"),
                                 mode = c("adjusted", "unadjusted", "both"),
                                 title_adjusted   = NULL,
                                 title_unadjusted = NULL) {
  
  mode <- match.arg(mode)
  
  df <- as.data.frame(data)
  
  # ---- basic column checks ----
  needed <- c(time_col, status_col, var_cols)
  missing_cols <- setdiff(needed, names(df))
  if (length(missing_cols) > 0) {
    stop("These columns are missing from 'data': ",
         paste(missing_cols, collapse = ", "))
  }
  
  # ---- time & status ----
  time_raw   <- as.numeric(df[[time_col]])
  status_raw <- df[[status_col]]
  
  # drop rows with missing time/status up front
  keep_ts <- !(is.na(time_raw) | is.na(status_raw))
  if (!any(keep_ts)) stop("No rows with non-missing time and status.")
  df        <- df[keep_ts, , drop = FALSE]
  time_raw   <- time_raw[keep_ts]
  status_raw <- status_raw[keep_ts]
  
  # ---- decide input time unit ----
  time_unit_in  <- match.arg(time_unit_in)
  time_unit_out <- match.arg(time_unit_out)
  
  if (time_unit_in == "auto") {
    max_t <- max(time_raw, na.rm = TRUE)
    if (max_t > 1000) {
      time_unit_in <- "days"
    } else if (max_t > 40) {
      time_unit_in <- "months"
    } else {
      time_unit_in <- "years"
    }
  }
  
  # conversion factors via days
  to_days <- c(
    days   = 1,
    months = 30.4375,
    years  = 365.25
  )
  factor_time <- to_days[time_unit_in] / to_days[time_unit_out]
  time <- time_raw * factor_time
  
  # event indicator: status == event_status
  status <- as.numeric(status_raw == event_status)
  
  df$time   <- time
  df$status <- status
  
  # ---- prepare each variable as factor with Negative baseline ----
  for (v in var_cols) {
    x <- df[[v]]
    x_chr <- as.character(x)
    
    vals <- sort(unique(x_chr[!is.na(x_chr)]))
    
    if (!"Negative" %in% vals) {
      warning("Variable '", v, "' has no 'Negative' level. Actual levels: ",
              paste(vals, collapse = ", "))
    }
    
    # Preferred order: Negative, Positive, Unknown, then anything else
    base_order <- c("Negative", "Positive", unknown_levels)
    lev_order <- unique(c(base_order[base_order %in% vals],
                          setdiff(vals, base_order)))
    
    f <- factor(x_chr, levels = lev_order)
    
    if ("Negative" %in% lev_order) {
      f <- stats::relevel(f, ref = "Negative")
    }
    
    df[[v]] <- f
  }
  
  results <- list(
    data_used      = df,
    time_unit_in   = time_unit_in,
    time_unit_out  = time_unit_out
  )
  
  # -----------------------------------------------------------------
  # Adjusted  (multivariable) Cox model: all vars in one model
  # -----------------------------------------------------------------
  if (mode %in% c("adjusted", "both")) {
    rhs <- paste(var_cols, collapse = " + ")
    cox_formula <- stats::as.formula(paste("survival::Surv(time, status) ~", rhs))
    
    fit_adj <- survival::coxph(cox_formula, data = df)
    
    if (is.null(title_adjusted)) {
      title_adjusted <- "Adjusted hazard ratios for overall survival"
    }
    
    forest_adj <- survminer::ggforest(
      fit_adj,
      data      = df,
      main      = title_adjusted,
      cpositions = c(0.02, 0.22, 0.4),
      fontsize   = 1.2
    )
    
    results$adjusted_fit     <- fit_adj
    results$adjusted_plot    <- forest_adj
    results$adjusted_formula <- cox_formula
  }
  
  # -----------------------------------------------------------------
  # Unadjusted: separate Cox model for each variable, merged forest
  #            using forestmodel::forest_model (like their code)
  # -----------------------------------------------------------------
  if (mode %in% c("unadjusted", "both")) {
    model_list <- list()
    stats_list <- list()
    
    for (v in var_cols) {
      # drop rows with NA for this specific variable
      df_v <- df[!is.na(df[[v]]), , drop = FALSE]
      if (!nrow(df_v)) next
      
      # build Surv(time, status) ~ `v` (backticks handle spaces in names)
      form_v <- stats::as.formula(
        paste("survival::Surv(time, status) ~", sprintf("`%s`", v))
      )
      
      fit_v <- survival::coxph(
        form_v,
        data = df_v,
        x = TRUE
      )
      model_list[[v]] <- fit_v
      
      s <- summary(fit_v)
      coef_tbl <- s$coefficients
      ci_tbl   <- s$conf.int
      
      if (!is.matrix(coef_tbl)) {
        coef_tbl <- t(as.matrix(coef_tbl))
        ci_tbl   <- t(as.matrix(ci_tbl))
      }
      
      term_names <- rownames(coef_tbl)
      hr    <- ci_tbl[, "exp(coef)"]
      lower <- ci_tbl[, "lower .95"]
      upper <- ci_tbl[, "upper .95"]
      pval  <- coef_tbl[, "Pr(>|z|)"]
      
      level <- gsub(paste0("^", v), "", term_names)
      
      df_rows <- data.frame(
        variable = v,
        level    = level,
        term     = term_names,
        HR       = hr,
        lower    = lower,
        upper    = upper,
        p.value  = pval,
        n        = nrow(df_v),
        stringsAsFactors = FALSE
      )
      
      stats_list[[v]] <- df_rows
    }
    
    if (length(model_list) == 0) {
      warning("No unadjusted Cox models could be fit.")
      results$unadjusted_fit_list  <- NULL
      results$unadjusted_stats     <- NULL
      results$unadjusted_plot      <- NULL
    } else {
      hr_df <- do.call(rbind, stats_list)
      
      if (is.null(title_unadjusted)) {
        title_unadjusted <- "Unadjusted hazard ratios (each variable separately)"
      }
      
      panels <- forestmodel::default_forest_panels()
      
      # Locate the existing combined HR/CI column, immediately after the forest
      forest_index <- which(vapply(
        panels,
        function(panel) identical(panel$item, "forest"),
        logical(1)
      ))
      
      estimate_index <- which(
        seq_along(panels) > forest_index &
          vapply(panels, function(panel) !is.null(panel$display), logical(1))
      )[1]
      
      # Replace that column with separate HR and confidence-interval columns
      panels <- append(
        panels[-estimate_index],
        list(
          list(
            width = 0.10, heading = "HR",
            display = ~ ifelse(
              reference, "1.00", sprintf("%.2f", trans(estimate))
            ),
            display_na = NA
          ),
          list(
            width = 0.18, heading = "95% CI",
            display = ~ ifelse(
              reference, "Reference",
              sprintf("%.2f–%.2f", trans(conf.low), trans(conf.high))
            ),
            display_na = NA
          )
        ),
        after = estimate_index - 1
      )
      
      forest_unadj <- forestmodel::forest_model(
        model_list = model_list,
        merge_models = TRUE,
        panels = panels,
        format_options = forestmodel::forest_model_format_options(
          text_size = 6
        )
      )
      
      forest_unadj <- forest_unadj +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(
            size = 20,
            colour = "black"
          )
        )
      
      
      
      results$unadjusted_fit_list  <- model_list
      results$unadjusted_stats     <- hr_df
      results$unadjusted_plot      <- forest_unadj
    }
  }
  
  results
}

HPV_call_forest_deduplicated <- multivar_forest_plot(
  data = penile_data_abbreviated_deduplicated,
  time_col = "OS.time",
  status_col = "OS",
  var_cols = c("HPV_ISH", "p16_IHC", "HPV_EM", "HPV_Signature"),
  mode = "unadjusted"
)

#### Testing for reviewer comment assessing the assumptions of the coxph

# Check proportional hazards for each fitted model
ph_checks_deduplicated <- lapply(
  HPV_call_forest_deduplicated$unadjusted_fit_list,
  survival::cox.zph
)

# Display test results
for (method in names(ph_checks_deduplicated)) {
  cat("\n", method, "\n")
  print(ph_checks_deduplicated[[method]])
}

# Display diagnostic plots
for (method in names(ph_checks_deduplicated)) {
  plot(ph_checks_deduplicated[[method]])
}

do.call(
  rbind,
  lapply(names(HPV_call_forest_deduplicated$unadjusted_fit_list), function(method) {
    fit <- HPV_call_forest_deduplicated$unadjusted_fit_list[[method]]
    
    data.frame(
      Method = method,
      Patients = fit$n,
      Deaths = fit$nevent
    )
  })
)

svg(here::here(export_folder_name, "HPV_Call_forest_deduplicated.svg"), width = 14)
print(HPV_call_forest_deduplicated$unadjusted_plot)
dev.off()


