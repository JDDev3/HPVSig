# app.R

library(shiny)
library(shinyjs)
library(bslib)
library(PCAtools)
library(ggplot2)
library(readr)
library(dplyr)
library(tibble)
library(gt)
library(purrr)
library(pROC)
library(mclust)
library(DT)
library(DESeq2)
library(biomaRt)
library(zip)
library(svglite)


options(shiny.maxRequestSize = 9999*1024^2)

## ------------------------------------------------------------------
## Internal gene filter list
## ------------------------------------------------------------------
FILTER_GENE_LIST <- as.vector(
  readr::read_tsv("App_files/HPV_BioSig_Genelist.tsv", show_col_types = FALSE)[[1]]
)

LOG_PSEUDOCOUNT <- 0.1


## ------------------------------------------------------------------
## Built-in example data
## ------------------------------------------------------------------
set.seed(123)

# Use at least 20 genes from the filter list + 80 random genes
example_sig_genes <- if (length(FILTER_GENE_LIST) >= 20) FILTER_GENE_LIST[1:20] else FILTER_GENE_LIST
example_other_genes <- paste0("GENE_", sprintf("%03d", 1:80))
example_genes <- c(example_sig_genes, example_other_genes)

example_rnaseq <- tibble::tibble(
  gene     = example_genes,
  sample_A = rpois(length(example_genes), lambda = 100),
  sample_B = rpois(length(example_genes), lambda = 120),
  sample_C = rpois(length(example_genes), lambda = 80),
  sample_D = rpois(length(example_genes), lambda = 150)
)

example_clin <- tibble::tibble(
  RNAseq.SampleID = c("sample_A", "sample_B", "sample_C", "sample_D"),
  HPV_status      = c("negative", "positive", "negative", "positive"),
  age             = c(65, 58, 70, 55)
)

## ------------------------------------------------------------------
## HPV helpers
## ------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (is.character(x) && !nzchar(x))) y else x
}

canon_norm_code <- function(code) {
  if (is.null(code) || !nzchar(as.character(code))) return("unknown")
  code <- as.character(code)
  if (identical(code, "upper_quartile")) return("upperquartile")
  code
}

canon_trans_code <- function(code) {
  if (is.null(code) || !nzchar(as.character(code))) return("unknown")
  as.character(code)
}


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
  unknown_mask[is.na(unknown_mask)] <- TRUE
  
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



# apply thresholds learned on (x_train, group_train) to all samples in x_all
threshold_predict_all <- function(x_all, x_train, group_train, res_list) {
  stopifnot(is.factor(group_train), nlevels(group_train) == 2L)
  
  g_levels <- levels(group_train)
  neg <- g_levels[1]
  pos <- g_levels[2]
  
  if (is.null(names(x_all))) {
    stop("x_all must have names = sample IDs for call mapping.")
  }
  
  calls_mat <- list()
  
  for (nm in names(res_list)) {
    res <- res_list[[nm]]
    ev  <- res$evaluation
    if (is.null(ev)) next
    
    thr <- ev$threshold
    if (is.na(thr)) {
      calls_mat[[paste0("call_", nm)]] <- rep(NA_character_, length(x_all))
      next
    }
    
    mean_pos <- mean(x_train[group_train == pos], na.rm = TRUE)
    mean_neg <- mean(x_train[group_train == neg], na.rm = TRUE)
    
    if (mean_pos >= mean_neg) {
      pred_pos_all <- x_all >= thr
    } else {
      pred_pos_all <- x_all <= thr
    }
    
    call_vec <- ifelse(pred_pos_all, pos, neg)
    calls_mat[[paste0("call_", nm)]] <- call_vec
  }
  
  if (length(calls_mat) == 0) {
    return(tibble::tibble(sample = names(x_all)))
  }
  
  calls_df <- as.data.frame(calls_mat, stringsAsFactors = FALSE)
  calls_df$sample <- names(x_all)
  calls_df
}

## ------------------------------------------------------------------
## PCA theme and plotting helpers
## ------------------------------------------------------------------
pca_theme <- function() {
  theme_classic() +
    theme(
      axis.title   = element_text(size = 34),
      axis.text    = element_text(size = 30),
      legend.title = element_text(size = 30),
      legend.text  = element_text(size = 28),
      plot.title   = element_text(size = 34, face = "bold")
    )
}

# Generic PCA plotting with count (%) in legend labels
# Generic PCA plotting with count (%) in legend labels + optional point labels
# Generic PCA plotting with count (%) in legend labels + optional point labels
plot_pca_colored <- function(df, col_var, label, label_var = NULL) {
  df2 <- df
  
  # Ensure the color var is discrete and includes NA as a level for counting
  df2[[col_var]] <- as.factor(df2[[col_var]])
  col_for_tab <- df2[[col_var]]
  tab <- table(col_for_tab, useNA = "ifany")
  total <- sum(tab)
  
  levs <- names(tab)
  lab_map <- paste0(
    levs,
    " (n=", as.integer(tab), ", ",
    round(100 * as.numeric(tab) / total, 1), "%)"
  )
  names(lab_map) <- levs
  
  p <- ggplot(df2, aes(x = PC1, y = PC2, color = .data[[col_var]])) +
    geom_point(size = 4, alpha = 0.9) +
    labs(
      title = paste("PCA (colored by", label, ")"),
      x     = "PC1",
      y     = "PC2",
      color = label
    ) +
    scale_color_discrete(
      labels = function(breaks) {
        out <- unname(lab_map[breaks])
        out[is.na(out)] <- breaks[is.na(out)]
        out
      }
    ) +
    pca_theme()
  
  # Optional labels (repelled)
  if (!is.null(label_var) && nzchar(label_var) && label_var %in% names(df2)) {
    
    lab_chr <- as.character(df2[[label_var]])
    lab_chr[is.na(lab_chr)] <- ""
    lab_chr <- trimws(lab_chr)
    
    df_lab <- df2
    df_lab$pca_label <- lab_chr
    df_lab <- df_lab[df_lab$pca_label != "", , drop = FALSE]
    
    if (nrow(df_lab) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = df_lab,
        aes(label = pca_label),
        size = 3.5,
        show.legend = FALSE,
        max.overlaps = Inf,      # IMPORTANT: don't drop labels
        box.padding = 0.35,
        point.padding = 0.25,
        min.segment.length = 0
      )
    }
  }
  
  p
}




## ------------------------------------------------------------------
## UI helpers for preprocessing cards
## ------------------------------------------------------------------
pill_tag <- function(text,
                     color = c("secondary","success","danger","warning","info","primary","dark","light"),
                     tooltip = NULL) {
  color <- match.arg(color)
  tags$span(
    class = paste0("pp-pill pp-pill-", color),
    title = if (!is.null(tooltip) && nzchar(tooltip)) tooltip else NULL,
    text
  )
}


pp_attr_pill <- function(text,
                         attr = c("norm", "trans", "scale", "compat"),
                         state = c("neutral", "yes", "no", "unknown"),
                         tooltip = NULL,
                         tier = c("default", "preferred", "acceptable"),
                         method_key = NULL) {
  attr  <- match.arg(attr)
  state <- match.arg(state)
  tier  <- match.arg(tier)
  
  mk <- NULL
  if (!is.null(method_key) && nzchar(method_key)) {
    mk <- paste0("pp-pill--m-", gsub("[^a-zA-Z0-9]+", "-", tolower(method_key)))
  }
  
  tags$span(
    class = paste(
      "pp-pill",
      paste0("pp-pill--", attr),
      paste0("pp-pill--", state),
      mk,
      if (tier != "default") paste0("pp-pill--tier-", tier) else NULL
    ),
    title = if (!is.null(tooltip) && nzchar(tooltip)) tooltip else NULL,
    text
  )
}

pp_pick_btn <- function(kind, value, pill,
                        selected = FALSE,
                        disabled = FALSE,
                        dim = FALSE,
                        tooltip = NULL) {
  tags$button(
    type = "button",
    class = paste(
      "pp-pick",
      if (selected) "is-selected" else NULL,
      if (disabled) "is-disabled" else NULL,
      if (dim) "is-dim" else NULL
    ),
    `data-kind`  = kind,
    `data-value` = value,
    title = if (!is.null(tooltip) && nzchar(tooltip)) tooltip else NULL,
    disabled = if (isTRUE(disabled)) "disabled" else NULL,
    pill
  )
}


pp_attr_item <- function(label, pill) {
  tags$div(
    class = "pp-attr-line",
    tags$div(class = "pp-attr-label", label),
    tags$div(class = "pp-attr-main", pill)
  )
}


row_item <- function(label, pill, src_chip = NULL, note = NULL) {
  tags$div(
    class = "pp-row",
    tags$div(class = "pp-label", label),
    tags$div(
      class = "pp-value",
      tags$div(
        class = "pp-inline",
        pill,
        if (!is.null(src_chip)) src_chip
      ),
      if (!is.null(note) && nzchar(note)) tags$div(class = "pp-note", note)
    )
  )
}

## ------------------------------------------------------------------
## UI
## ------------------------------------------------------------------
ui <- fluidPage(
  theme = bs_theme(version = 5, bg = "#FFFFFF", fg = "#212529"),
  useShinyjs(),
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "hpv_toggle.css"),
    tags$script(src = "hpv_toggle.js")
  ),
  
  titlePanel("Two-Group Threshold Explorer"),
  
  # Control card
  fluidRow(
    column(
      12,
      card(
        card_header("Report controls"),
        card_body(
          div(
            style = "display:flex; justify-content:flex-start; align-items:center; flex-wrap:wrap; gap:10px;",
            div(
              class = "btn-group",
              actionButton("edit_setup",    "Edit setup",       class = "btn-outline-secondary"),
              actionButton("show_meta_btn", "Preview metadata", class = "btn-outline-info")
            )
          )
        )
      )
    )
  ),
  
  # (Preprocessing overview card removed from main report area; preprocessing lives in the wizard
  #  and a textual summary is written into the ZIP export.)
  
  # PCA + Threshold plots side by side
  fluidRow(
    column(
      6,
      card(
        card_header("PCA plot"),
        card_body(
          uiOutput("pca_color_ui"),
          plotOutput("pca_plot", height = "400px")
        )
      )
    ),
    column(
      6,
      card(
        card_header("Threshold plot"),
        card_body(
          plotOutput("threshold_plot", height = "400px")
        )
      )
    )
  ),
  
  fluidRow(
    column(
      6,
      card(
        card_header("PCA Loadings Plot"),
        card_body(
          plotOutput("pca_loadings_plot", height = "400px")
        )
      )
    ),
    column(
      6,
      card(
        card_header("ROC diagnostics (Youden)"),
        card_body(
          shiny::fluidRow(
            shiny::column(
              4,
              uiOutput("youden_diag_ui")
            ),
            shiny::column(
              8,
              plotOutput("youden_roc_plot", height = "400px")
            )
          )
        )
      )
    )
  ),
  
  # Row 3: performance table full width
  fluidRow(
    column(
      12,
      card(
        card_header("Threshold performance table"),
        card_body(
          gt_output("threshold_gt")
        )
      )
    )
  ),
  
  # Export button
  fluidRow(
    column(
      12,
      card(
        card_body(
          div(
            style = "text-align:left;",
            downloadButton("download_zip", "Export report (ZIP)", class = "btn-success")
          )
        )
      )
    )
  )
)

## ------------------------------------------------------------------
## Server
## ------------------------------------------------------------------
server <- function(input, output, session) {
  
  wizard_step            <- reactiveVal(1)
  hpv_flip_state         <- reactiveVal(FALSE)
  example_defaults_applied <- reactiveVal(FALSE)
  
  
  pp_input <- reactiveValues(
    norm_status    = "raw",  # "raw" | "normalized"
    norm_method    = "unknown",     # only meaningful if norm_status == "normalized"
    trans_applied  = "no",          # "yes" | "no"
    trans_method   = "log2p0.1",      # only meaningful if trans_applied == "yes"
    scale_applied  = "no"           # "yes" | "no"
  )
  
  
  pp_plan <- reactiveValues(
    # what the app should APPLY (only when allowed by pp_input)
    norm  = "deseq2_sf",          # plan normalization if pp_input$norm_status=="raw"
    trans = "log2p0.1",             # plan transformation if pp_input$trans_applied=="none"
    scale = "no"                  # "no" | "yes"
  )
  
  get_pp_state <- shiny::reactive({
    list(
      norm_status    = pp_input$norm_status,
      norm_method    = pp_input$norm_method,
      log_applied    = pp_input$trans_applied,   # "yes"/"no"
      log_method     = pp_input$trans_method,    # "log2p0.1"/"vst"/"rlog"/"unknown" (only meaningful if log_applied=="yes")
      scaled_applied = pp_input$scale_applied,   # "yes"/"no"
      
      raw_norm_plan  = pp_plan$norm,             # plan normalization (if input is raw)
      plan_trans     = pp_plan$trans,            # plan transform (if input not transformed)
      plan_scale     = pp_plan$scale             # "yes"/"no"
    )
  })
  
  
  bucket_step <- reactiveVal(NULL)
  
  matrix_features <- reactive({
    mat <- tryCatch(rna_matrix_full(), error = function(e) NULL)
    if (is.null(mat)) return(NULL)
    if (ncol(mat) < 2 || nrow(mat) < 2) return(NULL)
    
    # sample values WITHOUT copying full matrix to a giant numeric vector
    n_take <- min(length(mat), 50000)
    x <- mat[sample.int(length(mat), n_take)]
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NULL)
    
    has_negative <- any(x < 0)
    frac_integer <- mean(abs(x - round(x)) < 1e-6)
    
    q01  <- as.numeric(stats::quantile(x, 0.01, na.rm = TRUE))
    q99  <- as.numeric(stats::quantile(x, 0.99, na.rm = TRUE))
    maxv <- max(x, na.rm = TRUE)
    
    likely_loglike <- (q99 <= 50) && (q01 >= -20) && (frac_integer < 0.8)
    
    likely_scaled <- FALSE
    if (nrow(mat) >= 30 && ncol(mat) >= 3) {
      rr <- mat[sample(seq_len(nrow(mat)), min(200, nrow(mat))), , drop = FALSE]
      row_means <- rowMeans(rr, na.rm = TRUE)
      row_sds   <- apply(rr, 1, sd, na.rm = TRUE)
      likely_scaled <- has_negative &&
        mean(abs(row_means), na.rm = TRUE) < 0.15 &&
        mean(abs(row_sds - 1), na.rm = TRUE) < 0.25
    }
    
    cs <- colSums(mat, na.rm = TRUE)
    med_cs <- median(cs, na.rm = TRUE)
    
    tpm_like <- !has_negative && frac_integer < 0.9 && med_cs > 5e5 && med_cs < 2e6
    
    list(
      has_negative   = has_negative,
      frac_integer   = frac_integer,
      q01            = q01,
      q99            = q99,
      maxv           = maxv,
      med_colsum     = med_cs,
      likely_loglike = likely_loglike,
      likely_scaled  = likely_scaled,
      tpm_like       = tpm_like
    )
  })
  
  harmonize_gene_ids_once <- function(mat, ref_genes, dds = NULL,
                                      min_match = 20,
                                      assume_counts = TRUE) {
    mat <- as.matrix(mat)
    ids0 <- rownames(mat)
    if (is.null(ids0)) return(mat)
    
    # normalize strings for matching
    norm <- function(x) {
      x <- toupper(trimws(as.character(x)))
      x[is.na(x)] <- ""
      x
    }
    
    ref_up <- unique(norm(ref_genes))
    
    overlap_score <- function(ids) {
      z <- unique(norm(ids))
      z <- z[nzchar(z)]
      length(intersect(z, ref_up))
    }
    
    collapse_dups <- function(m) {
      if (!anyDuplicated(rownames(m))) {
        return(m)
        }
      if (isTRUE(assume_counts)) {
        rowsum(m, group = rownames(m), reorder = FALSE)  # sum for counts
      } else {
        sm <- rowsum(m, group = rownames(m), reorder = FALSE)
        n  <- as.numeric(table(rownames(m)))
        n  <- n[rownames(sm)]
        sweep(sm, 1, n, "/")                             # mean otherwise
      }
    }
    
    # apply a candidate ID vector (same length as rows of mat), replacing where non-empty
    apply_ids <- function(candidate_ids) {
      candidate_ids <- trimws(as.character(candidate_ids))
      new_ids <- ids0
      ok <- !is.na(candidate_ids) & nzchar(candidate_ids)
      new_ids[ok] <- candidate_ids[ok]
      rownames(mat) <- new_ids
      mat <- collapse_dups(mat)
    }
    
    base <- overlap_score(ids0)
    if (base >= min_match) return(collapse_dups(mat))
    
    # ------------------------------------------------------------
    # 1) RowData scan: stop early when a column hits min_match
    # ------------------------------------------------------------
    if (!is.null(dds)) {
      rd <- as.data.frame(SummarizedExperiment::rowData(dds))
      if (nrow(rd) > 0 && ncol(rd) > 0) {
        
        # Align rowData rows to mat rows using rownames (strip Ensembl versions both sides)
        rd_rn <- rownames(rd)
        ids_key <- sub("\\.\\d+$", "", trimws(as.character(ids0)))
        
        idx <- NULL
        if (!is.null(rd_rn)) {
          rd_key <- sub("\\.\\d+$", "", trimws(as.character(rd_rn)))
          idx <- match(ids_key, rd_key)
        }
        
        preferred <- c(
          "symbol","SYMBOL","hgnc_symbol","HGNC",
          "gene_name","external_gene_name","gene","Gene"
        )
        
        # Preferred first, then other character/factor columns
        cand <- c(intersect(preferred, colnames(rd)),
                  setdiff(colnames(rd), preferred))
        
        for (cc in cand) {
          v <- rd[[cc]]
          if (!(is.character(v) || is.factor(v))) next
          v <- as.character(v)
          
          v_aligned <- NULL
          if (!is.null(idx) && length(idx) == length(ids0)) {
            v_aligned <- v[idx]               # aligned by match()
          } else if (length(v) == length(ids0)) {
            v_aligned <- v                    # fallback: assume already aligned
          } else {
            next
          }
          
          sc <- overlap_score(v_aligned)
          if (sc >= min_match && sc > base) {
            return(apply_ids(v_aligned))
          }
        }
      }
    }
    
    # ------------------------------------------------------------
    # 2) BioMart: ENSEMBL -> HGNC symbol (only if Ensembl-like)
    # ------------------------------------------------------------
    ids_nover <- sub("\\.\\d+$", "", trimws(as.character(ids0)))
    ens_like  <- grepl("^ENSG\\d+$", ids_nover)
    
    if (any(ens_like)) {
      mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
      bm <- biomaRt::getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol"),
        filters    = "ensembl_gene_id",
        values     = unique(ids_nover[ens_like]),
        mart       = mart
      )
      
      bm <- bm[!is.na(bm$hgnc_symbol) & nzchar(bm$hgnc_symbol), , drop = FALSE]
      bm <- bm[!duplicated(bm$ensembl_gene_id), , drop = FALSE]
      map <- setNames(bm$hgnc_symbol, bm$ensembl_gene_id)
      
      mapped <- unname(map[ids_nover])  # NA where not mapped
      # only take this if it improves overlap
      sc <- overlap_score(ifelse(!is.na(mapped) & nzchar(mapped), mapped, ids0))
      if (sc > base) {
        return(apply_ids(mapped))
      }
    }
    # Nothing improved; just collapse duplicates if any
    collapse_dups(mat)
  }
  
  apply_preproc_state <- function(norm_status = NULL,
                                  norm_method = NULL,
                                  raw_norm_plan = NULL,
                                  log_applied = NULL,
                                  log_method = NULL,
                                  plan_trans = NULL,
                                  scaled_applied = NULL,
                                  plan_scale = NULL) {
    
    if (!is.null(norm_status))    pp_input$norm_status <- norm_status
    if (!is.null(norm_method))    pp_input$norm_method <- norm_method
    if (!is.null(log_applied))    pp_input$trans_applied <- log_applied
    if (!is.null(log_method))     pp_input$trans_method <- log_method
    if (!is.null(scaled_applied)) pp_input$scale_applied <- scaled_applied
    if (!is.null(raw_norm_plan))  pp_plan$norm <- raw_norm_plan
    if (!is.null(plan_trans))     pp_plan$trans <- plan_trans
    if (!is.null(plan_scale))     pp_plan$scale <- plan_scale
    
    # Hard constraints:
    if (identical(pp_input$scale_applied, "yes")) {
      pp_plan$scale <- "no"
      pp_plan$trans <- "none"
    }
    if (identical(pp_input$trans_applied, "yes")) {
      pp_plan$trans <- "none"
    }
  }
  
  rds_obj <- reactiveVal(NULL)
  
  rds_selected_assay <- reactive({
    obj <- rds_obj()
    if (is.null(obj) || !inherits(obj, "SummarizedExperiment")) return(NULL)
    
    an <- SummarizedExperiment::assayNames(obj)
    if (length(an) == 0) return(NULL)
    
    sel <- input$rna_assay
    if (!is.null(sel) && sel %in% an) sel else an[[1]]
  })
  
  
  observeEvent(input$rna_file, {
    req(input$rna_file)
    
    ext <- tolower(tools::file_ext(input$rna_file$name))
    
    if (ext != "rds") {
      rds_obj(NULL)
      return()
    }
    
    obj <- tryCatch(readRDS(input$rna_file$datapath), error = function(e) NULL)
    if (is.null(obj)) return()
    
    rds_obj(obj)
    
    # Only keep your existing DESeq defaulting behavior (raw) if you want it.
    # SummarizedExperiment: do nothing (rely on user inputs, as requested).
    if (inherits(obj, "DESeqDataSet")) {
      apply_preproc_state(
        norm_status    = "raw",
        norm_method    = "unknown",
        log_applied    = "no",
        log_method     = "log2p0.1",
        scaled_applied = "no",
        raw_norm_plan  = "deseq2_sf",
        plan_trans     = "log2p0.1",
        plan_scale     = "no"
      )
    }
  }, ignoreInit = TRUE)
  
  output$rna_assay_ui <- renderUI({
    req(input$rna_file)
    ext <- tolower(tools::file_ext(input$rna_file$name))
    if (ext != "rds") return(NULL)
    
    obj <- rds_obj()
    if (!inherits(obj, "SummarizedExperiment")) return(NULL)
    
    an <- SummarizedExperiment::assayNames(obj)
    validate(need(length(an) > 0, "This RDS object has no assays."))
    
    if (length(an) == 1) return(NULL)  # 1 assay -> use it, no dropdown
    
    selectInput(
      "rna_assay",
      "Assay to use:",
      choices  = an,
      selected = an[[1]]  # default to first
    )
  })
  
  
  observeEvent(input$auto_preproc, {
    feat <- matrix_features()
    req(!is.null(feat))
    
    st <- get_pp_state()
    
    # Locks based ONLY on what the user said the INPUT already is
    lock_norm  <- identical(st$scaled_applied, "yes") || !identical(st$norm_status, "raw")
    lock_trans <- identical(st$scaled_applied, "yes") || identical(st$log_applied, "yes")
    lock_scale <- identical(st$scaled_applied, "yes")
    
    # Hard constraints for what the app can actually apply
    hard_no_counts <- isTRUE(feat$has_negative)   # count-based normalization needs non-negative
    hard_no_log    <- isTRUE(feat$has_negative)   # your preprocess_for_pca() blocks log2p0.1 on negatives
    
    # ----------------------------
    # 1) Plan: scaling (default no)
    # ----------------------------
    if (!lock_scale) {
      apply_preproc_state(plan_scale = "no")
    }
    
    # ----------------------------
    # 2) Plan: normalization (only meaningful if user said raw)
    # ----------------------------
    if (!lock_norm) {
      # Pick a sensible default; avoid changing if already set
      if (is.null(st$raw_norm_plan) || !nzchar(st$raw_norm_plan)) {
        if (!hard_no_counts) {
          apply_preproc_state(raw_norm_plan = "deseq2_sf")
        }
      } else {
        # If current choice is invalid given constraints, fall back
        if (hard_no_counts) {
          # Can't apply any count-based normalization safely; leave as-is and warn
          showNotification(
            "Auto-select: count-based normalization cannot be applied because negative values are present. Plan normalization left unchanged.",
            type = "warning",
            duration = 8
          )
        }
      }
    }
    
    # ----------------------------
    # 3) Plan: transformation
    # ----------------------------
    if (!lock_trans) {
      # If negatives present, we cannot safely apply log2p0.1 in-app per your preprocess_for_pca()
      if (hard_no_log) {
        apply_preproc_state(plan_trans = "none")
      } else {
        # If matrix already looks log-like, avoid double-logging by default
        if (isTRUE(feat$likely_loglike)) {
          apply_preproc_state(plan_trans = "none")
          
          # IMPORTANT: don't change input, but warn user their INPUT description might be off
          showNotification(
            sprintf(
              "Auto-select: matrix looks already log-scale (q99=%.1f). Plan sets Transformation = None to avoid double-logging. If you know it was already logged, update Input → Transformation accordingly.",
              feat$q99 %||% NA_real_
            ),
            type = "warning",
            duration = 10
          )
        } else {
          apply_preproc_state(plan_trans = "log2p0.1")
        }
      }
    }
    
    showNotification(
      "Auto-selected a processing plan (plan only). Input descriptions were not changed.",
      type = "message",
      duration = 6
    )
  })
  
  
  # Initial welcome
  observeEvent(TRUE, {
    showModal(modalDialog(
      title = "Welcome",
      size  = "l",
      easyClose = FALSE,
      footer = NULL,
      uiOutput("wizard_body")
    ))
  }, once = TRUE)
  
  # Re-open setup without resetting state
  observeEvent(input$edit_setup, {
    showModal(modalDialog(
      title = "Setup",
      size  = "l",
      easyClose = TRUE,
      footer = NULL,
      uiOutput("wizard_body")
    ))
  })
  
  # If RNA example is used, default to example clinical when entering clinical step (only once)
  observe({
    if (wizard_step() == 3 &&
        isTRUE(input$use_example_rna) &&
        !example_defaults_applied()) {
      updateRadioButtons(session, "include_clin", selected = "Yes")
      updateCheckboxInput(session, "use_example_clin", value = TRUE)
      example_defaults_applied(TRUE)
    }
  })
  
  observeEvent(input$use_example_rna, {
    if (!isTRUE(input$use_example_rna)) {
      example_defaults_applied(FALSE)
    }
  })
  
  rna_raw <- reactive({
    if (isTRUE(input$use_example_rna)) {
      return(example_rnaseq)
    }
    
    if (is.null(input$rna_file)) return(NULL)
    
    ext <- tolower(tools::file_ext(input$rna_file$name))
    
    if (ext == "rds") {
      obj <- rds_obj()
      validate(need(!is.null(obj), "Could not read the .rds file."))
      
      validate(need(
        inherits(obj, "DESeqDataSet") || inherits(obj, "SummarizedExperiment"),
        "The .rds must be a DESeqDataSet or SummarizedExperiment."
      ))
      
      return(obj)
    }
    
    
    # default TSV/TXT
    readr::read_delim(
      input$rna_file$datapath,
      delim = "\t",
      show_col_types = FALSE
    )
  })
  
  rna_matrix_full <- reactive({
    x <- rna_raw()
    req(!is.null(x))
    
    assume_counts <- identical(pp_input$norm_status, "raw")
    # A) DESeqDataSet
    if (inherits(x, "DESeqDataSet")) {
      an <- SummarizedExperiment::assayNames(x)
      validate(need(length(an) > 0, "DESeqDataSet has no assays."))
      
      first_assay <- an[[1]]
      assay_name  <- rds_selected_assay() %||% first_assay
      
      mat <- as.matrix(SummarizedExperiment::assay(x, assay_name))
      
      mat <- harmonize_gene_ids_once(
        mat          = mat,
        ref_genes     = FILTER_GENE_LIST,
        dds           = x,
        min_match     = 20,
        assume_counts = assume_counts
      )
      return(mat)
    }
    
    
    # B) SummarizedExperiment
    if (inherits(x, "SummarizedExperiment")) {
      an <- SummarizedExperiment::assayNames(x)
      validate(need(length(an) > 0, "SummarizedExperiment has no assays."))
      
      assay_name <- rds_selected_assay() %||% an[[1]]
      mat <- as.matrix(SummarizedExperiment::assay(x, assay_name))
      
      # Your rule: always rely on user inputs
      assume_counts <- identical(pp_input$norm_status, "raw")
      
      mat <- harmonize_gene_ids_once(
        mat          = mat,
        ref_genes     = FILTER_GENE_LIST,
        dds           = x,
        min_match     = 20,
        assume_counts = assume_counts
      )
      return(mat)
    }
    
    
    # C) TSV/TXT (data.frame / tibble)
    df <- x
    validate(need(ncol(df) >= 2, "RNA-seq file must have at least 2 columns (gene + one sample)."))
    
    gene_col <- df[[1]]
    mat_df   <- df[, -1, drop = FALSE]
    
    num_cols <- which(vapply(mat_df, is.numeric, logical(1)))
    validate(need(length(num_cols) >= 1, "Need at least one numeric sample column."))
    mat_df <- mat_df[, num_cols, drop = FALSE]
    
    mat <- as.matrix(mat_df)
    rownames(mat) <- trimws(as.character(gene_col))
    
    mat <- harmonize_gene_ids_once(
      mat       = mat,
      ref_genes = FILTER_GENE_LIST,
      dds       = NULL,              # no rowData for TSV
      min_match = 20,
      assume_counts = assume_counts
    )
    
    mat
  })
  
  
  
  
  
  rna_matrix_sig_unprocessed <- reactive({
    mat <- rna_matrix_full()
    if (!is.null(FILTER_GENE_LIST) && length(FILTER_GENE_LIST) > 0) {
      keep <- rownames(mat) %in% FILTER_GENE_LIST
      mat <- mat[keep, , drop = FALSE]
    }
    mat
  })
  
  
  rna_ok <- reactive({
    x <- rna_raw()
    if (is.null(x)) {
      return(list(ok = FALSE, messages = "Please upload an RNA-seq file or use the built-in example."))
    }
    
    msgs <- character()
    
    mat <- tryCatch(rna_matrix_sig_unprocessed(), error = function(e) e)
    if (inherits(mat, "error")) {
      msgs <- c(msgs, paste("Error processing matrix:", mat$message))
      return(list(ok = FALSE, messages = msgs))
    }
    
    if (nrow(mat) < 3) msgs <- c(msgs, "Need at least 3 filtered genes for PCA.")
    if (ncol(mat) < 2) msgs <- c(msgs, "Need at least 2 samples for PCA.")
    
    # TSV-specific light check (only if it's a data.frame/tibble)
    if (inherits(x, c("data.frame", "tbl_df"))) {
      if (is.numeric(x[[1]])) {
        msgs <- c(msgs, "First column appears numeric; it should contain gene IDs.")
      }
    }
    
    if (length(msgs) == 0) {
      list(
        ok       = TRUE,
        messages = c(
          "All checks passed.",
          paste("Filtered genes:", nrow(mat)),
          paste("Samples:", ncol(mat))
        )
      )
    } else {
      list(ok = FALSE, messages = msgs)
    }
  })
  
  
  output$rna_checks_text <- renderText({
    ck <- rna_ok()
    paste("-", ck$messages, collapse = "\n")
  })
  
  output$rna_example <- renderTable({
    head(example_rnaseq, 10)
  }, striped = TRUE, bordered = TRUE, spacing = "xs")
  
  observe({
    ck <- rna_ok()
    if (ck$ok) shinyjs::enable("step1_next") else shinyjs::disable("step1_next")
  })
  
  
  
  
  ## -----------------------------
  ## Gene overlap summary (reference list vs matched)
  ## -----------------------------
  gene_overlap <- reactive({
    mat <- tryCatch(rna_matrix_full(), error = function(e) NULL)
    if (is.null(mat)) return(NULL)
    
    genes_input <- rownames(mat)
    genes_input <- genes_input[!is.na(genes_input)]
    genes_input <- unique(as.character(genes_input))
    
    matched <- intersect(genes_input, FILTER_GENE_LIST)
    
    list(
      n_reference     = length(unique(FILTER_GENE_LIST)),
      n_input         = length(genes_input),
      n_matched       = length(matched),
      pct_matched_ref = if (length(unique(FILTER_GENE_LIST)) > 0) {
        100 * length(matched) / length(unique(FILTER_GENE_LIST))
      } else {
        NA_real_
      }
    )
  })
  
  
  # effective pipeline for PCA (resulting preprocessing before PCAtools::pca)
  preproc_effective <- reactive({
    # Represents the state of the matrix *entering PCA* after any app-applied transforms.
    
    norm_status    <- pp_input$norm_status
    norm_method    <- pp_input$norm_method
    log_applied    <- pp_input$trans_applied
    log_method     <- pp_input$trans_method
    scaled_applied <- pp_input$scale_applied
    raw_plan       <- pp_plan$norm
    plan_scale     <- pp_plan$scale
    
    # Canonical codes for matching against reference
    norm_code <- if (scaled_applied == "yes") {
      "as_is"
    } else if (norm_status == "raw") {
      raw_plan
    } else {  # normalized
      norm_method
    }
    
    log_code <- if (scaled_applied == "yes") {
      "as_is"
    } else if (log_applied == "yes") {
      log_method
    } else {
      if (!is.null(pp_plan$trans) && nzchar(pp_plan$trans)) {
        pp_plan$trans
      } else {
        "log2p0.1"
      }
    }
    
    
    scale_code <- if (scaled_applied == "yes") {
      "as_is"
    } else if (plan_scale == "yes") {
      "zscore"
    } else {
      "none"
    }
    
    # Human-readable labels for UI
    norm_label <- function(code) {
      switch(
        code,
        "none"           = "None",
        "deseq2_sf"      = "DESeq2",
        "tpm"            = "TPM",
        "fpkm"           = "FPKM/RPKM",
        "cpm"            = "CPM",
        "tmm"            = "TMM→CPM",
        "upperquartile"  = "Upper-quartile",
        "quantile"       = "Quantile",
        "as_is"          = "As provided",
        "unknown"        = "Unknown",
        code
      )
    }
    
    log_label <- function(code) {
      switch(
        code,
        "none"   = "No",
        "log2p0.1" = "Yes",
        "vst"    = "Yes",
        "rlog"   = "Yes",
        "unknown"= "Yes",
        "as_is"  = "As provided",
        code
      )
    }
    
    scale_label <- function(code) {
      switch(
        code,
        "none"  = "No",
        "zscore"= "Yes",
        "as_is" = "As provided",
        code
      )
    }
    
    list(
      norm_code     = norm_code,
      log_code      = log_code,
      scale_code    = scale_code,
      normalization = norm_label(norm_code),
      log_transform = log_label(log_code),
      scaled        = scale_label(scale_code)
    )
  })
  
  preproc_compatibility <- reactive({
    # Reference PCA expectations (hard-coded)
    ref_norm_code  <- "deseq2_sf"
    ref_log_code   <- "log2p0.1"
    ref_scale_code <- "none"
    
    eff <- preproc_effective()
    
    if (pp_input$scale_applied == "yes") {
      return(list(
        status         = "Not compatible",
        compatible     = FALSE,
        color          = "danger",
        reason         = "Input is already scaled/z-scored; this cannot be made comparable to an unscaled reference without reversing the scaling.",
        recommendation = "Provide unscaled gene expression (counts/TPM/CPM/etc.), or re-export without z-scoring, then rerun."
      ))
    }
    
    issues <- character()
    
    if (!identical(eff$norm_code, ref_norm_code)) {
      issues <- c(issues, "Normalization does not match the reference (reference uses DESeq2 size-factor normalized counts).")
    }
    if (!identical(eff$log_code, ref_log_code)) {
      issues <- c(issues, "Transformation does not match the reference (reference uses log2(x+0.1)).")
    }
    if (!identical(eff$scale_code, ref_scale_code)) {
      issues <- c(issues, "Scaling does not match the reference (reference is NOT gene-wise z-scored).")
    }
    
    if (length(issues) == 0) {
      return(list(
        status         = "Compatible",
        compatible     = TRUE,
        color          = "success",
        reason         = "Planned preprocessing matches the reference PCA preprocessing.",
        recommendation = "No action needed."
      ))
    }
    
    rec <- character()
    if (!identical(eff$norm_code, ref_norm_code)) {
      rec <- c(rec, "Use DESeq2 size-factor normalized counts (or provide raw counts and choose 'DESeq2 size factors' as the plan).")
    }
    if (!identical(eff$log_code, ref_log_code)) {
      rec <- c(rec, "Use log2(x+0.1) (avoid VST/rlog if you need strict comparability to the reference model).")
    }
    if (!identical(eff$scale_code, ref_scale_code)) {
      rec <- c(rec, "Turn OFF gene-wise scaling to match the reference model.")
    }
    
    list(
      status         = "Not compatible",
      compatible     = FALSE,
      color          = "warning",
      reason         = paste(issues, collapse = " "),
      recommendation = paste(unique(rec), collapse = " ")
    )
  })
  
  bucket_shell <- function(id, title, body_ui) {
    tags$div(
      id = id,
      class = "pp-bucket",
      tags$div(
        class = "pp-bucket-header",
        tags$div(
          tags$div(class = "pp-bucket-title", title),
          tags$div(class = "pp-bucket-subtitle", "Selections update immediately.")
        ),
        tags$button(type = "button", class = "pp-bucket-close btn-close", `aria-label` = "Close")
      ),
      body_ui
    )
  }
  
  
  ## -----------------------------
  ## Cards: Reference / Input / Acceptable
  ## -----------------------------
  output$preproc_cards_ui <- renderUI({
    
    # Reference pipeline (hard-coded to the stored PCA model)
    ref_norm_code  <- "deseq2_sf"
    ref_trans_code <- "log2p0.1"
    ref_scale_code <- "none"
    
    # Tooltips (centralized so all cards stay consistent)
    tt_norm <- list(
      deseq2_sf     = "DESeq2 size factors: dds <- DESeq2::DESeqDataSetFromMatrix(countData=counts_mat, ...) |> DESeq2::estimateSizeFactors(); norm_counts <- DESeq2::counts(dds, normalized=TRUE)",
      tmm           = "edgeR TMM → CPM: dge <- edgeR::DGEList(counts=counts_mat) |> edgeR::calcNormFactors(method='TMM'); x <- edgeR::cpm(dge, log=FALSE)",
      cpm           = "CPM: dge <- edgeR::DGEList(counts=counts_mat); x <- edgeR::cpm(dge, log=FALSE)",
      tpm           = "External (Salmon/kallisto/RSEM): TPM matrix; often imported via tximport::tximport().",
      fpkm          = "External (StringTie/Cufflinks/RSEM): FPKM/RPKM matrix.",
      upperquartile = "Upper-quartile: dge <- edgeR::DGEList(counts=counts_mat) |> edgeR::calcNormFactors(method='upperquartile'); x <- edgeR::cpm(dge, log=FALSE)",
      quantile      = "Quantile: preprocessCore::normalize.quantiles(x) or limma::normalizeBetweenArrays(x, method='quantile')",
      none          = "No upstream normalization (raw counts).",
      unknown       = "Normalization applied upstream; method unknown."
    )
    
    tt_trans <- list(
      none    = "No transformation applied.",
      log2p0.1= "base::log2(x + 0.1)",
      vst     = "VST: dds <- DESeq2::DESeqDataSetFromMatrix(countData=counts_mat, ...) |> DESeq2::estimateSizeFactors(); x <- DESeq2::varianceStabilizingTransformation(dds) |> SummarizedExperiment::assay()",
      rlog    = "rlog: dds <- DESeq2::DESeqDataSetFromMatrix(countData=counts_mat, ...) |> DESeq2::estimateSizeFactors(); x <- DESeq2::rlog(dds, blind=TRUE) |> SummarizedExperiment::assay()",
      unknown = "Transformation applied upstream; method unknown."
    )
    
    tt_scale <- list(
      yes     = "Per-gene z-score across samples: x <- t(base::scale(t(x)))",
      no      = "No per-gene z-scoring.",
      unknown = "Scaling/z-scoring status unknown (upstream preprocessing unclear)."
    )
    
    # Current user-selected preprocessing state (no auto/user provenance)
    st <- get_pp_state()
    
    # ---- Canonical tri-state inputs (YOUR requested logic) -------------------
    # norm_status:  raw -> no, unknown -> unknown, anything else -> yes
    norm_in <- st$norm_status %||% "unknown"
    norm_in_state <- if (identical(norm_in, "raw")) {
      "no"
    } else if (identical(norm_in, "unknown")) {
      "unknown"
    } else {
      "yes"
    }
    
    # trans_applied/log_applied: none/no -> no, unknown -> unknown, anything else -> yes
    trans_in <- st$log_applied %||% "unknown"
    trans_in_state <- if (trans_in %in% c("no", "none")) {
      "no"
    } else if (identical(trans_in, "unknown")) {
      "unknown"
    } else {
      "yes"
    }
    
    # scaled_applied: no/none -> no, unknown -> unknown, anything else -> yes
    scale_in <- st$scaled_applied %||% "unknown"
    scale_in_state <- if (scale_in %in% c("no", "none")) {
      "no"
    } else if (identical(scale_in, "unknown")) {
      "unknown"
    } else {
      "yes"
    }
    
    # ---- Canonical method codes for INPUT -----------------------------------
    norm_method_code  <- canon_norm_code(st$norm_method %||% "unknown")
    trans_method_code <- canon_trans_code(st$log_method  %||% "unknown")
    
    input_norm_code <- if (identical(norm_in_state, "no")) {
      "none"
    } else if (identical(norm_in_state, "unknown")) {
      "unknown"
    } else {
      # "yes": already normalized; method may still be unknown
      norm_method_code %||% "unknown"
    }
    
    input_trans_code <- if (identical(trans_in_state, "no")) {
      "none"
    } else if (identical(trans_in_state, "unknown")) {
      "unknown"
    } else {
      trans_method_code %||% "unknown"
    }
    
    input_scale_code <- if (identical(scale_in_state, "no")) {
      "none"
    } else if (identical(scale_in_state, "unknown")) {
      "unknown"
    } else {
      "zscore"
    }
    
    # ---- User-facing labels --------------------------------------------------
    norm_label <- function(code, status = NULL) {
      if (!is.null(status) && identical(status, "raw"))     return("None")
      if (!is.null(status) && identical(status, "unknown")) return("Unknown")
      switch(
        code,
        deseq2_sf     = "DESeq2",
        tmm           = "TMM→CPM",
        cpm           = "CPM",
        tpm           = "TPM",
        fpkm          = "FPKM/RPKM",
        upperquartile = "Upper-quartile",
        quantile      = "Quantile",
        none          = "None",
        unknown       = "Normalized (Method Unknown)",
        code
      )
    }
    
    trans_label <- function(code) {
      switch(
        code,
        log2p0.1 = "log2p0.1",
        vst      = "VST",
        rlog     = "rlog",
        none     = "None",
        unknown  = "Unknown",
        code
      )
    }
    
    scale_label <- function(code) {
      switch(
        code,
        zscore  = "Yes",
        none    = "No",
        unknown = "Unknown",
        code
      )
    }
    
    # ---- Pill right-stripe state (YES/NO/UNKNOWN) ----------------------------
    # For norm/trans, if upstream says "yes" but method unknown => show unknown stripe.
    norm_state <- if (identical(norm_in_state, "no")) {
      "no"
    } else if (identical(norm_in_state, "unknown")) {
      "unknown"
    } else {
      if (identical(norm_method_code, "unknown")) "unknown" else "yes"
    }
    
    trans_state <- if (identical(trans_in_state, "no")) {
      "no"
    } else if (identical(trans_in_state, "unknown")) {
      "unknown"
    } else {
      if (identical(trans_method_code, "unknown")) "unknown" else "yes"
    }
    
    scale_state <- if (identical(scale_in_state, "no")) {
      "no"
    } else if (identical(scale_in_state, "unknown")) {
      "unknown"
    } else {
      "yes"
    }
    
    comp <- preproc_compatibility()
    
    # ---- Build pills (reference) ---------------------------------------------
    ref_norm_pill  <- pp_attr_pill(
      norm_label(ref_norm_code),
      attr = "norm", state = "yes",
      tooltip = tt_norm[[ref_norm_code]] %||% tt_norm$unknown,
      method_key = ref_norm_code
    )
    
    ref_trans_pill <- pp_attr_pill(
      trans_label(ref_trans_code),
      attr = "trans", state = "yes",
      tooltip = tt_trans[[ref_trans_code]] %||% tt_trans$unknown,
      method_key = ref_trans_code
    )
    
    ref_scale_pill <- pp_attr_pill(
      "No",
      attr = "scale", state = "no",
      tooltip = tt_scale$no,
      method_key = "none"
    )
    
    # ---- Build pills (input) -------------------------------------------------
    in_norm_pill <- pp_attr_pill(
      norm_label(input_norm_code, norm_in),
      attr       = "norm",
      state      = norm_state,
      method_key = input_norm_code,
      tooltip    = if (identical(norm_in_state, "no")) {
        "Raw counts (no normalization applied upstream)."
      } else {
        tt_norm[[input_norm_code]] %||% tt_norm$unknown
      }
    )
    
    in_trans_pill <- pp_attr_pill(
      trans_label(input_trans_code),
      attr       = "trans",
      state      = trans_state,
      method_key = input_trans_code,
      tooltip    = tt_trans[[input_trans_code]] %||% tt_trans$unknown
    )
    
    in_scale_pill <- pp_attr_pill(
      scale_label(input_scale_code),
      attr       = "scale",
      state      = scale_state,
      method_key = input_scale_code,
      tooltip    = if (identical(input_scale_code, "zscore")) {
        tt_scale$yes
      } else if (identical(input_scale_code, "none")) {
        tt_scale$no
      } else {
        tt_scale$unknown
      }
    )
    
    tagList(
      fluidRow(
        column(
          6,
          card(
            card_header("Reference preprocessing"),
            card_body(
              pp_attr_item("Normalization",   ref_norm_pill),
              pp_attr_item("Transformation",  ref_trans_pill),
              pp_attr_item("Scaled/z-scored", ref_scale_pill)
            )
          )
        ),
        column(
          6,
          card(
            card_header(
              tags$div(
                style = "display:flex; align-items:center; justify-content:space-between; gap:8px;",
                tags$span("Input preprocessing"),
                bslib::popover(
                  trigger = tags$button(
                    type = "button",
                    class = "btn btn-outline-secondary btn-sm pp-help-btn",
                    `aria-label` = "Acceptable preprocessing before PCA",
                    "?"
                  ),
                  tags$div(
                    style = "min-width: 520px;",
                    tags$p(
                      class = "pp-muted",
                      "Preferred methods are most common for RNA-seq PCA. Acceptable methods are used in practice, but require clear documentation."
                    ),
                    tags$div(
                      class = "pp-acceptable-grid",
                      tags$div(
                        class = "pp-acceptable-block",
                        tags$div(class = "pp-acceptable-title", "Normalization"),
                        tags$div(
                          class = "pp-acceptable-tier",
                          tags$div(class = "pp-acceptable-tier-label", "Preferred"),
                          tags$div(
                            class = "pp-acceptable-pills",
                            pp_attr_pill("DESeq2",  attr="norm", state="neutral", tooltip=tt_norm$deseq2_sf, tier="preferred", method_key="deseq2_sf"),
                            pp_attr_pill("TMM→CPM", attr="norm", state="neutral", tooltip=tt_norm$tmm,      tier="preferred", method_key="tmm"),
                            pp_attr_pill("TPM",     attr="norm", state="neutral", tooltip=tt_norm$tpm,      tier="preferred", method_key="tpm")
                          )
                        ),
                        tags$div(
                          class = "pp-acceptable-tier",
                          tags$div(class = "pp-acceptable-tier-label", "Acceptable"),
                          tags$div(
                            class = "pp-acceptable-pills",
                            pp_attr_pill("CPM",            attr="norm", state="neutral", tooltip=tt_norm$cpm,           tier="acceptable", method_key="cpm"),
                            pp_attr_pill("FPKM/RPKM",      attr="norm", state="neutral", tooltip=tt_norm$fpkm,          tier="acceptable", method_key="fpkm"),
                            pp_attr_pill("Upper-quartile", attr="norm", state="neutral", tooltip=tt_norm$upperquartile, tier="acceptable", method_key="upperquartile"),
                            pp_attr_pill("Quantile",       attr="norm", state="neutral", tooltip=tt_norm$quantile,      tier="acceptable", method_key="quantile")
                          )
                        )
                      ),
                      tags$div(
                        class = "pp-acceptable-block",
                        tags$div(class = "pp-acceptable-title", "Transformation"),
                        tags$div(
                          class = "pp-acceptable-tier",
                          tags$div(class = "pp-acceptable-tier-label", "Preferred"),
                          tags$div(
                            class = "pp-acceptable-pills",
                            pp_attr_pill("log2p0.1", attr="trans", state="neutral", tooltip=tt_trans$log2p0.1, tier="preferred", method_key="log2p0.1"),
                            pp_attr_pill("VST",      attr="trans", state="neutral", tooltip=tt_trans$vst,      tier="preferred", method_key="vst"),
                            pp_attr_pill("rlog",     attr="trans", state="neutral", tooltip=tt_trans$rlog,     tier="preferred", method_key="rlog")
                          )
                        ),
                        tags$div(
                          class = "pp-acceptable-tier",
                          tags$div(class = "pp-acceptable-tier-label", "Acceptable"),
                          tags$div(
                            class = "pp-acceptable-pills",
                            pp_attr_pill("None", attr="trans", state="neutral", tooltip=tt_trans$none, tier="acceptable", method_key="none")
                          )
                        )
                      ),
                      tags$div(
                        class = "pp-acceptable-block",
                        tags$div(class = "pp-acceptable-title", "Scaling"),
                        tags$div(
                          class = "pp-acceptable-tier",
                          tags$div(class = "pp-acceptable-tier-label", "Preferred"),
                          tags$div(
                            class = "pp-acceptable-pills",
                            pp_attr_pill("No", attr="scale", state="no", tooltip=tt_scale$no, tier="preferred", method_key="none")
                          )
                        ),
                        tags$div(
                          class = "pp-acceptable-tier",
                          tags$div(class = "pp-acceptable-tier-label", "Acceptable"),
                          tags$div(
                            class = "pp-acceptable-pills",
                            pp_attr_pill("Yes", attr="scale", state="yes", tooltip=tt_scale$yes, tier="acceptable", method_key="zscore")
                          )
                        )
                      )
                    )
                  ),
                  title = "Acceptable preprocessing before PCA",
                  placement = "left",
                  options = list(
                    customClass = "pp-help-popover",
                    container = "body"
                  )
                )
              )
            ),
            card_body(
              pp_attr_item(
                "Normalization",
                tags$button(
                  type="button",
                  class="pp-card-pill-btn",
                  `data-step`="input_norm",
                  in_norm_pill
                )
              ),
              uiOutput("bucket_input_norm"),
              
              pp_attr_item(
                "Transformation",
                tags$button(
                  type="button",
                  class="pp-card-pill-btn",
                  `data-step`="input_trans",
                  in_trans_pill
                )
              ),
              uiOutput("bucket_input_trans"),
              
              pp_attr_item(
                "Scaled/z-scored",
                tags$button(
                  type="button",
                  class="pp-card-pill-btn",
                  `data-step`="input_scale",
                  in_scale_pill
                )
              ),
              uiOutput("bucket_input_scale")
            )
          )
        )
      )
    )
  })
  
  output$bucket_input_norm <- renderUI({
    if (!identical(bucket_step(), "input_norm")) return(NULL)
    
    st <- get_pp_state()
    
    pick_pill <- function(label, attr, selected, selected_state, method_key, tooltip = NULL) {
      pp_attr_pill(
        label,
        attr = attr,
        state = if (isTRUE(selected)) selected_state else "neutral",
        method_key = method_key,
        tooltip = tooltip
      )
    }
    
    bucket_shell(
      id    = "pp-bucket-input-norm",
      title = "Input: Normalization",
      body_ui = tagList(
        tags$div(class="pp-bucket-note",
                 "Describe what the uploaded matrix ALREADY is. The in-app normalization method is chosen in the Processing plan (below)."),
        
        tags$div(
          class="pp-pick-grid",
          
          local({
            sel <- identical(st$norm_status, "normalized")
            pp_pick_btn(
              "input_norm_status", "normalized",
              pick_pill("Already normalized", attr="norm", selected = sel, selected_state = "yes", method_key = st$norm_method),
              selected = sel
            )
          }),
          
          local({
            sel <- identical(st$norm_status, "raw")
            pp_pick_btn(
              "input_norm_status", "raw",
              pick_pill("Raw counts", attr="norm", selected = sel, selected_state = "no", method_key="none"),
              selected = sel,
              disabled = FALSE
            )
          })
        ),
        
        # If normalized: allow descriptive upstream method selection (optional).
        if (identical(st$norm_status, "normalized")) tagList(
          tags$div(class="pp-bucket-note",
                   "Optional: if you know the upstream normalization, select it (used for comparability checks)."),
          tags$div(
            class="pp-pick-grid",
            
            local({
              sel <- identical(st$norm_method, "unknown")
              pp_pick_btn(
                "input_norm_method","unknown",
                pick_pill("Unknown", attr="norm", selected = sel, selected_state = "unknown", method_key="unknown"),
                selected = sel
              )
            }),
            
            local({
              sel <- identical(st$norm_method, "deseq2_sf")
              pp_pick_btn(
                "input_norm_method","deseq2_sf",
                pick_pill("DESeq2", attr="norm", selected = sel, selected_state = "yes", method_key="deseq2_sf"),
                selected = sel
              )
            }),
            
            local({
              sel <- identical(st$norm_method, "tpm")
              pp_pick_btn(
                "input_norm_method","tpm",
                pick_pill("TPM", attr="norm", selected = sel, selected_state = "yes", method_key="tpm"),
                selected = sel
              )
            }),
            
            local({
              sel <- identical(st$norm_method, "cpm")
              pp_pick_btn(
                "input_norm_method","cpm",
                pick_pill("CPM", attr="norm", selected = sel, selected_state = "yes", method_key="cpm"),
                selected = sel
              )
            }),
            
            local({
              sel <- identical(st$norm_method, "fpkm")
              pp_pick_btn(
                "input_norm_method","fpkm",
                pick_pill("FPKM/RPKM", attr="norm", selected = sel, selected_state = "yes", method_key="fpkm"),
                selected = sel
              )
            })
          )
        )
      )
    )
  })
  
  
  
  output$bucket_input_trans <- renderUI({
    if (!identical(bucket_step(), "input_trans")) return(NULL)
    
    feat <- matrix_features()
    has_neg <- !is.null(feat) && isTRUE(feat$has_negative)
    st <- get_pp_state()
    
    pick_pill <- function(label, attr, selected, selected_state, method_key, tooltip = NULL) {
      pp_attr_pill(
        label,
        attr = attr,
        state = if (isTRUE(selected)) selected_state else "neutral",
        method_key = method_key,
        tooltip = tooltip
      )
    }
    
    bucket_shell(
      id    = "pp-bucket-input-trans",
      title = "Input: Transformation",
      body_ui = tagList(
        tags$div(class="pp-bucket-note",
                 "Describe whether the uploaded matrix already has a transformation applied (best match)."),
        if (has_neg) tags$div(
          class="pp-bucket-note",
          "Note: negative values are common for VST/rlog or centered/scaled matrices. ",
          "This section is descriptive, so methods are still selectable."
        ),
        
        tags$div(
          class="pp-pick-grid",
          
          local({
            sel <- identical(st$log_applied, "no")
            pp_pick_btn(
              "input_trans","none",
              pick_pill("None (not transformed)", attr="trans", selected = sel, selected_state = "no", method_key="none"),
              selected = sel
            )
          }),
          
          local({
            sel <- identical(st$log_applied, "yes") && identical(st$log_method, "log2p0.1")
            pp_pick_btn(
              "input_trans","log2p0.1",
              pick_pill("log2p0.1", attr="trans", selected = sel, selected_state = "yes", method_key="log2p0.1"),
              selected = sel
            )
          }),
          
          local({
            sel <- identical(st$log_applied, "yes") && identical(st$log_method, "vst")
            pp_pick_btn(
              "input_trans","vst",
              pick_pill("VST", attr="trans", selected = sel, selected_state = "yes", method_key="vst"),
              selected = sel
            )
          }),
          
          local({
            sel <- identical(st$log_applied, "yes") && identical(st$log_method, "rlog")
            pp_pick_btn(
              "input_trans","rlog",
              pick_pill("rlog", attr="trans", selected = sel, selected_state = "yes", method_key="rlog"),
              selected = sel
            )
          }),
          
          local({
            sel <- identical(st$log_applied, "yes") && identical(st$log_method, "unknown")
            pp_pick_btn(
              "input_trans","unknown",
              pick_pill("Unknown", attr="trans", selected = sel, selected_state = "unknown", method_key="unknown"),
              selected = sel
            )
          })
        )
      )
    )
  })
  
  
  
  output$bucket_input_scale <- renderUI({
    if (!identical(bucket_step(), "input_scale")) return(NULL)
    
    st <- get_pp_state()
    
    pick_pill <- function(label, attr, selected, selected_state, method_key, tooltip = NULL) {
      pp_attr_pill(
        label,
        attr = attr,
        state = if (isTRUE(selected)) selected_state else "neutral",
        method_key = method_key,
        tooltip = tooltip
      )
    }
    
    bucket_shell(
      id    = "pp-bucket-input-scale",
      title = "Input: Scaling",
      body_ui = tagList(
        tags$div(class="pp-bucket-note",
                 "Describe whether the uploaded matrix is already gene-wise z-scored/scaled."),
        
        tags$div(
          class="pp-pick-grid",
          
          local({
            sel <- identical(st$scaled_applied, "no")
            pp_pick_btn(
              "input_scale","none",
              pick_pill("Not scaled", attr="scale", selected = sel, selected_state = "no", method_key="none"),
              selected = sel
            )
          }),
          
          local({
            sel <- identical(st$scaled_applied, "yes")
            pp_pick_btn(
              "input_scale","zscore",
              pick_pill("Already scaled (z-score)", attr="scale", selected = sel, selected_state = "yes", method_key="zscore"),
              selected = sel
            )
          })
        )
      )
    )
  })
  
  
  ## -----------------------------
  ## Processing plan diagram
  ## -----------------------------
  output$preproc_plan_ui <- renderUI({
    req(rna_raw())
    
    # -----------------------------------------------------------------------
    # 1) INPUT tri-state (what the uploaded matrix ALREADY is)
    #    norm_status: "raw" -> no ; "unknown" -> unknown ; else -> yes
    #    trans_applied: c("none","no") -> no ; "unknown" -> unknown ; else -> yes
    #    scale_applied: c("none","no") -> no ; "unknown" -> unknown ; else -> yes
    # -----------------------------------------------------------------------
    
    norm_in  <- pp_input$norm_status   %||% "unknown"
    trans_in <- pp_input$trans_applied %||% "unknown"
    scale_in <- pp_input$scale_applied %||% "unknown"
    
    norm_in_state <- if (identical(norm_in, "raw")) {
      "no"
    } else if (identical(norm_in, "unknown")) {
      "unknown"
    } else {
      "yes"
    }
    
    trans_in_state <- if (identical(trans_in, "unknown")) {
      "unknown"
    } else if (trans_in %in% c("no", "none")) {
      "no"
    } else {
      "yes"
    }
    
    scale_in_state <- if (identical(scale_in, "unknown")) {
      "unknown"
    } else if (scale_in %in% c("no", "none")) {
      "no"
    } else {
      "yes"
    }
    
    # -----------------------------------------------------------------------
    # 2) LOCAL plan (what the app will do), derived from pp_plan + constraints
    # -----------------------------------------------------------------------
    
    plan_norm  <- pp_plan$norm  %||% "deseq2_sf"
    plan_trans <- pp_plan$trans %||% "log2p0.1"
    plan_scale <- pp_plan$scale %||% "no"         # expected: "yes" or "no"
    
    # Hard constraints:
    # - If input already scaled: do NOT do anything in-app (normalize/log/scale)
    if (identical(scale_in_state, "yes")) {
      plan_norm  <- "none"
      plan_trans <- "none"
      plan_scale <- "no"
    } else {
      # - If input already transformed: don't transform again
      if (identical(trans_in_state, "yes")) {
        plan_trans <- "none"
      }
      
      # - If input already normalized: don't normalize again
      if (identical(norm_in_state, "yes")) {
        plan_norm <- "none"
      }
      
      # - If user says input is raw (norm_in_state=="no"), allow plan_norm;
      #   if unknown, keep method but mark plan state as unknown (conditional).
      if (identical(norm_in_state, "unknown")) {
        # keep plan_norm as chosen, but state will become "unknown" below
      }
    }
    
    # -----------------------------------------------------------------------
    # 3) PLAN states (does the app apply the step?)
    #    yes/no/unknown reflect APPLY-IN-APP certainty.
    # -----------------------------------------------------------------------
    
    # Normalization is only meaningful if input is raw & not already scaled.
    norm_plan_state <- if (identical(plan_norm, "none")) {
      "no"
    } else if (identical(scale_in_state, "yes")) {
      "no"
    } else if (identical(norm_in_state, "no") && identical(scale_in_state, "no")) {
      "yes"
    } else {
      # norm_in_state unknown OR scale_in_state unknown
      "unknown"
    }
    
    trans_plan_state <- if (identical(plan_trans, "none")) {
      "no"
    } else if (identical(scale_in_state, "yes") || identical(trans_in_state, "yes")) {
      "no"
    } else if (identical(scale_in_state, "unknown") || identical(trans_in_state, "unknown")) {
      "unknown"
    } else {
      "yes"
    }
    
    scale_plan_state <- if (!identical(plan_scale, "yes")) {
      "no"
    } else if (identical(scale_in_state, "yes")) {
      "no"
    } else if (identical(scale_in_state, "unknown")) {
      "unknown"
    } else {
      "yes"
    }
    
    # -----------------------------------------------------------------------
    # 4) Label maps
    # -----------------------------------------------------------------------
    
    norm_map  <- c(
      deseq2_sf     = "DESeq2",
      cpm           = "CPM",
      tmm           = "TMM→CPM",
      upperquartile = "Upper-quartile→CPM",
      none          = "None"
    )
    trans_map <- c(
      log2p0.1 = "log2p0.1",
      vst      = "VST",
      rlog     = "rlog",
      none     = "None"
    )
    scale_map <- c(
      zscore = "Yes (z-score)",
      none   = "No"
    )
    
    # -----------------------------------------------------------------------
    # 5) Build pills
    # -----------------------------------------------------------------------
    
    # Normalization pill
    norm_pill <- if (!identical(plan_norm, "none")) {
      label <- norm_map[[plan_norm]] %||% plan_norm
      pp_attr_pill(
        label,
        attr = "norm",
        state = norm_plan_state,
        tooltip = if (identical(norm_plan_state, "unknown")) {
          "Normalization will be applied only if the uploaded matrix is raw counts and unscaled."
        } else {
          "Normalization applied in-app before PCA (raw counts only)."
        },
        method_key = plan_norm
      )
    } else {
      pp_attr_pill(
        "None",
        attr = "norm",
        state = "no",
        tooltip = "No in-app normalization will be applied.",
        method_key = "none"
      )
    }
    
    # Transformation pill
    trans_pill <- if (!identical(plan_trans, "none")) {
      label <- trans_map[[plan_trans]] %||% plan_trans
      pp_attr_pill(
        label,
        attr = "trans",
        state = trans_plan_state,
        tooltip = if (identical(trans_plan_state, "unknown")) {
          "Transformation will be applied only if the input is not already transformed/scaled."
        } else {
          "Transformation applied in-app before PCA (if any)."
        },
        method_key = plan_trans
      )
    } else {
      pp_attr_pill(
        "None",
        attr = "trans",
        state = "no",
        tooltip = "No in-app transformation will be applied.",
        method_key = "none"
      )
    }
    
    # Scaling pill
    scale_code  <- if (identical(plan_scale, "yes")) "zscore" else "none"
    scale_label <- scale_map[[scale_code]] %||% scale_code
    
    scale_pill <- pp_attr_pill(
      scale_label,
      attr = "scale",
      state = scale_plan_state,
      tooltip = if (identical(scale_plan_state, "unknown")) {
        "Scaling will be applied only if the input is not already scaled."
      } else {
        "Gene-wise scaling applied in-app before PCA (if any)."
      },
      method_key = scale_code
    )
    
    # -----------------------------------------------------------------------
    # 6) UI
    # -----------------------------------------------------------------------
    
    diagram_btn <- function(id, title, pill, attr) {
      actionButton(
        inputId = id,
        label   = tagList(tags$div(class = "pp-diag-title", title), pill),
        class   = paste("pp-diag-btn", paste0("pp-diag-btn--", attr))
      )
    }
    
    card(
      card_header(
        tags$div(
          style = "display:flex; align-items:center; justify-content:space-between; gap:10px; flex-wrap:wrap;",
          tags$span("Processing plan")
        )
      ),
      card_body(
        tags$div(style = "font-weight:600; margin-bottom:6px;",
                 "Plan steps (click to adjust):"),
        tags$div(
          class = "pp-diagram",
          diagram_btn("diag_norm",  "Normalization",  norm_pill,  "norm"),
          tags$span(class = "pp-diag-arrow", "→"),
          diagram_btn("diag_log",   "Transformation", trans_pill, "trans"),
          tags$span(class = "pp-diag-arrow", "→"),
          diagram_btn("diag_scale", "Scaling",        scale_pill, "scale")
        )
      )
    )
  })
  
  
  output$bucket_plan <- renderUI({
    req(wizard_step() == 2)
    step <- bucket_step()
    if (is.null(step) || !step %in% c("plan_norm","plan_trans","plan_scale")) return(NULL)
    
    feat <- matrix_features()
    
    has_neg        <- !is.null(feat) && isTRUE(feat$has_negative)
    likely_loglike <- !is.null(feat) && isTRUE(feat$likely_loglike)
    frac_integer   <- if (!is.null(feat)) feat$frac_integer else NA_real_
    
    has_edgeR <- requireNamespace("edgeR", quietly = TRUE)
    
    st <- get_pp_state()
    
    # Locks are based ONLY on what the user said the INPUT already is
    # (treating unknown as "not locked")
    lock_norm  <- identical(st$scaled_applied, "yes") || identical(st$norm_status, "normalized")
    lock_trans <- identical(st$scaled_applied, "yes") || identical(st$log_applied, "yes")
    lock_scale <- identical(st$scaled_applied, "yes")
    
    warn_box <- function(type = c("warning","info"), title = NULL, text = NULL) {
      type <- match.arg(type)
      tags$div(
        class = paste0("alert alert-", type, " py-2 px-3 mb-2"),
        role  = "alert",
        if (!is.null(title) && nzchar(title)) tags$div(tags$strong(title)),
        if (!is.null(text)  && nzchar(text))  tags$div(text)
      )
    }
    
    # Helper: unselected options stay neutral; selected option gets semantic color
    pick_pill <- function(label, attr, selected, selected_state, method_key, tooltip = NULL) {
      pp_attr_pill(
        label,
        attr = attr,
        state = if (isTRUE(selected)) selected_state else "neutral",
        method_key = method_key,
        tooltip = tooltip
      )
    }
    
    title <- switch(
      step,
      plan_norm  = "Plan: Normalization to apply",
      plan_trans = "Plan: Transformation to apply",
      plan_scale = "Plan: Scaling to apply",
      "Plan options"
    )
    
    body <- switch(
      step,
      
      plan_norm = {
        tagList(
          if (identical(st$scaled_applied, "yes")) {
            tags$div(class="pp-bucket-note",
                     "Plan normalization is disabled because input is already scaled/z-scored.")
          } else if (identical(st$norm_status, "normalized")) {
            tags$div(class="pp-bucket-note",
                     "Plan normalization is disabled because you marked the input as already normalized.")
          } else {
            tags$div(class="pp-bucket-note",
                     "Choose the normalization the app will apply (raw counts only).")
          },
          
          # Detection-based warnings ONLY (no disabling)
          if (!lock_norm && has_neg) {
            warn_box(
              "warning",
              "Negative values detected",
              "Raw counts are typically non-negative. If you run count-based normalization, it will likely error unless the input is truly raw counts."
            )
          },
          if (!lock_norm && !is.na(frac_integer) && frac_integer < 0.90) {
            warn_box(
              "info",
              "Matrix does not look integer-like",
              sprintf("Only %.2f of values look integer-like. Raw counts are often near-integer (heuristic only).", frac_integer)
            )
          },
          if (!lock_norm && likely_loglike) {
            warn_box(
              "info",
              "Matrix looks log-like",
              "If the matrix is already on a log scale, choosing raw-count normalization is likely incorrect (heuristic only)."
            )
          },
          
          tags$div(
            class="pp-pick-grid",
            
            local({
              sel <- identical(st$raw_norm_plan, "deseq2_sf")
              pp_pick_btn(
                "raw_norm_plan","deseq2_sf",
                pick_pill("DESeq2", attr="norm", selected = sel, selected_state = "yes", method_key="deseq2_sf"),
                selected = sel,
                disabled = lock_norm,
                dim      = lock_norm
              )
            }),
            
            local({
              sel <- identical(st$raw_norm_plan, "cpm")
              pp_pick_btn(
                "raw_norm_plan","cpm",
                pick_pill("CPM", attr="norm", selected = sel, selected_state = "yes", method_key="cpm"),
                selected = sel,
                disabled = lock_norm,
                dim      = lock_norm
              )
            }),
            
            local({
              sel <- identical(st$raw_norm_plan, "tmm")
              pp_pick_btn(
                "raw_norm_plan","tmm",
                pick_pill("TMM→CPM", attr="norm", selected = sel, selected_state = "yes", method_key="tmm"),
                selected = sel,
                disabled = lock_norm || !has_edgeR,
                dim      = lock_norm,
                tooltip  = if (!has_edgeR) "Disabled: edgeR not installed." else NULL
              )
            }),
            
            local({
              sel <- identical(st$raw_norm_plan, "upperquartile")
              pp_pick_btn(
                "raw_norm_plan","upperquartile",
                pick_pill("Upper-quartile→CPM", attr="norm", selected = sel, selected_state = "yes", method_key="upperquartile"),
                selected = sel,
                disabled = lock_norm || !has_edgeR,
                dim      = lock_norm,
                tooltip  = if (!has_edgeR) "Disabled: edgeR not installed." else NULL
              )
            })
          )
        )
      },
      
      plan_trans = {
        tagList(
          if (identical(st$scaled_applied, "yes")) {
            tags$div(class="pp-bucket-note",
                     "Plan transformation is disabled because input is already scaled/z-scored.")
          } else if (identical(st$log_applied, "yes")) {
            tags$div(class="pp-bucket-note",
                     "Plan transformation is disabled because you marked the input as already transformed.")
          } else {
            tags$div(class="pp-bucket-note",
                     "Choose what the app will apply before PCA (if anything).")
          },
          
          # Detection-based warnings ONLY (no disabling)
          if (!lock_trans && has_neg) {
            warn_box(
              "info",
              "Negatives are present",
              "If you choose log2(x+0.1) in-app, the matrix must be non-negative at that step. Negatives can happen for VST/rlog, centering, or log2(x+ε) with ε<1."
            )
          },
          if (!lock_trans && likely_loglike) {
            warn_box(
              "info",
              "Matrix looks log-like",
              "Applying log2(x+0.1) may be redundant (heuristic only)."
            )
          },
          
          tags$div(
            class="pp-pick-grid",
            
            local({
              sel <- identical(st$plan_trans, "log2p0.1")
              pp_pick_btn(
                "plan_trans","log2p0.1",
                pick_pill("log2p0.1", attr="trans", selected = sel, selected_state = "yes", method_key="log2p0.1"),
                selected = sel,
                disabled = lock_trans,
                dim      = lock_trans || likely_loglike || has_neg,
                tooltip  = if (has_neg) {
                  "If the matrix contains negatives at the log2(x+0.1) step, the app will error. Negatives may indicate prior transforms/centering."
                } else if (likely_loglike) {
                  "Matrix looks log-like; log2(x+0.1) may be redundant."
                } else NULL
              )
            }),
            
            local({
              sel <- identical(st$plan_trans, "vst")
              pp_pick_btn(
                "plan_trans","vst",
                pick_pill("VST", attr="trans", selected = sel, selected_state = "yes", method_key="vst"),
                selected = sel,
                disabled = lock_trans || !identical(st$norm_status, "raw"),
                dim      = lock_trans,
                tooltip  = if (!identical(st$norm_status, "raw")) "Disabled: VST requires raw counts input mode." else NULL
              )
            }),
            
            local({
              sel <- identical(st$plan_trans, "rlog")
              pp_pick_btn(
                "plan_trans","rlog",
                pick_pill("rlog", attr="trans", selected = sel, selected_state = "yes", method_key="rlog"),
                selected = sel,
                disabled = lock_trans || !identical(st$norm_status, "raw"),
                dim      = lock_trans,
                tooltip  = if (!identical(st$norm_status, "raw")) "Disabled: rlog requires raw counts input mode." else NULL
              )
            }),
            
            local({
              sel <- identical(st$plan_trans, "none")
              pp_pick_btn(
                "plan_trans","none",
                pick_pill("None", attr="trans", selected = sel, selected_state = "no", method_key="none"),
                selected = sel,
                disabled = lock_trans,
                dim      = lock_trans
              )
            })
          )
        )
      },
      
      plan_scale = {
        tagList(
          if (identical(st$scaled_applied, "yes")) {
            tags$div(class="pp-bucket-note",
                     "Plan scaling is disabled because you marked the input as already scaled/z-scored.")
          } else {
            tags$div(class="pp-bucket-note",
                     "Choose whether the app will apply per-gene z-scoring before PCA.")
          },
          
          tags$div(
            class="pp-pick-grid",
            
            local({
              sel <- identical(st$plan_scale, "no")
              pp_pick_btn(
                "plan_scale","no",
                pick_pill("No", attr="scale", selected = sel, selected_state = "no", method_key="none"),
                selected = sel,
                disabled = lock_scale,
                dim      = lock_scale
              )
            }),
            
            local({
              sel <- identical(st$plan_scale, "yes")
              pp_pick_btn(
                "plan_scale","yes",
                pick_pill("Yes (z-score)", attr="scale", selected = sel, selected_state = "yes", method_key="zscore"),
                selected = sel,
                disabled = lock_scale,
                dim      = lock_scale
              )
            })
          )
        )
      }
    )
    
    bucket_shell("pp-bucket-plan", title, body)
  })
  
  
  
  
  
  
  
  output$gene_overlap_ui <- renderUI({
    go <- gene_overlap()
    if (is.null(go)) return(NULL)
    
    match_color <- if (go$n_matched >= 50) {
      "success"
    } else if (go$n_matched > 0) {
      "warning"
    } else {
      "danger"
    }
    
    match_note <- if (!is.na(go$pct_matched_ref)) {
      base <- paste0(sprintf("%.1f", go$pct_matched_ref), "% of reference genes present in input.")
      if (go$n_matched < 50) {
        paste0(base, " Warning: fewer than 50 matched genes; PCA similarity may be unstable.")
      } else {
        base
      }
    } else {
      NULL
    }
    
    card(
      card_header("Gene overlap"),
      card_body(
        row_item("Reference genes",          pill_tag(as.character(go$n_reference), "dark")),
        row_item("Genes in uploaded file",   pill_tag(as.character(go$n_input), "secondary")),
        row_item("Matched to reference list", pill_tag(as.character(go$n_matched), match_color),
                 note = match_note)
      )
    )
  })
  
  output$compatibility_ui <- renderUI({
    compat <- preproc_compatibility()
    
    status  <- compat$status
    reason  <- compat$reason
    reco    <- compat$recommendation
    
    state <- if (identical(status, "Compatible")) "yes" else "no"
    method_key <- if (state == "yes") "match" else "diff"
    
    pill <- pp_attr_pill(
      status,                      # first arg = text
      attr       = "compat",
      state      = state,
      tooltip    = reason,
      method_key = method_key
    )
    
    card(
      class = "pp-card-compat",
      card_header("Compatibility check"),
      card_body(
        tags$div(
          class = "pp-compat-main",
          tags$div(class = "pp-compat-pill-wrapper", pill),
          tags$div(
            class = "pp-compat-text small",
            if (!is.null(reason) && nzchar(reason)) {
              tags$div(tags$strong("Reason: "), reason)
            },
            if (!is.null(reco) && nzchar(reco)) {
              tags$div(tags$strong("Recommendation: "), reco)
            }
          )
        )
      )
    )
  })
  
  
  observeEvent(input$pp_pick, {
    info <- input$pp_pick
    req(info$kind, info$value)
    
    kind0 <- as.character(info$kind)
    value <- as.character(info$value)
    
    # ----- INPUT side -----
    if (kind0 == "input_norm_status") {
      pp_input$norm_status <- value
      if (identical(value, "normalized") && !nzchar(pp_input$norm_method)) {
        pp_input$norm_method <- "unknown"
      }
    }
    
    if (kind0 == "input_norm_method") {
      pp_input$norm_method <- value
      pp_input$norm_status <- "normalized"
    }
    
    # FIX: input_trans sends method/none; store as (applied yes/no) + method
    if (kind0 == "input_trans") {
      if (identical(value, "none")) {
        pp_input$trans_applied <- "no"
        # pp_input$trans_method can stay as-is (it’s irrelevant when trans_applied=="no")
      } else {
        pp_input$trans_applied <- "yes"
        pp_input$trans_method  <- value  # log2p0.1/vst/rlog/unknown
      }
    }
    
    # FIX: input_scale sends none/zscore; store as yes/no
    if (kind0 == "input_scale") {
      pp_input$scale_applied <- if (identical(value, "zscore")) "yes" else "no"
    }
    
    # ----- PLAN side -----
    if (kind0 == "raw_norm_plan") pp_plan$norm  <- value
    if (kind0 == "plan_trans")    pp_plan$trans <- value
    if (kind0 == "plan_scale")    pp_plan$scale <- value
    
    # hard constraints after any pick
    if (identical(pp_input$scale_applied, "yes")) {
      pp_plan$scale <- "no"
      pp_plan$trans <- "none"
    }
    if (identical(pp_input$trans_applied, "yes")) {
      pp_plan$trans <- "none"
    }
    
    pp_warn_pick(kind0, value, matrix_features())
  })
  
  
  pp_warn_pick <- function(kind, value, feat) {
    if (is.null(feat)) return(invisible(NULL))
    
    has_neg       <- isTRUE(feat$has_negative)
    frac_integer  <- feat$frac_integer %||% NA_real_
    maxv          <- feat$maxv %||% NA_real_
    likely_loglike<- isTRUE(feat$likely_loglike)
    likely_scaled <- isTRUE(feat$likely_scaled)
    
    # convenience
    warn <- function(txt) showNotification(txt, type="warning", duration = 7)
    info <- function(txt) showNotification(txt, type="message", duration = 7)
    
    # INPUT: norm status
    if (kind == "input_norm_status" && value == "raw") {
      if (has_neg) {
        warn("You selected 'Raw counts', but negative values were detected. Raw counts are typically non-negative. This choice is probably incorrect.")
      }
      if (!is.na(frac_integer) && frac_integer < 0.90) {
        warn(sprintf("You selected 'Raw counts', but only %.2f of values look integer-like. Raw counts are usually near-integer. This choice is probably incorrect.", frac_integer))
      }
      if (likely_loglike) {
        info("You selected 'Raw counts', but the matrix looks log-like (value range/structure). If you proceed, count-based steps may fail later.")
      }
    }
    
    if (kind == "input_norm_status" && value == "normalized") {
      # strong count-like signal
      if (!has_neg && !is.na(frac_integer) && !is.na(maxv) && frac_integer >= 0.95 && maxv > 50) {
        info("You selected 'Already normalized', but the matrix looks count-like (high integer fraction and large values). If it is truly normalized, you can keep this—just be aware this contradicts the heuristic.")
      }
    }
    
    # INPUT: transformation
    if (kind == "input_trans" && value == "none") {
      if (likely_loglike) {
        info("You selected 'Not transformed', but the matrix looks log-like. If you later apply log2(x+0.1), it may be redundant.")
      }
    }
    
    if (!is.null(feat) && kind == "input_trans" && !(value %in% c("none", "unknown"))) {
      likely_loglike <- isTRUE(feat$likely_loglike)
      q99 <- feat$q99 %||% NA_real_
      
      # only warn on strong evidence of non-log scale
      strong_not_log <- isFALSE(likely_loglike) && is.finite(q99) && q99 >= 200
      
      if (strong_not_log) {
        info(sprintf(
          "You selected a log transformation, but values look non-log-scale (q99=%.1f). Please double-check upstream preprocessing.",
          q99
        ))
      }
    }
    
    
    
    # INPUT: scaling
    if (kind == "input_scale" && value == "zscore") {
      if (!likely_scaled) {
        warn("You indicated the matrix is already z-scored/scaled, but the heuristic did not detect strong scaling signatures. This is allowed, but compatibility may be affected.")
      }
    }
    if (kind == "input_scale" && value == "none") {
      if (likely_scaled) {
        info("You indicated 'Not scaled', but the heuristic suggests scaling may already be present. This is allowed if you're sure.")
      }
    }
    
    invisible(NULL)
  }

  observeEvent(input$diag_norm,  { bucket_step("plan_norm");  session$sendCustomMessage("pp_scroll_to", list(id="pp-bucket-plan")) })
  observeEvent(input$diag_log,   { bucket_step("plan_trans"); session$sendCustomMessage("pp_scroll_to", list(id="pp-bucket-plan")) })
  observeEvent(input$diag_scale, { bucket_step("plan_scale"); session$sendCustomMessage("pp_scroll_to", list(id="pp-bucket-plan")) })
  
  bucket_dom_id <- function(step) {
    switch(step,
           input_norm  = "pp-bucket-input-norm",
           input_trans = "pp-bucket-input-trans",
           input_scale = "pp-bucket-input-scale",
           plan_norm   = "pp-bucket-plan",
           plan_trans  = "pp-bucket-plan",
           plan_scale  = "pp-bucket-plan",
           "pp-bucket-plan")
  }
  
  toggle_bucket <- function(step) {
    step <- as.character(step)
    if (!nzchar(step)) return()
    
    if (identical(bucket_step(), step)) {
      bucket_step(NULL)
    } else {
      bucket_step(step)
      session$sendCustomMessage("pp_scroll_to", list(id = bucket_dom_id(step)))
    }
  }
  
  
  observeEvent(input$pp_open, {
    req(input$pp_open)
    toggle_bucket(input$pp_open)
  }, ignoreInit = TRUE)
  
  observeEvent(input$pp_bucket_close, {
    bucket_step(NULL)
  })
  
  
  observeEvent(wizard_step(), {
    if (wizard_step() != 2) bucket_step(NULL)
  })
  
  ## -----------------------------
  ## Reference PCA model
  ## -----------------------------
  hnscc_pca_reactive <- shiny::reactive({
    readRDS("App_files/hpv_signature_tcga_hnscc_pca_model.rds")
  })
  
  rna_dds_obj <- reactive({
    x <- rna_raw()
    if (inherits(x, "DESeqDataSet")) x else NULL
  })
  
  
  ## -----------------------------
  ## Preprocessing application before PCA
  ## -----------------------------
  apply_raw_normalization <- function(counts_mat, method, dds_obj = NULL) {
    stopifnot(all(counts_mat >= 0, na.rm = TRUE))
    if (method == "deseq2_sf") {
      if (!is.null(dds_obj) && inherits(dds_obj, "DESeqDataSet")) {
        dds <- dds_obj
        if (is.null(DESeq2::sizeFactors(dds)) && is.null(DESeq2::normalizationFactors(dds))) {
          dds <- DESeq2::estimateSizeFactors(dds)
        }
        
        norm_mat <- DESeq2::counts(dds, normalized = TRUE)
        # IMPORTANT: bring gene IDs back to HGNC symbols so filtering works
        norm_mat <- harmonize_gene_ids_once(
          mat       = norm_mat,
          ref_genes = FILTER_GENE_LIST,
          dds       = dds,          # lets rowData columns win if available
          min_match = 20,
          assume_counts = TRUE     # normalized -> average duplicates (per your rule)
        )
        return(norm_mat)
      }
      
      # fallback TSV-only counts
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts_mat,
        colData   = data.frame(row.names = colnames(counts_mat)),
        design    = ~ 1
      )
      dds <- DESeq2::estimateSizeFactors(dds)
      return(DESeq2::counts(dds, normalized = TRUE))
    }
    
    if (method == "cpm") {
      lib <- colSums(counts_mat, na.rm = TRUE)
      lib[lib == 0] <- NA_real_
      return(sweep(counts_mat, 2, lib / 1e6, "/"))
    }
    
    if (method == "tmm") {
      if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR is required for TMM normalization. Please install edgeR or choose another method.")
      }
      dge <- edgeR::DGEList(counts = counts_mat)
      dge <- edgeR::calcNormFactors(dge, method = "TMM")
      # return CPM after TMM scaling
      return(edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = FALSE, prior.count = 0))
    }
    
    if (method == "upperquartile") {
      if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR is required for upper-quartile normalization. Please install edgeR or choose another method.")
      }
      dge <- edgeR::DGEList(counts = counts_mat)
      dge <- edgeR::calcNormFactors(dge, method = "upperquartile")
      return(edgeR::cpm(dge, normalized.lib.sizes = TRUE, log = FALSE, prior.count = 0))
    }
    
    
    stop("Unknown raw normalization plan method.")
  }
  
  preprocess_for_pca <- function(mat_full, dds_obj = NULL) {
    # mat_full: numeric matrix, rows=genes, cols=samples (FULL gene set)
    
    mat_num <- mat_full
    
    # If input already scaled, do not modify; still filter to signature genes for PCA
    if (identical(pp_input$scale_applied, "yes")) {
      mat_proc <- mat_num
      
      if (!is.null(FILTER_GENE_LIST) && length(FILTER_GENE_LIST) > 0) {
        keep <- rownames(mat_proc) %in% FILTER_GENE_LIST
        mat_proc <- mat_proc[keep, , drop = FALSE]
      }
      
      return(list(
        mat_proc   = mat_proc,
        scale_unit = FALSE,
        comparable = FALSE
      ))
    }
    
    # Decide planned in-app transformation up front (used to gate normalization)
    tr_plan <- if (!identical(pp_input$trans_applied, "yes")) {
      pp_plan$trans %||% "log2p0.1"
    } else {
      "none"
    }
    
    mat_proc <- mat_num
    
    # 1) Normalize (ONLY if user says raw) — do this on FULL gene set
    # IMPORTANT: if we plan VST/rlog, skip separate normalization; DESeq2 handles size factors internally.
    if (identical(pp_input$norm_status, "raw") && !(tr_plan %in% c("vst", "rlog"))) {
      validate(need(!any(mat_num < 0, na.rm = TRUE),
                    "Raw counts-based normalization selected but matrix contains negative values."))
      method <- pp_plan$norm %||% "deseq2_sf"
      mat_proc <- apply_raw_normalization(mat_num, method, dds_obj = dds_obj)
      mat_proc[!is.finite(mat_proc)] <- 0
      mat_proc[is.na(mat_proc)] <- 0
    }
    
    # 2) Transform unless user says already transformed — do this on FULL gene set
    if (!identical(pp_input$trans_applied, "yes")) {
      
      if (tr_plan %in% c("vst", "rlog")) {
        
        validate(need(identical(pp_input$norm_status, "raw"),
                      "VST/rlog can only be applied when the input matrix is raw counts."))
        
        validate(need(!any(mat_num < 0, na.rm = TRUE),
                      "Cannot apply VST/rlog: matrix contains negative values."))
        
        frac_int <- mean(abs(mat_num - round(mat_num)) < 1e-6, na.rm = TRUE)
        validate(need(frac_int >= 0.90,
                      "Matrix does not appear to be raw integer counts; VST/rlog require raw counts."))
        
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = round(mat_num),
          colData   = data.frame(row.names = colnames(mat_num)),
          design    = ~ 1
        )
        dds <- DESeq2::estimateSizeFactors(dds)
        
        mat_proc <- if (tr_plan == "vst") {
          SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = TRUE))
        } else {
          SummarizedExperiment::assay(DESeq2::rlog(dds, blind = TRUE))
        }
        
        mat_proc[!is.finite(mat_proc)] <- 0
        mat_proc[is.na(mat_proc)] <- 0
        
      } else if (tr_plan == "log2p0.1") {
        
        validate(need(!any(mat_proc < 0, na.rm = TRUE),
                      "Cannot apply log2(x+0.1): negatives present."))
        mat_proc <- log2(mat_proc + 0.1)
        
      } else if (tr_plan == "none") {
        # no-op
        
      } else {
        # fallback: treat unknown as log2(x+pseudocount)
        validate(need(!any(mat_proc < 0, na.rm = TRUE),
                      "Cannot apply log2(x+pseudocount): matrix contains negative values."))
        mat_proc <- log2(mat_proc + LOG_PSEUDOCOUNT)
      }
    }
    
    mat_proc[!is.finite(mat_proc)] <- 0
    mat_proc[is.na(mat_proc)] <- 0
    
    # 3) Filter AFTER normalization/transform (size factors come from FULL matrix)
    if (!is.null(FILTER_GENE_LIST) && length(FILTER_GENE_LIST) > 0) {
      keep <- rownames(mat_proc) %in% FILTER_GENE_LIST
      mat_proc <- mat_proc[keep, , drop = FALSE]
    }
    
    # 4) Scaling plan (PCAtools will do it via scale=TRUE/FALSE)
    scale_unit <- identical(pp_plan$scale, "yes")
    
    list(
      mat_proc   = mat_proc,
      scale_unit = scale_unit,
      comparable = isTRUE(preproc_compatibility()$compatible)
    )
  }
  
  
  ## -----------------------------
  ## PCA
  ## -----------------------------
  pca_res <- reactive({
    mat <- rna_matrix_full()
    req(mat)
    
    pp <- preprocess_for_pca(mat, dds_obj = rna_dds_obj())
    p <- PCAtools::pca(
      pp$mat_proc,
      metadata  = NULL,
      center    = TRUE,
      scale     = pp$scale_unit,
      removeVar = 0,
      rank      = min(10, ncol(pp$mat_proc) - 1)
    )
    
    # Align sign only if preprocessing is compatible with reference
    if (isTRUE(pp$comparable)) {
      hnscc_pca <- hnscc_pca_reactive()
      aln <- align_pcatools_pc_sign(hnscc_pca, p, pc = "PC1")
      p   <- aln$new_pca_aligned
    }
    
    p
  })
  
  pca_scores <- reactive({
    p <- pca_res()
    scores <- as.data.frame(p$rotated)
    scores$sample <- rownames(scores)
    scores
  })
  
  output$pca_loadings_plot <- shiny::renderPlot({
    req(analysis())
    req(pca_res())
    comp <- preproc_compatibility()
    if (!isTRUE(comp$compatible)) {
      plot.new()
      text(0.5, 0.5, "Loadings correlation suppressed:\npreprocessing not compatible with reference.", cex = 1.2)
      return(invisible())
    }
    plot_pcatools_loadings_correlation(hnscc_pca_reactive(), pca_res(), pc = "PC1")
  })
  
  ## -----------------------------
  ## Clinical data
  ## -----------------------------
  
  clin_raw <- reactive({
    if (is.null(input$include_clin) || input$include_clin != "Yes") return(NULL)
    
    if (isTRUE(input$use_rds_meta)) {
      md <- rds_coldata()
      validate(need(!is.null(md), "No metadata available from the uploaded .rds (colData is missing or RNA input is not an .rds)."))
      return(md)
    }
    
    if (isTRUE(input$use_example_clin)) {
      return(example_clin)
    }
    
    if (!is.null(input$clin_file)) {
      return(readr::read_delim(input$clin_file$datapath, delim = "\t", show_col_types = FALSE))
    }
    
    NULL
  })
  
  
  rds_coldata <- reactive({
    x <- rna_raw()
    if (is.null(x)) return(NULL)
    if (!(inherits(x, "DESeqDataSet") || inherits(x, "SummarizedExperiment"))) return(NULL)
    
    md <- as.data.frame(SummarizedExperiment::colData(x))
    
    # Add an explicit sample ID column (first), without clobbering an existing "sample"
    sid <- if ("sample" %in% names(md)) ".sample_id" else "sample"
    md[[sid]] <- rownames(md)
    
    md <- md[, c(sid, setdiff(names(md), sid)), drop = FALSE]
    md
  })
  
  output$clin_rds_meta_ui <- renderUI({
    md <- rds_coldata()
    if (is.null(md)) return(NULL)
    
    checkboxInput(
      "use_rds_meta",
      "Use sample metadata from uploaded .rds (colData)",
      value = if (!is.null(input$use_rds_meta)) isTRUE(input$use_rds_meta) else TRUE
    )
  })
  
  
  
  output$clin_column_ui <- renderUI({
    df <- clin_raw()
    req(!is.null(df))
    selected_col <- if (!is.null(input$sample_id_col) && input$sample_id_col %in% names(df)) {
      input$sample_id_col
    } else {
      names(df)[1]
    }
    selectInput(
      "sample_id_col",
      "Column with sample IDs:",
      choices  = names(df),
      selected = selected_col
    )
  })
  
  output$hpv_column_ui <- renderUI({
    df <- clin_raw()
    req(!is.null(df))
    
    if (ncol(df) == 2) {
      selectInput(
        "hpv_col",
        "HPV status column:",
        choices  = names(df),
        selected = names(df)[2]
      )
    } else {
      selectInput(
        "hpv_col",
        "HPV status column:",
        choices  = c("Select HPV status column..." = "", names(df)),
        selected = ""
      )
    }
  })
  
  
  hpv_infer <- reactive({
    req(clin_raw(), input$hpv_col)
    col <- clin_raw()[[input$hpv_col]]
    
    infer_binary_group(
      x             = col,
      flip          = FALSE,              # "auto guess" (toggle applies flip)
      drop_unknown  = TRUE
    )
  })
  
  output$hpv_toggle_ui <- renderUI({
    req(hpv_infer())
    inf <- hpv_infer()
    
    neg0 <- as.character(inf$neg)
    pos0 <- as.character(inf$pos)
    
    flip <- isTRUE(hpv_flip_state())
    neg  <- if (flip) pos0 else neg0
    pos  <- if (flip) neg0 else pos0
    
    div(
      class = "hpv-toggle-container",
      div(
        class = "hpv-toggle-column hpv-toggle-negative",
        div(class = "hpv-toggle-label", "Negative"),
        actionButton("hpv_neg_btn", neg, class = "hpv-pill")
      ),
      div(
        class = "hpv-toggle-column hpv-toggle-positive",
        div(class = "hpv-toggle-label", "Positive"),
        actionButton("hpv_pos_btn", pos, class = "hpv-pill")
      )
    )
  })
  
  
  observeEvent(input$hpv_neg_btn, {
    hpv_flip_state(!hpv_flip_state())
  })
  observeEvent(input$hpv_pos_btn, {
    hpv_flip_state(!hpv_flip_state())
  })
  
  output$hpv_mapping_text <- renderText({
    req(clin_raw(), input$hpv_col)
    
    col <- clin_raw()[[input$hpv_col]]
    
    inf <- infer_binary_group(
      col,
      flip = isTRUE(hpv_flip_state()),
      drop_unknown = TRUE
    )
    
    paste0(
      "Detected non-unknown HPV levels: ", paste(inf$observed_levels, collapse = ", "), "\n",
      "Current mapping:\n",
      "  Negative: ", inf$neg, "\n",
      "  Positive: ", inf$pos, "\n",
      "Reason: ", inf$reason
    )
  })
  
  
  clin_ok <- reactive({
    if (is.null(input$include_clin)) {
      return(list(ok = FALSE, messages = "Please choose whether to include clinical data."))
    }
    
    if (input$include_clin == "No") {
      return(list(ok = TRUE, messages = "No clinical data: GMM-only mode will be used."))
    }
    
    msgs <- character()
    df <- clin_raw()
    if (is.null(df)) {
      msgs <- c(msgs, "Please upload clinical data or use the built-in example.")
      return(list(ok = FALSE, messages = msgs))
    }
    
    if (is.null(input$sample_id_col) || !input$sample_id_col %in% names(df)) {
      msgs <- c(msgs, "Select a valid sample ID column.")
    }
    if (is.null(input$hpv_col) || !input$hpv_col %in% names(df)) {
      msgs <- c(msgs, "Select a valid HPV status column.")
    }
    
    if (length(msgs) == 0) {
      hpv <- df[[input$hpv_col]]
      known <- !is_unknown_label(hpv)
      nlev <- length(unique(hpv[known]))
      if (nlev < 2) {
        msgs <- c(msgs, "HPV status column must have at least 2 distinct non-unknown values.")
      }
    }
    
    if (length(msgs) == 0) {
      list(ok = TRUE, messages = "Clinical checks passed.")
    } else {
      list(ok = FALSE, messages = msgs)
    }
  })
  
  observe({
    if (wizard_step() == 3) {
      ck <- clin_ok()
      if (ck$ok) shinyjs::enable("generate_report") else shinyjs::disable("generate_report")
    }
  })
  
  ## -----------------------------
  ## Wizard body UI (3 steps)
  ## -----------------------------
  output$wizard_body <- renderUI({
    if (wizard_step() == 1) {
      tagList(
        h3("Step 1: Upload RNA-seq data"),
        p("Expected format (tab-delimited TSV):"),
        tags$ul(
          tags$li("Rows = genes"),
          tags$li("Columns = samples"),
          tags$li("First column = gene identifiers"),
          tags$li("Remaining columns = numeric expression/count values")
        ),
        div(
          class = "pp-upload-row",
          checkboxInput(
            "use_example_rna",
            "Use built-in example RNA-seq dataset",
            value = isTRUE(input$use_example_rna)
          ),
          fileInput("rna_file", "Upload RNA-seq file or DeseqDataset", accept = c(".tsv", ".txt", ".rds")),
          uiOutput("rna_assay_ui")
        ),
        conditionalPanel(
          condition = "input.use_example_rna === true",
          hr(),
          h4("Example format (first 10 rows)"),
          tableOutput("rna_example")
        ),
        h4("File checks"),
        div(class = "wrap-pre", verbatimTextOutput("rna_checks_text")),
        tags$div(
          style = "margin-top: 15px; text-align: right;",
          actionButton("step1_next", "Next: Preprocessing plan", class = "btn-primary")
        )
      )
    } else if (wizard_step() == 2) {
      tagList(
        h3("Step 2: Preprocessing"),
        tags$p(class="text-muted",
               "Specify what the uploaded matrix already is (click the pills in the Input preprocessing card). ",
               "Then choose what the app should apply (click steps in the Processing plan diagram)."),
        
        
        uiOutput("preproc_cards_ui"),
        
        div(
          style = "display:flex; align-items:center; gap:10px; flex-wrap:wrap; margin-bottom: 6px;",
          actionButton("auto_preproc", "Auto-select plan", class = "btn-outline-primary btn-sm"),
          tags$span(class = "text-muted small",
                    "Uses matrix heuristics to recommend settings; you can override via the cards and plan diagram.")
        ),
        uiOutput("preproc_plan_ui"),
        uiOutput("bucket_plan"),
        
        tags$div(
          style = "margin-top:12px;",
          fluidRow(
            column(6, uiOutput("gene_overlap_ui")),
            column(6, uiOutput("compatibility_ui"))
          )
        ),
        
        tags$div(
          style = "margin-top: 15px; text-align: right;",
          actionButton("step2_back", "Back"),
          actionButton("step2_next", "Next: Clinical options", class = "btn-primary")
        )
      )
      
    } else {
      tagList(
        h3("Step 3: Clinical information (optional)"),
        p("You can choose whether to include clinical metadata in the analysis."),
        fluidRow(
          column(
            6,
            h4("If clinical data is provided"),
            tags$ul(
              tags$li("Runs ROC + Youden, Gaussian (auto), and GMM thresholding on PC1 (using non-unknown HPV values)."),
              tags$li("Computes sensitivity, specificity, and accuracy for each method."),
              tags$li("Appends threshold-based HPV calls to your clinical dataset; unknowns get NA in HPV status but still get call assignments.")
            )
          ),
          column(
            6,
            h4("If clinical data is not provided"),
            tags$ul(
              tags$li("Runs GMM thresholding only (unsupervised)."),
              tags$li("Derives a two-group call based solely on the PC1 distribution."),
              tags$li("Creates a new metadata table with sample IDs, PC1, and GMM call.")
            )
          )
        ),
        hr(),
        radioButtons(
          "include_clin",
          "Include clinical data?",
          choices  = c("No", "Yes"),
          selected = if (!is.null(input$include_clin)) input$include_clin else "No",
          inline   = TRUE
        ),
        conditionalPanel(
          condition = "input.include_clin == 'Yes'",
          
          uiOutput("clin_rds_meta_ui"),
          
          conditionalPanel(
            condition = "!input.use_rds_meta",
            
            checkboxInput(
              "use_example_clin",
              "Use built-in example clinical data",
              value = isTRUE(input$use_example_clin)
            ),
            fileInput("clin_file", "Upload clinical file (TSV)", accept = c(".tsv", ".txt"))
          ),
          
          uiOutput("clin_column_ui"),
          uiOutput("hpv_column_ui"),
          br(),
          h4("Negative / positive mapping (unknown/NA/blank ignored)"),
          uiOutput("hpv_toggle_ui"),
          div(class = "wrap-pre", verbatimTextOutput("hpv_mapping_text"))
        ),
        tags$div(
          style = "margin-top: 15px; text-align: right;",
          actionButton("step3_back", "Back"),
          actionButton("generate_report", "Generate report", class = "btn-primary")
        )
      )
    }
  })
  
  observeEvent(input$step1_next, {
    ck <- rna_ok()
    if (!ck$ok) {
      showNotification("Please fix RNA-seq issues before proceeding.", type = "error")
    } else {
      wizard_step(2)
    }
  })
  
  observeEvent(input$step2_back, {
    wizard_step(1)
  })
  observeEvent(input$step2_next, {
    wizard_step(3)
  })
  observeEvent(input$step3_back, {
    wizard_step(2)
  })
  
  ## -----------------------------
  ## Main analysis
  ## -----------------------------
  analysis <- eventReactive(input$generate_report, {
    scores <- pca_scores()
    pc_df  <- scores %>% dplyr::select(sample, dplyr::starts_with("PC"))
    
    include_clinical <- !is.null(input$include_clin) && input$include_clin == "Yes"
    
    if (include_clinical) {
      clin <- clin_raw()
      validate(need(!is.null(clin), "Clinical data not available."))
      
      sample_col <- input$sample_id_col
      hpv_col    <- input$hpv_col
      
      validate(
        need(sample_col %in% names(clin), "Sample ID column must be valid."),
        need(hpv_col %in% names(clin), "HPV status column must be valid.")
      )
      
      # Ensure sample IDs compare cleanly (avoid factor vs character mismatches)
      pc_df$sample <- as.character(pc_df$sample)
      clin[[sample_col]] <- as.character(clin[[sample_col]])
      
      # Guard against duplicates (common source of ordering mismatch)
      if (anyDuplicated(pc_df$sample)) {
        dup <- unique(pc_df$sample[duplicated(pc_df$sample)])
        validate(need(FALSE, paste0("Duplicate sample IDs in PCA scores: ", paste(dup, collapse = ", "))))
      }
      if (anyDuplicated(clin[[sample_col]])) {
        dup <- unique(clin[[sample_col]][duplicated(clin[[sample_col]])])
        validate(need(FALSE, paste0("Duplicate sample IDs in clinical sample ID column: ", paste(dup, collapse = ", "))))
      }
      
      common_ids <- intersect(pc_df$sample, clin[[sample_col]])
      validate(need(length(common_ids) >= 2, "Need overlapping samples between RNA-seq and clinical."))
      
      # Align deterministically by common_ids ordering
      pc_sub <- pc_df %>%
        dplyr::filter(sample %in% common_ids) %>%
        dplyr::arrange(sample)
      
      clin_sub <- clin %>%
        dplyr::filter(.data[[sample_col]] %in% common_ids) %>%
        dplyr::arrange(.data[[sample_col]])
      
      validate(need(all(pc_sub$sample == clin_sub[[sample_col]]),
                    "Sample ordering mismatch after alignment. Check duplicated sample IDs or hidden whitespace."))
      
      hpv_vec <- clin_sub[[hpv_col]]
      
      # Known HPV indices (for training)
      known_idx <- !is_unknown_label(hpv_vec)
      validate(need(any(known_idx), "No non-unknown HPV values available for thresholding."))
      
      # Infer neg/pos from knowns + apply user flip, return ready-to-use factor
      inf <- infer_binary_group(
        x            = hpv_vec,
        flip         = isTRUE(hpv_flip_state()),
        drop_unknown = TRUE
      )
      
      group_train <- inf$train_group
      validate(need(is.factor(group_train) && nlevels(group_train) == 2L,
                    "HPV status must resolve to exactly 2 levels after inference."))
      
      # Subset PC scores to the same known subset
      pc_thr <- pc_sub[known_idx, , drop = FALSE]
      validate(need(nrow(pc_thr) == length(group_train),
                    "Internal mismatch: training samples and training groups differ in length."))
      
      # Training vector: only known HPV
      x_pc1_train <- pc_thr$PC1
      names(x_pc1_train) <- pc_thr$sample
      # Full vector: all overlapping samples (including unknown)
      x_pc1_all <- pc_sub$PC1
      names(x_pc1_all) <- pc_sub$sample
      
      # Thresholds learned on known subset
      thr_res <- determine_threshold_two_groups(
        x              = x_pc1_train,
        group          = group_train,  # ordered as (neg, pos)
        methods        = c("youden", "gaussian_auto", "gmm"),
        var_label      = "PC1 score",
        digits         = 3,
        group_colors   = NULL,
        method_palette = "Dark2"
      )
      
      # Calls for ALL overlapping samples (including unknowns) using learned thresholds
      call_all_df <- threshold_predict_all(
        x_all       = x_pc1_all,
        x_train     = x_pc1_train,
        group_train = group_train,
        res_list    = thr_res$results
      )
      
      pc_sub_for_join <- pc_sub %>%
        { if (!"PC2" %in% names(.)) dplyr::mutate(., PC2 = NA_real_) else . } %>%
        dplyr::select(sample, PC1, PC2)
      
      # Build augmented clinical table
      clin_aug <- clin_sub %>%
        dplyr::left_join(
          pc_sub_for_join,
          by = setNames("sample", sample_col)
        ) %>%
        dplyr::left_join(
          call_all_df,
          by = setNames("sample", sample_col)
        )
      
      # HPV status for plotting, with unknowns explicitly labeled "Unknown"
      hpv_vec_full <- clin_aug[[hpv_col]]
      hpv_plot <- as.character(hpv_vec_full)
      hpv_plot[is_unknown_label(hpv_vec_full)] <- "Unknown"
      
      plot_df <- clin_aug
      plot_df$hpv_plot <- factor(hpv_plot)  # forces discrete palette (prevents 0/1 from becoming continuous/black)
      
      # Color options for PCA:
      color_opts <- c("HPV status" = "hpv_plot")
      call_cols <- grep("^call_", names(plot_df), value = TRUE)
      for (cc in call_cols) {
        method_name <- sub("^call_", "", cc)
        label <- paste0("Call: ", method_name)
        color_opts[label] <- cc
      }
      
      list(
        mode              = "supervised",
        pca_scores        = pc_df,
        threshold_res     = thr_res,
        threshold_plot    = thr_res$plot,
        pca_data          = plot_df,
        pca_color_options = color_opts,
        clin_augmented    = clin_aug,
        meta_unsupervised = NULL,
        gmm_threshold     = NULL,
        gmm_summary_table = NULL,
        gmm_summary_gt    = NULL,
        hpv_col_name      = hpv_col
      )
      
    } else {
      # Unsupervised GMM-only branch
      pc_all <- pc_df %>% dplyr::arrange(sample)
      x_pc1  <- pc_all$PC1
      names(x_pc1) <- pc_all$sample
      
      gmm_res <- threshold_gmm(x_pc1, group = NULL)
      gmm_thr <- gmm_res$threshold
      fit     <- gmm_res$mclust_fit
      
      means    <- gmm_res$params$means
      pos_comp <- which.max(means)
      
      comp_idx <- apply(fit$z, 1, which.max)
      
      gmm_call <- ifelse(comp_idx == pos_comp, "GMM_high", "GMM_low")
      gmm_call <- factor(gmm_call, levels = c("GMM_low", "GMM_high"))
      
      meta_df <- tibble::tibble(
        sample   = names(x_pc1),
        PC1      = as.numeric(x_pc1),
        gmm_call = gmm_call
      )
      
      summary_raw <- tibble::tibble(
        method      = "GMM (unsupervised)",
        threshold   = gmm_thr,
        accuracy    = NA_real_,
        sensitivity = NA_real_,
        specificity = NA_real_
      )
      
      summary_table <- format_threshold_summary(
        summary_raw,
        var_label         = "PC1 score",
        digits            = 3,
        include_confusion = FALSE
      )
      
      dens_plot <- plot_thresholds(
        x              = x_pc1,
        group          = gmm_call,
        summary_table  = summary_table,
        var_label      = "PC1 score",
        group_colors   = NULL,
        method_palette = NULL
      )
      
      summary_gt <- threshold_gt_table(
        summary_table,
        digits   = 3,
        subtitle = NULL
      )
      
      plot_df <- pc_all %>%
        dplyr::left_join(
          meta_df %>% dplyr::select(sample, gmm_call),
          by = "sample"
        )
      
      color_opts <- c("GMM call" = "gmm_call")
      
      list(
        mode              = "unsupervised",
        pca_scores        = pc_all,
        threshold_res     = NULL,
        threshold_plot    = dens_plot,
        pca_data          = plot_df,
        pca_color_options = color_opts,
        clin_augmented    = NULL,
        meta_unsupervised = meta_df,
        gmm_threshold     = gmm_thr,
        gmm_summary_table = summary_table,
        gmm_summary_gt    = summary_gt,
        hpv_col_name      = NULL
      )
    }
  })
  
  observeEvent(input$generate_report, {
    removeModal()
  })
  
  ## -----------------------------
  ## PCA color selector + plot
  ## -----------------------------
  output$pca_color_ui <- renderUI({
    req(analysis())
    res  <- analysis()
    opts <- res$pca_color_options
    validate(need(length(opts) > 0, "No coloring options available."))
    
    # Use the metadata table the user sees (supervised: clin_augmented; unsupervised: meta_unsupervised)
    meta_df <- if (identical(res$mode, "supervised")) res$clin_augmented else res$meta_unsupervised
    if (is.null(meta_df)) meta_df <- res$pca_data
    
    # Label choices from metadata columns
    label_choices <- c("None" = "", stats::setNames(names(meta_df), names(meta_df)))
    
    sel <- isolate(input$pca_label_by)
    if (is.null(sel) || !(sel %in% c("", names(meta_df)))) sel <- ""
    
    bslib::accordion(
      id = "pca_controls_acc",
      bslib::accordion_panel(
        title = "PCA aesthetics",
        radioButtons(
          "pca_color_by",
          "Color PCA by:",
          choices  = names(opts),
          selected = if (!is.null(input$pca_color_by) && input$pca_color_by %in% names(opts)) input$pca_color_by else names(opts)[1],
          inline   = TRUE
        ),
        selectInput(
          "pca_label_by",
          "Label points by:",
          choices  = label_choices,
          selected = sel
        ),
        tags$p(
          class = "text-muted small",
          "Labels can get crowded for large cohorts—set to 'None' to hide labels."
        )
      ),
      open = "PCA aesthetics"
    )
  })
  
  
  output$pca_plot <- renderPlot({
    req(analysis())
    res  <- analysis()
    df   <- res$pca_data
    opts <- res$pca_color_options
    validate(need(!is.null(df), "PCA data not available."))
    validate(need(!is.null(opts) && length(opts) > 0, "No PCA color options defined."))
    
    label <- if (!is.null(input$pca_color_by) && input$pca_color_by %in% names(opts)) {
      input$pca_color_by
    } else {
      names(opts)[1]
    }
    col_var <- opts[[label]]
    
    label_var <- NULL
    if (!is.null(input$pca_label_by) && nzchar(input$pca_label_by) && input$pca_label_by %in% names(df)) {
      label_var <- input$pca_label_by
    }
    
    plot_pca_colored(df, col_var, label, label_var = label_var)
  })
  
  
  ## -----------------------------
  ## Threshold plot and table
  ## -----------------------------
  output$threshold_plot <- renderPlot({
    req(analysis())
    analysis()$threshold_plot
  })
  
  output$threshold_gt <- render_gt({
    req(analysis())
    res <- analysis()
    if (res$mode == "supervised") {
      req(res$threshold_res$summary_gt)
      res$threshold_res$summary_gt
    } else {
      if (!is.null(res$gmm_summary_gt)) {
        res$gmm_summary_gt
      } else {
        gt::gt(
          tibble(
            Method    = "GMM (unsupervised)",
            Threshold = signif(res$gmm_threshold, 3)
          )
        ) |>
          gt::fmt_number(columns = "Threshold", decimals = 3)
      }
    }
  })
  
  youden_res <- reactive({
    req(analysis())
    res <- analysis()
    if (!identical(res$mode, "supervised")) return(NULL)
    
    y <- res$threshold_res$results[["youden"]]
    if (is.null(y)) return(NULL)
    
    y
  })
  
  output$youden_diag_ui <- renderUI({
    y <- youden_res()
    if (is.null(y)) {
      return(tags$div(class = "text-muted", "ROC diagnostics available only in supervised mode (clinical labels required)."))
    }
    
    ci <- y$auc_ci
    ci_txt <- if (!is.null(ci) && length(ci) == 3 && all(is.finite(ci))) {
      sprintf("[%.3f, %.3f]", ci[1], ci[3])
    } else {
      "NA"
    }
    
    tags$div(
      style = "margin-bottom: 8px;",
      tags$div(tags$strong("AUC: "), sprintf("%.3f", y$auc)),
      tags$div(tags$strong("AUC 95% CI: "), ci_txt),
      tags$div(tags$strong("DeLong p-value (AUC>0.5): "), if (is.finite(y$p_value)) sprintf("%.3g", y$p_value) else "NA"),
      tags$div(tags$strong("Youden threshold: "), sprintf("%.3f", y$threshold)),
      tags$div(tags$strong("Sensitivity / Specificity: "),
               sprintf("%.3f / %.3f", y$sensitivity, y$specificity)),
      tags$div(tags$strong("Youden J: "), sprintf("%.3f", y$youden_J))
    )
  })
  
  output$youden_roc_plot <- renderPlot({
    y <- youden_res()
    req(y)
    print(y$roc_plot)
  })
  
  
  ## -----------------------------
  ## Metadata preview modal
  ## -----------------------------
  output$meta_preview_dt <- renderDT({
    req(analysis())
    res <- analysis()
    df <- if (res$mode == "supervised") res$clin_augmented else res$meta_unsupervised
    datatable(
      df,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  observeEvent(input$show_meta_btn, {
    showModal(modalDialog(
      title = "Preview of augmented metadata",
      size  = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),
      DTOutput("meta_preview_dt")
    ))
  })
  
  ## -----------------------------
  ## Download ZIP report (with preprocessing summary)
  ## -----------------------------
  output$download_zip <- downloadHandler(
    filename = function() {
      paste0("threshold_report_", Sys.Date(), ".zip")
    },
    content = function(file) {
      req(analysis())
      res <- analysis()
      
      rep_dir <- file.path(tempdir(), paste0("threshold_report_", format(Sys.time(), "%Y%m%d_%H%M%S")))
      dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
      
      
      # 1. PCA scores
      readr::write_csv(res$pca_scores, file.path(rep_dir, "pca_scores.csv"))
      
      # 2. Metadata + threshold summary + confusion tables
      if (res$mode == "supervised") {
        readr::write_csv(res$clin_augmented, file.path(rep_dir, "clinical_with_calls.csv"))
        readr::write_csv(res$threshold_res$summary_table, file.path(rep_dir, "threshold_summary.csv"))
        
        # Confusion tables per method
        conf_list <- lapply(res$threshold_res$results, function(r) {
          ev <- r$evaluation
          if (is.null(ev) || is.null(ev$confusion)) return(NULL)
          
          tab <- ev$confusion
          method_name <- if (!is.null(ev$method) && !is.na(ev$method)) ev$method else r$method
          
          # Convert to a plain data.frame/tibble
          conf <- if (is.data.frame(tab)) {
            tibble::as_tibble(tab)
          } else {
            tibble::as_tibble(as.data.frame(tab, stringsAsFactors = FALSE))
          }
          
          # Standardize column names:
          # First two columns are the dimension labels (whatever they're currently called)
          if (ncol(conf) < 3) return(NULL)
          
          names(conf)[1:2] <- c("Observed", "Predicted")
          
          # Third column is the count; might be Freq/n/Count/etc.
          cnt_col <- setdiff(names(conf), c("Observed", "Predicted"))[1]
          names(conf)[names(conf) == cnt_col] <- "Count"
          
          conf %>%
            dplyr::mutate(Method = method_name, .before = 1)
        })
        
        conf_df <- bind_rows(conf_list)
        if (nrow(conf_df) > 0) {
          readr::write_csv(conf_df, file.path(rep_dir, "confusion_tables.csv"))
        }
        
        if (!is.null(res$threshold_res$summary_gt)) {
          gt::gtsave(
            res$threshold_res$summary_gt,
            filename = file.path(rep_dir, "threshold_summary.html")
          )
        }
      } else {
        readr::write_csv(res$meta_unsupervised, file.path(rep_dir, "meta_gmm_calls.csv"))
        if (!is.null(res$gmm_summary_table)) {
          readr::write_csv(res$gmm_summary_table, file.path(rep_dir, "gmm_threshold_summary.csv"))
        }
        if (!is.null(res$gmm_summary_gt)) {
          gt::gtsave(
            res$gmm_summary_gt,
            filename = file.path(rep_dir, "gmm_threshold_summary.html")
          )
        }
      }
      
      # 3. Preprocessing summary text (instead of showing overview in the report UI)
      eff    <- preproc_effective()
      compat <- preproc_compatibility()
      go     <- gene_overlap()
      
      preproc_lines <- c(
        "Preprocessing summary for PCA",
        "====================================",
        "",
        sprintf("Normalization entering PCA: %s", eff$normalization),
        sprintf("Transformation entering PCA: %s", eff$log_transform),
        sprintf("Scaled/z-scored entering PCA: %s", eff$scaled),
        "",
        "Compatibility with reference PCA:",
        sprintf("Status: %s", compat$status),
        sprintf("Reason: %s", compat$reason),
        sprintf("Recommendation: %s", compat$recommendation),
        ""
      )
      
      if (!is.null(go)) {
        preproc_lines <- c(
          preproc_lines,
          "Gene overlap with reference signature:",
          sprintf("Reference genes: %d", go$n_reference),
          sprintf("Genes in uploaded matrix: %d", go$n_input),
          sprintf("Matched genes: %d", go$n_matched),
          if (!is.na(go$pct_matched_ref)) sprintf("Percent of reference present: %.1f%%", go$pct_matched_ref) else NULL
        )
      }
      
      writeLines(preproc_lines, file.path(rep_dir, "preprocessing_summary.txt"))
      
      # 4. PCA plots (all coloring options), as SVG, larger size
      df   <- res$pca_data
      opts <- res$pca_color_options
      
      if (!is.null(df) && !is.null(opts) && length(opts) > 0) {
        
        lvar <- if (!is.null(input$pca_label_by) && nzchar(input$pca_label_by)) input$pca_label_by else NULL
        
        for (lab in names(opts)) {
          var <- opts[[lab]]
          p   <- plot_pca_colored(df, var, lab, label_var = lvar)
          safe <- gsub("[^A-Za-z0-9_]+", "_", lab)
          ggplot2::ggsave(
            filename = file.path(rep_dir, paste0("pca_", safe, ".svg")),
            plot     = p,
            width    = 9, height = 6, dpi = 300,
            device   = "svg"
          )
        }
      }
      
      # 5. Threshold plot, as SVG, larger size
      ggplot2::ggsave(
        filename = file.path(rep_dir, "threshold_plot.svg"),
        plot     = res$threshold_plot,
        width    = 9, height = 6, dpi = 300,
        device   = "svg"
      )
      
      # 5b. PCA loadings plot (only if comparable; mirrors UI logic)
      comp <- preproc_compatibility()
      if (isTRUE(comp$compatible)) {
        try({
          grDevices::svg(file.path(rep_dir, "pca_loadings_PC1.svg"), width = 9, height = 6)
          print(plot_pcatools_loadings_correlation(hnscc_pca_reactive(), pca_res(), pc = "PC1"))
          grDevices::dev.off()
        }, silent = TRUE)
      }
      
      # 5c. ROC plot (Youden) in supervised mode
      if (identical(res$mode, "supervised")) {
        y <- res$threshold_res$results[["youden"]]
        if (!is.null(y) && !is.null(y$roc_plot)) {
          try({
            ggplot2::ggsave(
              filename = file.path(rep_dir, "youden_roc.svg"),
              plot     = y$roc_plot,
              width    = 9, height = 6, dpi = 300,
              device   = "svg"
            )
          }, silent = TRUE)
        }
      }
      
      
      # 6. Zip everything
      owd <- setwd(rep_dir)
      on.exit(setwd(owd), add = TRUE)
      files_to_zip <- list.files(rep_dir, full.names = TRUE, recursive = FALSE)
      zip::zipr(
        zipfile = file,
        files   = files_to_zip,
        root    = rep_dir
      )
    }
  )
}

shinyApp(ui, server)
