library(SummarizedExperiment)
library(dplyr)
library(mclust)

#load app functions
source("HPV_BioSig_Functions.R")

#load tcga data
tcga <- readRDS("data/tcga_primary_hnc.rds")

output_dir <- file.path(
  "Results",
  format(Sys.time(), "%Y%m%d_%H%M%S")
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

folds <- 5L
repetitions <- 10L
random_seed <- 20261001L

labels <- tolower(trimws(as.character(colData(tcga)$hpv_status)))
keep <- labels %in% c("negative", "positive")
x <- as.matrix(assay(tcga))[, keep, drop = FALSE]
y <- factor(labels[keep], levels = c("negative", "positive"))

write_result <- function(data, filename) {
  write.csv(data, file.path(output_dir, filename), row.names = FALSE)
}
write_result(as.data.frame(table(HPV_status = labels, useNA = "ifany")),
             "label_counts.csv")

fit_signature <- function(expression, status) {
  genes <- matrixTests::row_t_welch(
    expression[, status == "positive", drop = FALSE],
    expression[, status == "negative", drop = FALSE]
  )
  genes$gene <- rownames(expression)
  genes <- genes |>
    filter(is.finite(pvalue)) |>
    mutate(fdr = p.adjust(pvalue, method = "BH")) |>
    filter(fdr < 0.05, abs(mean.diff) > 1, pmax(mean.x, mean.y) > 2) |>
    arrange(pvalue)
  genes <- genes[seq_len(100), , drop = FALSE]
  
  expression <- expression[genes$gene, , drop = FALSE]
  pca <- PCAtools::pca(expression, center = TRUE, scale = FALSE, rank = 2)
  gene_order <- rownames(pca$loadings)
  centers <- rowMeans(expression[gene_order, , drop = FALSE])
  weights <- pca$loadings[gene_order, "PC1"]
  scores <- as.numeric(pca$rotated[colnames(expression), "PC1"])
  
  if (mean(scores[status == "positive"]) < mean(scores[status == "negative"])) {
    weights <- -weights
    scores <- -scores
  }
  
  roc <- pROC::roc(status, scores, levels = levels(status),
                   direction = "<", quiet = TRUE)
  youden <- as.numeric(unlist(pROC::coords(
    roc, "best", best.method = "youden", ret = "threshold", transpose = FALSE
  )))
  cutoffs <- c(
    Youden = youden[is.finite(youden)][1],
    Gaussian = threshold_gaussian_auto(scores, status)$threshold,
    GMM = threshold_gmm(scores, group = NULL,
                        method = "equal_component", fallback = "none")$threshold
  )
  list(genes = gene_order, centers = centers, weights = weights,
       cutoffs = cutoffs, selected_genes = genes)
}

predict_signature <- function(model, expression) {
  centered <- sweep(expression[model$genes, , drop = FALSE],
                    1, model$centers, "-")
  scores <- as.numeric(crossprod(model$weights, centered))
  predictions <- data.frame(Sample = colnames(expression), PC1 = scores)
  for (method in names(model$cutoffs)) {
    predictions[[method]] <- ifelse(
      scores >= model$cutoffs[[method]], "positive", "negative"
    )
  }
  predictions
}

performance <- function(observed, predictions) {
  auc <- as.numeric(pROC::auc(pROC::roc(
    observed, predictions$PC1, levels = c("negative", "positive"),
    direction = "<", quiet = TRUE
  )))
  bind_rows(lapply(c("Youden", "Gaussian", "GMM"), function(method) {
    called <- predictions[[method]]
    tp <- sum(observed == "positive" & called == "positive")
    tn <- sum(observed == "negative" & called == "negative")
    fp <- sum(observed == "negative" & called == "positive")
    fn <- sum(observed == "positive" & called == "negative")
    data.frame(Method = method, N = length(observed),
               Positive = tp + fn, Negative = tn + fp,
               TP = tp, TN = tn, FP = fp, FN = fn, AUC = auc,
               Sensitivity = tp / (tp + fn), Specificity = tn / (tn + fp),
               Accuracy = (tp + tn) / length(observed))
  }))
}

set.seed(random_seed)
full_model <- fit_signature(x, y)

#Validation that the models match with the functions used. 

reference_pca <- readRDS(
  "data/hpv_signature_tcga_hnscc_pca_model.rds"
)

full_predictions <- predict_signature(full_model, x)

# Match the original scores to the same sample order
reference_scores <- reference_pca$rotated$PC1[
  match(full_predictions$Sample, rownames(reference_pca$rotated))
]

# Allow either PCA sign orientation
max_score_difference <- min(
  max(abs(full_predictions$PC1 - reference_scores)),
  max(abs(full_predictions$PC1 + reference_scores))
)

pca_check_passed <-
  setequal(full_model$genes, rownames(reference_pca$loadings)) &&
  setequal(full_predictions$Sample, rownames(reference_pca$rotated)) &&
  is.finite(max_score_difference) &&
  max_score_difference <= 1e-6

cat(
  "\nFull-cohort PCA reproduction check:",
  if (pca_check_passed) "PASS" else "FAIL",
  "\nMaximum score difference after allowing sign reversal:",
  format(max_score_difference, scientific = TRUE),
  "\n\n"
)

apparent <- performance(y, predict_signature(full_model, x))
write_result(apparent, "apparent_performance.csv")
write_result(full_model$selected_genes, "full_cohort_selected_genes.csv")


# Full-cohort HNSCC signature boxplot
fig.tcga_hnscc_hpv_pca_boxplot <- ggpubr::ggboxplot(
  data = data.frame(
    HPVSIG = full_predictions$PC1,
    hpv_status = y
  ),
  x = "hpv_status",
  y = "HPVSIG",
  add = "jitter",
  fill = "hpv_status",
  bxp.errorbar = TRUE,
  xlab = "HPV Status",
  ylab = "HPV Signature",
  alpha = 0.4
) +
  ggplot2::scale_fill_manual(
    values = c(
      "negative" = "#0072B5FF",
      "positive" = "#BC3C29FF"
    )
  ) +
  ggpubr::stat_compare_means(
    method = "wilcox.test"
  )

print(fig.tcga_hnscc_hpv_pca_boxplot)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "fig-tcga_hnscc_hpv_boxplot.pdf"
  ),
  plot = fig.tcga_hnscc_hpv_pca_boxplot,
  width = 5,
  height = 6,
  units = "in"
)

set.seed(random_seed)
partitions <- lapply(seq_len(repetitions), function(i) {
  assignment <- integer(length(y))
  for (group in levels(y)) {
    indices <- which(y == group)
    shuffled <- sample(indices, length(indices), replace = FALSE)
    assignment[shuffled] <- rep(seq_len(folds), length.out = length(shuffled))
  }
  assignment
})

fold_results <- predictions <- selected_genes <- list()
for (r in seq_len(repetitions)) {
  for (f in seq_len(folds)) {
    key <- paste(r, f, sep = "_")
    test <- which(partitions[[r]] == f)
    train <- which(partitions[[r]] != f)
    set.seed(random_seed + r * folds + f)
    
    model <- fit_signature(x[, train, drop = FALSE], y[train])
    predicted <- predict_signature(model, x[, test, drop = FALSE])
    metrics <- performance(y[test], predicted)
    metrics$Threshold <- unname(model$cutoffs[metrics$Method])
    
    fold_results[[key]] <- cbind(Repetition = r, Fold = f, metrics)
    predictions[[key]] <- cbind(
      Repetition = r, Fold = f, Observed = as.character(y[test]), predicted
    )
    selected_genes[[key]] <- cbind(Repetition = r, Fold = f, model$selected_genes)
  }
}

fold_metrics <- bind_rows(fold_results)
write_result(fold_metrics, "fold_performance.csv")
write_result(bind_rows(predictions), "held_out_predictions.csv")
write_result(bind_rows(selected_genes), "fold_selected_genes.csv")

repeat_metrics <- fold_metrics |>
  group_by(Repetition, Method) |>
  summarise(
    N = sum(N),
    AUC = weighted.mean(AUC, Positive * Negative),
    Sensitivity = sum(TP) / (sum(TP) + sum(FN)),
    Specificity = sum(TN) / (sum(TN) + sum(FP)),
    Accuracy = (sum(TP) + sum(TN)) / (sum(TP) + sum(TN) + sum(FP) + sum(FN)),
    .groups = "drop"
  )
write_result(repeat_metrics, "repeat_performance.csv")

performance_summary <- bind_rows(lapply(c("Youden", "Gaussian", "GMM"), function(method) {
  bind_rows(lapply(c("AUC", "Sensitivity", "Specificity", "Accuracy"), function(metric) {
    values <- repeat_metrics[repeat_metrics$Method == method, ][[metric]]
    data.frame(Method = method, Metric = metric,
               Apparent = apparent[[metric]][apparent$Method == method],
               CV_mean = mean(values), SD_across_repetitions = sd(values))
  }))
}))
write_result(performance_summary, "performance_summary.csv")
print(performance_summary, row.names = FALSE)
