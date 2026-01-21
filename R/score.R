#' Calculate metabolic pathway scores for single cells
#'
#' @description
#' scMetaTraj_score() maps gene-level expression to pathway-level
#' metabolic activity scores. The resulting matrix defines the
#' metabolic feature space used by all downstream scMetaTraj modules.
#'
#' IMPORTANT:
#' - Scores represent relative metabolic states, NOT metabolic flux.
#' - Designed to be robust to dropout in scRNA-seq data.
#'
#' @param x A Seurat object or a gene x cell expression matrix.
#' @param gene_sets A named list: pathway -> character vector of genes.
#' @param assay Character. Seurat assay to use. Default "RNA".
#' @param slot Character. Expression slot. Default "data".
#' @param method Character. Scoring method: "mean" or "zscore".
#' @param min_genes Integer. Minimal number of detected genes per pathway.
#' @param scale Logical. Whether to z-score pathway scores across cells.
#'
#' @return A numeric matrix: cells x pathways.
#'
scMetaTraj_score <- function(
    x,
    gene_sets,
    assay = "RNA",
    slot = "data",
    method = c("mean", "zscore"),
    min_genes = 3,
    scale = TRUE
) {
  method <- match.arg(method)
  
  # -----------------------------
  # Step 1: Extract expression
  # -----------------------------
  if (inherits(x, "Seurat")) {
    expr <- Seurat::GetAssayData(
      x,
      assay = assay,
      slot = slot
    )
  } else if (is.matrix(x)) {
    expr <- x
  } else {
    stop("x must be a Seurat object or a gene x cell matrix.")
  }
  
  genes_available <- rownames(expr)
  cells <- colnames(expr)
  
  # -----------------------------
  # Step 2: Score pathways
  # -----------------------------
  scores <- lapply(names(gene_sets), function(pw) {
    genes <- gene_sets[[pw]]
    genes_use <- intersect(genes, genes_available)
    
    # Drop unstable pathways
    if (length(genes_use) < min_genes) {
      return(rep(NA_real_, length(cells)))
    }
    
    sub_expr <- expr[genes_use, , drop = FALSE]
    
    if (method == "mean") {
      colMeans(sub_expr)
    } else {
      z <- t(scale(t(sub_expr)))
      colMeans(z, na.rm = TRUE)
    }
  })
  
  score_mat <- do.call(rbind, scores)
  rownames(score_mat) <- names(gene_sets)
  colnames(score_mat) <- cells
  
  # Transpose to cells x pathways
  score_mat <- t(score_mat)
  
  # -----------------------------
  # Step 3: Remove empty pathways
  # -----------------------------
  valid <- colSums(!is.na(score_mat)) > 0
  score_mat <- score_mat[, valid, drop = FALSE]
  
  # -----------------------------
  # Step 4: Optional scaling
  # -----------------------------
  if (scale) {
    score_mat <- scale(score_mat)
  }
  
  score_mat
}
