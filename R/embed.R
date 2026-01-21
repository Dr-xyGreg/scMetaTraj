#' Embed cells in metabolic feature space
#'
#' @description
#' scMetaTraj_embed() constructs a low-dimensional representation
#' of cells based on pathway-level metabolic scores.
#'
#' DESIGN PRINCIPLES:
#' - PCA is the true analysis space (for graph construction).
#' - UMAP is ONLY for visualization.
#'
#' @param scores Numeric matrix, cells x pathways.
#' @param method Character. "PCA" (default) or "UMAP".
#' @param n_pcs Integer. Number of PCs to return / use.
#' @param umap_n_neighbors Integer. UMAP n_neighbors.
#' @param umap_min_dist Numeric. UMAP min_dist.
#' @param seed Integer. Random seed.
#'
#' @return A numeric matrix:
#' - PCA: cells x PCs
#' - UMAP: cells x 2
#'
scMetaTraj_embed <- function(
    scores,
    method = c("PCA", "UMAP"),
    n_pcs = 10,
    umap_n_neighbors = 30,
    umap_min_dist = 0.3,
    seed = 123
) {
  method <- match.arg(method)
  
  # -----------------------------
  # Step 0: Input validation
  # -----------------------------
  if (!is.matrix(scores)) {
    stop("scores must be a numeric matrix (cells x pathways).")
  }
  if (nrow(scores) < 2 || ncol(scores) < 2) {
    stop("scores must have at least 2 cells and 2 pathways.")
  }
  
  max_pcs <- min(ncol(scores), nrow(scores) - 1)
  if (n_pcs > max_pcs) {
    warning("n_pcs too large, reduced to ", max_pcs)
    n_pcs <- max_pcs
  }
  
  # -----------------------------
  # Step 1: PCA (analysis space)
  # -----------------------------
  set.seed(seed)
  pca <- stats::prcomp(
    scores,
    center = FALSE,
    scale. = FALSE
  )
  
  pca_embed <- pca$x[, seq_len(n_pcs), drop = FALSE]
  colnames(pca_embed) <- paste0("PC_", seq_len(n_pcs))
  rownames(pca_embed) <- rownames(scores)
  
  if (method == "PCA") {
    return(pca_embed)
  }
  
  # -----------------------------
  # Step 2: UMAP (visualization)
  # -----------------------------
  if (!requireNamespace("uwot", quietly = TRUE)) {
    stop("Package 'uwot' is required for UMAP.")
  }
  
  set.seed(seed)
  umap_embed <- uwot::umap(
    X = pca_embed,
    n_neighbors = umap_n_neighbors,
    min_dist = umap_min_dist,
    n_components = 2,
    metric = "euclidean",
    verbose = FALSE
  )
  
  colnames(umap_embed) <- c("UMAP_1", "UMAP_2")
  rownames(umap_embed) <- rownames(scores)
  
  umap_embed
}
