#' Cluster cells in metabolic PCA space
#'
#' @description
#' scMetaTraj_cluster() identifies metabolic subclusters by constructing
#' a kNN graph in metabolic PCA space and applying community detection.
#'
#' IMPORTANT DESIGN PRINCIPLES:
#' - Clustering is performed ONLY in metabolic PCA space.
#' - UMAP coordinates must NEVER be used for clustering.
#' - Results are independent of transcriptomic clustering.
#'
#' @param embedding Numeric matrix (cells x PCs).
#'   Output of scMetaTraj_embed(method = "PCA").
#' @param k Integer. Number of nearest neighbors for kNN graph.
#' @param resolution Numeric. Resolution parameter for clustering
#'   (used for Leiden only).
#' @param method Character. "leiden" (default) or "louvain".
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A factor of length equal to number of cells,
#'   giving metabolic cluster labels per cell.
#'
#' @export
scMetaTraj_cluster <- function(
    embedding,
    k = 20,
    resolution = 0.5,
    method = c("leiden", "louvain"),
    seed = 123
) {
  
  method <- match.arg(method)
  
  # -----------------------------
  # Step 0: Input validation
  # -----------------------------
  if (!is.matrix(embedding)) {
    stop("embedding must be a numeric matrix (cells x PCs).")
  }
  
  if (nrow(embedding) <= k) {
    stop("Number of cells must be greater than k.")
  }
  
  if (is.null(rownames(embedding))) {
    stop("embedding must have rownames (cell IDs).")
  }
  
  set.seed(seed)
  
  # -----------------------------
  # Step 1: Build kNN graph
  # -----------------------------
  nn <- RANN::nn2(
    data = embedding,
    k = k + 1   # include self, removed below
  )
  
  nn_idx <- nn$nn.idx[, -1, drop = FALSE]
  
  # -----------------------------
  # Step 2: Convert kNN to igraph
  # -----------------------------
  edges <- cbind(
    rep(seq_len(nrow(nn_idx)), each = ncol(nn_idx)),
    as.vector(nn_idx)
  )
  
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  
  # -----------------------------
  # Step 3: Community detection
  # -----------------------------
  if (method == "leiden") {
    
    # CRAN-safe Leiden implementation (igraph >= 1.3)
    cl <- igraph::cluster_leiden(
      g,
      resolution_parameter = resolution
    )$membership
    
  } else {
    
    cl <- igraph::cluster_louvain(g)$membership
    
  }
  
  clusters <- factor(cl)
  names(clusters) <- rownames(embedding)
  
  clusters
}
