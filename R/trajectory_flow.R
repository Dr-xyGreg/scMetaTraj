#' Compute local directional consistency along mPT
#'
#' @description
#' Estimate local direction vectors pointing toward neighbors with higher
#' metabolic pseudotime. Intended for visualization only.
#'
#' @param emb_pca Matrix (cells x PCs) used for neighborhood definition.
#' @param emb_umap Matrix (cells x 2) used only for visualization.
#' @param pseudotime Numeric vector of mPT.
#' @param k Integer. Number of nearest neighbors.
#' @param min_delta Minimum mPT difference to consider a neighbor "forward".
#'
#' @return Data.frame with UMAP coordinates and dx, dy vectors.
#'
#' @export
scMetaTraj_flow <- function(
    emb_pca,
    emb_umap,
    pseudotime,
    k = 15,
    min_delta = 0.02
) {
  emb_pca <- as.matrix(emb_pca)
  emb_umap <- as.matrix(emb_umap)
  
  n <- nrow(emb_pca)
  D <- as.matrix(stats::dist(emb_pca))
  
  dx <- dy <- numeric(n)
  
  for (i in seq_len(n)) {
    nn <- order(D[i, ])[2:(k + 1)]
    forward <- nn[pseudotime[nn] - pseudotime[i] > min_delta]
    
    if (length(forward) < 2) {
      dx[i] <- dy[i] <- NA
      next
    }
    
    w <- pseudotime[forward] - pseudotime[i]
    center <- colSums(emb_umap[forward, , drop = FALSE] * w) / sum(w)
    
    dx[i] <- center[1] - emb_umap[i, 1]
    dy[i] <- center[2] - emb_umap[i, 2]
  }
  
  data.frame(
    UMAP_1 = emb_umap[, 1],
    UMAP_2 = emb_umap[, 2],
    dx = dx,
    dy = dy
  )
}
