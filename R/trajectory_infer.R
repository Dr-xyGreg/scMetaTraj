#' Infer metabolic pseudotime from a k-nearest-neighbor graph
#'
#' @description
#' \code{scMetaTraj_infer()} builds a weighted k-nearest-neighbor graph in
#' metabolic PCA space and computes graph distances from a selected root cell.
#' The rescaled distances define metabolic pseudotime (mPT).
#'
#' @param embedding Numeric matrix of cells x PCs, usually returned by
#'   \code{\link{scMetaTraj_embed}(method = "PCA")}.
#' @param k Integer. Number of nearest neighbors used to define the graph.
#' @param root_mode Character. Strategy used to choose the root cell. One of
#'   \code{"pc1_min"}, \code{"pc1_max"}, \code{"axis_min"},
#'   \code{"axis_max"}, or \code{"manual"}.
#' @param axis_score Optional numeric vector used when \code{root_mode} is
#'   \code{"axis_min"} or \code{"axis_max"}.
#' @param root_cell Optional character scalar giving the row name of the root
#'   cell when \code{root_mode = "manual"}.
#' @param scale Logical. Whether to rescale graph distances to the interval
#'   \code{[0, 1]}.
#'
#' @return A named list with elements:
#' \itemize{
#'   \item \code{mPT}: numeric vector of metabolic pseudotime values.
#'   \item \code{root}: selected root cell name.
#'   \item \code{dist}: raw graph distances from the root cell.
#' }
#'
#' @examples
#' set.seed(123)
#' embedding <- matrix(rnorm(120 * 5), nrow = 120, ncol = 5)
#' rownames(embedding) <- paste0("Cell", seq_len(nrow(embedding)))
#' mpt <- scMetaTraj_infer(embedding, k = 15, root_mode = "pc1_min")
#' head(mpt$mPT)
#'
#' @export
scMetaTraj_infer <- function(
    embedding,
    k = 20,
    root_mode = c("pc1_min", "pc1_max", "axis_min", "axis_max", "manual"),
    axis_score = NULL,
    root_cell = NULL,
    scale = TRUE
) {
  root_mode <- match.arg(root_mode)
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required.")
  }
  
  emb <- as.matrix(embedding)
  if (is.null(rownames(emb))) stop("embedding must have rownames (cell IDs).")
  
  n <- nrow(emb)
  k <- min(k, n - 1)
  
  ## ---- FIX HERE ----
  D <- as.matrix(stats::dist(emb))
  
  nn_idx <- base::apply(D, 1, function(x) order(x)[2:(k + 1)])
  edges <- data.frame(
    from = rep(seq_len(n), each = k),
    to = as.vector(nn_idx)
  )
  edges$weight <- D[cbind(edges$from, edges$to)]
  edges$from_u <- pmin(edges$from, edges$to)
  edges$to_u <- pmax(edges$from, edges$to)
  edges <- stats::aggregate(
    weight ~ from_u + to_u,
    data = edges,
    FUN = min
  )
  
  g <- igraph::graph_from_data_frame(
    edges[, c("from_u", "to_u", "weight")],
    directed = FALSE,
    vertices = data.frame(name = seq_len(n))
  )
  igraph::E(g)$weight <- edges$weight
  
  root_idx <- switch(
    root_mode,
    pc1_min  = which.min(emb[, 1]),
    pc1_max  = which.max(emb[, 1]),
    axis_min = which.min(axis_score),
    axis_max = which.max(axis_score),
    manual   = which(rownames(emb) == root_cell)
  )
  
  root_name <- rownames(emb)[root_idx]
  
  d <- igraph::distances(g, v = root_idx, weights = igraph::E(g)$weight)
  d <- as.numeric(d)
  d[is.infinite(d)] <- NA
  
  mPT <- d
  if (scale) {
    rng <- range(mPT, na.rm = TRUE)
    if (diff(rng) > 0) {
      mPT <- (mPT - rng[1]) / diff(rng)
    }
  }
  
  names(mPT) <- rownames(emb)
  names(d)   <- rownames(emb)
  
  list(mPT = mPT, root = root_name, dist = d)
}
