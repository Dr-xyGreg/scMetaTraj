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
  edges <- cbind(
    from = rep(seq_len(n), each = k),
    to   = as.vector(nn_idx)
  )
  weights <- D[edges]
  
  g <- igraph::graph_from_edgelist(edges, directed = TRUE)
  igraph::E(g)$weight <- weights
  g <- igraph::as_undirected(
    g, mode = "collapse",
    edge_attr_comb = list(weight = "min")
  )
  
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
