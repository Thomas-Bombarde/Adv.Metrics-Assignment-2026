# helpers.R
# DGP, grouping, and variance estimators for Bester, Conley & Hansen (2011)
# simulations on a K x K spatial lattice.

# --- DGP: spatial MA(2) field, Chebyshev (max-coordinate) metric ---
# eps_s = sum_{||h|| <= 2} gamma^||h|| * u_{s+h}, with ||h|| = max(|h1|, |h2|).
gen_spatial_field <- function(K, gamma, pad = 2L) {
  M <- K + 2L * pad
  U <- matrix(rnorm(M * M), M, M)
  out <- matrix(0, K, K)
  for (hi in -pad:pad) {
    for (hj in -pad:pad) {
      h_norm <- max(abs(hi), abs(hj))
      if (h_norm > pad) next
      w <- gamma^h_norm
      ri <- (pad + 1L + hi):(pad + K + hi)
      cj <- (pad + 1L + hj):(pad + K + hj)
      out <- out + w * U[ri, cj]
    }
  }
  out
}

# --- Group construction: partition K x K lattice into G contiguous square blocks ---
# Requires sqrt(G) integer and K divisible by sqrt(G).
make_groups <- function(K, G) {
  g_side <- sqrt(G)
  stopifnot(abs(g_side - round(g_side)) < 1e-9)
  g_side <- as.integer(round(g_side))
  stopifnot(K %% g_side == 0L)
  block_size <- K %/% g_side
  i_idx <- rep(seq_len(K), times = K)
  j_idx <- rep(seq_len(K), each  = K)
  bi <- (i_idx - 1L) %/% block_size
  bj <- (j_idx - 1L) %/% block_size
  as.integer(bi * g_side + bj + 1L)
}

# --- Cluster covariance estimator (BCH eq. on p. 140) ---
# hat V = (1/N) sum_g (X_g' e_g)(X_g' e_g)'
compute_cce <- function(X, resid, groups) {
  k <- ncol(X)
  G <- max(groups)
  V <- matrix(0, k, k)
  for (g in seq_len(G)) {
    idx <- which(groups == g)
    Xe <- crossprod(X[idx, , drop = FALSE], resid[idx])
    V <- V + tcrossprod(Xe)
  }
  V / length(resid)
}

# --- HC0 (heteroskedasticity-robust, ignores dependence) ---
compute_hc0 <- function(X, resid) {
  N <- length(resid)
  crossprod(X * as.vector(resid)) / N
}

# --- 2D Bartlett HAC (product kernel on Chebyshev metric, cutoff H) ---
# Loops over lag offsets, so cost is O(H^2 * N).
compute_hac_bartlett <- function(X, resid, K, H) {
  N <- K^2
  k <- ncol(X)
  V <- matrix(0, k, k)
  rows <- rep(seq_len(K), times = K)
  cols <- rep(seq_len(K), each  = K)
  for (dx in -H:H) {
    for (dy in -H:H) {
      w <- (1 - abs(dx) / H) * (1 - abs(dy) / H)
      if (w <= 0) next
      r_j <- rows + dy
      c_j <- cols + dx
      valid <- r_j >= 1L & r_j <= K & c_j >= 1L & c_j <= K
      if (!any(valid)) next
      i_idx <- which(valid)
      j_idx <- (c_j[valid] - 1L) * K + r_j[valid]
      Xe_i <- X[i_idx, , drop = FALSE] * resid[i_idx]
      Xe_j <- X[j_idx, , drop = FALSE] * resid[j_idx]
      V <- V + w * crossprod(Xe_i, Xe_j)
    }
  }
  V / N
}

# --- One Monte Carlo replication ---
# Returns t-statistics for H0: beta_slope = 1 under various variance estimators,
# and the slope-row of the CCE for tracking the V_hat distribution.
run_one_replication <- function(K, gamma, G_grid, H_hac) {
  x_mat <- gen_spatial_field(K, gamma)
  e_mat <- gen_spatial_field(K, gamma)
  N <- K^2
  X <- cbind(1, as.vector(x_mat))
  y <- as.vector(x_mat) + as.vector(e_mat)  # beta_0 = 0, beta_1 = 1

  XtX <- crossprod(X)
  XtX_inv <- solve(XtX)
  beta_hat <- as.vector(XtX_inv %*% crossprod(X, y))
  resid <- as.vector(y - X %*% beta_hat)
  Q_hat <- XtX / N
  Q_inv <- solve(Q_hat)

  slope_se <- function(V_hat) sqrt((Q_inv %*% V_hat %*% Q_inv)[2, 2] / N)
  slope_t  <- function(V_hat) (beta_hat[2] - 1) / slope_se(V_hat)

  # HC0
  V_hc0 <- compute_hc0(X, resid)
  t_hc0 <- slope_t(V_hc0)

  # HAC (Bartlett)
  V_hac <- compute_hac_bartlett(X, resid, K, H_hac)
  t_hac <- slope_t(V_hac)
  V_hac_asymp <- (Q_inv %*% V_hac %*% Q_inv)[2, 2]  # for tracking concentration of HAC

  # CCE for each G in grid
  cce_rows <- lapply(G_grid, function(G) {
    groups <- make_groups(K, G)
    V_cce <- compute_cce(X, resid, groups)
    t_raw <- slope_t(V_cce)
    V_asymp <- (Q_inv %*% V_cce %*% Q_inv)[2, 2]  # asymptotic-variance scale: N * Var(beta_hat)
    data.frame(
      G = G,
      V_slope_asymp = V_asymp,        # converges in distribution to a non-degenerate random variable
      V_slope = V_asymp / N,          # Var(beta_hat): goes to zero like 1/N
      t_cce = t_raw,
      t_cce_rescaled = t_raw * sqrt((G - 1) / G)
    )
  })
  cce_df <- do.call(rbind, cce_rows)

  list(
    beta_slope = beta_hat[2],
    t_hc0 = t_hc0,
    t_hac = t_hac,
    V_hac_asymp = V_hac_asymp,
    cce = cce_df
  )
}
