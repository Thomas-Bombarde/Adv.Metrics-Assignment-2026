# run_sim.R
# Monte Carlo driver. Sweeps gamma (spatial MA decay), K (lattice size),
# and G (number of CCE groups), and saves per-replication t-statistics
# and slope-variance estimates for downstream plotting.
#
# Usage from the report project root:
#   Rscript code_and_data/run_sim.R [n_sim]

source("code_and_data/paths.R")

suppressPackageStartupMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
n_sim <- if (length(args) >= 1) as.integer(args[1]) else 500L

source(file.path(code_dir, "helpers.R"))

set.seed(20260512)

# --- Design ---
gamma_grid <- c(0, 0.3, 0.6)
K_grid     <- c(12, 24, 36, 48)       # all divisible by sqrt(G) for G in {4, 9, 16}
G_grid     <- c(4, 9, 16)             # values of G for CCE
H_hac      <- 5L                      # HAC kernel cutoff (matches MA range + buffer)

# Sanity: each (K, G) must be compatible
compat <- function(K, G) {
  g_side <- sqrt(G)
  abs(g_side - round(g_side)) < 1e-9 && K %% round(g_side) == 0
}
for (K in K_grid) for (G in G_grid) stopifnot(compat(K, G))

# --- Run MC ---
cat(sprintf("Running %d reps per (gamma, K) cell; total cells = %d\n",
            n_sim, length(gamma_grid) * length(K_grid)))

records <- list()
i_row <- 1L
t0 <- Sys.time()
for (gamma in gamma_grid) {
  for (K in K_grid) {
    cat(sprintf("  gamma=%.2f, K=%d (L_N for G=4: %d) ...\n",
                gamma, K, (K / 2)^2))
    for (rep_id in seq_len(n_sim)) {
      r <- run_one_replication(K, gamma, G_grid, H_hac)
      base <- data.frame(
        rep_id = rep_id, gamma = gamma, K = K, N = K^2,
        beta_slope = r$beta_slope, t_hc0 = r$t_hc0,
        t_hac = r$t_hac, V_hac_asymp = r$V_hac_asymp
      )
      records[[i_row]] <- cbind(base, r$cce)
      i_row <- i_row + 1L
    }
  }
}
cat(sprintf("Total elapsed: %.1f s\n", as.numeric(Sys.time() - t0, units = "secs")))

df <- bind_rows(records)
saveRDS(df, file.path(results_dir, "sim_results.rds"))
cat(sprintf("Saved %d rows to results/sim_results.rds\n", nrow(df)))
