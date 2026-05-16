# Simulate the fixed-G trade-off for the same group counts used later
# in the HiClimR empirical application.
#
# Usage from the report project root:
#   Rscript code_and_data/simulate_group_tradeoff.R [n_sim]

source("code_and_data/paths.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

source(file.path(code_dir, "helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
n_sim <- if (length(args) >= 1) as.integer(args[1]) else 400L

set.seed(20260516)

K <- 60L
gamma <- 0.6
N <- K^2

group_specs <- tibble::tribble(
  ~G, ~row_blocks, ~col_blocks,
   4L,          2L,          2L,
   8L,          2L,          4L,
  12L,          3L,          4L,
  16L,          4L,          4L,
  20L,          4L,          5L
)

make_groups_rect <- function(K, row_blocks, col_blocks) {
  stopifnot(K %% row_blocks == 0L, K %% col_blocks == 0L)
  block_h <- K %/% row_blocks
  block_w <- K %/% col_blocks
  rows <- rep(seq_len(K), times = K)
  cols <- rep(seq_len(K), each = K)
  br <- (rows - 1L) %/% block_h
  bc <- (cols - 1L) %/% block_w
  as.integer(br * col_blocks + bc + 1L)
}

groups_by_g <- setNames(
  lapply(seq_len(nrow(group_specs)), function(i) {
    make_groups_rect(K, group_specs$row_blocks[i], group_specs$col_blocks[i])
  }),
  group_specs$G
)

one_rep <- function(rep_id) {
  x_mat <- gen_spatial_field(K, gamma)
  e_mat <- gen_spatial_field(K, gamma)
  X <- cbind(1, as.vector(x_mat))
  y <- as.vector(x_mat) + as.vector(e_mat)

  XtX <- crossprod(X)
  beta_hat <- as.vector(solve(XtX, crossprod(X, y)))
  resid <- as.vector(y - X %*% beta_hat)
  Q_inv <- solve(XtX / N)

  bind_rows(lapply(seq_len(nrow(group_specs)), function(i) {
    G <- group_specs$G[i]
    V_cce <- compute_cce(X, resid, groups_by_g[[as.character(G)]])
    V_asymp <- (Q_inv %*% V_cce %*% Q_inv)[2, 2]
    se <- sqrt(V_asymp / N)
    t_raw <- (beta_hat[2] - 1) / se
    tibble::tibble(
      rep_id = rep_id,
      G = G,
      group_size = N / G,
      V_slope_asymp = V_asymp,
      t_cce = t_raw,
      t_bch = t_raw * sqrt((G - 1) / G),
      reject_normal = abs(t_raw) > qnorm(0.975),
      reject_bch = abs(t_raw * sqrt((G - 1) / G)) > qt(0.975, df = G - 1)
    )
  }))
}

cat(sprintf("Running %d trade-off replications on K=%d, gamma=%.1f\n",
            n_sim, K, gamma))
sim <- bind_rows(lapply(seq_len(n_sim), one_rep))

out_dir <- figures_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(sim, file.path(out_dir, "figD_group_tradeoff_sim.csv"))

variance_summary <- sim |>
  group_by(G, group_size) |>
  summarise(
    q25 = quantile(V_slope_asymp, 0.25),
    median = median(V_slope_asymp),
    q75 = quantile(V_slope_asymp, 0.75),
    .groups = "drop"
  )

rejection_summary <- sim |>
  group_by(G) |>
  summarise(
    `CCE + N(0,1)` = mean(reject_normal),
    `CCE + BCH t_{G-1}` = mean(reject_bch),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(-G, names_to = "method", values_to = "rate")

p_var <- ggplot(variance_summary, aes(x = G, y = median)) +
  geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 0.8,
                 colour = "#6b7280") +
  geom_point(size = 2.4, colour = "#1b6ca8") +
  scale_x_continuous(breaks = group_specs$G) +
  labs(
    x = "Number of groups, G",
    y = expression("Asymptotic slope variance scale"),
    caption = "Points are medians; bars are interquartile ranges."
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

p_rej <- ggplot(rejection_summary,
                aes(x = G, y = rate, colour = method, group = method)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey35") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_x_continuous(breaks = group_specs$G) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Number of groups, G",
    y = "Null rejection rate",
    colour = NULL
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

fig <- p_var + p_rej + plot_layout(widths = c(1, 1))

ggsave(file.path(out_dir, "figD_group_tradeoff.png"), fig,
       width = 9.2, height = 4.1, dpi = 180)

cat("Saved code_and_data/output/figures/figD_group_tradeoff.png\n")
