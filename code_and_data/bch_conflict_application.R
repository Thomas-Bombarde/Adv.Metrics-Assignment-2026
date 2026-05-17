# bch_conflict_application.R
# Monthly empirical BCH application for the Africa weather-conflict panel.
#
# This script reproduces the report's main non-descriptive results table:
#
#   conflict_0_1_it = beta * temp_anom_it + cell FE + year-month FE + e_it
#
# The estimate is obtained with fixest::feols and the BCH scalar correction is
# applied through fixest::summary() with the small-sample settings reported in
# the essay. The raw CCE standard error is recovered from the corrected one.

source("code_and_data/paths.R")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(fixest)
  library(readr)
})

data_path <- file.path(data_dir, "weather_conflict_panel.rdata")
cluster_path <- file.path(cluster_dir, "cluster_assignments.csv")
out_cluster_dir <- cluster_dir
out_table_dir <- tables_dir
out_figure_dir <- figures_dir
dir.create(out_cluster_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)

load(data_path)
if (!exists("conf")) stop("Expected object `conf` in ", data_path)

dt <- as.data.table(conf)
dt[, `:=`(
  year = as.integer(year),
  month = as.integer(month),
  cell_id = as.integer(cell_id),
  conflict_0_1 = as.numeric(conflict_0_1),
  temp_anom = as.numeric(temp - climatology)
)]
dt[, ym := sprintf("%04d-%02d", year, month)]
setorder(dt, cell_id, ym)

n <- nrow(dt)

assignments <- read_csv(cluster_path, show_col_types = FALSE) |>
  transmute(
    cell_id = as.integer(cell_id),
    proposal_id,
    k = as.integer(k),
    method,
    hybrid,
    contigConst,
    nPC,
    cluster = as.integer(cluster)
)

one_group_result <- function(a) {
  d <- dt[a, on = "cell_id"]
  if (anyNA(d$cluster)) stop("Missing cluster assignment for ", unique(a$proposal_id))

  d[, bch_region := factor(cluster)]
  est <- feols(conflict_0_1 ~ temp_anom | cell_id + ym, data = d)
  est_sum <- summary(
    est,
    vcov = ~ bch_region,
    ssc = ssc(
      K.adj = FALSE,
      K.fixef = "nonnested",
      G.adj = TRUE,
      t.df = "min"
    )
  )

  g <- nlevels(d$bch_region)
  beta_hat <- coef(est)[["temp_anom"]]
  r2 <- summary(est)$sq.cor
  bch_se <- se(est_sum)[["temp_anom"]]
  raw_se <- bch_se / sqrt(g / (g - 1))
  raw_var <- raw_se^2
  t_value <- tstat(est_sum)[["temp_anom"]]
  p_value <- pvalue(est_sum)[["temp_anom"]]
  crit <- qt(0.975, df = g - 1)

  tibble(
    proposal_id = unique(a$proposal_id),
    estimate = beta_hat,
    std_error = bch_se,
    t_value = t_value,
    p_value = p_value,
    ci_low = beta_hat - crit * bch_se,
    ci_high = beta_hat + crit * bch_se,
    raw_cce_std_error = raw_se,
    raw_cce_variance = raw_var,
    asymptotic_variance = n * raw_var,
    nobs = n,
    r2 = r2,
    k = unique(a$k),
    method = unique(a$method),
    hybrid = unique(a$hybrid),
    contigConst = unique(a$contigConst),
    nPC = unique(a$nPC)
  )
}

results <- bind_rows(
  lapply(split(assignments, 
               assignments$proposal_id), 
         one_group_result)) |>
  arrange(k)

write_csv(results, file.path(out_cluster_dir, "fixest_bch_results.csv"))
write_csv(results, file.path(out_table_dir, "monthly_bch_results.csv"))
write_csv(results, file.path(out_table_dir, "bch_application_variance_results.csv"))

compact_table <- results |>
  transmute(
    Method = paste0("BCH CCE: ", proposal_id),
    G = k,
    Estimate = sprintf("%.5f", estimate),
    `Raw CCE variance` = sprintf("%.6g", raw_cce_variance),
    `Asymptotic variance` = sprintf("%.4f", asymptotic_variance),
    SE = sprintf("%.5f", std_error),
    Statistic = sprintf("%.2f", t_value),
    Reference = paste0("t_", k - 1),
    `p-value` = sprintf("%.3f", p_value)
  )
write_csv(compact_table, file.path(out_table_dir, "bch_application_variance_table.csv"))
saveRDS(list(results = results, compact_table = compact_table), file.path(out_table_dir, "bch_application_variance_results.rds"))

variance_plot_data <- results |>
  mutate(k = as.integer(k))

p_var <- ggplot(variance_plot_data, aes(x = k, y = asymptotic_variance)) +
  geom_line(colour = "#1b6ca8", linewidth = 0.8) +
  geom_point(colour = "#1b6ca8", size = 2.3) +
  scale_x_continuous(breaks = variance_plot_data$k) +
  labs(
    x = "Number of BCH spatial groups (G)",
    y = expression("Asymptotic variance scale: " * N * widehat(Var)(widehat(beta)))
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(out_figure_dir, "bch_application_variance_by_g.png"),
  p_var,
  width = 7.2,
  height = 4.2,
  dpi = 300
)

print(results |>
  transmute(
    proposal_id,
    G = k,
    estimate,
    std_error,
    t_value,
    p_value,
    ci_low,
    ci_high,
    nobs
  ))

cat("Wrote monthly BCH application outputs to:\n")
cat(" - ", file.path(out_cluster_dir, "fixest_bch_results.csv"), "\n", sep = "")
cat(" - ", file.path(out_table_dir, "monthly_bch_results.csv"), "\n", sep = "")
cat(" - ", file.path(out_table_dir, "bch_application_variance_results.csv"), "\n", sep = "")
cat(" - ", file.path(out_figure_dir, "bch_application_variance_by_g.png"), "\n", sep = "")
