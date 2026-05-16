# Score-independence diagnostics for the monthly BCH application.
#
# The diagnostic follows the logic of BCH: validity concerns group-level
# regression scores, not climate similarity alone. We residualize the outcome
# and treatment on the same fixed effects as the empirical model, form the
# scalar coefficient score, aggregate it by group and month, and summarize
# cross-group correlations.

source("code_and_data/paths.R")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})

data_path <- file.path(data_dir, "weather_conflict_panel.rdata")
assign_path <- file.path(cluster_dir, "cluster_assignments.csv")
out_dir <- tables_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(data_path)
if (!exists("conf")) stop("Expected object `conf` in ", data_path)

dt <- as.data.table(conf)
dt[, `:=`(
  year = as.integer(year),
  month = as.integer(month),
  cell_id = as.integer(cell_id),
  conflict_0_1 = as.numeric(conflict_0_1),
  temp_anom = as.numeric(temp - climatology),
  ADM0_ISO = as.character(ADM0_ISO)
)]
dt[, ym := sprintf("%04d-%02d", year, month)]
setorder(dt, cell_id, ym)

demean_two_way <- function(x, id, time) {
  x - ave(x, id, FUN = mean) - ave(x, time, FUN = mean) + mean(x)
}

dt[, y_tilde := demean_two_way(conflict_0_1, cell_id, ym)]
dt[, x_tilde := demean_two_way(temp_anom, cell_id, ym)]

den <- dt[, sum(x_tilde^2, na.rm = TRUE)]
beta_hat <- dt[, sum(x_tilde * y_tilde, na.rm = TRUE) / den]
dt[, resid := y_tilde - beta_hat * x_tilde]
dt[, score := x_tilde * resid]

score_diagnostics <- function(base, groups, scheme, proposal_id = NA_character_) {
  d <- copy(base)[groups, on = "cell_id"]
  if (anyNA(d$group)) stop("Missing group assignments for ", scheme)

  sizes <- unique(d[, .(cell_id, group)])[, .N, by = group]$N
  g_scores <- d[, .(S = sum(score, na.rm = TRUE)), by = .(group, ym)]
  wide <- dcast(g_scores, ym ~ group, value.var = "S", fill = 0)
  score_mat <- as.matrix(wide[, -"ym"])
  cmat <- suppressWarnings(cor(score_mat))
  off <- cmat[upper.tri(cmat)]
  abs_off <- abs(off)
  energy <- colSums(score_mat^2)

  total_score <- d[, sum(score, na.rm = TRUE), by = group]
  v_cce <- sum(total_score$V1^2) / den^2
  se <- sqrt(v_cce)
  raw_t <- beta_hat / se
  g <- length(unique(d$group))
  bch_t <- raw_t * sqrt((g - 1) / g)
  p_value <- 2 * pt(abs(bch_t), df = g - 1, lower.tail = FALSE)

  tibble(
    scheme = scheme,
    proposal_id = proposal_id,
    G = g,
    cell_size_min = min(sizes),
    cell_size_median = median(sizes),
    cell_size_max = max(sizes),
    mean_abs_score_corr = mean(abs_off, na.rm = TRUE),
    median_abs_score_corr = median(abs_off, na.rm = TRUE),
    max_abs_score_corr = max(abs_off, na.rm = TRUE),
    share_abs_corr_gt_025 = mean(abs_off > 0.25, na.rm = TRUE),
    max_score_energy_share = max(energy / sum(energy)),
    estimate = beta_hat,
    se_cce = se,
    t_bch = bch_t,
    p_bch = p_value
  )
}

assignments <- read_csv(assign_path, show_col_types = FALSE) |>
  transmute(
    cell_id = as.integer(cell_id),
    proposal_id,
    group = paste0("c", cluster)
  )

hiclmr_diag <- bind_rows(lapply(split(assignments, assignments$proposal_id), function(a) {
  k <- unique(sub(".*k([0-9]+).*", "\\1", a$proposal_id))
  score_diagnostics(
    dt,
    a |> select(cell_id, group),
    scheme = paste0("HiClimR G=", as.integer(k)),
    proposal_id = unique(a$proposal_id)
  )
})) |>
  arrange(G)

region_lookup <- tibble::tribble(
  ~ADM0_ISO, ~macro_region,
  "B28", "North Africa", "DZA", "North Africa", "EGY", "North Africa",
  "LBY", "North Africa", "MAR", "North Africa", "MRT", "North Africa",
  "SDZ", "North Africa", "TUN", "North Africa",
  "BEN", "West Africa", "BFA", "West Africa", "CIV", "West Africa",
  "CPV", "West Africa", "GHA", "West Africa", "GIN", "West Africa",
  "GNB", "West Africa", "LBR", "West Africa", "MLI", "West Africa",
  "NER", "West Africa", "NGA", "West Africa", "SEN", "West Africa",
  "SLE", "West Africa", "TGO", "West Africa",
  "CAF", "Central Africa", "CMR", "Central Africa", "TCD", "Central Africa",
  "COG", "Central Africa", "COD", "Central Africa", "GNQ", "Central Africa",
  "GAB", "Central Africa",
  "BDI", "East Africa", "COM", "East Africa", "DJI", "East Africa",
  "ERI", "East Africa", "ETH", "East Africa", "KEN", "East Africa",
  "MDG", "East Africa", "RWA", "East Africa", "SOM", "East Africa",
  "SSD", "East Africa", "TZA", "East Africa", "UGA", "East Africa",
  "AGO", "Southern Africa", "BWA", "Southern Africa", "LSO", "Southern Africa",
  "MOZ", "Southern Africa", "MWI", "Southern Africa", "NAM", "Southern Africa",
  "SWZ", "Southern Africa", "ZAF", "Southern Africa", "ZMB", "Southern Africa",
  "ZWE", "Southern Africa"
)

missing_region <- setdiff(unique(dt$ADM0_ISO), region_lookup$ADM0_ISO)
if (length(missing_region) > 0) {
  stop("Missing macro-region mapping for: ", paste(missing_region, collapse = ", "))
}

macro_groups <- unique(dt[, .(cell_id, ADM0_ISO)]) |>
  left_join(region_lookup, by = "ADM0_ISO") |>
  transmute(cell_id, group = macro_region)

macro_diag <- score_diagnostics(
  dt,
  macro_groups,
  scheme = "Administrative macroregions",
  proposal_id = "AU-style macroregions"
)

country_groups <- unique(dt[, .(cell_id, group = ADM0_ISO)])
country_diag <- score_diagnostics(
  dt,
  country_groups,
  scheme = "Administrative countries",
  proposal_id = "ADM0_ISO"
)

diag <- bind_rows(hiclmr_diag, macro_diag, country_diag) |>
  arrange(case_when(
    scheme == "Administrative macroregions" ~ 998L,
    scheme == "Administrative countries" ~ 999L,
    TRUE ~ G
  ))

write_csv(diag, file.path(out_dir, "score_independence_diagnostics.csv"))

print(diag, n = Inf, width = Inf)
