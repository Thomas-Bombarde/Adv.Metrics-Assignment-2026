# propose_region_clusterings.R
# Build candidate HiClimR regionalizations for the weather-conflict panel.
#
# Usage from the report project root:
#   Rscript code_and_data/propose_region_clusterings.R

source("code_and_data/paths.R")

required_packages <- c(
  "HiClimR", "data.table", "dplyr", "ggplot2", "ggrepel",
  "RColorBrewer", "patchwork", "readr", "rnaturalearth", "sf"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them first, e.g. install.packages(..., lib = '.r-lib')."
  )
}

suppressPackageStartupMessages({
  library(HiClimR)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(readr)
  library(rnaturalearth)
  library(sf)
})

set.seed(20260514)

out_dir <- cluster_dir
map_dir <- cluster_map_dir
dir.create(map_dir, showWarnings = FALSE, recursive = TRUE)

load(file.path(data_dir, "weather_conflict_panel.rdata"))
if (!exists("conf")) stop("Expected object `conf` in weather_conflict_panel.rdata.")

dt <- as.data.table(conf)
dt[, time_id := sprintf("%04d-%02d", year, month)]
dt[, temp_anom := temp - climatology]

coords <- unique(dt[, .(cell_id, longitude, latitude, ADM0_ISO)])
setorder(coords, cell_id)

wide_temp <- dcast(dt, cell_id ~ time_id, value.var = "temp_anom")
setorder(wide_temp, cell_id)
stopifnot(identical(wide_temp$cell_id, coords$cell_id))

x_temp <- as.matrix(wide_temp[, -"cell_id"])
storage.mode(x_temp) <- "double"
colnames(x_temp) <- seq_len(ncol(x_temp))

cell_summary <- dt[, .(
  mean_temp = mean(temp),
  mean_temp_anom = mean(temp_anom),
  sd_temp_anom = sd(temp_anom),
  conflict_month_share = mean(conflict_0_1),
  total_conflicts = sum(conflicts),
  conflict_months = sum(conflict_0_1)
), by = .(cell_id)]

proposal_specs <- tibble::tribble(
  ~proposal_id,              ~k,  ~method,     ~hybrid, ~kH, ~contigConst, ~nPC, ~description,
  "temp_ward_k04",            4L, "ward",      FALSE,   NA,       0,       NA,  "Broad Ward clusters on detrended, standardized monthly temperature anomalies.",
  "temp_regional_k08",        8L, "regional",  FALSE,   NA,       0,       NA,  "Regional-linkage clusters emphasizing separation of regional mean time series.",
  "temp_ward_contig_k12",    12L, "ward",      FALSE,   NA,       0.35,    NA,  "Ward clusters with a moderate geographic contiguity weight.",
  "temp_ward_pca_k16",       16L, "ward",      FALSE,   NA,       0,       12L, "Ward clusters after PCA filtering of the anomaly matrix.",
  "temp_hybrid_k20",         20L, "ward",      TRUE,    60L,      0,       NA,  "Hybrid Ward-regional clusters: Ward tree with regional reconstruction above 60 clusters."
)

run_hiclimr <- function(spec) {
  n_pc <- if (is.na(spec$nPC)) NULL else as.integer(spec$nPC)
  k_h <- if (is.na(spec$kH)) NULL else as.integer(spec$kH)

  HiClimR(
    x = x_temp,
    lon = coords$longitude,
    lat = coords$latitude,
    lonStep = 1,
    latStep = 1,
    geogMask = FALSE,
    contigConst = spec$contigConst,
    varThresh = 0,
    detrend = TRUE,
    standardize = TRUE,
    nPC = n_pc,
    method = spec$method,
    hybrid = spec$hybrid,
    kH = k_h,
    nSplit = 1,
    upperTri = TRUE,
    verbose = FALSE,
    validClimR = TRUE,
    rawStats = TRUE,
    k = spec$k,
    minSize = 1,
    alpha = 0.05,
    plot = FALSE,
    dendrogram = FALSE,
    labels = FALSE
  )
}

palette_for_k <- function(k) {
  if (k <= 12) {
    RColorBrewer::brewer.pal(max(3, k), "Set3")[seq_len(k)]
  } else {
    grDevices::hcl.colors(k, palette = "Spectral", rev = TRUE)
  }
}

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
bbox_pad <- 0.5
xlim <- range(coords$longitude) + c(-bbox_pad, bbox_pad)
ylim <- range(coords$latitude) + c(-bbox_pad, bbox_pad)

make_cluster_map <- function(assignments, spec) {
  n_clusters <- as.integer(spec$k)
  plot_df <- assignments %>%
    mutate(cluster = factor(cluster, levels = seq_len(n_clusters)))

  labels_df <- plot_df %>%
    group_by(cluster) %>%
    summarise(
      longitude = weighted.mean(longitude, sd_temp_anom + 1e-6),
      latitude = weighted.mean(latitude, sd_temp_anom + 1e-6),
      .groups = "drop"
    )

  ggplot() +
    geom_sf(data = world, fill = NA, color = "grey40", linewidth = 0.25) +
    geom_tile(
      data = plot_df,
      aes(x = longitude, y = latitude, fill = cluster),
      width = 1,
      height = 1,
      alpha = 0.92
    ) +
    ggrepel::geom_text_repel(
      data = labels_df,
      aes(x = longitude, y = latitude, label = cluster),
      seed = 20260514,
      min.segment.length = 0.35,
      max.overlaps = Inf,
      size = 3.2,
      fontface = "bold",
      color = "grey10",
      box.padding = 0.6,
      point.padding = 0.4,
      segment.color = "grey35",
      segment.size = 0.25
    ) +
    scale_fill_manual(values = palette_for_k(n_clusters), drop = FALSE) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(x = NULL, y = NULL, fill = "Cluster") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.position = "right",
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, margin = margin(b = 6))
    ) +
    ggtitle(
      sprintf("%s: %d clusters", spec$proposal_id, n_clusters),
      subtitle = spec$description
    )
}

objects <- list()
all_assignments <- list()
all_summaries <- list()
all_plots <- list()

for (i in seq_len(nrow(proposal_specs))) {
  spec <- as.list(proposal_specs[i, ])
  message("Running ", spec$proposal_id, " ...")
  t0 <- proc.time()
  fit <- run_hiclimr(spec)
  elapsed <- unname((proc.time() - t0)[["elapsed"]])
  objects[[spec$proposal_id]] <- fit

  assignments <- coords %>%
    mutate(cluster = as.integer(fit$region)) %>%
    left_join(as_tibble(cell_summary), by = "cell_id") %>%
    mutate(
      proposal_id = spec$proposal_id,
      k = spec$k,
      method = spec$method,
      hybrid = spec$hybrid,
      kH = ifelse(is.na(spec$kH), NA_integer_, spec$kH),
      contigConst = spec$contigConst,
      nPC = ifelse(is.na(spec$nPC), NA_integer_, spec$nPC)
    )
  all_assignments[[spec$proposal_id]] <- assignments

  summary_df <- assignments %>%
    group_by(proposal_id, cluster) %>%
    summarise(
      cells = n(),
      countries = paste(sort(unique(ADM0_ISO)), collapse = ","),
      lon_min = min(longitude),
      lon_max = max(longitude),
      lat_min = min(latitude),
      lat_max = max(latitude),
      mean_temp = mean(mean_temp),
      sd_temp_anom = mean(sd_temp_anom),
      conflict_month_share = mean(conflict_month_share),
      total_conflicts = sum(total_conflicts),
      .groups = "drop"
    ) %>%
    mutate(
      elapsed_seconds = elapsed,
      interCor_mean = fit$statSum["Mean", "interCor"],
      intraCor_mean = fit$statSum["Mean", "intraCor"],
      diffCor_mean = fit$statSum["Mean", "diffCor"]
    )
  all_summaries[[spec$proposal_id]] <- summary_df

  p <- make_cluster_map(assignments, spec)
  all_plots[[spec$proposal_id]] <- p
  ggsave(
    filename = file.path(map_dir, paste0(spec$proposal_id, ".png")),
    plot = p,
    width = 8.5,
    height = 7,
    dpi = 180
  )
}

assignments_out <- bind_rows(all_assignments)
summary_out <- bind_rows(all_summaries)

readr::write_csv(assignments_out, file.path(out_dir, "cluster_assignments.csv"))
readr::write_csv(summary_out, file.path(out_dir, "cluster_summaries.csv"))
readr::write_csv(proposal_specs, file.path(out_dir, "proposal_specs.csv"))
saveRDS(objects, file.path(out_dir, "hiclimr_objects.rds"))

overview <- wrap_plots(all_plots, ncol = 2) +
  plot_annotation(title = "Candidate HiClimR Temperature-Anomaly Regionalizations")
ggsave(
  filename = file.path(map_dir, "all_clusterings_overview.png"),
  plot = overview,
  width = 15,
  height = 18,
  dpi = 160
)

summary_lines <- c(
  "# Candidate HiClimR Region Clusterings",
  "",
  "Input: `code_and_data/analysis_data/weather_conflict_panel.rdata`, object `conf`.",
  "",
  "Clustering variable: monthly temperature anomaly, `temp - climatology`, arranged as cells by months.",
  "HiClimR preprocessing: linear detrending and row standardization for all proposals.",
  "Conflict variables are not used in the correlation distance because many cells have no conflict variation; they are reported as diagnostics in `cluster_summaries.csv`.",
  "",
  "## Proposals",
  ""
)

for (i in seq_len(nrow(proposal_specs))) {
  spec <- proposal_specs[i, ]
  ss <- summary_out %>% filter(proposal_id == spec$proposal_id)
  summary_lines <- c(
    summary_lines,
    sprintf(
      "- `%s`: %d groups, method `%s`, hybrid `%s`, contiguity `%s`, nPC `%s`. %s Mean intra-cluster correlation: %.3f; mean inter-cluster correlation: %.3f.",
      spec$proposal_id,
      spec$k,
      spec$method,
      spec$hybrid,
      spec$contigConst,
      ifelse(is.na(spec$nPC), "none", as.character(spec$nPC)),
      spec$description,
      unique(ss$intraCor_mean),
      unique(ss$interCor_mean)
    )
  )
}

summary_lines <- c(
  summary_lines,
  "",
  "## Outputs",
  "",
  "- `code_and_data/output/clusterings/proposal_specs.csv`: technique and parameter grid.",
  "- `code_and_data/output/clusterings/cluster_assignments.csv`: one row per grid cell and proposal.",
  "- `code_and_data/output/clusterings/cluster_summaries.csv`: cluster-level geography, temperature, validation, and conflict diagnostics.",
  "- `code_and_data/output/clusterings/hiclimr_objects.rds`: fitted HiClimR objects for reuse.",
  "- `code_and_data/output/clusterings/maps/*.png`: one map per proposal plus `all_clusterings_overview.png`."
)

writeLines(summary_lines, file.path(out_dir, "README.md"))

print(proposal_specs)
cat("\nWrote outputs to:\n", normalizePath(out_dir), "\n", sep = "")
