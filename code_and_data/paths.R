# Shared path convention for the reproducibility scripts.
#
# Open `report.Rproj` first, or otherwise set the working directory to the
# report project root before running these scripts.

project_root <- normalizePath(getwd(), mustWork = TRUE)

if (!file.exists(file.path(project_root, "report.Rproj"))) {
  stop(
    "These scripts expect the working directory to be the report project root ",
    "(the folder containing report.Rproj). Current working directory: ",
    project_root,
    call. = FALSE
  )
}

code_dir <- file.path(project_root, "code_and_data")
data_dir <- file.path(code_dir, "analysis_data")
results_dir <- file.path(code_dir, "results")
output_dir <- file.path(code_dir, "output")
figures_dir <- file.path(output_dir, "figures")
cluster_dir <- file.path(output_dir, "clusterings")
cluster_map_dir <- file.path(cluster_dir, "maps")
tables_dir <- file.path(output_dir, "tables")

local_lib <- file.path(code_dir, ".r-lib")
if (dir.exists(local_lib)) {
  .libPaths(c(normalizePath(local_lib), .libPaths()))
}

