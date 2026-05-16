packages <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "ggrepel",
  "HiClimR",
  "knitr",
  "patchwork",
  "readr",
  "rmarkdown",
  "rnaturalearth",
  "RColorBrewer",
  "scales",
  "sf",
  "tibble",
  "tidyr"
)

project_root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(project_root, "report.Rproj"))) {
  stop(
    "Run this from the report project root, the folder containing report.Rproj.",
    call. = FALSE
  )
}

local_lib <- file.path(project_root, "code_and_data", ".r-lib")
dir.create(local_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(normalizePath(local_lib), .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"))

missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) == 0) {
  message("All required packages are already installed.")
} else {
  install.packages(missing, lib = local_lib)
}
