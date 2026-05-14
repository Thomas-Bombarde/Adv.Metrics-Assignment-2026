# =======================
# Debut
# =======================
library(tidyverse)
library(data.table)
library(fixest)
library(modelsummary)
library(sf)

# Keep tidyverse verbs stable even if another package masks them in the session
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
summarize <- dplyr::summarize
transmute <- dplyr::transmute
left_join <- dplyr::left_join
arrange <- dplyr::arrange
group_by <- dplyr::group_by
slice_head <- dplyr::slice_head
bind_rows <- dplyr::bind_rows
first <- dplyr::first
desc <- dplyr::desc

# define directories
project_dir <- getwd()
output_dir <- file.path(project_dir, "output")
data_dir <- file.path(project_dir, "analysis_data")
raw_data_dir <- file.path(project_dir, "raw_data")

table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_minimal(base_size = 11))
setFixest_notes(FALSE)

# LaTeX setup
latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([&_#$%])", "\\\\\\1", x)
  x
}

format_value <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "",
    formatC(x, digits = digits, format = "f", big.mark = ",")
  )
}

write_booktabs <- function(df, path, caption, label, digits = 3) {
  out <- df
  numeric_cols <- vapply(out, is.numeric, logical(1))
  out[numeric_cols] <- lapply(out[numeric_cols], format_value, digits = digits)
  out[] <- lapply(out, latex_escape)

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{", paste(c("l", rep("r", ncol(out) - 1)), collapse = ""), "}"),
    "\\toprule",
    paste(names(out), collapse = " & "),
    "\\\\",
    "\\midrule",
    apply(out, 1, function(row) paste0(paste(row, collapse = " & "), "\\\\")),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
  writeLines(lines, path, useBytes = TRUE)
}

# =======================
# 0. Data loading and preparation
# =======================

# 1. Load conflict panel (conf): cell-month panel with conflict and temperature data
load_conflict_panel <- function() {
  rdata_path <- file.path(data_dir, "weather_conflict_panel.rdata")
  csv_path <- file.path(data_dir, "weather_conflict_panel.csv")

  if (file.exists(rdata_path)) {
    loaded_name <- load(rdata_path)
    get(loaded_name[1])
  } else if (file.exists(csv_path)) {
    read_csv(csv_path, show_col_types = FALSE)
  } else {
    stop("Cannot find weather_conflict_panel.rdata or weather_conflict_panel.csv in analysis_data/.")
  }
}

conf <- load_conflict_panel() %>%
  mutate(
    year = as.integer(year),
    month = as.integer(month),
    cell_id = as.integer(cell_id),
    conflict_0_1 = as.integer(conflict_0_1),
    conflicts = replace_na(as.numeric(conflicts), 0),
    temp = as.numeric(temp),
    climatology = as.numeric(climatology),
    temp_anomaly = as.numeric(temp - climatology, na.rm = TRUE)
  )

# 2. Aggregate by year (conf_year): annual grid-cell level
conf_year <- conf %>%
  group_by(cell_id, year) %>%
  summarise(
    temp = mean(temp, na.rm = TRUE),
    temp_anomaly = mean(temp_anomaly, na.rm = TRUE),
    climatology = mean(climatology, na.rm = TRUE),
    conflict = as.integer(any(conflict_0_1 == 1, na.rm = TRUE)),
    conflict_events = sum(conflicts, na.rm = TRUE),
    conflict_months = sum(conflict_0_1, na.rm = TRUE),
    iso3c = first(ADM0_ISO),
    latitude = first(latitude),
    longitude = first(longitude),
    .groups = "drop"
  ) %>%
  arrange(cell_id, year)

# set up lags for temperature and conflict
setDT(conf_year)
conf_year[, temp_lag := shift(temp, 1), by = cell_id]
conf_year[, conflict_lag := shift(conflict, 1), by = cell_id]

# 3. Add country-level covariates from WDI (conf_year_cov): GDP per capita, poverty, region, income group
wdi_path <- file.path(raw_data_dir, "WDI_1_11_2024.csv")
if (file.exists(wdi_path)) {
  wdi <- read_csv(wdi_path, show_col_types = FALSE) %>%
    transmute(
      iso3c,
      year = as.integer(year),
      gdpcap = suppressWarnings(as.numeric(`NY.GDP.PCAP.KD`)),
      poverty = suppressWarnings(as.numeric(`SI.POV.DDAY`)),
      region,
      income
    ) %>%
    filter(!is.na(iso3c), iso3c != "", !is.na(year))

  conf_year_cov <- conf_year %>%
    left_join(wdi, by = c("iso3c", "year")) %>%
    mutate(
      loggdpcap = log(gdpcap),
      loggdpcapdm = loggdpcap - mean(loggdpcap, na.rm = TRUE),
      povertydm = poverty - mean(poverty, na.rm = TRUE)
    )
} else {
  conf_year_cov <- conf_year %>%
    mutate(
      gdpcap = NA_real_,
      poverty = NA_real_,
      region = NA_character_,
      income = NA_character_,
      loggdpcap = NA_real_,
      loggdpcapdm = NA_real_,
      povertydm = NA_real_
    )
}
setDT(conf_year_cov)

# =======================
# 2. Descriptive analysis
# =======================

# 1. Descriptive statistics for panel coverage and main variables
panel_summary <- tibble(
  Statistic = c(
    "Cell-year observations",
    "Grid cells",
    "Countries",
    "First year",
    "Last year",
    "Conflict incidence",
    "Mean annual temperature",
    "Mean annual temperature anomaly",
    "Mean annual conflict events"
  ),
  Value = c(
    as.character(nrow(conf_year)),
    as.character(n_distinct(conf_year$cell_id)),
    as.character(n_distinct(conf_year$iso3c)),
    as.character(min(conf_year$year, na.rm = TRUE)),
    as.character(max(conf_year$year, na.rm = TRUE)),
    mean(conf_year$conflict, na.rm = TRUE),
    mean(conf_year$temp, na.rm = TRUE),
    mean(conf_year$temp_anomaly, na.rm = TRUE),
    mean(conf_year$conflict_events, na.rm = TRUE)
  )
)
write_booktabs(
  panel_summary,
  file.path(table_dir, "descriptive_panel_summary.tex"),
  "Panel coverage and main variables.",
  "tab:descriptive-panel-summary",
  digits = 3
)

# 2. Descriptive statistics for key variables
variable_summary <- conf_year_cov %>%
  as_tibble() %>%
  select(
    `Conflict indicator` = conflict,
    `Conflict events` = conflict_events,
    `Conflict months` = conflict_months,
    `Temperature` = temp,
    `Climateology` = climatology,
    `Temperature anomaly` = temp_anomaly,
    # `Lagged temperature` = temp_lag,
    # `GDP per capita` = gdpcap,
    # `Log GDP per capita` = loggdpcap
  ) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "value") %>%
  group_by(Variable) %>%
  summarise(
    N = sum(!is.na(value)),
    Mean = mean(value, na.rm = TRUE),
    SD = sd(value, na.rm = TRUE),
    Min = min(value, na.rm = TRUE),
    Median = median(value, na.rm = TRUE),
    Max = max(value, na.rm = TRUE),
    .groups = "drop"
  )
write_booktabs(
  variable_summary,
  file.path(table_dir, "descriptive_variable_summary.tex"),
  "Descriptive statistics for the annual grid-cell panel.",
  "tab:descriptive-variable-summary",
  digits = 3
)

# 3. Top 15 countries by grid-cell total conflict incidence
country_summary <- conf_year %>%
  as_tibble() %>%
  group_by(iso3c) %>%
  summarise(
    `Conflict incidence` = mean(conflict, na.rm = TRUE),
    `Conflict events` = as.character(sum(conflict_events, na.rm = TRUE)),
    `Mean temperature` = mean(temp, na.rm = TRUE),
    `Mean temperature anomaly` = mean(temp_anomaly, na.rm = TRUE),
    `Cells` = as.character(n_distinct(cell_id)),
    .groups = "drop"
  ) %>%
  arrange(desc(`Conflict incidence`), desc(`Conflict events`)) %>%
  slice_head(n = 15)
write_booktabs(
  country_summary,
  file.path(table_dir, "descriptive_top_conflict_countries.tex"),
  "Countries with the highest grid-cell conflict incidence.",
  "tab:descriptive-top-conflict-countries",
  digits = 3
)

# 4. Summary by year
annual_summary <- conf_year %>%
  as_tibble() %>%
  group_by(year) %>%
  summarise(
    conflict_incidence = mean(conflict, na.rm = TRUE),
    mean_temperature = mean(temp, na.rm = TRUE),
    event_count = sum(conflict_events, na.rm = TRUE),
    .groups = "drop"
  )

annual_plot <- ggplot(annual_summary, aes(x = year)) +
  # conflict bars
  geom_col(
    aes(y = conflict_incidence, fill = "Conflict incidence"),
    alpha = 0.80,
    width = 0.75
  ) +
  # temperature line
  geom_line(
    aes(
      y = scales::rescale(
        mean_temperature,
        to = range(conflict_incidence)
      ),
      color = "Mean temperature"
    ),
    linewidth = 0.9
  ) +
  scale_fill_manual(
    values = c("Conflict incidence" = "#e4839e"),
    name = NULL
  ) +
  scale_color_manual(
    values = c("Mean temperature" = "#1b6ca8"),
    name = NULL
  ) +
  scale_y_continuous(
    name = "Conflict incidence",
    sec.axis = sec_axis(
      ~ scales::rescale(
        .,
        to = range(annual_summary$mean_temperature)
      ),
      name = "Mean temperature (°C)"
    )
  ) +
  scale_x_continuous(
    breaks = seq(
      min(annual_summary$year),
      max(annual_summary$year),
      by = 5
    )
  ) +
  labs(
    x = "Year",
    title = "Conflict Incidence and Mean Temperature Over Time"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    plot.title = element_text(face = "bold")
  )
ggsave(
  file.path(
    figure_dir,
    "annual_conflict_temperature.png"
  ),
  annual_plot,
  width = 7.2,
  height = 4.4,
  dpi = 300
)

# 5. Conflict map
cell_summary <- conf_year %>%
  as_tibble() %>%
  group_by(cell_id, latitude, longitude, iso3c) %>%
  summarise(
    conflict_incidence = mean(conflict, na.rm = TRUE),
    conflict_events = sum(conflict_events, na.rm = TRUE),
    mean_temp = mean(temp, na.rm = TRUE),
    mean_temp_anomaly = mean(temp_anomaly, na.rm = TRUE),
    .groups = "drop"
  )

map_xlim <- range(conf$longitude, na.rm = TRUE) + c(-1, 1)
map_ylim <- range(conf$latitude, na.rm = TRUE) + c(-1, 1)
map_bbox <- st_bbox(
  c(xmin = map_xlim[1], xmax = map_xlim[2], ymin = map_ylim[1], ymax = map_ylim[2]),
  crs = st_crs(4326)
)

countries_path <- file.path(raw_data_dir, "ne_10m_admin_0_countries", "ne_10m_admin_0_countries.shp")
world <- NULL
if (file.exists(countries_path)) {
  world <- st_read(countries_path, quiet = TRUE) %>%
    st_transform(4326) %>%
    st_make_valid() %>%
    suppressWarnings(st_crop(map_bbox))
}

base_map <- ggplot()
if (!is.null(world)) {
  base_map <- base_map +
    geom_sf(data = world, fill = "grey95", color = NA)
}

conflict_map <- base_map +
  geom_tile(
    data = cell_summary,
    aes(x = longitude, y = latitude, fill = conflict_incidence),
    width = 1,
    height = 1,
    alpha = 0.90
  ) +
  coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
  scale_fill_gradient(low = "#f7f7f2", high = "#7c2d3a", name = "Incidence") +
  labs(x = NULL, y = NULL, title = "Mean conflict incidence") +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))
temperature_map <- base_map +
  geom_tile(
    data = cell_summary,
    aes(x = longitude, y = latitude, fill = mean_temp),
    width = 1,
    height = 1,
    alpha = 0.90
  ) +
  coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
  scale_fill_gradient(low = "#4f8fc0", high = "#c65f36", name = "Deg. C") +
  labs(x = NULL, y = NULL, title = "Mean annual temperature") +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))

events_map <- base_map +
  geom_tile(
    data = cell_summary,
    aes(x = longitude, y = latitude, fill = log1p(conflict_events)),
    width = 1,
    height = 1,
    alpha = 0.90
  ) +
  coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
  scale_fill_gradient(low = "#f5f5f5", high = "#3b2d5c", name = "log(1 + events)") +
  labs(x = NULL, y = NULL, title = "Total recorded conflict events") +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))


ggsave(file.path(figure_dir, "map_conflict_incidence.png"), conflict_map, width = 7.5, height = 5.2, dpi = 300)
ggsave(file.path(figure_dir, "map_mean_temperature.png"), temperature_map, width = 7.5, height = 5.2, dpi = 300)
ggsave(file.path(figure_dir, "map_conflict_events.png"), events_map, width = 7.5, height = 5.2, dpi = 300)

# =======================
# 3. Regression analysis
# =======================

# data prep for regression: drop rows with missing values in key variables
reg_data <- conf_year_cov %>%
  as_tibble() %>%
  filter(!is.na(temp), !is.na(temp_lag), !is.na(conflict))

# run models with feols
m_base <- feols(
  conflict ~ temp | cell_id + year,
  data = reg_data,
  cluster = ~cell_id
)

m_lag <- feols(
  conflict ~ temp + temp_lag | cell_id + year,
  data = reg_data,
  cluster = ~cell_id
)

m_events <- feols(
  log1p(conflict_events) ~ temp | cell_id + year,
  data = reg_data,
  cluster = ~cell_id
)

m_events_lag <- feols(
  log1p(conflict_events) ~ temp + temp_lag | cell_id + year,
  data = reg_data,
  cluster = ~cell_id
)

m_gdp <- feols(
  conflict ~ temp + temp_lag + temp:loggdpcapdm + temp_lag:loggdpcapdm | cell_id + year,
  data = filter(reg_data, !is.na(loggdpcapdm)),
  cluster = ~cell_id
)

m_events_gdp <- feols(
  log1p(conflict_events) ~ temp + temp_lag + temp:loggdpcapdm + temp_lag:loggdpcapdm | cell_id + year,
  data = filter(reg_data, !is.na(loggdpcapdm)),
  cluster = ~cell_id
)

coef_map <- c(
  "temp" = "Temperature",
  "temp_lag" = "Lagged temperature",
  "temp:loggdpcapdm" = "Temperature x log GDP pc",
  "temp_lag:loggdpcapdm" = "Lagged temperature x log GDP pc"
)

modelsummary(
  list(
    "Incidence" = m_base,
    "Incidence (with lag)" = m_lag,
    # "Conflict + GDP" = m_gdp,
    # "Events + GDP" = m_events_gdp,
    "Events" = m_events,
    "Events (with lag)" = m_events_lag
   
  ),
  coef_map = coef_map,
  output = file.path(table_dir, "regression_results.tex"),
  stars = TRUE,
  statistic = "std.error",
  title = "Temperature and Conflict",
  notes = "All specifications include grid-cell and year fixed effects. Standard errors are clustered by grid cell.",
  gof_omit = "IC|Log|Adj|Within|RMSE"
)

# Function to calculate combined effects of temp and temp_lag
combined_effects <- function(model, interaction_var = NULL) {
  b <- coef(model)
  vc <- vcov(model)
  terms <- c("temp", "temp_lag")
  estimate <- sum(b[terms], na.rm = TRUE)
  se <- sqrt(sum(vc[terms, terms], na.rm = TRUE))

  out <- tibble(
    Effect = "Temperature + lagged temperature",
    Estimate = estimate,
    SE = se,
    `t statistic` = estimate / se
  )

  if (!is.null(interaction_var)) {
    int_terms <- paste0(terms, ":", interaction_var)
    int_estimate <- sum(b[int_terms], na.rm = TRUE)
    int_se <- sqrt(sum(vc[int_terms, int_terms], na.rm = TRUE))
    out <- bind_rows(
      out,
      tibble(
        Effect = paste("Interaction sum:", interaction_var),
        Estimate = int_estimate,
        SE = int_se,
        `t statistic` = int_estimate / int_se
      )
    )
  }
  out
}

# Create a table of combined effects for the main models
effect_table <- bind_rows(
  combined_effects(m_lag) %>% mutate(Model = "Incidence (with lag)"),
  # combined_effects(m_gdp, "loggdpcapdm") %>% mutate(Model = "Conflict + GDP"),
  # combined_effects(m_events_gdp, "loggdpcapdm") %>% mutate(Model = "Events + GDP"),
  combined_effects(m_events_lag) %>% mutate(Model = "Events (with lag)")
) %>%
  select(Model, Effect, Estimate, SE, `t statistic`)
write_booktabs(
  effect_table,
  file.path(table_dir, "regression_combined_effects.tex"),
  "Combined contemporaneous and lagged temperature effects.",
  "tab:regression-combined-effects",
  digits = 3
)

baseline_incidence <- mean(reg_data$conflict, na.rm = TRUE)
base_combined <- sum(coef(m_lag)[c("temp", "temp_lag")], na.rm = TRUE)
diagnostics <- tibble(
  Statistic = c(
    "Regression sample mean conflict incidence",
    "Combined temperature effect",
    "Percent of baseline incidence per 1 degree C",
    "Cells in regression sample",
    "Countries in regression sample"
  ),
  Value = c(
    baseline_incidence,
    base_combined,
    100 * base_combined / baseline_incidence,
    n_distinct(reg_data$cell_id),
    n_distinct(reg_data$iso3c)
  )
)
write_booktabs(
  diagnostics,
  file.path(table_dir, "regression_diagnostics.tex"),
  "Regression diagnostics and scale of the baseline estimate.",
  "tab:regression-diagnostics",
  digits = 3
)

saveRDS(
  list(
    conf_year = conf_year,
    panel_summary = panel_summary,
    variable_summary = variable_summary,
    country_summary = country_summary,
    annual_summary = annual_summary,
    models = list(
      conflict_base = m_base,
      conflict_lag = m_lag,
      # conflict_gdp = m_gdp,
      # events_gdp = m_events_gdp,
      events = m_events,
      events_lag = m_events_lag
    )
  ),
  file.path(output_dir, "completed_analysis_objects.rds")
)

message("Completed descriptive tables, maps, and regressions in: ", output_dir)
