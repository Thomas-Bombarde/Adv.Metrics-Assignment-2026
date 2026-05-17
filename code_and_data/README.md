# Code and Data for the BCH/HiClimR Report

This folder is the reproducible entry point for the report:

`Fixed-G Inference with Climate-Based Clusters`

It contains the report source, the derived analysis data, the scripts, and the generated outputs needed to rerun the analysis. It is intentionally smaller than the full project directory. It does not duplicate the raw replication archive or the literature PDFs.

All commands below assume that the working directory is the report R project root, namely the folder containing `report.Rproj`. Open `report.Rproj` first, or set the working directory to that folder before running the scripts.

## Folder Structure

```text
code_and_data/
├── README.md
├── install_packages.R
├── Report_HiClimR_BCH.Rmd
├── Report_HiClimR_BCH.pdf
├── references.bib
├── helpers.R
├── run_sim.R
├── make_figures.R
├── simulate_group_tradeoff.R
├── propose_region_clusterings.R
├── bch_conflict_application.R
├── score_independence_diagnostics.R
├── analysis_data/
│   └── weather_conflict_panel.rdata
├── results/
│   └── sim_results.rds
└── output/
    ├── figures/
    ├── clusterings/
    └── tables/
```

## Data Provenance

The file `code_and_data/analysis_data/weather_conflict_panel.rdata` is the
derived weather-conflict panel used in the report. It contains the object `conf`.

The dataset was obtained from the replication package 
for Burke et al. (2024), `Will Wealth Weaken Weather Wars?`.
In the parent project, the raw package is stored under:

```text
analysis_data/5W_replication/
```

The raw replication archive is not duplicated here because
it is large and not needed to reproduce the report from the derived panel. 
The included panel keeps the public weather-conflict fields 
used for this assignment: grid-cell identifiers, coordinates,
country identifiers, monthly temperature, climatology,
conflict counts, and conflict indicators.

## Software

The analysis is written in R. Install the required packages with:
```bash
Rscript code_and_data/install_packages.R
```
This creates and uses a project-local library,
`code_and_data/.r-lib`. The shared path helper adds that local library to `.libPaths()` when scripts are run from the project root.

Main packages used:
```text
data.table, dplyr, ggplot2, ggrepel, HiClimR, knitr, patchwork,
readr, rmarkdown, rnaturalearth, scales, sf, tibble, tidyr
```

Rendering the PDF also requires a working LaTeX installation.

## Reproduction Order

From the report project root, run:
```bash
Rscript code_and_data/run_sim.R 500
Rscript code_and_data/make_figures.R
Rscript code_and_data/simulate_group_tradeoff.R 400
Rscript code_and_data/propose_region_clusterings.R
Rscript code_and_data/bch_conflict_application.R
Rscript code_and_data/score_independence_diagnostics.R
Rscript -e 'rmarkdown::render("code_and_data/Report_HiClimR_BCH.Rmd", output_file = "Report_HiClimR_BCH.pdf", clean = TRUE)'
```

The first and fourth commands are the slowest:
- `run_sim.R` regenerates the main Monte Carlo file `results/sim_results.rds`.
- `propose_region_clusterings.R` reruns HiClimR and regenerates cluster assignments, cluster summaries, and maps.

The folder already includes the current outputs used in the report, so the report can be rendered directly if those outputs do not need to be regenerated.

If you only want the report PDF and the already-generated figures and tables, the minimal command is:

```bash
Rscript -e 'rmarkdown::render("code_and_data/Report_HiClimR_BCH.Rmd", output_file = "Report_HiClimR_BCH.pdf", clean = TRUE)'
```

## Script Roles

`helpers.R`  
Shared Monte Carlo functions: 
spatial moving-average data generation,
BCH group construction, 
CCE, HC0, and Bartlett HAC estimators.

`run_sim.R`  
Runs the main fixed-\(G\) 
Monte Carlo simulation and saves `results/sim_results.rds`.

`make_figures.R`  
Creates the simulation figures `figA_t_density.png`, 
`figB_rejection_rates.png`, and `figC_Vhat_density.png`.

`simulate_group_tradeoff.R`  
Creates the \(G=4,8,12,16,20\) Monte Carlo trade-off 
figure used to motivate the empirical application.

`propose_region_clusterings.R`  
Uses HiClimR to build candidate climate groupings 
from monthly temperature anomalies. 
Outputs are saved under `code_and_data/output/clusterings/`.

`bch_conflict_application.R`  
Computes the empirical fixed-effects
coefficient and BCH-style inference across HiClimR groupings. 
Outputs are saved under `code_and_data/output/tables/` and `code_and_data/output/figures/`.

`score_independence_diagnostics.R`  
Computes score-correlation diagnostics for HiClimR groups,
administrative macroregions, and country clusters.

`Report_HiClimR_BCH.Rmd`  
R Markdown source for the final report.

## Current Report Outputs

Key generated files included here:

```text
code_and_data/output/figures/figA_t_density.png
code_and_data/output/figures/figB_rejection_rates.png
code_and_data/output/figures/figC_Vhat_density.png
code_and_data/output/figures/figD_group_tradeoff.png
code_and_data/output/figures/bch_application_variance_by_g.png
code_and_data/output/clusterings/cluster_assignments.csv
code_and_data/output/clusterings/cluster_summaries.csv
code_and_data/output/clusterings/proposal_specs.csv
code_and_data/output/clusterings/maps/all_clusterings_overview.png
code_and_data/output/clusterings/maps/temp_ward_k04.png
code_and_data/output/tables/bch_application_variance_results.csv
code_and_data/output/tables/score_independence_diagnostics.csv
```

## Notes on Interpretation

The report uses HiClimR to define candidate spatial groups 
from temperature-anomaly co-movement.
Conflict outcomes are not used to choose the groups. 
BCH validity still concerns regression scores, 
so the score diagnostics are part of the empirical 
argument rather than an optional robustness table.

## Practical Notes

- Keep working from the report project root so the relative paths in the scripts resolve correctly.
- `code_and_data/.r-lib` is a project-local library path used by the install script and should not be committed.
- The committed bundle is designed to be sufficient for rerunning the report without the raw replication archive.
