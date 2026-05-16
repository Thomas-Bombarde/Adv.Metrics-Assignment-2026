# Candidate HiClimR Region Clusterings

Input: `code_and_data/analysis_data/weather_conflict_panel.rdata`, object `conf`.

Clustering variable: monthly temperature anomaly, `temp - climatology`, arranged as cells by months.
HiClimR preprocessing: linear detrending and row standardization for all proposals.
Conflict variables are not used in the correlation distance because many cells have no conflict variation; they are reported as diagnostics in `cluster_summaries.csv`.

## Proposals

- `temp_ward_k04`: 4 groups, method `ward`, hybrid `FALSE`, contiguity `0`, nPC `none`. Broad Ward clusters on detrended, standardized monthly temperature anomalies. Mean intra-cluster correlation: 0.683; mean inter-cluster correlation: 0.239.
- `temp_regional_k08`: 8 groups, method `regional`, hybrid `FALSE`, contiguity `0`, nPC `none`. Regional-linkage clusters emphasizing separation of regional mean time series. Mean intra-cluster correlation: 0.763; mean inter-cluster correlation: 0.213.
- `temp_ward_contig_k12`: 12 groups, method `ward`, hybrid `FALSE`, contiguity `0.35`, nPC `none`. Ward clusters with a moderate geographic contiguity weight. Mean intra-cluster correlation: 0.816; mean inter-cluster correlation: 0.210.
- `temp_ward_pca_k16`: 16 groups, method `ward`, hybrid `FALSE`, contiguity `0`, nPC `12`. Ward clusters after PCA filtering of the anomaly matrix. Mean intra-cluster correlation: 0.849; mean inter-cluster correlation: 0.229.
- `temp_hybrid_k20`: 20 groups, method `ward`, hybrid `TRUE`, contiguity `0`, nPC `none`. Hybrid Ward-regional clusters: Ward tree with regional reconstruction above 60 clusters. Mean intra-cluster correlation: 0.855; mean inter-cluster correlation: 0.208.

## Outputs

- `code_and_data/output/clusterings/proposal_specs.csv`: technique and parameter grid.
- `code_and_data/output/clusterings/cluster_assignments.csv`: one row per grid cell and proposal.
- `code_and_data/output/clusterings/cluster_summaries.csv`: cluster-level geography, temperature, validation, and conflict diagnostics.
- `code_and_data/output/clusterings/hiclimr_objects.rds`: fitted HiClimR objects for reuse.
- `code_and_data/output/clusterings/maps/*.png`: one map per proposal plus `all_clusterings_overview.png`.
