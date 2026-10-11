# Data dictionary

UTF-8 CSV, comma separator, period decimal separator. Integers are unformatted counts. Age strings use en dashes and ≥; preserve UTF-8.

| File in data/processed | Unit | Columns |
|---|---|---|
| recuentos_por_ano_sexo_edad_diagnostico.csv | Country × year × sex × age × diagnosis; 1,848 rows | country (Chile/Brazil), year (2015–2025), sex (Women/Men), age (7 adult bands), diagnostic_group (PE/I80_I82/Varices_I83/Coronary/Cerebrovascular/Other_arterial), count (discharge events, zero cells explicit) |
| denominadores_INE_2015_2025.csv | Chile year × sex × age; 154 rows | year, sex, age, population (June 30 population, persons) |
| chile_detailed_diagnoses.csv | Chile year × sex × age × principal ICD code | year, sex, age, code (uppercase, no decimal), count; includes all eligible I00–I99 diagnoses, not just selected contrast |
| chile_logistic_cells.csv | Chile year × sex × age × region × insurance × outcome | year, sex, age, region (source code/string), insurance (source code), venous (1 for I26/I80–I82, 0 for selected arterial), count |
| chile_masked_counts.csv | Source year × principal ICD three-character code | year (source file year, not imputed patient field), code, count (jointly masked age/sex/year records within circulatory diagnoses) |
| auditoria_fuentes.json | Source provenance and flow | Exact filenames, uncompressed source SHA-256 hashes, byte counts, Brazil headers/mapping, Chile exclusions by year |
| population_source.json | INE provenance | Original name, SHA-256, reference date, series base |

No patient identifiers or synthetic records are included. The logistic-cell table is a frequency table, not an individual-level patient dataset. No claim of statistical anonymization beyond the supplied public source data is made.

## Model outputs

- `OR`, `MH_OR`, `OR_age_adjusted_MH`: female-to-male case–case odds ratios, venous relative to the specified arterial comparator.
- `RR_year_adjusted`: female-to-male population-denominated discharge-rate ratio adjusted for categorical year within age.
- `RR_pooled`: corresponding ratio of pooled absolute rates; not year-adjusted.
- `quasi_low`, `quasi_high`: main dispersion-scaled t-based 95% interval.
- `low`, `high` in the RR output: unscaled Poisson normal intervals kept for audit, not the manuscript's main intervals.
- `rate_100k`: discharge events per 100,000 population-years.
- `pearson_dispersion`: unbounded Pearson estimate; model covariance actually uses max(1, value).
- `residual_df`: 10 for 2015–2025, 9 for 2015–2024.

Counts in `Varices_I83` remain visible in the compositional file even on `without_varices` rows; those counts are NOT included in `selected_total` or `venous_count` for that definition.
