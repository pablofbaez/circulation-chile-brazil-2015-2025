# Sex differences in the arterial and venous composition of hospitalizations across adulthood
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23306108.svg)](https://doi.org/10.5281/zenodo.23306108)

Reproducibility package for the article:

> Carrasco-Wong I, Bertini A, Baez P, Amador C, Gutiérrez J, Rotarou E, Díaz-Navarro R, Illanes SE, Chiarello DI. **Sex Differences in the Arterial and Venous Composition of Hospitalizations Across Adulthood in Chile and Brazil, 2015–2025.** [TODO: journal, year, DOI once available]

The reference tables used for regression checks are recorded in `docs/manuscript_tables_reference.json`, together with the SHA-256 of the manuscript version they were taken from.

## Authors and affiliations

Ivo Carrasco-Wong, Ayleen Bertini, Pablo Baez, Cristian Amador, Jaime Gutiérrez, Elena Rotarou, Rienzi Díaz-Navarro, Sebastián E. Illanes, Delia I. Chiarello.

- Facultad de Ciencias, Universidad San Sebastián, Santiago, Chile.
- School of Medical Technology, Faculty of Sciences, Universidad San Sebastián, Santiago, Chile.
- PhD Program in Chronic Disease, Faculty of Sciences, Universidad San Sebastián, Santiago, Chile.
- Faculty of Dentistry, Universidad San Sebastián, Santiago, Chile.
- Laboratorio de Fisiología Cardiovascular, Departamento de Medicina Interna, Escuela de Medicina, Facultad de Medicina, MEDING, Universidad de Valparaíso, Chile.
- Reproductive Biology Program, Center for Biomedical Research and Innovation (CiiB), Universidad de los Andes, Santiago, Chile.
- IMPACT, Center of Interventional Medicine for Precision and Advanced Cellular Therapy, Faculty of Medicine, Universidad de los Andes, Santiago, Chile.
- Department of Obstetrics and Gynecology, Faculty of Medicine, Universidad de los Andes, Santiago, Chile.

Contact: Delia I. Chiarello, delia.chiarello@uss.cl

## Scope

Adults aged ≥20 years, 2015–2025 inclusive. Main venous definition: I26 + I80–I82. Arterial comparators: I20–I25, I60–I69 and I70–I79, separately and pooled. I83 is included only in sensitivity analyses. Chile: 834,780 selected discharges. Brazil: 7,409,594 SUS-financed hospitalizations. Units are discharge events, not unique patients or first events. Brazil is not a census of all hospitalizations.

## Quick reproduction from included data

Requires Python 3.11 or newer. The exact package versions tested are in `requirements.txt`.

```bash
python -m venv .venv
# Linux/macOS:
source .venv/bin/activate
# Windows PowerShell instead: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python run_pipeline.py
```

This command uses the included aggregate counts and INE population denominators. It refits the models, regenerates numeric table sources, creates Figures 1, 2 and S1 with confidence intervals and the graphical abstract, and compares the estimates with the reference outputs. Results are written to `results/`. See `results/validation.json` for computational checks. No external account, API or network is required after dependency installation.

Tested on: Ubuntu 24.04 LTS and Python 3.12.10. Approximate run time: 2 minutes.

## Rebuild from original files

1. Put the eleven Chilean source datasets in one input folder using the packaging and names listed in `docs/SOURCES.csv`: `CHILE(1).zip` contains 2015–2023; the 2024 and 2025 CSV files are separate. Do not substitute an older compiled dataset.
2. Add `BRASIL TODO(1).zip` and the two 2015 sex-specific exports to that folder. Copies of these aggregate Brazilian files are included in `data/raw/`.
3. Use the included original INE workbook (base 2024). Its original name and checksum are in `docs/SOURCES.csv`.
4. Execute, replacing the input folder with its actual path:

```bash
python run_pipeline.py --raw-dir /path/to/originals --ine data/raw/INE_base2024_original.xlsx
```

The raw mode rebuilds all analytical aggregates from the original sources and then runs the same models. Large Chilean source files are not bundled. Public acquisition portals and every source/member filename and checksum are listed in `docs/SOURCES.csv`. Sources were downloaded on March 2026. Exact historical download URLs were not retained; the checksums identify the actual source snapshots used. A future source revision may not reproduce the same counts.

## Files and order of execution

`run_pipeline.py` is the single entry point and defines the execution order. Script numbers are identifiers, not the execution order (for example, `00` runs after `01`, and `06_validate.py` runs last); the table below lists the scripts in the order they are executed.

| Step | Script | Inputs | Main output |
|---|---|---|---|
| Optional raw 1 | `scripts/01_extract.py` | Chile originals and 22 Brazil exports | Harmonized counts; detailed Chile ICD counts; logistic sufficient statistics; masked counts; source audit |
| Optional raw 2 | `scripts/00_population.py` | Original INE XLSX | June 30 populations, year × sex × age |
| 1 | `scripts/02_composition.py` | Main counts | Age-specific/MH ORs, annual MH ORs, standardized proportions, Brazil excluding 2023 |
| 2 | `scripts/03_rates.py` | Chile counts and populations | Absolute rates, year-adjusted RRs and dispersion-scaled intervals; excluding 2025; subtype/I83 sensitivities |
| 3 | `scripts/04_secondary.py` | Detailed diagnostic/covariate counts | Logistic model, LR interactions, narrower ICD sensitivities, full circulatory composition, masking bounds |
| 4 | `scripts/05_figures.py` | Numeric model outputs | Figures 1, 2 and S1, vector PDF and PNG |
| 5 | `scripts/07_tables.py` | Numeric results and curated reference text | CSV sources for all 3 main and 20 supplemental tables (A1–A20); Wilson proportions |
| 6 | `scripts/08_other_arterial_posthoc.py` | Detailed Chile ICD counts and INE denominators | Post hoc I70–I79 decomposition; Tables A21 and A22, counts, rates, metadata |
| 7 | `scripts/09_graphical_abstract.py` | `results/tables/A20.csv`; Natural Earth country outlines | Graphical abstract, PDF and PNG |
| 8 | `scripts/06_validate.py` | New outputs and independent reference results | `validation.json` (12 validation groups) |

A12 (risk models), A14 (diagnostic mapping), and A15 (data availability) contain curated text transcribed from the manuscript. They are not statistical models refitted by this pipeline. Table CSVs preserve numeric precision and use long format where useful; their layout is not a reproduction of the typeset tables.

## What is included

- `data/processed/`: actual aggregate analytical data, not simulated data. Data dictionary: `docs/DATA_DICTIONARY.md`.
- `data/raw/`: 22 Brazilian TabNet exports within their original containers, and the original INE workbook under a portable filename.
- `data/reference/`: retained numeric estimates used for regression checks, not the inputs to the fitted models.
- `data/external/`: Chile and Brazil outlines from Natural Earth (1:110m Admin 0 Countries, v5.1.2; public domain), used only in the graphical abstract.
- `scripts/`: executable Python code. NumPy/SciPy implement the grouped binomial and Poisson likelihoods explicitly.
- `results/`: regenerated model outputs, table sources, figures and validation report.
- `docs/`: equations, preprocessing, exact filenames, table map and audit record.

## Statistical interpretation

Case–case OR = (venous women / arterial women) / (venous men / arterial men). A larger venous proportion need not imply a higher population-denominated venous hospitalization rate. The Chile rate analysis distinguishes these quantities. Readmissions remain in the counts. Nominal confidence intervals do not correct for correlated readmissions, coding error, missingness or differences in coverage. Population projections are treated as fixed. No Brazilian population rate model is fitted.

## Post hoc analysis: other arterial diseases in Chile

`scripts/08_other_arterial_posthoc.py` is part of `python run_pipeline.py` and can also be run alone with `python scripts/08_other_arterial_posthoc.py`. It reads the included detailed diagnostic counts and INE denominators; no raw Chilean records or network connection are needed. It partitions I70–I79 into I70/I73, I71, I72, I74 and I77–I79, and also examines I70–I79 excluding I71. All comparisons use I26 plus I80–I82 as the venous numerator. This exploratory analysis was added during manuscript revision and was not prespecified.

- `results/tables/A21.csv`: wide table of pooled female-to-male discharge rate ratios.
- `results/tables/A22.csv`: long table of age-specific and age-adjusted case–case ORs and 95% CIs.
- `results/chile_other_arterial_posthoc_rates.csv`: counts, population-years, rates and pooled RRs.
- `results/chile_other_arterial_posthoc_counts.csv`: year × sex × age counts for all components.
- `results/chile_other_arterial_posthoc_odds.csv`: full-precision ORs, CIs and four-cell counts.
- `results/chile_other_arterial_posthoc_metadata.json`: definitions, inputs and interpretation.

The post hoc validation checks reconcile annual diagnostic counts with the main analytical input, recover all coronary and other-arterial A3 ORs, CIs and MH summaries, check the pooled rate-ratio identity, and compare A21 and A22 with the values reported in the manuscript. The manuscript reference values are kept under `data/reference/` together with their source file and checksum. See `docs/POSTHOC_OTHER_ARTERIAL.md`.

## Graphical abstract

`scripts/09_graphical_abstract.py` draws the graphical abstract from Table A20 (`results/tables/A20.csv`). Every number shown is computed from that table, and the statements in the figure are checked against the data; the script stops with an error if a check fails. The country outlines are read from `data/external/`, so no network connection is needed.

## Data availability and sources

All inputs are public, aggregate or de-identified administrative data. Source custodians:

- Chile, hospital discharges: Departamento de Estadísticas e Información de Salud (DEIS), Ministerio de Salud.
- Brazil, SUS-financed hospitalizations: DATASUS (SIH/SUS, TabNet), Ministério da Saúde.
- Chile, population projections (base 2024): Instituto Nacional de Estadísticas (INE).
- Country outlines in the graphical abstract: Natural Earth (public domain).

Source data remain subject to the terms of their custodians; this repository asserts no new license over third-party data.

## License

Code is released under the MIT License (see `LICENSE`). The license covers the code and documentation written by the authors, not the third-party source data in `data/raw/` or the aggregates derived from them.

## How to cite

Please cite the article above and this repository release (https://doi.org/10.5281/zenodo.23306108), as well as DEIS, DATASUS and INE as source custodians. Citation metadata are in `CITATION.cff`.

## Use of AI tools

AI-assisted tools were used in the preparation and verification of the code. The authors reviewed all code and outputs and take full responsibility for them.
