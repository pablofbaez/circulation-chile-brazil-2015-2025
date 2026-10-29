# Analytical specification

## Source processing

Chile originals are semicolon-delimited CSV, decoded as Latin-1, read in chunks as strings. Fields used: `SEXO`, `GRUPO_EDAD`, `ANO_EGRESO`, `DIAG1`; supplementary logistic processing additionally uses `REGION_RESIDENCIA` and `PREVISION`. No individual identifier is used or generated. No deduplication is performed: duplicate-looking rows may represent different discharge events.

Normalize DIAG1: strip whitespace, uppercase, remove decimal point. Retain principal diagnoses beginning I plus two digits (I00–I99). Exclude records with unusable adult age or unrecorded male/female sex. Records jointly masked in age, sex and year (`*`) are counted separately. An age label in days/months is pediatric; otherwise use its lower bound. Collapse adult labels to 20–29, 30–39, 40–49, 50–59, 60–69, 70–79, ≥80. This harmonizes five-year 2021 bands with the other years. Sex 1/HOMBRE = Men and 2/MUJER = Women. Eligible records' ANO_EGRESO must match the source year. The 2025 ERROR field is not an additional selection filter.

Each eligible diagnosis is classified into PE (I26), I80–I82, varices (I83), coronary (I20–I25), cerebrovascular (I60–I69), other arterial (I70–I79), or outside selected. Detailed ICD counts retain subcodes to support narrower sensitivity definitions. I63.6/I67.6 are not moved from cerebrovascular to venous in the main analysis.

Brazil uses 22 supplied TabNet exports, one year and sex per file. Read `Período` and `Sexo` from file headers rather than guessing from filenames. Seven age columns are checked, as are row totals and the sum of disjoint morbidity categories. The full exact category mapping is in `data/processed/auditoria_fuentes.json`. The reported download date 2026-09-20 was supplied by the author. The exact administrative meaning of the selected TabNet period must be confirmed from extraction documentation; do not silently equate it with Chile discharge year.

## Case–case odds ratios

For each age: a = venous women, b = arterial women, c = venous men, d = arterial men. OR = ad/(bc). Log-Wald variance = 1/a + 1/b + 1/c + 1/d. CI = exp(log OR ± 1.95996398454 × SE). No zero-cell correction is needed for these observed contrasts.

For Mantel–Haenszel, within stratum i define n=a+b+c+d, R=sum(ad/n), S=sum(bc/n). OR_MH=R/S. Robins–Breslow–Greenland variance of log OR_MH:

`sum((a+d)ad/n²)/(2R²) + sum(((a+d)bc+(b+c)ad)/n²)/(2RS) + sum((b+c)bc/n²)/(2S²)`.

Use the same seven strata in each country and in each annual analysis. Weights are contrast/country specific. All eleven years are in the main period; Brazil excluding 2023 is sensitivity. I83 is added only in the sensitivity definition.

## Proportions and direct standardization

Venous fraction = V/(V+A), within sex/age/country. Figure S1 uses nominal binomial Wilson 95% intervals. Directly standardized fractions use each country's combined female+male selected sample age distribution; women and men have identical weights within a country, but countries have different weights. Absolute percentage-point differences are calculated before rounding.

## Chile rates

INE base 2024 original series: retain June 30 populations for 2015–2025, map H=Men, M=Women, sum single ages into the same adult bands. Absolute annual rate = discharges/population × 100,000. Period rate = sum annual discharges/sum annual populations × 100,000. Nominal unscaled count-based intervals for absolute rates use chi-square/Garwood limits; these additional output columns are not readmission-corrected intervals.

For each diagnostic outcome and age band separately:

`Y_sex,year ~ Poisson(mu)`

`log(mu) = log(population) + intercept + beta × female + categorical_year`.

Full period: 22 cells, 12 parameters, 10 residual df. RR = exp(beta), common over years within that age/outcome model. Pearson phi = sum((Y-mu)²/mu)/df. Main covariance = inverse information × max(1,phi). Main CI = exp(beta ± t_(0.975,df) × scaled SE). Report `quasi_low` and `quasi_high`; `low/high` are unscaled Poisson audit columns. Excluding 2025 leaves 20 cells, 11 parameters and 9 residual df.

The sex coefficient is independently checked with the conditional binomial score in year strata, with offset log(female population/male population). Separate outcomes: main venous, PE, I80–I82, I83, venous+I83, pooled arterial and three arterial territories. Fixed internal Chile 2015 adult population weights produce an additional age-standardized annual-rate output; this is an ancillary result, not an additional main-paper model.

## Additional Chile logistic model

`logit(P(venous | selected discharge)) = intercept + female + C(age) + C(year) + C(region) + C(insurance)`.

All 834,780 discharges contribute. Region codes 01–16 retain their identities; 99 and Extranjero are grouped as Unknown (17 levels). Insurance retains codes 1,2,3,4,5,96,99 (7 levels). Grouping identical predictor/outcome rows with frequency counts preserves the individual Bernoulli log likelihood, coefficient estimates and information matrix. It does not identify patients or correct readmission dependence.

Each interaction (female×age/year/region/insurance) is added separately to the same base model. LR = 2(log likelihood interaction − log likelihood base), with df equal to the added parameter count. Binomial combination constants are omitted consistently and cancel in LR differences. The reconstructed model reproduces all four manuscript LR tests at their displayed precision.

## Diagnostic sensitivities and masking

A19 code sets are explicit in `04_secondary.py`: exclude I69 and I80.0; narrower venous I26 + I80.1–I80.9 + I82 versus I21–I22 + I60–I64 + I74; I26 versus I21; each arterial territory; add I83. These are diagnostic definitions, not proof of urgency or acute clinical status.

Masked-record bounds use 3,325 unallocated selected venous records. For each age separately, add all to men for the lower bound, or all to women for the upper bound. Other strata are not simultaneously assigned the same records. Bounds are hypothetical marginal scenarios, not imputations or confidence intervals.

## Limits

No incidence, first-event identification, individual risk prediction, causal hormonal mechanism, score-performance validation or population rate for Brazil is estimated. Age strata are cross-sectional. All P values and CIs are nominal. Dispersion scaling does not replace clustered patient inference. Source coding and completeness were not externally validated by this pipeline.

## Post hoc decomposition of other arterial diagnoses

The exploratory extension requested on 7 October 2026 uses
`chile_detailed_diagnoses.csv` and the same pooled 2015–2025 INE population
exposures. It partitions I70–I79 into I70/I73, I71, I72, I74, and I77–I79.
Every annual sex/age cell must reconcile to the original other-arterial
count; nonzero codes outside this partition cause an error. Comparators
are I71, I70–I79 excluding I71, all I70–I79, and coronary I20–I25. Venous
counts always include I26 and I80–I82 and exclude I83.

A21 uses pooled discharge rates and female-to-male rate ratios; these RRs
are not year-adjusted, and no confidence intervals are assigned to them.
A22 uses age-specific log-Wald case–case OR intervals and the same seven-age
Mantel–Haenszel/RBG estimator defined above. At a fixed age, pooled
RR_venous / RR_comparator equals the case–case OR because the sex-specific
population exposures cancel. The age-adjusted MH OR must not be obtained
by dividing separately standardized rate ratios. No zero-cell correction
is applied; an observed zero causes explicit failure. All these results
describe hospital-discharge episodes, not disease incidence or mechanisms.
