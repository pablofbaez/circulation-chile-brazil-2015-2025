# Reproduction audit — 2026-09-27

Anchor: `Manuscrito_Rienzi_Carrasco_hibrido_v2(1).docx`; checksum in `manuscript_tables_reference.json`.

## Executed in this preparation

- Full extraction from the nine annual Chile files inside CHILE(1).zip and the separate 2024/2025 files; no cached-count shortcut.
- Parsing of all 22 definitive Brazilian sex/year exports.
- Re-export of denominators from the original INE workbook.
- Refit of 126 age/outcome/period Poisson mean models and independent conditional-binomial coefficient checks.
- Recalculation of 128 age-specific/age-adjusted OR rows with/without I83; 176 annual MH rows across four comparators and two definitions.
- Reconstruction and fit of the Chile multivariable logistic model and four separate interaction models from frequency counts of identical predictor/outcome combinations.
- Recalculation of eight Chile diagnostic sensitivity contrasts and masked-record marginal bounds.
- Regeneration of the three statistical figures and machine-readable sources for the 23 manuscript tables.

## Numerical findings

Main counts: Chile 834,780; Brazil 7,409,594. I83 additions: Chile 69,243; Brazil 786,875.

Additional Chile logistic OR: 1.8169946675 (95% CI 1.7932781572–1.8410248341), which rounds to 1.82 (1.79–1.84).

Interaction LR tests: age 261.9487707 (6 df); year 23.7051602 (10 df); region 70.3811130 (16 df); insurance 11.9069883 (6 df). All round to the manuscript values.

`results/validation.json` contains the executed checks. Comparisons with retained CSV results use tolerances 1e-9 relative/1e-10 absolute; manuscript A18 and A19 comparisons use their displayed precision. Twenty-eight pooled rate-ratio identities were checked: RR_venous/RR_arterial equals the corresponding case–case OR.

## Code provenance

Core extraction, composition, rate models and figure code were adapted from the analytical scripts retained with the current manuscript work. Paths were made portable and extraction was extended to export the sufficient statistics needed for supplementary models.

A complete original implementation of the region/insurance model and the narrow diagnostic sensitivities was not recovered. `04_secondary.py` is an independently reconstructed implementation of the manuscript methods; it reproduces the displayed results. Do not describe it as the historical source script. The old early-stage pipeline was inspected but was not bundled as the final pipeline because it omitted these covariates and used older data structures.

The figures reproduce the current estimands and data; they need not be binary-identical to images embedded in Word. Figure S1 uses purple/gray for women/men.

## What this audit does not establish

- Clinical validity of ICD coding, national reporting completeness or independence of readmissions.
- Independent re-execution of the online TabNet queries against the live source. This package validates the supplied export contents.
- The exact administrative period option originally selected in TabNet beyond the headers retained in the exports. Confirm that option in the extraction log before publication.
- Independent literature validation of the curated risk-model table A12. Its text is a manuscript reference, not a computed output.
- A license grant over third-party data or a public repository URL. Neither has been assigned here.

The raw Chile data are not bundled due to their size. Their exact filenames/hashes and processing code are provided; all aggregate numerical inputs needed for the included models are bundled.

## Portability check

The complete aggregate-based pipeline was executed in a separate directory after deleting all generated results. Every implemented validation check passed. No paths to the original chat workspace are needed for this mode.

## Post hoc extension verification — 2026-10-18

The I70–I79 component extension requested on 7 October was added as
`08_other_arterial_posthoc.py` and integrated into `run_pipeline.py` before
final validation. The full aggregate-based pipeline was executed successfully
in a separate directory with no pre-existing results. All twelve validation
groups passed, including the original eight groups and four new post hoc
groups. No live queries or raw Chilean extraction were performed in this
extension verification.

The new checks reproduce 87,899 other-arterial discharges (33,221 women,
54,678 men) and 109,820 venous discharges. All 154 annual sex/age cells for
venous, coronary, and other-arterial diagnoses reconcile with the original
analytical input. The other-arterial MH OR is 1.8437668253 (95% CI
1.81019746–1.87795872), or 1.84 (1.81–1.88) after rounding. The coronary
and other-arterial estimates and intervals match the existing main outputs
at full computational precision. All 56 displayed A21 rate ratios and all
32 displayed A22 OR/CI triplets match the author-supplied manuscript.
Twenty-eight additional pooled rate-ratio identities pass.

The original source snapshots and historical source audit remain unchanged.
The package checksum manifest is regenerated for this updated delivery.
The new reference JSON records the author-supplied A21/A22 manuscript file
and its SHA-256; it is an editorial numerical reference, not a new data input.
