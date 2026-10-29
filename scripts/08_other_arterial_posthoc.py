"""Post hoc Chile I70–I79 decomposition supporting manuscript Tables A21/A22.

Uses aggregate discharge episodes and pooled INE population-years, 2015–2025.
No annual adjustment, incidence estimation, or patient deduplication is performed.
Run directly or through run_pipeline.py. Results retain unrounded precision.
"""
from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
AGES = ['20–29', '30–39', '40–49', '50–59', '60–69', '70–79', '≥80']
SEXES = ['Women', 'Men']
YEARS = list(range(2015, 2026))
COMPONENTS = {
    'I70_I73': {'I70', 'I73'},
    'I71': {'I71'},
    'I72': {'I72'},
    'I74': {'I74'},
    'I77_I79': {'I77', 'I78', 'I79'},
}
COMPARATORS = ['I71', 'Other_arterial_excluding_I71', 'Other_arterial', 'Coronary']
Z = norm.ppf(0.975)


def require(condition, message):
    """Fail even if Python is invoked with -O."""
    if not bool(condition):
        raise ValueError(message)


def crude_or(a, b, c, d):
    cells = np.asarray([a, b, c, d], dtype=float)
    require(np.isfinite(cells).all() and (cells > 0).all(),
            'Positive finite cells required; no automatic zero-cell correction is applied.')
    estimate = a * d / (b * c)
    limits = np.exp(np.log(estimate) + np.array([-1, 1]) * Z * np.sqrt(sum(1 / cells)))
    return float(estimate), float(limits[0]), float(limits[1])


def mantel_haenszel(a, b, c, d):
    a, b, c, d = [np.asarray(x, dtype=float) for x in (a, b, c, d)]
    cells = np.array([a, b, c, d])
    require(cells.shape == (4, len(AGES)) and np.isfinite(cells).all()
            and (cells > 0).all(), 'MH requires seven positive age-stratum tables.')
    n = a + b + c + d
    r, s = np.sum(a * d / n), np.sum(b * c / n)
    estimate = r / s
    variance = (
        np.sum((a + d) * a * d / n**2) / (2 * r**2)
        + np.sum(((a + d) * b * c + (b + c) * a * d) / n**2) / (2 * r * s)
        + np.sum((b + c) * b * c / n**2) / (2 * s**2)
    )
    limits = np.exp(np.log(estimate) + np.array([-1, 1]) * Z * np.sqrt(variance))
    return float(estimate), float(limits[0]), float(limits[1])


def calculate(root=ROOT):
    data = root / 'data' / 'processed'
    detail = pd.read_csv(data / 'chile_detailed_diagnoses.csv', dtype={'code': str})
    population = pd.read_csv(data / 'denominadores_INE_2015_2025.csv')
    main = pd.read_csv(data / 'recuentos_por_ano_sexo_edad_diagnostico.csv')
    keys = ['year', 'sex', 'age']
    expected = pd.MultiIndex.from_product([YEARS, SEXES, AGES], names=keys)
    for name, frame, value, unique in [
        ('diagnoses', detail, 'count', keys + ['code']),
        ('population', population, 'population', keys),
    ]:
        require(not frame[unique + [value]].isna().any().any(), f'Missing {name} values.')
        require(not frame.duplicated(unique).any(), f'Duplicate {name} cells.')
        require(set(frame.year) == set(YEARS), f'Unexpected {name} years.')
        require(set(frame.sex) == set(SEXES) and set(frame.age) == set(AGES),
                f'Unexpected {name} sex/age labels.')
        require(pd.api.types.is_numeric_dtype(frame[value])
                and np.isfinite(frame[value]).all() and (frame[value] >= 0).all(),
                f'Invalid {name} values.')
        if value == 'count':
            require((frame[value] == np.floor(frame[value])).all(), 'Counts must be integers.')
    require(population.set_index(keys).index.sort_values().equals(expected.sort_values()),
            'Population must have exactly one positive cell for every year/sex/age.')
    require((population.population > 0).all(), 'Population denominators must be positive.')
    detail['code'] = detail.code.str.strip().str.upper().str.replace('.', '', regex=False)
    require(not detail.duplicated(keys + ['code']).any(), 'Duplicate normalized diagnostic codes.')
    require(detail.code.str.match(r'^I\d{2}').all(), 'Expected Chapter IX diagnostic codes.')
    block = detail.code.str[:3]
    number = pd.to_numeric(block.str[1:])
    masks = {name: block.isin(codes) for name, codes in COMPONENTS.items()}
    masks.update({
        'Venous': (block == 'I26') | number.between(80, 82),
        'Other_arterial_excluding_I71': number.between(70, 79) & (block != 'I71'),
        'Other_arterial': number.between(70, 79),
        'Coronary': number.between(20, 25),
    })
    covered = block.isin(set().union(*COMPONENTS.values()))
    require(detail.loc[number.between(70, 79) & ~covered, 'count'].sum() == 0,
            'Nonzero I70–I79 codes outside the specified component partition.')
    annual = pd.DataFrame({
        name: detail.loc[mask].groupby(keys)['count'].sum().reindex(expected, fill_value=0)
        for name, mask in masks.items()
    })
    require((annual[list(COMPONENTS)].sum(axis=1) == annual.Other_arterial).all(),
            'Components do not sum to I70–I79 within year/sex/age.')
    require((annual.Other_arterial - annual.I71 == annual.Other_arterial_excluding_I71).all(),
            'Non-aortic counts do not reconcile.')
    # Check against a separately supplied analytical input, at all 154 cells.
    main = main[main.country == 'Chile']
    require(not main.duplicated(keys + ['diagnostic_group']).any(), 'Duplicate main Chile cells.')
    for name, groups in [('Venous', ['PE', 'I80_I82']),
                         ('Other_arterial', ['Other_arterial']), ('Coronary', ['Coronary'])]:
        counts = main[main.diagnostic_group.isin(groups)].groupby(keys)['count'].sum()
        require(np.array_equal(annual[name].to_numpy(), counts.reindex(expected).to_numpy()),
                f'Detailed and main-input {name} counts disagree by year/sex/age.')
    counts = annual.groupby(['sex', 'age']).sum()
    populations = population.groupby(['sex', 'age']).population.sum()
    rows = []
    for age in AGES:
        for name in masks:
            women = int(counts.loc[('Women', age), name])
            men = int(counts.loc[('Men', age), name])
            pw = float(populations.loc[('Women', age)])
            pm = float(populations.loc[('Men', age)])
            require(men > 0, f'Zero male denominator count in {name}/{age}.')
            rows.append([age, name, women, men, pw, pm, 100000 * women / pw,
                         100000 * men / pm, (women / pw) / (men / pm)])
    rates = pd.DataFrame(rows, columns=[
        'age', 'diagnostic_group', 'women_count', 'men_count', 'women_population_years',
        'men_population_years', 'women_rate_100k', 'men_rate_100k', 'pooled_RR',
    ])
    odds = []
    for comparator in COMPARATORS:
        arrays = [np.array([counts.loc[(sex, age), group] for age in AGES], float)
                  for sex, group in [('Women', 'Venous'), ('Women', comparator),
                                     ('Men', 'Venous'), ('Men', comparator)]]
        for i, age in enumerate(AGES):
            cells = [x[i] for x in arrays]
            odds.append([comparator, age, *map(int, cells), *crude_or(*cells)])
        odds.append([comparator, 'Age_adjusted_MH', *[int(x.sum()) for x in arrays],
                     *mantel_haenszel(*arrays)])
    odds = pd.DataFrame(odds, columns=[
        'arterial_comparator', 'age_or_estimator', 'venous_women', 'arterial_women',
        'venous_men', 'arterial_men', 'OR', 'CI95_lower', 'CI95_upper',
    ])
    a21_order = [*COMPONENTS, 'Other_arterial_excluding_I71', 'Other_arterial', 'Coronary']
    a21 = rates.pivot(index='age', columns='diagnostic_group', values='pooled_RR').reindex(
        index=AGES, columns=a21_order)
    a21.columns.name = None
    return rates, odds, a21, annual


def main():
    rates, odds, a21, annual = calculate()
    out = ROOT / 'results'
    tables = out / 'tables'
    tables.mkdir(parents=True, exist_ok=True)
    rates.to_csv(out / 'chile_other_arterial_posthoc_rates.csv', index=False)
    odds.to_csv(out / 'chile_other_arterial_posthoc_odds.csv', index=False)
    annual.reset_index().to_csv(out / 'chile_other_arterial_posthoc_counts.csv', index=False)
    a21.to_csv(tables / 'A21.csv')
    odds.to_csv(tables / 'A22.csv', index=False)
    totals = annual.sum().astype(int).to_dict()
    metadata = {
        'analysis': 'Post hoc Chile I70–I79 decomposition; not prespecified',
        'period': '2015–2025', 'ages': AGES,
        'unit': 'Discharge episodes; repeat admissions retained',
        'venous_codes': ['I26', 'I80', 'I81', 'I82'],
        'component_codes': {k: sorted(v) for k, v in COMPONENTS.items()},
        'rate_estimator': 'Sum discharges / sum annual INE midyear populations; per 100000 population-years',
        'rate_ratio': 'Pooled female-to-male rate ratio; not year-adjusted; no RR CIs assigned',
        'odds_estimator': 'Female-to-male case–case OR; age-specific log-Wald 95% CI; age-adjusted MH with RBG variance',
        'zero_cell_policy': 'Fail if an OR cell is zero; no continuity correction',
        'totals': totals,
        'inputs': ['data/processed/chile_detailed_diagnoses.csv',
                   'data/processed/denominadores_INE_2015_2025.csv',
                   'data/processed/recuentos_por_ano_sexo_edad_diagnostico.csv'],
        'limitations': 'No incidence, causal inference, patient linkage, or correction for correlated readmissions; nominal OR intervals.',
    }
    (out / 'chile_other_arterial_posthoc_metadata.json').write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    pooled = odds[(odds.arterial_comparator == 'Other_arterial')
                  & (odds.age_or_estimator == 'Age_adjusted_MH')].iloc[0]
    print(f'Post hoc I70–I79: n={totals["Other_arterial"]:,}; '
          f'MH OR={pooled.OR:.8f} ({pooled.CI95_lower:.8f}–{pooled.CI95_upper:.8f})')
    print('Generated results/tables/A21.csv and A22.csv')


if __name__ == '__main__':
    main()
