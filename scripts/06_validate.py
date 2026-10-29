"""Compare outputs with retained numeric estimates and manuscript tables, not with rounded graphs."""
from pathlib import Path
import json
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[1];O=R/'results';reports=[]
for file in ['odds_ratios_con_y_sin_varices.csv','composicion_por_edad_sexo.csv','OR_anuales_ajustados_por_edad.csv','razones_tasas_ajustadas_ano.csv']:
 a=pd.read_csv(O/file);b=pd.read_csv(R/'data/reference'/file)
 assert list(a.columns)==list(b.columns) and a.shape==b.shape,file
 for c in a:
  if pd.api.types.is_numeric_dtype(a[c]):assert np.allclose(a[c],b[c],rtol=1e-9,atol=1e-10),f'{file}:{c}'
  else:assert a[c].equals(b[c]),f'{file}:{c}'
 reports.append({'check':file,'status':'PASS','rows':len(a)})
# Compare regenerated A19 and A18 to the actual manuscript at displayed precision.
reference=json.loads((R/'docs/manuscript_tables_reference.json').read_text())['tables']
a=pd.read_csv(O/'A19_chile_diagnostic_sensitivities.csv')
for i,row in a.iterrows():
 old=reference[21][i+1]
 assert int(old[3].replace(',',''))==row['n'],('A19 n',i,old,row.to_dict())
 import re
 numbers=[float(x) for x in re.findall(r'\d+\.\d+',old[2])]
 assert [round(row[k],2) for k in ['MH_OR','MH_lower','MH_upper']]==numbers,('A19 MH',i,numbers,row.to_dict())
reports.append({'check':'A19 Chile eight diagnostic contrasts','status':'PASS'})
a=pd.read_csv(O/'A18_logistic_interactions.csv');mismatches=[]
for i,row in a.iterrows():
 old=reference[20][i+1]
 if round(row.LR,2)!=float(old[1]) or row.df!=int(old[2]):mismatches.append({'interaction':row.interaction,'manuscript_LR':old[1],'computed_LR':row.LR,'manuscript_df':old[2],'computed_df':row.df})
reports.append({'check':'A18 logistic interactions','status':'PASS' if not mismatches else 'REQUIRES_RECONCILIATION','differences':mismatches})
m=json.loads((O/'logistic_main.json').read_text());reports.append({'check':'Additional adjusted OR rounds to 1.82 (1.79–1.84)','status':'PASS' if [round(m[k],2) for k in ['OR','lower','upper']]==[1.82,1.79,1.84] else 'REQUIRES_RECONCILIATION','computed':m})
counts=pd.read_csv(R/'data/processed/recuentos_por_ano_sexo_edad_diagnostico.csv');main=counts[counts.diagnostic_group!='Varices_I83']
assert main.groupby('country')['count'].sum().to_dict()=={'Chile':834780,'Brazil':7409594}
# Independent algebraic check: pooled RR_venous / RR_arterial = case-case OR.
r=pd.read_csv(O/'tasas_periodo_por_edad_sexo.csv');o=pd.read_csv(O/'odds_ratios_con_y_sin_varices.csv')
for age in r.age.unique():
 def rr(k):
  z=r[(r.age==age)&(r.outcome==k)].set_index('sex');return z.loc['Women','rate_100k']/z.loc['Men','rate_100k']
 for c,k in [('Pooled_arterial','Arterial'),('Coronary','Coronary'),('Cerebrovascular','Cerebrovascular'),('Other_arterial','Other_arterial')]:
  val=o[(o.country=='Chile')&(o.venous_definition=='without_varices')&(o.arterial_comparator==c)&(o.age_or_estimator==age)].OR.iloc[0]
  assert abs(rr('Venous_without_I83')/rr(k)-val)<1e-10
reports.append({'check':'28 pooled rate-ratio identities','status':'PASS'})
# Post hoc A21/A22: reconcile with independently supplied main counts/results and
# author-supplied manuscript tables, rather than comparing outputs to themselves.
post_rates = pd.read_csv(O/'chile_other_arterial_posthoc_rates.csv')
post_odds = pd.read_csv(O/'chile_other_arterial_posthoc_odds.csv')
post_counts = pd.read_csv(O/'chile_other_arterial_posthoc_counts.csv')
post_ref = json.loads((R/'data/reference/chile_other_arterial_posthoc_manuscript.json').read_text())
assert int(post_counts.Other_arterial.sum()) == 87899
assert int(post_counts.Venous.sum()) == 109820
component_totals = post_rates.groupby('diagnostic_group')[['women_count','men_count']].sum()
expected_totals = {'I70_I73': [6207,11397], 'I71': [5160,12276],
                   'I72': [3359,3853], 'I74': [6853,9906], 'I77_I79': [11642,17246]}
for name, total in expected_totals.items():
    assert component_totals.loc[name].tolist() == total, ('Post hoc component totals', name)
reports.append({'check':'Post hoc annual cell reconciliation and component totals', 'status':'PASS',
                'other_arterial_n':87899,'venous_n':109820})
base = pd.read_csv(O/'odds_ratios_con_y_sin_varices.csv')
base = base[(base.country=='Chile') & (base.venous_definition=='without_varices')]
for name in ['Other_arterial', 'Coronary']:
    left = post_odds[post_odds.arterial_comparator==name].set_index('age_or_estimator')
    right = base[base.arterial_comparator==name].set_index('age_or_estimator').reindex(left.index)
    assert np.allclose(left[['OR','CI95_lower','CI95_upper']],
                       right[['OR','CI95_lower','CI95_upper']],rtol=1e-12,atol=1e-12), name
reports.append({'check':'Post hoc reproduces all A3 coronary/other-arterial ORs and CIs plus MH summaries',
                'status':'PASS','estimates':16})
# Check all displayed A21 ratios and all A22 age-specific/MH ORs and intervals.
for row in post_ref['A21']:
    for name, expected in row.items():
        if name=='age': continue
        estimate = post_rates[(post_rates.age==row['age']) &
                              (post_rates.diagnostic_group==name)].pooled_RR.iloc[0]
        assert round(estimate,2)==expected, ('A21', row['age'],name,estimate,expected)
for row in post_ref['A22']:
    estimate = post_odds[(post_odds.age_or_estimator==row['age_or_estimator']) &
                        (post_odds.arterial_comparator==row['arterial_comparator'])].iloc[0]
    for name in ['OR','CI95_lower','CI95_upper']:
        assert round(estimate[name],2)==row[name], ('A22',row,name,estimate[name])
reports.append({'check':'A21 56 displayed RRs and A22 32 displayed OR/CI triplets', 'status':'PASS'})
for age in post_rates.age.unique():
    values = post_rates[post_rates.age==age].set_index('diagnostic_group').pooled_RR
    for name in ['I71','Other_arterial_excluding_I71','Other_arterial','Coronary']:
        estimate = post_odds[(post_odds.age_or_estimator==age) &
                             (post_odds.arterial_comparator==name)].OR.iloc[0]
        assert np.isclose(values['Venous']/values[name],estimate,rtol=1e-12,atol=1e-12)
reports.append({'check':'Post hoc 28 pooled RR_venous/RR_comparator identities', 'status':'PASS'})

(O/'validation.json').write_text(json.dumps(reports,indent=2,ensure_ascii=False));print(json.dumps(reports,indent=2,ensure_ascii=False))
if any(x['status']!='PASS' for x in reports):print('Secondary discrepancies remain; see validation.json.')
