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
(O/'validation.json').write_text(json.dumps(reports,indent=2,ensure_ascii=False));print(json.dumps(reports,indent=2,ensure_ascii=False))
if any(x['status']!='PASS' for x in reports):print('Secondary discrepancies remain; see validation.json.')
