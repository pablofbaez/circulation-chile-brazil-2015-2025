"""Machine-readable table sources. Numeric outputs are regenerated; A12/A14/A15 are curated text."""
from pathlib import Path
import pandas as pd,numpy as np,json
R=Path(__file__).resolve().parents[1];O=R/'results';T=O/'tables';T.mkdir(exist_ok=True)
c=pd.read_csv(O/'composicion_por_edad_sexo.csv');c=c[c.venous_definition=='without_varices']
e=pd.read_csv(O/'odds_ratios_con_y_sin_varices.csv');e=e[e.venous_definition=='without_varices']
n=pd.read_csv(R/'data/processed/recuentos_por_ano_sexo_edad_diagnostico.csv')
r=pd.read_csv(O/'tasas_periodo_por_edad_sexo.csv');rr=pd.read_csv(O/'razones_tasas_ajustadas_ano.csv');full=rr[rr.period=='2015–2025']
rows=[]
for (ct,s),g in c.groupby(['country','sex']):
 total=g.selected_total.sum()
 rows.append([ct,s,'Total',int(total),100.])
 for row in g.itertuples():rows.append([ct,s,row.age,int(row.selected_total),100*row.selected_total/total])
 for field in ['venous_count','PE','I80_I82','Coronary','Cerebrovascular','Other_arterial']:
  val=g[field].sum();rows.append([ct,s,field,int(val),100*val/total])
pd.DataFrame(rows,columns=['country','sex','category','n','percent']).to_csv(T/'Table_1.csv',index=False)
e[e.age_or_estimator=='Age_adjusted_MH'].to_csv(T/'Table_2.csv',index=False)
r[r.outcome.isin(['Venous_without_I83','Coronary'])].to_csv(T/'Table_3_rates.csv',index=False)
full[full.outcome.isin(['Venous_without_I83','Coronary'])][['outcome','age','RR_pooled']].to_csv(T/'Table_3_pooled_RR.csv',index=False)
for ct,a,b in [('Chile','A1','A3'),('Brazil','A2','A4')]:
 q=e[(e.country==ct)&(e.age_or_estimator!='Age_adjusted_MH')]
 z=q[q.arterial_comparator=='Pooled_arterial'].merge(c[c.country==ct],left_on=['country','age_or_estimator'],right_on=['country','age'],suffixes=('_OR','_counts'))
 z.to_csv(T/f'{a}.csv',index=False);q[q.arterial_comparator!='Pooled_arterial'].to_csv(T/f'{b}.csv',index=False)
a=pd.read_csv(O/'OR_anuales_ajustados_por_edad.csv')
for ct,label in [('Chile','A5'),('Brazil','A6')]:
 q=a[(a.country==ct)&(a.venous_definition=='without_varices')&(a.arterial_comparator=='Pooled_arterial')]
 totals=n[(n.country==ct)&(n.diagnostic_group!='Varices_I83')].groupby('year')['count'].sum().rename('selected_n')
 q.merge(totals,on='year').to_csv(T/f'{label}.csv',index=False)
full[full.outcome.isin(['PE','I80_I82','Varices_I83','Venous_with_I83'])].to_csv(T/'A7.csv',index=False)
audit=json.loads((R/'data/processed/auditoria_fuentes.json').read_text())
flow=pd.DataFrame(audit['chile_flow']).drop(columns=['age_labels','sex_labels']);flow['masked_percent']=100*flow.jointly_masked/flow.chapter_IX;flow.sort_values('year').to_csv(T/'A8.csv',index=False)
rr[rr.outcome=='Venous_without_I83'].to_csv(T/'A9.csv',index=False)
x=e[e.age_or_estimator=='Age_adjusted_MH'].copy();x['n']=[int(c[c.country==row.country].venous_count.sum()+ (c[c.country==row.country][['Coronary','Cerebrovascular','Other_arterial']].sum().sum() if row.arterial_comparator=='Pooled_arterial' else c[c.country==row.country][row.arterial_comparator].sum())) for row in x.itertuples()];x.to_csv(T/'A10.csv',index=False)
full[full.outcome.isin(['Venous_without_I83','Arterial'])].to_csv(T/'A11_RR.csv',index=False);r[r.outcome.isin(['Venous_without_I83','Arterial'])].to_csv(T/'A11_rates.csv',index=False)
ref=json.loads((R/'docs/manuscript_tables_reference.json').read_text())['tables']
for ix,label in [(14,'A12'),(16,'A14'),(17,'A15')]:pd.DataFrame(ref[ix][1:],columns=ref[ix][0]).to_csv(T/f'{label}_curated.csv',index=False)
f=flow.sum(numeric_only=True);eligible=int(f.eligible_chapter_IX);v=int(f.varices_I83);selected=int(f.selected_without_varices)
pd.DataFrame([['Chapter IX',int(f.chapter_IX)],['Age missing or ineligible',int(f.jointly_masked+f.other_ineligible_age)],['Other adult sex missing',int(f.adult_unknown_sex)],['Eligible adults',eligible],['Varices',v],['Other outside main contrast',eligible-v-selected],['Main contrast',selected]],columns=['stage','n']).to_csv(T/'A16.csv',index=False)
for filename,label in [('A13_masked_bounds.csv','A13'),('A17_chapter_IX.csv','A17'),('A18_logistic_interactions.csv','A18'),('A19_chile_diagnostic_sensitivities.csv','A19_Chile')]:pd.read_csv(O/filename).to_csv(T/f'{label}.csv',index=False)
# Brazil A19 expanded/main crude and MH results.
allc=pd.read_csv(O/'composicion_por_edad_sexo.csv');alle=pd.read_csv(O/'odds_ratios_con_y_sin_varices.csv');rows=[]
for definition,g in allc[allc.country=='Brazil'].groupby('venous_definition'):
 w=g[g.sex=='Women'];m=g[g.sex=='Men'];a=w.venous_count.sum();b=w.selected_total.sum()-a;c0=m.venous_count.sum();d=m.selected_total.sum()-c0;est=a*d/(b*c0);se=np.sqrt(1/a+1/b+1/c0+1/d)
 mh=alle[(alle.country=='Brazil')&(alle.venous_definition==definition)&(alle.arterial_comparator=='Pooled_arterial')&(alle.age_or_estimator=='Age_adjusted_MH')].iloc[0]
 rows.append([definition,int(g.selected_total.sum()),est,np.exp(np.log(est)-1.95996398454*se),np.exp(np.log(est)+1.95996398454*se),mh.OR,mh.CI95_lower,mh.CI95_upper])
pd.DataFrame(rows,columns=['definition','n','crude_OR','crude_lower','crude_upper','MH_OR','MH_lower','MH_upper']).to_csv(T/'A19_Brazil.csv',index=False)
n.groupby(['country','age','sex','diagnostic_group'])['count'].sum().unstack('diagnostic_group').to_csv(T/'A20.csv')
# Wilson CIs used in Figure S1, exported rather than only displayed.
c=pd.read_csv(O/'composicion_por_edad_sexo.csv');c=c[c.venous_definition=='without_varices'].copy();p=c.venous_count/c.selected_total;z=1.959964;nn=c.selected_total;mid=(p+z*z/(2*nn))/(1+z*z/nn);half=z*np.sqrt(p*(1-p)/nn+z*z/(4*nn*nn))/(1+z*z/nn)
c['lower_percent']=100*(mid-half);c['upper_percent']=100*(mid+half);c.to_csv(O/'S1_proportions_Wilson.csv',index=False)
