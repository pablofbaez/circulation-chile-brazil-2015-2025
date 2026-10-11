from pathlib import Path
import numpy as np,pandas as pd,json
R=Path(__file__).resolve().parents[1];OUT=R/'results';OUT.mkdir(exist_ok=True)
df=pd.read_csv(R/'data/processed/recuentos_por_ano_sexo_edad_diagnostico.csv')
ages=['20–29','30–39','40–49','50–59','60–69','70–79','≥80'];groups=['PE','I80_I82','Varices_I83','Coronary','Cerebrovascular','Other_arterial']
assert len(df)==1848 and not df.duplicated(['country','year','sex','age','diagnostic_group']).any()
assert set(df.year)==set(range(2015,2026)) and (df['count']>=0).all()
def mh(vw,aw,vm,am):
 vw,aw,vm,am=[np.array(x,dtype=float) for x in [vw,aw,vm,am]];n=vw+aw+vm+am
 R=np.sum(vw*am/n);S=np.sum(aw*vm/n);o=R/S
 var=np.sum((vw+am)*vw*am/n**2)/(2*R**2)+np.sum(((vw+am)*aw*vm+(aw+vm)*vw*am)/n**2)/(2*R*S)+np.sum((aw+vm)*aw*vm/n**2)/(2*S**2)
 lo,hi=np.exp(np.log(o)+np.array([-1,1])*1.95996398454*np.sqrt(var));return float(o),float(lo),float(hi)
est=[];composition=[]
for country in ['Chile','Brazil']:
 sub=df[df.country==country];period=f'{sub.year.min()}–{sub.year.max()}'
 pivot=sub.groupby(['age','sex','diagnostic_group'])['count'].sum().unstack('diagnostic_group').reindex(columns=groups)
 for definition,vg in [('with_varices',['PE','I80_I82','Varices_I83']),('without_varices',['PE','I80_I82'])]:
  for age in ages:
   for sex in ['Women','Men']:
    row=pivot.loc[(age,sex)];total=int(row[vg+groups[3:]].sum());ven=int(row[vg].sum())
    composition.append([country,period,definition,age,sex,ven,total,100*ven/total,*[int(row[g]) for g in groups]])
  for label,ag in [('Pooled_arterial',groups[3:])]+[(g,[g]) for g in groups[3:]]:
   vw=np.array([pivot.loc[(a,'Women'),vg].sum() for a in ages],float);vm=np.array([pivot.loc[(a,'Men'),vg].sum() for a in ages],float)
   aw=np.array([pivot.loc[(a,'Women'),ag].sum() for a in ages],float);am=np.array([pivot.loc[(a,'Men'),ag].sum() for a in ages],float)
   est.append([country,period,definition,label,'Age_adjusted_MH',*mh(vw,aw,vm,am)])
   odds=vw*am/(aw*vm);se=np.sqrt(1/vw+1/am+1/aw+1/vm)
   for age,o,e in zip(ages,odds,se):est.append([country,period,definition,label,age,float(o),float(np.exp(np.log(o)-1.95996398454*e)),float(np.exp(np.log(o)+1.95996398454*e))])
pd.DataFrame(est,columns=['country','period','venous_definition','arterial_comparator','age_or_estimator','OR','CI95_lower','CI95_upper']).to_csv(OUT/'odds_ratios_con_y_sin_varices.csv',index=False)
pd.DataFrame(composition,columns=['country','period','venous_definition','age','sex','venous_count','selected_total','venous_percent',*groups]).to_csv(OUT/'composicion_por_edad_sexo.csv',index=False)
annual=[]
for (country,year),s in df.groupby(['country','year']):
 p=s.pivot(index='age',columns=['sex','diagnostic_group'],values='count').reindex(ages)
 for definition,vg in [('with_varices',groups[:3]),('without_varices',groups[:2])]:
  for label,ag in [('Pooled_arterial',groups[3:])]+[(g,[g]) for g in groups[3:]]:
   vals=[p.loc[:,[(sex,g) for g in gs]].sum(axis=1).values for sex,gs in [('Women',vg),('Women',ag),('Men',vg),('Men',ag)]]
   annual.append([country,year,definition,label,*mh(*vals)])
pd.DataFrame(annual,columns=['country','year','venous_definition','arterial_comparator','OR_age_adjusted_MH','CI95_lower','CI95_upper']).to_csv(OUT/'OR_anuales_ajustados_por_edad.csv',index=False)


# Direct standardization within country, same weights for women and men.
c=pd.read_csv(OUT/'composicion_por_edad_sexo.csv');c=c[c.venous_definition=='without_varices']
rows=[]
for country,g in c.groupby('country'):
 w=g[g.sex=='Women'].set_index('age');m=g[g.sex=='Men'].set_index('age')
 weights=(w.selected_total+m.selected_total)/(w.selected_total.sum()+m.selected_total.sum())
 for sex,x in [('Women',w),('Men',m)]:
  rows.append([country,sex,100*x.venous_count.sum()/x.selected_total.sum(),sum(weights*x.venous_percent)])
pd.DataFrame(rows,columns=['country','sex','crude_percent','age_standardized_percent']).to_csv(OUT/'standardized_venous_fractions.csv',index=False)
# Brazil sensitivity excluding 2023, using identical MH estimator.
z=df[(df.country=='Brazil')&(df.year!=2023)].groupby(['age','sex','diagnostic_group'])['count'].sum().unstack('diagnostic_group')
values=[np.array([z.loc[(age,sex),gs].sum() for age in ages]) for sex,gs in [('Women',groups[:2]),('Women',groups[3:]),('Men',groups[:2]),('Men',groups[3:])]]
pd.DataFrame([['Brazil',2023,*mh(*values)]],columns=['country','excluded_year','OR','lower','upper']).to_csv(OUT/'brazil_exclude_2023.csv',index=False)
