from pathlib import Path
import pandas as pd, numpy as np, json, hashlib, sys
from scipy.stats import chi2, norm, t
from scipy.optimize import brentq
from scipy.special import expit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
src=ROOT/'data/processed/denominadores_INE_2015_2025.csv'
csrc=ROOT/'data/processed/recuentos_por_ano_sexo_edad_diagnostico.csv'
pop=pd.read_csv(src)
assert len(pop)==154 and not pop.duplicated(['year','sex','age']).any() and (pop.population>0).all()
pop.to_csv(OUT/'denominadores_INE_2015_2025.csv',index=False)
c=pd.read_csv(csrc);c=c[c.country=='Chile'].copy()
groups={'Venous_with_I83':['PE','I80_I82','Varices_I83'],'Venous_without_I83':['PE','I80_I82'],'PE':['PE'],'I80_I82':['I80_I82'],'Varices_I83':['Varices_I83'],'Arterial':['Coronary','Cerebrovascular','Other_arterial'],'Coronary':['Coronary'],'Cerebrovascular':['Cerebrovascular'],'Other_arterial':['Other_arterial']}
assert set(c.diagnostic_group)==set(sum(groups.values(),[])),set(c.diagnostic_group)
records=[]
for label,diagnoses in groups.items():
 d=c[c.diagnostic_group.isin(diagnoses)].groupby(['year','sex','age'],as_index=False)['count'].sum().merge(pop,validate='one_to_one',on=['year','sex','age'])
 assert len(d)==154
 d['outcome']=label;records.append(d)
d=pd.concat(records,ignore_index=True)
def add_rates(x):
 x=x.copy();n=x['count'];den=x.population
 x['rate_100k']=n/den*1e5
 x['rate_low']=np.where(n>0,chi2.ppf(.025,2*n)/2/den*1e5,0)
 x['rate_high']=chi2.ppf(.975,2*(n+1))/2/den*1e5
 return x
annual=add_rates(d);annual.to_csv(OUT/'tasas_anuales_por_edad_sexo.csv',index=False)
pooled=add_rates(d.groupby(['outcome','sex','age'],as_index=False)[['count','population']].sum())
pooled.to_csv(OUT/'tasas_periodo_por_edad_sexo.csv',index=False)
# Poisson log-linear, separate outcome/age, offset log(population), sex + categorical year.
# Main inference uses Pearson dispersion floored at 1 and t(df); df=10 or 9.
models=[]
for (outcome,age,period),g in pd.concat([d.assign(period='2015–2025'),d[d.year<=2024].assign(period='2015–2024')]).groupby(['outcome','age','period']):
 g=g.sort_values(['year','sex']);y=g['count'].to_numpy(float);exposure=g.population.to_numpy(float)
 X=np.column_stack([np.ones(len(g)),(g.sex=='Women').astype(float),*[(g.year==yr).astype(float) for yr in sorted(g.year.unique())[1:]]])
 b=np.zeros(X.shape[1]);b[0]=np.log(y.sum()/exposure.sum())
 for iteration in range(100):
  mu=exposure*np.exp(X@b);H=X.T@(mu[:,None]*X);step=np.linalg.solve(H,X.T@(y-mu));b+=step
  if np.max(np.abs(step))<1e-10:break
 else:raise RuntimeError('No convergence')
 mu=exposure*np.exp(X@b);H=X.T@(mu[:,None]*X);cov=np.linalg.inv(H)
 assert np.max(np.abs(X.T@(y-mu)))<1e-6
 se=np.sqrt(cov[1,1]);df=len(y)-X.shape[1];phi=np.sum((y-mu)**2/mu)/df
 qse=se*np.sqrt(max(1,phi));crit=t.ppf(.975,df)
 F=g.sex=='Women';raw=(y[F].sum()/exposure[F].sum())/(y[~F].sum()/exposure[~F].sum())
 # Independent conditional-binomial score yields the same sex coefficient.
 pairs=g.pivot(index='year',columns='sex',values=['count','population'])
 nf=pairs['count']['Women'].to_numpy();nt=(pairs['count']['Women']+pairs['count']['Men']).to_numpy();off=np.log(pairs['population']['Women']/pairs['population']['Men']).to_numpy()
 independent=brentq(lambda z:np.sum(nf-nt*expit(off+z)),-10,10)
 assert abs(independent-b[1])<1e-8
 models.append(dict(outcome=outcome,age=age,period=period,RR_year_adjusted=np.exp(b[1]),low=np.exp(b[1]-1.959964*se),high=np.exp(b[1]+1.959964*se),RR_pooled=raw,pearson_dispersion=phi,quasi_low=np.exp(b[1]-crit*qse),quasi_high=np.exp(b[1]+crit*qse),residual_df=df))
allmodels=pd.DataFrame(models);allmodels.to_csv(OUT/'razones_tasas_ajustadas_ano.csv',index=False)
m=allmodels[allmodels.period=='2015–2025'].copy()
# Fixed internal standard: national adult population, both sexes, 2015 midyear, seven age bands.
w=pop[pop.year==2015].groupby('age').population.sum();w=w/w.sum();w.rename('weight').to_csv(OUT/'pesos_estandar_Chile_2015.csv')
asr=[]
for (outcome,sex,year),g in d.groupby(['outcome','sex','year']):
 g=g.set_index('age').loc[w.index];r=np.sum(w*g['count']/g.population)*1e5;v=np.sum(w*w*g['count']/g.population**2)*1e10
 asr.append(dict(outcome=outcome,sex=sex,year=year,age_standardized_rate=r,low=max(0,r-1.959964*np.sqrt(v)),high=r+1.959964*np.sqrt(v)))
a=pd.DataFrame(asr);a.to_csv(OUT/'tasas_anuales_estandarizadas_edad.csv',index=False)
audit={'population_source':src.name,'population_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'counts_sha256':hashlib.sha256(csrc.read_bytes()).hexdigest(),'reference_date':'30 June each year','years':list(range(2015,2026)),'population_cells':len(pop),'model_count':len(allmodels),'method':'Poisson sex + categorical year, separate model per age and outcome; log population offset; quasi-Poisson Pearson dispersion floored at 1 with t CI; unscaled Poisson CI retained for audit; exclude-2025 sensitivity','standard':'Chile adult both-sex population, 30 June 2015, seven age bands','scope':'Discharges, not incident persons. No correction for suppressed records or repeated admissions; population projection treated fixed.'}
(OUT/'auditoria_tasas.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2))
print(m[m.outcome.isin(['Venous_with_I83','Venous_without_I83','Arterial'])].to_string(index=False))
print('Pooled expanded rates');print(pooled[pooled.outcome=='Venous_with_I83'].to_string(index=False))
