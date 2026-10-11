"""Reconstructed supplementary analyses, checked against the manuscript's published precision.
Grouped binomial likelihood is algebraically equivalent to individual logistic likelihood.
No patient identifiers; count multiplicity is retained. Regions 99/Extranjero collapse to Unknown.
"""
from pathlib import Path
import json,math
import numpy as np,pandas as pd
from scipy.special import expit
from scipy.stats import chi2
R=Path(__file__).resolve().parents[1];O=R/'results';D=R/'data/processed'
ages=['20–29','30–39','40–49','50–59','60–69','70–79','≥80']
def mh(a,b,c,d):
 a,b,c,d=[np.asarray(x,float) for x in [a,b,c,d]];n=a+b+c+d
 u=sum(a*d/n);v=sum(b*c/n);est=u/v
 var=sum((a+d)*a*d/n**2)/(2*u*u)+sum(((a+d)*b*c+(b+c)*a*d)/n**2)/(2*u*v)+sum((b+c)*b*c/n**2)/(2*v*v)
 return est,*np.exp(np.log(est)+np.array([-1,1])*1.95996398454*np.sqrt(var))
def crude(a,b,c,d):
 est=a*d/(b*c);se=np.sqrt(1/a+1/b+1/c+1/d)
 return est,*np.exp(np.log(est)+np.array([-1,1])*1.95996398454*se)
# Detailed diagnostic sensitivities and full Chapter IX descriptive composition.
d=pd.read_csv(D/'chile_detailed_diagnoses.csv',dtype={'code':str})
c3=d.code.str[:3];num=pd.to_numeric(c3.str[1:]);v=(c3=='I26')|num.between(80,82);a=num.between(20,25)|num.between(60,79)
restricted_v=(c3=='I26')|(d.code.str.match(r'^I80[1-9]'))|(c3=='I82')
restricted_a=num.between(21,22)|num.between(60,64)|(c3=='I74')
scenarios=[('Main',v,a),('Exclude I69 and I800',v&~d.code.str.startswith('I800'),a&(c3!='I69')),('Restricted acute',restricted_v,restricted_a),('PE vs I21',c3=='I26',c3=='I21'),('Venous vs coronary',v,num.between(20,25)),('Venous vs cerebrovascular',v,num.between(60,69)),('Venous vs other arterial',v,num.between(70,79)),('Include I83',v|(c3=='I83'),a)]
rows=[]
for name,vmask,amask in scenarios:
 vals=[]
 for sex,mask in [('Women',vmask),('Women',amask),('Men',vmask),('Men',amask)]:
  vals.append(d[mask&(d.sex==sex)].groupby('age')['count'].sum().reindex(ages,fill_value=0).to_numpy(float))
 rows.append([name,int(sum(x.sum() for x in vals)),*crude(*[x.sum() for x in vals]),*mh(*vals)])
pd.DataFrame(rows,columns=['contrast','n','crude_OR','crude_lower','crude_upper','MH_OR','MH_lower','MH_upper']).to_csv(O/'A19_chile_diagnostic_sensitivities.csv',index=False)
d['group']=np.select([num.between(20,25),num.between(60,69),num.between(70,79),v,c3=='I83'],['Coronary','Cerebrovascular','Other_arterial','Venous','Varices_I83'],default='Other_circulatory')
t=d.groupby(['sex','group'],as_index=False)['count'].sum();t['denominator']=t.groupby('sex')['count'].transform('sum');t['percent']=100*t['count']/t.denominator;t.to_csv(O/'A17_chapter_IX.csv',index=False)
# Logistic main effects and four separate interaction likelihood-ratio tests.
g=pd.read_csv(D/'chile_logistic_cells.csv',dtype={'region':str,'insurance':str})
g['region']=g.region.where(g.region.isin([f'{i:02d}' for i in range(1,17)]),'Unknown')
keys=['year','sex','age','region','insurance']
g=g.groupby(keys+['venous'])['count'].sum().unstack('venous',fill_value=0).reset_index()
y=g[1].to_numpy(float);n=(g[0]+g[1]).to_numpy(float);sex=(g.sex=='Women').to_numpy(float)
assert int(n.sum())==834780
parts={k:pd.get_dummies(g[k].astype(str),drop_first=True,dtype=float).to_numpy() for k in ['age','year','region','insurance']}
X=np.column_stack([np.ones(len(g)),sex,*parts.values()])
def fit(X):
 b=np.zeros(X.shape[1]);b[0]=np.log(y.sum()/(n.sum()-y.sum()))
 for it in range(100):
  eta=X@b;p=expit(eta);H=X.T@((n*p*(1-p))[:,None]*X);score=X.T@(y-n*p)
  step=np.linalg.solve(H,score);b+=step
  if np.max(np.abs(step))<1e-9:break
 else:raise RuntimeError('Logistic model did not converge')
 eta=X@b;p=expit(eta);H=X.T@((n*p*(1-p))[:,None]*X)
 ll=float(np.sum(y*eta-n*np.logaddexp(0,eta)))
 assert np.max(np.abs(X.T@(y-n*p)))<1e-5
 return b,np.linalg.inv(H),ll
b,cov,ll=fit(X);se=np.sqrt(cov[1,1])
main={'OR':float(np.exp(b[1])),'lower':float(np.exp(b[1]-1.95996398454*se)),'upper':float(np.exp(b[1]+1.95996398454*se)),'n':int(n.sum()),'parameters':X.shape[1],'loglik_without_binomial_constant':ll,'region_levels':sorted(g.region.unique()),'insurance_levels':sorted(g.insurance.unique()),'region_rule':'01 through 16; 99 and Extranjero -> Unknown'}
(O/'logistic_main.json').write_text(json.dumps(main,indent=2))
rows=[]
for name,part in parts.items():
 xx=np.column_stack([X,sex[:,None]*part]);_,_,ll1=fit(xx);lr=2*(ll1-ll);df=part.shape[1]
 rows.append([name,lr,df,chi2.sf(lr,df)])
pd.DataFrame(rows,columns=['interaction','LR','df','p']).to_csv(O/'A18_logistic_interactions.csv',index=False)
# Extreme marginal bounds: each age considered independently.
m=pd.read_csv(D/'chile_masked_counts.csv');M=int(m[m.code.isin(['I26','I80','I81','I82'])]['count'].sum());assert M==3325
rate=pd.read_csv(O/'tasas_periodo_por_edad_sexo.csv');rate=rate[rate.outcome=='Venous_without_I83'];rows=[]
for age in ages:
 z=rate[rate.age==age].set_index('sex');w=z.loc['Women'];m=z.loc['Men'];rr=(w['count']/w.population)/(m['count']/m.population)
 lo=(w['count']/w.population)/((m['count']+M)/m.population);hi=((w['count']+M)/w.population)/(m['count']/m.population)
 needed=math.ceil(w['count']*m.population/w.population-m['count']) if rr>1 else math.ceil(m['count']*w.population/m.population-w['count'])
 rows.append([age,rr,lo,hi,needed,'Men' if rr>1 else 'Women',M])
pd.DataFrame(rows,columns=['age','pooled_RR','extreme_lower','extreme_upper','added_to_cross_1','assigned_sex','available_masked']).to_csv(O/'A13_masked_bounds.csv',index=False)
print('Logistic main',main)
print(pd.DataFrame(rows).to_string(index=False))
