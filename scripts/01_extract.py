"""Extract and audit original exports. Rate models are fitted separately in 03_rates.py.
Run: python recalcular_composicion.py --inputs UPLOAD_FOLDER --output OUTPUT_FOLDER
Requires pandas and numpy. Does not modify the manuscript.
"""
import argparse, csv, hashlib, io, json, re
from collections import Counter
from pathlib import Path
from zipfile import ZipFile
import numpy as np
import pandas as pd

ap=argparse.ArgumentParser();ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--reuse',type=Path,help='Previously verified Chile aggregates and manifest; optional');args=ap.parse_args()
args.output.mkdir(parents=True,exist_ok=True)
ages=['20–29','30–39','40–49','50–59','60–69','70–79','≥80']
groups=['PE','I80_I82','Varices_I83','Coronary','Cerebrovascular','Other_arterial']
mapping={
 'PE':['Embolia pulmonar'],
 'I80_I82':['Flebite tromboflebite embolia e trombose venosa'],
 'Varices_I83':['Veias varicosas das extremidades inferiores'],
 'Coronary':['Infarto agudo do miocárdio','Outras doenças isquêmicas do coração'],
 'Cerebrovascular':['Hemorragia intracraniana','Infarto cerebral','Acid vascular cerebr não espec hemorrág ou isq','Outras doenças cerebrovasculares'],
 'Other_arterial':['Arteroesclerose','Outras doenças vasculares periféricas','Embolia e trombose arteriais','Outras doenças das artérias arteríolas e capil']}
out=[];manifest=[];audit=[];rawcounts=Counter(); detail=Counter(); covariates=Counter(); masked=Counter()

def age_map(s):
 if re.search(r'\b(D[IÍ]AS?|MESES?|MES)\b',s.upper()):return None
 m=re.match(r'^(\d+)\s',s.strip())
 if not m:return None
 n=int(m[1])
 if n<20:return None
 return ages[min(n//10-2,6)]

reuse_years=set()
if args.reuse:
 cached=pd.read_csv(args.reuse/'recuentos_por_ano_sexo_edad_diagnostico.csv')
 cached=cached[cached.country=='Chile']
 ca=json.loads((args.reuse/'auditoria_fuentes.json').read_text())
 audit=ca['chile_flow'];manifest=[x for x in ca['files'] if x['country']=='Chile']
 reuse_years=set(cached.year.unique());assert reuse_years=={x['year'] for x in audit}
 for r in cached.itertuples(index=False):rawcounts[(r.year,r.sex,r.age,r.diagnostic_group)]=int(r.count)
 print('Reusing previously verified Chile years:',sorted(reuse_years),flush=True)

def process_chile(name,opener,size):
  year=int(re.search(r'20\d\d',name)[0]);metrics=Counter();age_labels=Counter();sex_labels=Counter();digest=hashlib.sha256()
  if year in reuse_years:return
  # Hash the uncompressed original, separately from parsing for reproducibility.
  with opener() as f:
   for block in iter(lambda:f.read(1024*1024),b''):digest.update(block)
  with opener() as f:
   for df in pd.read_csv(f,sep=';',encoding='latin1',dtype=str,usecols=['SEXO','GRUPO_EDAD','ANO_EGRESO','DIAG1','REGION_RESIDENCIA','PREVISION'],keep_default_na=False,chunksize=200000):
    metrics['all_source_rows']+=len(df)
    norm=df.DIAG1.str.strip().str.upper().str.replace('.','',regex=False)
    chapter=norm.str.match(r'^I\d{2}')
    d=df.loc[chapter].copy();d['code']=norm.loc[chapter].str[:3]
    metrics['chapter_IX']+=len(d);age_labels.update(d.GRUPO_EDAD);sex_labels.update(d.SEXO)
    joint=(d.SEXO=='*')&(d.GRUPO_EDAD=='*')&(d.ANO_EGRESO=='*')
    metrics['jointly_masked']+=int(joint.sum())
    for code,cnt in d.loc[joint,'code'].value_counts().items():masked[(year,code)]+=int(cnt)
    d['age']=d.GRUPO_EDAD.map(age_map)
    d['sex']=d.SEXO.str.strip().str.upper().map({'MUJER':'Women','HOMBRE':'Men','2':'Women','1':'Men'})
    metrics['other_ineligible_age']+=int((~joint & d.age.isna()).sum())
    metrics['adult_unknown_sex']+=int((d.age.notna() & d.sex.isna()).sum())
    eligible=d.age.notna() & d.sex.notna()
    assert (d.loc[eligible,'ANO_EGRESO']==str(year)).all()
    d=d.loc[eligible].copy();metrics['eligible_chapter_IX']+=len(d)
    n=pd.to_numeric(d.code.str[1:],errors='raise')
    d['group']=np.select([d.code.eq('I26'),n.between(80,82),d.code.eq('I83'),n.between(20,25),n.between(60,69),n.between(70,79)],groups,default='Outside_selected')
    for (sex,age,group),cnt in d.groupby(['sex','age','group']).size().items():rawcounts[(year,sex,age,group)]+=int(cnt)
    for (sex,age,code),cnt in d.assign(full_code=norm.loc[d.index]).groupby(['sex','age','full_code']).size().items():detail[(year,sex,age,code)]+=int(cnt)
    selected=d[~d.group.isin(['Outside_selected','Varices_I83'])].copy()
    selected['venous']=selected.group.isin(['PE','I80_I82']).astype(int)
    for (sex,age,reg,ins,v),cnt in selected.groupby(['sex','age','REGION_RESIDENCIA','PREVISION','venous'],dropna=False).size().items():covariates[(year,sex,age,reg,ins,v)]+=int(cnt)
    metrics['selected_with_varices']+=int(d.group.ne('Outside_selected').sum())
    metrics['varices_I83']+=int(d.group.eq('Varices_I83').sum())
    metrics['selected_without_varices']+=int((~d.group.isin(['Outside_selected','Varices_I83'])).sum())
  assert metrics['chapter_IX']==metrics['jointly_masked']+metrics['other_ineligible_age']+metrics['adult_unknown_sex']+metrics['eligible_chapter_IX']
  audit.append({'year':year,**metrics,'age_labels':dict(age_labels),'sex_labels':dict(sex_labels)})
  manifest.append({'country':'Chile','file':name,'sha256':digest.hexdigest(),'source_bytes':size})
  print('Chile',year,dict(metrics),flush=True)
with ZipFile(args.inputs/'CHILE(1).zip') as z:
 for name in sorted(n for n in z.namelist() if n.endswith('.csv') and not n.startswith('__MACOSX')):
  process_chile(name,lambda:z.open(name),z.getinfo(name).file_size)
for year in [2024,2025]:
 p=args.inputs/f'EGR_DATOS_ABIERTO_{year}(2).csv'
 if p.exists():process_chile(p.name,lambda:p.open('rb'),p.stat().st_size)
assert len({x['year'] for x in audit})==len(audit)
for year in sorted({k[0] for k in rawcounts}):
 for sex in ['Women','Men']:
  for age in ages:
   for group in groups:out.append(['Chile',year,sex,age,group,rawcounts[(year,sex,age,group)]])

bfiles=[(p.name,p.read_bytes()) for p in args.inputs.glob('*.csv') if p.name.startswith(('FEM_SOLO_2015','Masc_SOLO2015'))]
with ZipFile(args.inputs/'BRASIL TODO(1).zip') as z:
 bfiles +=[(n,z.read(n)) for n in z.namelist() if n.startswith('BRASIL TODO/') and n.endswith('.csv')]
seen=set()
for name,bb in sorted(bfiles):
 s=bb.decode('latin1');year=int(re.search(r'Período:(\d{4})',s)[1]);sex={'Fem':'Women','Masc':'Men'}[re.search(r'Sexo: (\w+)',s)[1]]
 assert (year,sex) not in seen;seen.add((year,sex))
 ls=s.splitlines();ix=next(i for i,l in enumerate(ls) if l.startswith('"Lista Morb'))
 rows=list(csv.reader(ls[ix:],delimiter=';'));tab={r[0].removeprefix('.. '):[int(v) for v in r[1:]] for r in rows[1:] if len(r)==9}
 assert all(sum(v[:7])==v[7] for v in tab.values())
 assert np.array_equal(np.sum([v for k,v in tab.items() if k not in ['Total','09 Doenças do aparelho circulatório']],axis=0),tab['Total'])
 assert rows[0][1:-1]==['20 a 29 anos','30 a 39 anos','40 a 49 anos','50 a 59 anos','60 a 69 anos','70 a 79 anos','80 anos e mais']
 for group,labels in mapping.items():
  counts=np.sum([tab[l][:7] for l in labels],axis=0)
  for age,c in zip(ages,counts):out.append(['Brazil',year,sex,age,group,int(c)])
 manifest.append({'country':'Brazil','file':name,'sha256':hashlib.sha256(bb).hexdigest(),'source_bytes':len(bb),'declared_download_date':'2026-09-20','header':ls[:ix]})
assert seen=={(y,s) for y in range(2015,2026) for s in ['Women','Men']}
df=pd.DataFrame(out,columns=['country','year','sex','age','diagnostic_group','count']);df.to_csv(args.output/'recuentos_por_ano_sexo_edad_diagnostico.csv',index=False)

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
pd.DataFrame(est,columns=['country','period','venous_definition','arterial_comparator','age_or_estimator','OR','CI95_lower','CI95_upper']).to_csv(args.output/'odds_ratios_con_y_sin_varices.csv',index=False)
pd.DataFrame(composition,columns=['country','period','venous_definition','age','sex','venous_count','selected_total','venous_percent',*groups]).to_csv(args.output/'composicion_por_edad_sexo.csv',index=False)
(args.output/'auditoria_fuentes.json').write_text(json.dumps({'files':manifest,'chile_flow':audit,'brazil_group_mapping':mapping,'coverage':{ct:[int(df[df.country==ct].year.min()),int(df[df.country==ct].year.max())] for ct in ['Chile','Brazil']},'limitations':['This extraction stage does not fit rates; see 00_population.py and 03_rates.py.','Brazil download date supplied by author, not encoded in CSV.','IC95% nominal, based on discharge counts; does not account for repeat admissions or coding/coverage uncertainty.']},ensure_ascii=False,indent=2))
print('COUNTS',df.groupby(['country','diagnostic_group'])['count'].sum().to_string(),flush=True)
print('MH',pd.DataFrame(est,columns=['country','period','definition','comparator','stratum','OR','lower','upper']).query("stratum=='Age_adjusted_MH'").to_string(index=False),flush=True)

if not args.reuse:
 pd.DataFrame([(*k,n) for k,n in detail.items()],columns=['year','sex','age','code','count']).to_csv(args.output/'chile_detailed_diagnoses.csv',index=False)
 pd.DataFrame([(*k,n) for k,n in covariates.items()],columns=['year','sex','age','region','insurance','venous','count']).to_csv(args.output/'chile_logistic_cells.csv',index=False)
 pd.DataFrame([(*k,n) for k,n in masked.items()],columns=['year','code','count']).to_csv(args.output/'chile_masked_counts.csv',index=False)
