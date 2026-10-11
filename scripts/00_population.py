"""Export Chile midyear population denominators from the INE original workbook."""
from pathlib import Path
import sys,json,hashlib
import pandas as pd
R=Path(__file__).resolve().parents[1]
src=Path(sys.argv[1]);p=pd.read_excel(src);p.columns=p.columns.str.strip()
p['date']=pd.to_datetime(p.FECHA,dayfirst=True)
p=p[(p.date.dt.month==6)&(p.date.dt.day==30)&p.date.dt.year.between(2015,2025)].copy()
assert len(p)==2222 and not p.duplicated(['date','SEXO','EDAD']).any()
assert set(p.SEXO)=={'H','M'} and (p.POBLACION>0).all()
p['year']=p.date.dt.year;p['sex']=p.SEXO.map({'H':'Men','M':'Women'})
p=p[p.EDAD>=20].copy();p['age']=p.EDAD.map(lambda a:f'{a//10*10}–{a//10*10+9}' if a<80 else '≥80')
out=p.groupby(['year','sex','age'],as_index=False).POBLACION.sum().rename(columns={'POBLACION':'population'})
assert len(out)==154
out.to_csv(R/'data/processed/denominadores_INE_2015_2025.csv',index=False)
(R/'data/processed/population_source.json').write_text(json.dumps({'file':src.name,'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'reference_date':'30 June','base':2024},indent=2))
