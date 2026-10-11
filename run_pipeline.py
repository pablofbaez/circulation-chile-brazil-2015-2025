"""Run from supplied aggregates, or rebuild them from originals with --raw-dir and --ine."""
from pathlib import Path
import argparse,subprocess,sys
R=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--raw-dir',type=Path);p.add_argument('--ine',type=Path);a=p.parse_args()
def run(name,*args):
 print('Running',name,flush=True)
 subprocess.run([sys.executable,str(R/'scripts'/name),*map(str,args)],check=True)
if a.raw_dir:run('01_extract.py','--inputs',a.raw_dir.resolve(),'--output',R/'data/processed')
if a.ine:run('00_population.py',a.ine.resolve())
for name in ['02_composition.py','03_rates.py','04_secondary.py','05_figures.py','07_tables.py','06_validate.py']:run(name)
