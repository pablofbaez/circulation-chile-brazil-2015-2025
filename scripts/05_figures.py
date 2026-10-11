from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];O=R/'results/figures';O.mkdir(parents=True,exist_ok=True)
OR=pd.read_csv(R/'results/odds_ratios_con_y_sin_varices.csv')
C=pd.read_csv(R/'results/composicion_por_edad_sexo.csv');C=C[C.venous_definition=='without_varices']
RR=pd.read_csv(R/'results/razones_tasas_ajustadas_ano.csv')
ages=['20–29','30–39','40–49','50–59','60–69','70–79','≥80'];comps=['Pooled_arterial','Coronary','Cerebrovascular','Other_arterial']
# Figures are generated from numeric results, with confidence intervals.
plt.rcParams.update({'font.size':9})
fig,axs=plt.subplots(2,2,figsize=(7,6.0),sharex=True,sharey=True)
for ax,co,letter in zip(axs.flat,comps,'ABCD'):
 for ct,color,marker in [('Chile','#0072B2','o'),('Brazil','#D55E00','s')]:
  q=OR[(OR.country==ct)&(OR.venous_definition=='without_varices')&(OR.arterial_comparator==co)].set_index('age_or_estimator').loc[ages]
  x=np.arange(7)+(-.04 if ct=='Chile' else .04)
  ax.errorbar(x,q.OR,yerr=[q.OR-q.CI95_lower,q.CI95_upper-q.OR],fmt=marker+'-',color=color,capsize=2,markersize=3,label=ct,linewidth=1)
 ax.axhline(1,color='gray',linestyle='--',linewidth=.8);ax.set_yscale('log');ax.set_ylim(.8,9);ax.set_yticks([1,2,4,8],[1,2,4,8]);ax.grid(axis='y',alpha=.2);ax.spines[['top','right']].set_visible(False)
 ax.text(-.13,1.03,letter,transform=ax.transAxes,fontweight='bold',fontsize=12);ax.set_xticks(range(7),ages,rotation=45,ha='right')
for ax in axs[:,0]:ax.set_ylabel('Case–case OR (women vs men)')
for ax in axs[1]:ax.set_xlabel('Age group (years)')
axs[0,1].legend(frameon=False);fig.tight_layout()
fig.savefig(O/'Figura_1.png',dpi=300);fig.savefig(O/'Figura_1.pdf');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(7,3.7),sharey=True)
for ax,out,color,letter in zip(axs,['Venous_without_I83','Arterial'],['#0072B2','#D55E00'],'AB'):
 q=RR[(RR.period=='2015–2025')&(RR.outcome==out)].set_index('age').loc[ages]
 ax.errorbar(q.RR_year_adjusted,range(7),xerr=[q.RR_year_adjusted-q.quasi_low,q.quasi_high-q.RR_year_adjusted],fmt='o',color=color,capsize=3,markersize=4)
 ax.axvline(1,color='gray',linestyle='--',linewidth=1);ax.set_xscale('log');ax.set_xlim(.33,1.9);ax.set_xticks([.4,.6,1,1.5],[.4,.6,1,1.5]);ax.grid(axis='x',alpha=.2)
 ax.set_yticks(range(7),ages);ax.set_xlabel('Rate ratio (women vs men)');ax.spines[['top','right']].set_visible(False);ax.text(-.12,1.03,letter,transform=ax.transAxes,fontweight='bold',fontsize=12)
axs[0].invert_yaxis();axs[0].set_ylabel('Age group (years)');fig.tight_layout()
fig.savefig(O/'Figura_2.png',dpi=300);fig.savefig(O/'Figura_2.pdf');plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(7,3.7),sharey=True)
for ax,ct,letter in zip(axs,['Chile','Brazil'],'AB'):
 for sex,color,marker in [('Women','#6A3D9A','s'),('Men','#606060','o')]:
  q=C[(C.country==ct)&(C.sex==sex)].set_index('age').loc[ages];n=q.selected_total.to_numpy(float);p=q.venous_count.to_numpy(float)/n;z=1.959964;mid=(p+z*z/(2*n))/(1+z*z/n);half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
  ax.errorbar(range(7),100*p,yerr=[100*(p-mid+half),100*(mid+half-p)],fmt=marker+'-',color=color,label=sex,capsize=2,markersize=3)
 ax.set_xticks(range(7),ages,rotation=45,ha='right');ax.set_ylim(0,55);ax.set_xlabel('Age group (years)');ax.text(-.12,1.04,letter,transform=ax.transAxes,fontweight='bold',fontsize=12);ax.spines[['top','right']].set_visible(False)
axs[0].set_ylabel('Venous proportion (%)');axs[1].legend(frameon=False);fig.tight_layout()
fig.savefig(O/'Figura_S1.png',dpi=300);fig.savefig(O/'Figura_S1.pdf');plt.close(fig)
