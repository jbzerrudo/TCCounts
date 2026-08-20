"""Step 4. Redraw Figure 1. Usage: python src/figure.py data/par_annual_series.csv figures/"""
import sys
import numpy as np, pandas as pd, scipy.stats as st, warnings, json; warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],'font.size':8,
 'axes.linewidth':0.7,'xtick.major.width':0.7,'ytick.major.width':0.7,'axes.spines.top':False,
 'axes.spines.right':False,'legend.frameon':False,'axes.labelsize':8,'legend.fontsize':7,'axes.titlesize':8.5})
BLUE='#0072B2'; ORANGE='#D55E00'; GREY='#666666'; TEAL='#009E73'
S=pd.read_csv((sys.argv[1] if len(sys.argv)>1 else 'data/par_annual_series.csv')).set_index('SEASON')
sw=pd.read_csv('data/start_year_sweep.csv')
fig,ax=plt.subplots(2,2,figsize=(7.2,5.3))

a=ax[0,0]
a.fill_between(S.index,0,S.any_int,color=BLUE,alpha=.85,lw=0,label='Classified by an agency')
a.fill_between(S.index,S.any_int,S.total,color=ORANGE,alpha=.85,lw=0,label='No intensity from any agency')
a.plot(S.index,S.total,color='k',lw=.8,label='All track entries')
for yr,txt,ha in [(1945,'JTWC 1945','right'),(1951,'JMA 1951','left'),(1977,'WMO 1977','left')]:
    a.axvline(yr,color='k',lw=.6,ls=':')
    a.annotate(txt,(yr,1.2),xytext=(-2 if ha=='right' else 2,0),textcoords='offset points',fontsize=6.2,
               ha=ha,va='bottom',bbox=dict(fc='white',ec='none',pad=0.7),zorder=6)
a.set_xlabel('Year'); a.set_ylabel('Storms entering the PAR'); a.set_ylim(0,40); a.set_xlim(1884,2023)
a.set_title('(a)  The archive fills in, it does not just vary',loc='left'); a.legend(loc='upper left',fontsize=6.4,ncol=1)

b=ax[0,1]
b.fill_between(sw.start,sw.lo,sw.hi,color=GREY,alpha=.16,lw=0)
b.plot(sw.start,sw.slope,color='k',lw=1.9,label='All track entries')
up=sw[sw.sig=='+']; dn=sw[sw.sig=='-']
b.plot(up.start,up.slope,'o',ms=3.4,color=TEAL,label='significant increase')
b.plot(dn.start,dn.slope,'o',ms=3.4,color=ORANGE,label='significant decrease')
c2=sw.dropna(subset=['cslope'])
b.plot(c2.start,c2.cslope,color=BLUE,lw=1.7,label='Classified only')
b.axhline(0,color='k',lw=.7)
b.set_xlabel('First year of the fitted window (all end 2023)'); b.set_ylabel('OLS slope (storms yr$^{-1}$)')
b.set_title('(b)  The trend is whatever the start year makes it',loc='left'); b.legend(loc='upper left',fontsize=6.3,ncol=1)

c=ax[1,0]
eras=[(1884,1922),(1923,1944),(1945,1950),(1951,1976),(1977,2000),(2001,2023)]
lab=[f"{x}-\n{y}" for x,y in eras]; x=np.arange(len(eras)); w=.38
cl=[S.loc[i:j,'any_int'].mean() for i,j in eras]; un=[S.loc[i:j,'pos_only'].mean() for i,j in eras]
c.bar(x-w/2,cl,w,color=BLUE,label='Classified'); c.bar(x+w/2,un,w,color=ORANGE,label='Unclassified')
for xi,v in zip(x-w/2,cl): c.annotate(f'{v:.1f}',(xi,v),xytext=(0,2),textcoords='offset points',ha='center',fontsize=6.3,color=BLUE)
for xi,v in zip(x+w/2,un): c.annotate(f'{v:.1f}',(xi,v),xytext=(0,2),textcoords='offset points',ha='center',fontsize=6.3,color=ORANGE)
c.set_xticks(x); c.set_xticklabels(lab,fontsize=6.3); c.set_ylabel('Mean storms per year'); c.set_ylim(0,25)
c.set_title('(c)  One column empties as the other fills',loc='left'); c.legend(loc='upper left',fontsize=6.6)

e=ax[1,1]
T2=json.load(open('data/backtest.json'))
names=[r[0] for r in T2][::-1]; mae=[float(r[1]) for r in T2][::-1]
cols=[ORANGE if n.startswith('Linear') else (GREY if n.startswith('Last') else TEAL) for n in names]
e.barh(range(len(names)),mae,color=cols)
for i,v in enumerate(mae): e.annotate(f'{v:.2f}',(v,i),xytext=(3,0),textcoords='offset points',va='center',fontsize=7)
e.set_yticks(range(len(names))); e.set_yticklabels(names,fontsize=6.8)
e.set_xlabel('MAE predicting the next 20-year mean (storms)'); e.set_xlim(0,3.3)
e.set_title('(d)  Climatology beats trend extrapolation',loc='left')
plt.tight_layout(pad=.6,w_pad=1.9,h_pad=1.7)
plt.savefig(''+(sys.argv[2] if len(sys.argv)>2 else 'figures/')+'Fig1.png',dpi=600,bbox_inches='tight'); plt.savefig(''+(sys.argv[2] if len(sys.argv)>2 else 'figures/')+'Fig1.pdf',bbox_inches='tight')
print("saved")
