"""Step 4. Redraw Figure 1. Usage: python src/figure.py data/par_annual_series.csv figures/"""
import sys
import numpy as np, pandas as pd, scipy.stats as st, warnings, json; warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans'],'font.size':8,
 'axes.linewidth':0.7,'xtick.major.width':0.7,'ytick.major.width':0.7,'axes.spines.top':False,
 'axes.spines.right':False,'legend.frameon':False,'axes.labelsize':8,'legend.fontsize':7,'axes.titlesize':8.5})
BLUE='#0072B2'; ORANGE='#D55E00'; GREY='#666666'; TEAL='#009E73'
S=pd.read_csv((sys.argv[1] if len(sys.argv)>1 else 'data/par_annual_series.csv')).set_index('SEASON')
rng=np.random.default_rng(3)
def band(s,L=10,n=2500):
    v=s.values; N=len(v); out=[]
    for _ in range(n):
        idx=np.concatenate([np.arange(i,i+L)%N for i in rng.integers(0,N,size=int(np.ceil(N/L)))])[:N]
        out.append(st.linregress(s.index,v[idx]).slope)
    return np.percentile(out,[2.5,97.5])

fig,ax=plt.subplots(2,2,figsize=(7.2,5.2))
# (a) the three series
a=ax[0,0]
a.plot(S.index,S.total,color=GREY,lw=1.1,label='All track entries in the PAR')
a.plot(S.index,S.jma_class,color=BLUE,lw=1.6,label='Intensity-classified by JMA')
a.plot(S.index,S.pos_only,color=ORANGE,lw=1.4,ls='--',label='No intensity from any agency')
a.axvline(1945,color='k',lw=0.6,ls=':'); a.axvline(1951,color='k',lw=0.6,ls=':')
a.annotate('JTWC\n1945',(1945,38),fontsize=6.5,ha='right',va='top',color='k')
a.annotate('JMA\n1951',(1952,38),fontsize=6.5,ha='left',va='top',color='k')
a.set_xlabel('Year'); a.set_ylabel('Storms entering the PAR'); a.set_ylim(0,40)
a.set_title('(a)  Two components moving in opposite directions',loc='left'); a.legend(loc='upper right',ncol=1,fontsize=6.6)
# (b) trend vs start year
b=ax[0,1]
for col,c,lab,y0 in [('total',GREY,'All track entries',1923),('jma_class',BLUE,'JMA-classified',1951)]:
    xs=list(range(y0,2001,3)); sl=[]; lo=[]; hi=[]
    for st0 in xs:
        s=S.loc[st0:2023,col]; sl.append(st.linregress(s.index,s.values).slope); l,h=band(s); lo.append(l); hi.append(h)
    b.fill_between(xs,lo,hi,color=c,alpha=0.13,lw=0)
    b.plot(xs,sl,color=c,lw=1.8,label=lab)
b.axhline(0,color='k',lw=0.7)
b.set_xlabel('First year of the fitted window (all end 2023)'); b.set_ylabel('OLS slope (storms yr$^{-1}$)')
b.set_title('(b)  The decline exists only in the raw count',loc='left'); b.legend(loc='lower right')
b.annotate('shading = 95% range under\nblock-resampled no-trend null',(1962,0.115),fontsize=6.5,color=GREY)
# (c) era means
c=ax[1,0]
eras=[(1951,1976),(1977,2000),(2001,2023)]; lab=[f"{a0}-{b0}" for a0,b0 in eras]
tot=[S.loc[a0:b0,'total'].mean() for a0,b0 in eras]; cls=[S.loc[a0:b0,'jma_class'].mean() for a0,b0 in eras]
un=[S.loc[a0:b0,'pos_only'].mean() for a0,b0 in eras]
x=np.arange(3); w=0.27
c.bar(x-w,tot,w,color=GREY,label='All entries'); c.bar(x,cls,w,color=BLUE,label='JMA-classified'); c.bar(x+w,un,w,color=ORANGE,label='Unclassified')
for xi,v in zip(x-w,tot): c.annotate(f'{v:.1f}',(xi,v),xytext=(0,2),textcoords='offset points',ha='center',fontsize=6.5,color=GREY)
for xi,v in zip(x,cls): c.annotate(f'{v:.1f}',(xi,v),xytext=(0,2),textcoords='offset points',ha='center',fontsize=6.5,color=BLUE)
for xi,v in zip(x+w,un): c.annotate(f'{v:.1f}',(xi,v),xytext=(0,2),textcoords='offset points',ha='center',fontsize=6.5,color=ORANGE)
c.set_xticks(x); c.set_xticklabels(lab); c.set_ylabel('Mean storms per year'); c.set_ylim(0,30)
c.set_title('(c)  The raw decline is the unclassified column',loc='left'); c.legend(loc='upper right',ncol=1)
# (d) backtest
e=ax[1,1]
s=S.loc[1951:2023,'jma_class']; origins=list(range(1971,2004))
meth={"20-yr trailing mean":[], "Full-record mean":[], "Last observed value":[], "Linear trend extrapolated":[]}
for t in origins:
    tr=s[s.index<=t]; tgt=s[(s.index>t)&(s.index<=t+20)].mean(); lr=st.linregress(tr.index,tr.values)
    meth["20-yr trailing mean"].append(tr.iloc[-20:].mean()-tgt); meth["Full-record mean"].append(tr.mean()-tgt)
    meth["Last observed value"].append(tr.iloc[-1]-tgt); meth["Linear trend extrapolated"].append(lr.slope*(t+10)+lr.intercept-tgt)
names=list(meth); mae=[np.abs(np.array(meth[k])).mean() for k in names]
o=np.argsort(mae); names=[names[i] for i in o]; mae=[mae[i] for i in o]
cols=[TEAL if i<len(names)-1 else ORANGE for i in range(len(names))]
e.barh(range(len(names)),mae,color=cols)
for i,v in enumerate(mae): e.annotate(f'{v:.2f}',(v,i),xytext=(3,0),textcoords='offset points',va='center',fontsize=7)
e.set_yticks(range(len(names))); e.set_yticklabels(names,fontsize=7); e.invert_yaxis()
e.set_xlabel('MAE predicting the next 20-year mean (storms)'); e.set_xlim(0,4.0)
e.set_title('(d)  Trend extrapolation is the worst method',loc='left')
plt.tight_layout(pad=0.6,w_pad=1.9,h_pad=1.7)
plt.savefig(''+(sys.argv[2] if len(sys.argv)>2 else 'figures/')+'Fig1.png',dpi=600,bbox_inches='tight'); plt.savefig(''+(sys.argv[2] if len(sys.argv)>2 else 'figures/')+'Fig1.pdf',bbox_inches='tight')
print("saved. backtest MAE:", {n:round(m,3) for n,m in zip(names,mae)})
