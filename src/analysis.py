"""
Step 3. Recompute every value quoted in the manuscript and check each one.
Prints OK or MISMATCH beside every number. A clean run prints zero MISMATCH.

Usage   python src/analysis.py data/par_annual_series.csv par_clipped.csv
        (the second argument is optional; it enables the archive-provenance checks)
"""
import sys
import pandas as pd, numpy as np, scipy.stats as st, warnings; warnings.filterwarnings("ignore")
import statsmodels.api as sm, statsmodels.formula.api as smf
from statsmodels.tsa.stattools import acf
from statsmodels.stats.diagnostic import acorr_ljungbox
S=pd.read_csv((sys.argv[1] if len(sys.argv)>1 else 'data/par_annual_series.csv')).set_index('SEASON')
d=pd.read_csv((sys.argv[2] if len(sys.argv)>2 else "par_clipped.csv"), low_memory=False)
d=d[(d.SEASON>=1923)&(d.SEASON<=2023)]
ok=lambda n,a,b,tol: print(f"  {'OK ' if abs(a-b)<=tol else 'MISMATCH'}  {n:44s} draft={b}  recomputed={round(float(a),4)}")
ok("n",len(S),101,0); ok("mean total",S.total.mean(),21.95,.005); ok("max total",S.total.max(),37,0)
ok("argmax year",S.total.idxmax(),1971,0)
for col,sl,r2,p in [('total',-0.0070,0.0018,0.673),('any_int',0.2082,0.505,9e-17),('pos_only',-0.2152,0.708,3e-28)]:
    lr=st.linregress(S.index,S[col]); ok(f"slope {col} 1923-2023",lr.slope,sl,.0005)
    ok(f"R2 {col}",lr.rvalue**2,r2,.001)
a=st.linregress(S.index,S.any_int).slope; b=st.linregress(S.index,S.pos_only).slope
ok("components sum to total slope",a+b,-0.0070,.0005)
for col,a0,b0,sl in [('total',1951,2023,-0.1215),('total',1977,2023,-0.1220),('jma_class',1951,2023,0.0749),
                     ('jma_class',1977,2023,-0.0457),('pos_only',1951,2023,-0.1323)]:
    s=S.loc[a0:b0,col]; ok(f"slope {col} {a0}-{b0}",st.linregress(s.index,s.values).slope,sl,.0005)
ok("total decline 1951-2023 (storms)",st.linregress(S.loc[1951:2023].index,S.loc[1951:2023,'total']).slope*73,-8.9,.05)
ok("unclassified decline 1951-2023",st.linregress(S.loc[1951:2023].index,S.loc[1951:2023,'pos_only']).slope*73,-9.7,.05)
for a0,b0,t,c,u in [(1951,1976,25.8,12.4,7.0),(1977,2000,22.8,17.2,2.4),(2001,2023,19.3,16.2,0.4)]:
    w=S.loc[a0:b0]
    ok(f"era total {a0}-{b0}",w.total.mean(),t,.05); ok(f"era classified {a0}-{b0}",w.jma_class.mean(),c,.05)
    ok(f"era unclassified {a0}-{b0}",w.pos_only.mean(),u,.05)
ok("unclassified share 1951-76 (%)",100*S.loc[1951:1976,'pos_only'].mean()/S.loc[1951:1976,'total'].mean(),27,.6)
ok("unclassified share 2001-23 (%)",100*S.loc[2001:2023,'pos_only'].mean()/S.loc[2001:2023,'total'].mean(),2,.6)
per=(S.pos_only/S.total)
for a0,b0,v in [(1923,1944,100),(1945,1950,32),(1951,1960,26),(1961,1970,31),(1977,1986,17),(1987,1999,5),(2000,2023,2)]:
    ok(f"unclassified fraction {a0}-{b0} (%)",100*per.loc[a0:b0].mean(),v,0.6)
for c,y in [('USA_LAT',1945),('TOK_GRADE',1951),('WMO_WIND',1977),('TD6_LAT',1923)]:
    v=pd.to_numeric(d[c],errors='coerce'); m=v.notna()&(v!=0); ok(f"first PAR season with {c}",d.loc[m,'SEASON'].min(),y,0)
v=pd.to_numeric(d['TD6_LAT'],errors='coerce'); m=v.notna()&(v!=0); ok("last PAR season with TD6_LAT",d.loc[m,'SEASON'].max(),1989,0)
g=pd.to_numeric(d['TOK_LAT'],errors='coerce'); mm=g.notna()&(g!=0)
ok("JMA fix coverage 1950s (%)",100*mm[(d.SEASON>=1950)&(d.SEASON<=1959)].mean(),71,1.0)
ok("JMA fix coverage 2010s (%)",100*mm[(d.SEASON>=2010)&(d.SEASON<=2019)].mean(),88,1.0)
y=S.loc[1977:2023,'jma_class']; ok("mean classified 1977-2023",y.mean(),16.7,.05)
df2=pd.DataFrame({'y':y.values}); g0=smf.glm("y ~ 1",data=df2,family=sm.families.Poisson()).fit()
pr=(df2.y-g0.fittedvalues)/np.sqrt(g0.fittedvalues)
ok("dispersion 1977-2023",pr.var(ddof=1),0.644,.002); ok("ACF1 1977-2023",acf(pr,nlags=1)[1],0.269,.002)
ok("LB(10) p 1977-2023",acorr_ljungbox(pr,lags=[10],return_df=True)['lb_pvalue'].iloc[0],0.856,.003)
lo,hi=st.poisson.ppf([.025,.975],y.mean()); ok("Poisson lo",lo,9,0); ok("Poisson hi",hi,25,0)
s=S.loc[1951:2023,'jma_class']; origins=list(range(1971,2004)); ok("backtest origins",len(origins),33,0)
E={"trail":[],"full":[],"last":[],"lin":[]}
for t in origins:
    tr=s[s.index<=t]; tgt=s[(s.index>t)&(s.index<=t+20)].mean(); lr=st.linregress(tr.index,tr.values)
    E["trail"].append(tr.iloc[-20:].mean()-tgt); E["full"].append(tr.mean()-tgt)
    E["last"].append(tr.iloc[-1]-tgt); E["lin"].append(lr.slope*(t+10)+lr.intercept-tgt)
for k,v in [("trail",2.19),("full",2.86),("last",3.12),("lin",3.24)]:
    ok(f"backtest MAE {k}",np.abs(np.array(E[k])).mean(),v,.006)
ok("backtest bias linear",np.mean(E["lin"]),2.55,.006)
