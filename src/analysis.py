"""
Step 3. Recompute every value quoted in the manuscript and check each one.
Prints OK or MISMATCH beside every number. A clean run prints zero MISMATCH (58 checks).

Usage   python src/analysis.py data/par_annual_series.csv data/start_year_sweep.csv
"""
import sys, pandas as pd, numpy as np, scipy.stats as st, warnings; warnings.filterwarnings("ignore")
S=pd.read_csv(sys.argv[1] if len(sys.argv)>1 else 'data/par_annual_series.csv').set_index('SEASON')
sw=pd.read_csv(sys.argv[2] if len(sys.argv)>2 else 'data/start_year_sweep.csv')
ok=lambda n,a,b,t: print(("  OK  " if abs(float(a)-b)<=t else "  MISMATCH")+"  %-40s draft=%s recomputed=%s"%(n,b,round(float(a),4)))
ok("n",len(S),140,0); ok("mean",S.total.mean(),19.69,.006); ok("storms total",S.total.sum(),2757,0)
ok("max",S.total.max(),35,0); ok("argmax",S.total.idxmax(),1964,0)
ok("min",S.total.min(),8,0); ok("argmin",S.total.idxmin(),1885,0)
ok("storms 1884-1944",S.loc[1884:1944,'total'].sum(),1004,0); ok("classified 1884-1944",S.loc[1884:1944,'any_int'].sum(),0,0)
for a,b,v in [(1884,1922,14.87),(1951,1976,25.12),(2001,2023,19.26)]: ok("era total %d-%d"%(a,b),S.loc[a:b,'total'].mean(),v,.006)
for a,b,v in [(1945,1950,14.7),(1951,1976,18.6),(1977,2000,20.0),(2001,2023,18.9)]: ok("era classified %d-%d"%(a,b),S.loc[a:b,'any_int'].mean(),v,.05)
for a,b,v in [(1945,1950,6.3),(1951,1976,6.5),(1977,2000,2.1),(2001,2023,0.4)]: ok("era unclassified %d-%d"%(a,b),S.loc[a:b,'pos_only'].mean(),v,.05)
for a,b,v in [(1884,1944,100),(1945,1950,30),(1951,1976,26),(1977,2000,10),(2001,2023,2)]:
    ok("unclassified share %d-%d"%(a,b),100*S.loc[a:b,'pos_only'].mean()/S.loc[a:b,'total'].mean(),v,.7)
for c,a,b,v in [('total',1884,2023,0.0584),('any_int',1884,2023,0.2042),('pos_only',1884,2023,-0.1458),
                ('total',1951,2023,-0.1086),('any_int',1951,2023,0.0131),('pos_only',1951,2023,-0.1217),
                ('any_int',1977,2023,-0.0192)]:
    s=S.loc[a:b,c]; ok("slope %s %d-%d"%(c,a,b),st.linregress(s.index,s.values).slope,v,.0005)
for c,a,b,v in [('total',1884,2023,0.204),('pos_only',1884,2023,0.619),('total',1951,2023,0.237),('pos_only',1951,2023,0.551)]:
    s=S.loc[a:b,c]; ok("R2 %s %d-%d"%(c,a,b),st.linregress(s.index,s.values).rvalue**2,v,.001)
c=st.linregress(S.index,S.any_int).slope; u=st.linregress(S.index,S.pos_only).slope
ok("identity 1884-2023",c+u,0.0584,.0005)
c2=st.linregress(S.loc[1951:].index,S.loc[1951:,'any_int']).slope; u2=st.linregress(S.loc[1951:].index,S.loc[1951:,'pos_only']).slope
ok("identity 1951-2023",c2+u2,-0.1086,.0005)
ok("sweep start years",len(sw),59,0); ok("significant rises",(sw.sig=='+').sum(),4,0); ok("significant declines",(sw.sig=='-').sum(),16,0)
ok("sweep min slope",sw.slope.min(),-0.1455,.0005); ok("sweep max slope",sw.slope.max(),0.0584,.0005)
ok("last rise start",sw[sw.sig=='+'].start.max(),1890,0)
ok("first decline start",sw[sw.sig=='-'].start.min(),1946,0); ok("last decline start",sw[sw.sig=='-'].start.max(),1978,0)
ok("classified sweep n",sw.csig.notna().sum(),25,0); ok("classified significant",sw.csig.dropna().astype(str).isin(['+','-']).sum(),0,0)
m=S.loc[1977:2023,'any_int'].mean(); ok("classified mean 1977-2023",m,19.45,.006)
lo,hi=st.poisson.ppf([.025,.975],m); ok("Poisson lo",lo,11,0); ok("Poisson hi",hi,29,0)
s=S.loc[1951:2023,'any_int']; org=list(range(1971,2004)); ok("backtest origins",len(org),33,0)
E={"trail":[],"full":[],"lin":[],"last":[]}
for t in org:
    tr=s[s.index<=t]; tg=s[(s.index>t)&(s.index<=t+20)].mean(); lr=st.linregress(tr.index,tr.values)
    E["trail"].append(tr.iloc[-20:].mean()-tg); E["full"].append(tr.mean()-tg)
    E["lin"].append(lr.slope*(t+10)+lr.intercept-tg); E["last"].append(tr.iloc[-1]-tg)
for k,v in [("full",0.95),("trail",1.09),("lin",1.37),("last",2.79)]: ok("backtest MAE %s"%k,np.abs(np.array(E[k])).mean(),v,.006)
ok("backtest bias linear",np.mean(E["lin"]),0.61,.006)
ok("trend premium (%)",100*(np.abs(np.array(E['lin'])).mean()/np.abs(np.array(E['full'])).mean()-1),44,1.0)
