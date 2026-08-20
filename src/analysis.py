"""
Step 3. Recompute every value quoted in the manuscript and check each one.
Prints OK or MISMATCH beside every number. A clean run prints zero MISMATCH.

Usage   python src/analysis.py data/par_annual_series.csv data/start_year_sweep.csv
"""
import sys
import pandas as pd, numpy as np, scipy.stats as st, warnings; warnings.filterwarnings("ignore")
S=pd.read_csv((sys.argv[1] if len(sys.argv)>1 else 'data/par_annual_series.csv')).set_index('SEASON')
sw=pd.read_csv((sys.argv[2] if len(sys.argv)>2 else 'data/start_year_sweep.csv'))
ok=lambda n,a,b,t: print(f"  {'OK ' if abs(float(a)-b)<=t else 'MISMATCH'}  {n:42s} draft={b}  recomputed={round(float(a),4)}")
ok("n 1884-2023",len(S),140,0); ok("mean",S.total.mean(),20.00,.005)
ok("max",S.total.max(),37,0); ok("argmax",S.total.idxmax(),1971,0)
ok("min",S.total.min(),8,0); ok("argmin",S.total.idxmin(),1885,0)
ok("storms 1884-1944",S.loc[1884:1944,'total'].sum(),1009,0)
ok("classified 1884-1944",S.loc[1884:1944,'any_int'].sum(),0,0)
for a,b,v in [(1884,1922,14.95),(1923,1944,19.36),(1945,1950,21.67),(1951,1976,25.81),(1977,2000,22.79),(2001,2023,19.26)]:
    ok(f"era total {a}-{b}",S.loc[a:b,'total'].mean(),v,.006)
for a,b,v in [(1945,1950,14.7),(1951,1976,18.8),(1977,2000,20.4),(2001,2023,18.9)]:
    ok(f"era classified {a}-{b}",S.loc[a:b,'any_int'].mean(),v,.05)
for a,b,v in [(1945,1950,7.0),(1951,1976,7.0),(1977,2000,2.4),(2001,2023,0.4)]:
    ok(f"era unclassified {a}-{b}",S.loc[a:b,'pos_only'].mean(),v,.05)
for a,b,v in [(1884,1944,100),(1945,1950,32),(1951,1976,27),(1977,2000,10),(2001,2023,2)]:
    ok(f"unclassified share {a}-{b} (%)",100*S.loc[a:b,'pos_only'].mean()/S.loc[a:b,'total'].mean(),v,.7)
for col,a,b,v in [('total',1884,2023,0.0608),('any_int',1884,2023,0.2061),('pos_only',1884,2023,-0.1452),
                  ('total',1951,2023,-0.1215),('any_int',1951,2023,0.0108),('pos_only',1951,2023,-0.1323),
                  ('any_int',1977,2023,-0.0287)]:
    s=S.loc[a:b,col]; ok(f"slope {col} {a}-{b}",st.linregress(s.index,s.values).slope,v,.0005)
c=st.linregress(S.index,S.any_int).slope; u=st.linregress(S.index,S.pos_only).slope
ok("identity 1884-2023",c+u,0.0608,.0005)
c2=st.linregress(S.loc[1951:2023].index,S.loc[1951:2023,'any_int']).slope; u2=st.linregress(S.loc[1951:2023].index,S.loc[1951:2023,'pos_only']).slope
ok("identity 1951-2023",c2+u2,-0.1215,.0005)
for col,a,b,v in [('total',1884,2023,0.2010),('pos_only',1884,2023,0.6155),('total',1951,2023,0.2650),('pos_only',1951,2023,0.5764)]:
    s=S.loc[a:b,col]; ok(f"R2 {col} {a}-{b}",st.linregress(s.index,s.values).rvalue**2,v,.001)
ok("sweep start years sampled",len(sw),59,0)
ok("significant increases",(sw.sig=='+').sum(),4,0); ok("significant decreases",(sw.sig=='-').sum(),20,0)
ok("min slope in sweep",sw.slope.min(),-0.1595,.0005); ok("max slope in sweep",sw.slope.max(),0.0608,.0005)
ok("last significant-increase start",sw[sw.sig=='+'].start.max(),1890,0)
ok("classified sweep n",sw.csig.notna().sum(),25,0); ok("classified significant",(sw.csig.dropna()!=0).sum(),0,0)
ok("classified mean 1977-2023",S.loc[1977:2023,'any_int'].mean(),19.66,.006)
lo,hi=st.poisson.ppf([.025,.975],S.loc[1977:2023,'any_int'].mean()); ok("Poisson lo",lo,11,0); ok("Poisson hi",hi,29,0)
s=S.loc[1951:2023,'any_int']; org=list(range(1971,2004)); ok("backtest origins",len(org),33,0)
E={"trail":[],"full":[],"lin":[],"last":[]}
for t in org:
    tr=s[s.index<=t]; tg=s[(s.index>t)&(s.index<=t+20)].mean(); lr=st.linregress(tr.index,tr.values)
    E["trail"].append(tr.iloc[-20:].mean()-tg); E["full"].append(tr.mean()-tg)
    E["lin"].append(lr.slope*(t+10)+lr.intercept-tg); E["last"].append(tr.iloc[-1]-tg)
for k,v in [("full",1.12),("trail",1.32),("lin",1.62),("last",3.03)]: ok(f"backtest MAE {k}",np.abs(np.array(E[k])).mean(),v,.006)
ok("backtest bias linear",np.mean(E["lin"]),0.60,.006)
ok("trend-extrap error premium (%)",100*(np.abs(np.array(E['lin'])).mean()/np.abs(np.array(E['full'])).mean()-1),45,1.0)
