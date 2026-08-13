import os, sys, math
import numpy as np
import pandas as pd
from scipy.stats import linregress
ROOT=os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(ROOT,'src'))
from fidelity_model import *

alphas=[0.8,1.2,1.6]; rhos=[0.15,0.30]
Ns=[64,128,256,512,1024,2048,4096,8192,16384,32768]
rows_modes=[]; rows_disp=[]; rows_static=[]; rows_tail=[]
qgrid=np.linspace(0,np.pi,401)

# Publication dispersion uses the same infinite-lattice polylogarithmic symbol as the theory.
# This avoids a hidden real-space cutoff in Figure 2.
for a in alphas:
    lam_grid=np.array([lambda_symbol_float(a,q) for q in qgrid],dtype=float)
    for M in [16,32,64,128,256,512,1024,2048,4096,8192,16384,32768]:
        exact=float(tail(a,M)); asym=1/(a*zeta_norm(a))*M**(-a)
        rows_tail.append(dict(alpha=a,M=M,tail=exact,asymptotic=asym,relative_error=(asym-exact)/exact))
    for rho in rhos:
        qstar,gmax,qc,gc,P,m,ratio,lamc=select_commensurate_mode(a,rho)
        rows_modes.append(dict(alpha=a,rho=rho,q_star=qstar,g_max=gmax,q_c=qc,g_c=gc,period=P,mode_index=m,g_ratio=ratio,lambda_qc=lamc))
        gains=np.sqrt(np.clip(lam_grid*(2*rho-lam_grid),0,None))
        for q,lam,gain in zip(qgrid,lam_grid,gains):
            rows_disp.append(dict(alpha=a,rho=rho,q=q,lambda_inf=lam,mi_gain=gain))
    # Cache open operators once per alpha,N and evaluate both rho-selected q values.
    modes=[r for r in rows_modes if r['alpha']==a]
    for N in Ns:
        op=OpenOperator(a,N); n=np.arange(N); sl=central_slice(N,.25)
        for md in modes:
            rho=md['rho']; P=md['period']; qc=md['q_c']
            if N%P!=0: continue
            v=np.cos(qc*n); lam=lambda_exact(a,qc); Linf=-lam*v
            Lo=np.real_if_close(op.linear_open(v)); Lc=np.real_if_close(Lo-op.kout*v)
            den=np.sqrt(np.mean(v[sl]**2))
            eo=np.sqrt(np.mean(np.abs(Lo[sl]-Linf[sl])**2))/den
            ec=np.sqrt(np.mean(np.abs(Lc[sl]-Linf[sl])**2))/den
            J=kernel(a,N//2); r=np.arange(1,N//2)
            lamn=2*np.sum(J[:N//2-1]*(1-np.cos(r*qc)))
            if N%2==0: lamn+=J[N//2-1]*(1-np.cos((N//2)*qc))
            en=abs(lamn-lam)
            rows_static.append(dict(alpha=a,rho=rho,N=N,q_c=qc,period=P,open_error=eo,corrected_error=ec,naive_periodic_error=en,exact_periodic_error=0.0))

df_modes=pd.DataFrame(rows_modes); df_disp=pd.DataFrame(rows_disp); df_static=pd.DataFrame(rows_static); df_tail=pd.DataFrame(rows_tail)
df_modes.to_csv(os.path.join(ROOT,'data/theory/selected_modes.csv'),index=False)
df_disp.to_csv(os.path.join(ROOT,'data/theory/dispersion_mi.csv'),index=False)
df_tail.to_csv(os.path.join(ROOT,'data/theory/kernel_tail.csv'),index=False)
df_static.to_csv(os.path.join(ROOT,'data/static/static_closure_scaling.csv'),index=False)
slopes=[]; fits=[]
for (a,rho),g in df_static.groupby(['alpha','rho']):
    g=g.sort_values('N').copy()
    for col,name in [('open_error','open'),('corrected_error','corrected'),('naive_periodic_error','naive_periodic')]:
        vals=g[col].to_numpy(); ns=g.N.to_numpy()
        loc=np.r_[np.nan,-np.log(vals[1:]/vals[:-1])/np.log(ns[1:]/ns[:-1])]
        for N,v,b in zip(ns,vals,loc): slopes.append(dict(alpha=a,rho=rho,closure=name,N=N,error=v,beta_local=b))
        gg=g.tail(4) if name=='corrected' else g[g.N>=256]
        lr=linregress(np.log(gg.N),np.log(gg[col]))
        fits.append(dict(alpha=a,rho=rho,closure=name,beta=-lr.slope,logA=lr.intercept,r2=lr.rvalue**2,N_min=int(gg.N.min()),N_max=int(gg.N.max()),n_points=len(gg)))
pd.DataFrame(slopes).to_csv(os.path.join(ROOT,'data/static/static_local_slopes.csv'),index=False)
pd.DataFrame(fits).to_csv(os.path.join(ROOT,'data/static/static_fits.csv'),index=False)
print('MODES')
print(df_modes.to_string(index=False))
print('\nFITS')
print(pd.DataFrame(fits).to_string(index=False))
