import os,sys,math,time
import numpy as np
import pandas as pd
from scipy.optimize import brentq
ROOT=os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0,os.path.join(ROOT,'src'))
from fidelity_model import *

modes=pd.read_csv(os.path.join(ROOT,'data/theory/selected_modes.csv'))
static_fits=pd.read_csv(os.path.join(ROOT,'data/static/static_fits.csv'))
Ns_all=[128,256,512,1024]
h=0.02; eps=1e-6; save_dt=0.10
all_diag=[]; all_horiz=[]; all_fit=[]; all_lono=[]; validation=[]
map_rows=[]

def rk4_pair(w,h,op,rho):
    def rhs(x,rr):
        Lo=op.linear_open(x)
        L=Lo.copy()
        L[1]=L[1]+op.kout*(1-x[1])
        return 1j*(L+rr*(np.abs(x)**2-1)*x)
    return rk4_step(w,h,rhs,rho)

def simulate_case(alpha,rho,qc,gc,P,m,Zmax):
    Ns=[N for N in Ns_all if N%P==0]
    per=PeriodicExactOperator(alpha,int(P))
    ops={N:OpenOperator(alpha,N) for N in Ns}
    wr=1+eps*np.cos(qc*np.arange(P))
    states={N:np.vstack([1+eps*np.cos(qc*np.arange(N)),1+eps*np.cos(qc*np.arange(N))]).astype(complex) for N in Ns}
    p0_ref=float(np.sum(np.abs(wr)**2)); p0={N:float(np.sum(np.abs(states[N][0])**2)) for N in Ns}
    save_every=max(1,int(round(save_dt/h))); steps=int(round(Zmax/h))
    rows=[]
    maps=[]
    for k in range(steps+1):
        if k%save_every==0:
            Z=k*h; vv=wr-1; amp=float(abs(np.fft.fft(vv)[m])/P)
            pref=float(np.sum(np.abs(wr)**2)); mxref=float(np.max(np.abs(wr)**2))
            for N in Ns:
                ref=periodic_repeat(wr,N); sl=central_slice(N,.25); wo,wc=states[N]
                eo=float(np.sqrt(np.mean(np.abs(wo[sl]-ref[sl])**2)))
                ec=float(np.sqrt(np.mean(np.abs(wc[sl]-ref[sl])**2)))
                eio=float(np.sqrt(np.mean((np.abs(wo[sl])**2-np.abs(ref[sl])**2)**2)))
                eic=float(np.sqrt(np.mean((np.abs(wc[sl])**2-np.abs(ref[sl])**2)**2)))
                po=float(np.sum(np.abs(wo)**2))
                rows.append(dict(alpha=alpha,rho=rho,N=N,Z=Z,tau=gc*Z,reference_mode_amplitude=amp,reference_max_intensity=mxref,open_field_error=eo,corrected_field_error=ec,open_intensity_error=eio,corrected_intensity_error=eic,reference_power_rel=(pref-p0_ref)/p0_ref,open_power_rel=(po-p0[N])/p0[N]))
                if alpha==0.8 and rho==0.30 and N==256 and abs((Z/0.2)-round(Z/0.2))<1e-8:
                    center=np.arange(N)
                    # Map rows are intentionally long-form CSV for full reproducibility.
                    for n in range(N):
                        maps.append((Z,n,float(abs(ref[n])**2),float(abs(wo[n])**2),float(abs(wc[n])**2)))
        if k<steps:
            wr=rk4_step(wr,h,per.rhs,rho)
            for N in Ns:
                states[N]=rk4_pair(states[N],h,ops[N],rho)
    return pd.DataFrame(rows),maps

for _,md in modes.iterrows():
    a=float(md.alpha); rho=float(md.rho); qc=float(md.q_c); gc=float(md.g_c); P=int(md.period); m=int(md.mode_index)
    Zmax=90.0 if rho<0.2 else 48.0
    print(f'Running alpha={a}, rho={rho}, P={P}, Zmax={Zmax}',flush=True)
    t0=time.time(); df,maps=simulate_case(a,rho,qc,gc,P,m,Zmax); print(' seconds',time.time()-t0,flush=True)
    all_diag.append(df)
    map_rows.extend(maps)
    # Numerical MI fit from unique reference series.
    ref=df[df.N==df.N.min()][['Z','reference_mode_amplitude']].drop_duplicates().sort_values('Z')
    mask=(ref.reference_mode_amplitude>=2e-6)&(ref.reference_mode_amplitude<=5e-4)
    if mask.sum()<10: mask=(ref.Z>=0.15*Zmax)&(ref.Z<=0.55*Zmax)
    gnum=np.polyfit(ref.loc[mask,'Z'],np.log(ref.loc[mask,'reference_mode_amplitude']),1)[0]
    # Static beta for matching case.
    beta=float(static_fits[(static_fits.alpha==a)&(static_fits.rho==rho)&(static_fits.closure=='open')].beta.iloc[0])
    # Dynamic fitting window: same MI amplitude interval, exclude tiny error floor.
    fitdf=df.merge(ref.assign(in_mi=mask.to_numpy()),on=['Z','reference_mode_amplitude'],how='left')
    fitdf=fitdf[(fitdf.in_mi==True)&(fitdf.open_field_error>1e-10)&(fitdf.Z>0)]
    y=np.log(fitdf.open_field_error.to_numpy()) + beta*np.log(fitdf.N.to_numpy()) - gc*fitdf.Z.to_numpy()
    X=np.c_[np.ones(len(fitdf)),np.log(fitdf.Z.to_numpy())]
    coef=np.linalg.lstsq(X,y,rcond=None)[0]; logA,mu=coef
    resid=y-X@coef; r2=1-np.sum(resid**2)/np.sum((y-y.mean())**2)
    # Free diagnostic fit logE = c - beta_dyn logN + mu_free logZ + g_eff Z
    Y=np.log(fitdf.open_field_error.to_numpy())
    XF=np.c_[np.ones(len(fitdf)),-np.log(fitdf.N.to_numpy()),np.log(fitdf.Z.to_numpy()),fitdf.Z.to_numpy()]
    cf=np.linalg.lstsq(XF,Y,rcond=None)[0]
    all_fit.append(dict(alpha=a,rho=rho,g_theory=gc,g_num=gnum,g_ratio=gnum/gc,beta_static=beta,logA=logA,A=np.exp(logA),mu=mu,restricted_r2=r2,beta_dynamic=cf[1],mu_free=cf[2],g_eff=cf[3],g_eff_ratio=cf[3]/gc,n_fit=len(fitdf)))
    # Horizons and predictions.
    A=np.exp(logA)
    for N,gN in df.groupby('N'):
        gN=gN.sort_values('Z')
        for closure,col in [('open','open_field_error'),('corrected','corrected_field_error')]:
            for thr in [1e-3,1e-2]:
                obs=first_sustained_crossing(gN.Z.to_numpy(),gN[col].to_numpy(),thr,2)
                pred=np.nan
                if closure=='open':
                    f=lambda z: A*N**(-beta)*z**mu*np.exp(gc*z)-thr
                    try:
                        pred=brentq(f,max(h,1e-4),max(Zmax*2,10))
                    except Exception: pass
                all_horiz.append(dict(alpha=a,rho=rho,N=int(N),closure=closure,threshold=thr,Z_observed=obs,Z_predicted=pred,beta=beta,mu=mu,g=gc))
    # Leave-one-N-out restricted fits/predictions for 1e-3 open horizon.
    for Nhold in sorted(df.N.unique()):
        tr=fitdf[fitdf.N!=Nhold]
        yy=np.log(tr.open_field_error.to_numpy())+beta*np.log(tr.N.to_numpy())-gc*tr.Z.to_numpy()
        XX=np.c_[np.ones(len(tr)),np.log(tr.Z.to_numpy())]
        cc=np.linalg.lstsq(XX,yy,rcond=None)[0]; Ah=np.exp(cc[0]); muh=cc[1]
        f=lambda z: Ah*Nhold**(-beta)*z**muh*np.exp(gc*z)-1e-3
        try: pred=brentq(f,max(h,1e-4),max(Zmax*2,10))
        except Exception: pred=np.nan
        ghold=df[df.N==Nhold].sort_values('Z'); obs=first_sustained_crossing(ghold.Z.to_numpy(),ghold.open_field_error.to_numpy(),1e-3,2)
        all_lono.append(dict(alpha=a,rho=rho,N_holdout=int(Nhold),Z_observed=obs,Z_predicted=pred,relative_error=(pred-obs)/obs if np.isfinite(obs) and obs else np.nan,mu_fit=muh,A_fit=Ah))

D=pd.concat(all_diag,ignore_index=True); F=pd.DataFrame(all_fit); H=pd.DataFrame(all_horiz); L=pd.DataFrame(all_lono)
D.to_csv(os.path.join(ROOT,'data/dynamic/dynamic_diagnostics.csv'),index=False)
F.to_csv(os.path.join(ROOT,'data/dynamic/dynamic_fits.csv'),index=False)
H.to_csv(os.path.join(ROOT,'data/dynamic/fidelity_horizons.csv'),index=False)
L.to_csv(os.path.join(ROOT,'data/dynamic/heldout_predictions.csv'),index=False)
if map_rows:
    pd.DataFrame(map_rows,columns=['Z','n','reference_intensity','open_intensity','corrected_intensity']).to_csv(os.path.join(ROOT,'data/dynamic/propagation_map_a080_r030_N256.csv'),index=False)
print('\nDYNAMIC FITS\n',F.to_string(index=False))
print('\nLONO SUMMARY median abs rel err=',np.nanmedian(np.abs(L.relative_error)))
