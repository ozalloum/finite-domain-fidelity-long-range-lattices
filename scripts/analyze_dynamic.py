import os,sys
import numpy as np
import pandas as pd
from scipy.stats import linregress
ROOT=os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0,os.path.join(ROOT,'src'))
from fidelity_model import first_sustained_crossing
D=pd.read_csv(os.path.join(ROOT,'data/dynamic/dynamic_diagnostics.csv'))
F=pd.read_csv(os.path.join(ROOT,'data/dynamic/dynamic_fits.csv'))
rows_g=[]; rows_h=[]; rows_hfit=[]
for (a,rho),g in D.groupby(['alpha','rho']):
    gt=float(F[(F.alpha==a)&(F.rho==rho)].g_theory.iloc[0]); beta=float(F[(F.alpha==a)&(F.rho==rho)].beta_static.iloc[0])
    for N,h in g.groupby('N'):
        mask=(h.reference_mode_amplitude>=2e-6)&(h.reference_mode_amplitude<=2e-4)&(h.open_field_error>1e-8)&(h.open_field_error<3e-5)
        if mask.sum()>=10:
            lr=linregress(h.loc[mask,'Z'],np.log(h.loc[mask,'open_field_error']))
            ge=lr.slope; r2=lr.rvalue**2
        else: ge=np.nan; r2=np.nan
        rows_g.append(dict(alpha=a,rho=rho,N=int(N),g_theory=gt,g_error=ge,g_ratio=ge/gt if np.isfinite(ge) else np.nan,r2=r2,n_points=int(mask.sum())))
        for thr in [1e-5,3e-5,1e-4,3e-4,1e-3,1e-2]:
            for closure,col in [('open','open_field_error'),('corrected','corrected_field_error')]:
                z=first_sustained_crossing(h.Z.to_numpy(),h[col].to_numpy(),thr,2)
                rows_h.append(dict(alpha=a,rho=rho,N=int(N),closure=closure,threshold=thr,Z_f=z,g=gt,beta=beta))
    for thr in [1e-5,3e-5,1e-4,3e-4,1e-3]:
        hh=pd.DataFrame([r for r in rows_h if r['alpha']==a and r['rho']==rho and r['closure']=='open' and r['threshold']==thr]).dropna()
        if len(hh)>=3:
            lr=linregress(np.log(hh.N),hh.Z_f)
            pred=beta/gt
            rows_hfit.append(dict(alpha=a,rho=rho,threshold=thr,n_points=len(hh),slope_Z_vs_lnN=lr.slope,intercept=lr.intercept,r2=lr.rvalue**2,predicted_slope_beta_over_g=pred,slope_ratio=lr.slope/pred))
G=pd.DataFrame(rows_g); H=pd.DataFrame(rows_h); HF=pd.DataFrame(rows_hfit)
G.to_csv(os.path.join(ROOT,'data/dynamic/error_growth_rates.csv'),index=False)
H.to_csv(os.path.join(ROOT,'data/dynamic/fidelity_horizons_all_thresholds.csv'),index=False)
HF.to_csv(os.path.join(ROOT,'data/dynamic/horizon_scaling_fits.csv'),index=False)
print(G.groupby(['alpha','rho']).g_ratio.agg(['median','min','max']).to_string())
print('\nEarly-horizon fits (1e-5):')
print(HF[HF.threshold==1e-5].to_string(index=False))
