import os, numpy as np, pandas as pd
ROOT=os.path.dirname(os.path.dirname(__file__))
sf=pd.read_csv(os.path.join(ROOT,'data/static/static_fits.csv'))
loc=pd.read_csv(os.path.join(ROOT,'data/static/static_local_slopes.csv'))
g=pd.read_csv(os.path.join(ROOT,'data/dynamic/error_growth_rates.csv'))
hf=pd.read_csv(os.path.join(ROOT,'data/dynamic/horizon_scaling_fits.csv'))
st=pd.read_csv(os.path.join(ROOT,'data/stochastic/stochastic_pilot_summary.csv'))
val=pd.read_csv(os.path.join(ROOT,'data/validation/operator_validation.csv'))
step=pd.read_csv(os.path.join(ROOT,'data/validation/dynamic_step_sensitivity.csv'))
rows=[]
for a,gg in sf[sf.closure=='open'].groupby('alpha'):
 rows.append(dict(metric=f'beta_open_alpha_{a}',value=gg.beta.mean(),uncertainty=np.nan,range_low=gg.beta.min(),range_high=gg.beta.max(),notes='mean across rho=0.15,0.30; range_low/range_high are a descriptive two-case range, not a confidence interval or sampling uncertainty'))
last=loc[loc.closure=='corrected'].sort_values('N').groupby(['alpha','rho'],as_index=False).tail(1)
for _,r in last.iterrows(): rows.append(dict(metric=f'beta_corrected_local_alpha_{r.alpha}_rho_{r.rho}',value=r.beta_local,uncertainty=np.nan,range_low=np.nan,range_high=np.nan,notes=f'local exponent at N={int(r.N)}'))
rows.append(dict(metric='median_g_error_over_g_MI',value=np.nanmedian(g.g_ratio),uncertainty=np.nan,range_low=np.nan,range_high=np.nan,dispersion_measure=np.nanmedian(np.abs(g.g_ratio-1)),notes='median ratio across 24 N/case estimates; dispersion_measure is the median absolute deviation from unity and is descriptive, not inferential uncertainty'))
early=hf[hf.threshold==1e-5]
rows.append(dict(metric='median_horizon_slope_ratio',value=np.median(early.slope_ratio),uncertainty=np.nan,range_low=np.nan,range_high=np.nan,dispersion_measure=np.median(np.abs(early.slope_ratio-1)),notes='measured dZf/dlnN divided by beta/g across six cases; dispersion_measure is the median absolute deviation from unity and is descriptive'))
rows.append(dict(metric='median_horizon_slope_ratio_excluding_smallq_stress',value=np.median(early[~((early.alpha==.8)&(early.rho==.15))].slope_ratio),uncertainty=np.nan,range_low=np.nan,range_high=np.nan,dispersion_measure=np.median(np.abs(early[~((early.alpha==.8)&(early.rho==.15))].slope_ratio-1)),notes='excludes alpha=0.8,rho=0.15 small-q pre-asymptotic stress case; dispersion_measure is descriptive'))
for N,gg in st[st.Z==25].groupby('N'):
 p=gg.set_index('closure').D_W
 rows.append(dict(metric=f'pilot_DW_reduction_N{N}_Z25',value=1-p['corrected']/p['open'],uncertainty=np.nan,range_low=np.nan,range_high=np.nan,notes='R=64 paired stochastic pilot'))
rows.append(dict(metric='max_periodization_identity_error',value=val[val.test=='periodization_identity'].abs_error.max(),uncertainty=np.nan,range_low=np.nan,range_high=np.nan,notes='machine-precision identity check'))
rows.append(dict(metric='max_fft_vs_direct_operator_error',value=val[val.test=='fft_vs_direct'].abs_error.max(),uncertainty=np.nan,range_low=np.nan,range_high=np.nan,notes='N=64 random complex vector'))
pd.DataFrame(rows).to_csv(os.path.join(ROOT,'data/headline_results.csv'),index=False)
# Tables
params=pd.DataFrame([
 ['alpha','0.8, 1.2, 1.6','interaction-range exponent'],['rho','0.15, 0.30','dimensionless Kerr strength'],['N deterministic','128, 256, 512, 1024','dynamic domains'],['N static','64 to 32768','operator-scaling sweep'],['epsilon','1e-6','single-mode deterministic perturbation'],['h','0.02','production RK4 step'],['central window','N/4','error/statistical observation region'],['stochastic sigma','1e-3','pilot noise amplitude'],['stochastic R','64','paired pilot ensemble'],['stochastic Nref','2048','exact-periodized pilot reference']
],columns=['parameter','value','description'])
params.to_csv(os.path.join(ROOT,'data/table01_parameters.csv'),index=False)
# one row per case
cases=sf[sf.closure=='open'][['alpha','rho','beta']].rename(columns={'beta':'beta_static'})
gsum=g.groupby(['alpha','rho'],as_index=False).agg(g_error_ratio_median=('g_ratio','median'))
res=cases.merge(gsum,on=['alpha','rho']).merge(early[['alpha','rho','slope_Z_vs_lnN','predicted_slope_beta_over_g','slope_ratio','r2']],on=['alpha','rho'])
res.to_csv(os.path.join(ROOT,'data/table02_main_results.csv'),index=False)
print(pd.DataFrame(rows).to_string(index=False))
print(res.to_string(index=False))
