"""Fast independent reproducibility checks for Paper 1.

Writes logs/reproducibility_smoke_test.csv and exits nonzero if a required check fails.
This does not rerun the full parameter campaign.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from fidelity_model import tail, PeriodicExactOperator, OpenOperator, kernel, lambda_exact, static_defects, rk4_step, central_slice, periodic_repeat

rows=[]
def add(name,value,target,tolerance,mode='abs'):
    if mode=='abs':
        err=abs(value-target); passed=err<=tolerance
    elif mode=='rel':
        err=abs(value-target)/max(abs(target),1e-300); passed=err<=tolerance
    else:
        raise ValueError(mode)
    rows.append(dict(test=name,value=float(value),target=float(target),error=float(err),tolerance=float(tolerance),error_mode=mode,status='PASS' if passed else 'FAIL'))
    return passed

# 1) Kernel tail local exponent for alpha=0.8.
a=0.8
m1,m2=2048,4096
beta_tail=-np.log(float(tail(a,m2))/float(tail(a,m1)))/np.log(m2/m1)
add('kernel_tail_local_exponent_alpha_0.8',beta_tail,a,0.015,'abs')

# 2) Exact periodization identity for period 64.
P=64
per=PeriodicExactOperator(a,P)
errs=[]
for m in range(P//2+1):
    q=2*np.pi*m/P
    errs.append(abs(per.lam[m]-lambda_exact(a,q)))
add('exact_periodization_max_abs_error',max(errs),0.0,2e-13,'abs')

# 3) FFT open operator vs direct summation.
N=32
rng=np.random.default_rng(1907)
v=rng.normal(size=N)+1j*rng.normal(size=N)
op=OpenOperator(a,N)
Lf=op.linear_open(v)
J=kernel(a,N-1)
Ld=np.zeros(N,dtype=complex)
for n in range(N):
    for k in range(N):
        if k!=n:
            Ld[n]+=J[abs(k-n)-1]*(v[k]-v[n])
add('fft_open_vs_direct_max_abs_error',np.max(np.abs(Lf-Ld)),0.0,1e-12,'abs')

# 4) Static open scaling in discovery mode, independent from stored fit.
rho=0.30
m=1; qc=2*np.pi/64
Ns=np.array([512,1024,2048,4096])
errs=np.array([static_defects(a,int(N),qc)[0] for N in Ns])
beta_static=-np.polyfit(np.log(Ns),np.log(errs),1)[0]
add('static_open_beta_alpha_0.8',beta_static,a,0.025,'abs')

# 5) Short deterministic benchmark: reproduce frozen open error at Z=10 for N=128.
N=128; h=0.02; Z=10.0; eps=1e-6
per=PeriodicExactOperator(a,64); oo=OpenOperator(a,N)
wr=(1+eps*np.cos(qc*np.arange(64))).astype(complex)
wo=(1+eps*np.cos(qc*np.arange(N))).astype(complex)
for _ in range(round(Z/h)):
    wr=rk4_step(wr,h,per.rhs,rho)
    wo=rk4_step(wo,h,oo.rhs_open,rho)
ref=periodic_repeat(wr,N); sl=central_slice(N,.25)
err=float(np.sqrt(np.mean(np.abs(wo[sl]-ref[sl])**2)))
frozen=pd.read_csv(ROOT/'data/dynamic/dynamic_diagnostics.csv')
fr=frozen[(frozen.alpha==a)&(frozen.rho==rho)&(frozen.N==N)]
j=(fr.Z-Z).abs().idxmin(); target=float(fr.loc[j,'open_field_error'])
add('dynamic_open_error_a080_r030_N128_Z10',err,target,0.01,'rel')

out=pd.DataFrame(rows)
(ROOT/'logs').mkdir(exist_ok=True)
out.to_csv(ROOT/'logs/reproducibility_smoke_test.csv',index=False)
print(out.to_string(index=False))
if (out.status!='PASS').any():
    raise SystemExit(1)
print('\nAll reproducibility smoke tests PASS.')
