import os,sys
import numpy as np, pandas as pd
ROOT=os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0,os.path.join(ROOT,'src'))
from fidelity_model import *

a=.8;rho=.30;N=256;P=64;m=1;qc=2*np.pi/P
# Operator identity and FFT-vs-dense
rows=[]
per=PeriodicExactOperator(a,P)
for mm in range(0,P//2+1):
 q=2*np.pi*mm/P
 rows.append(dict(test='periodization_identity',index=mm,value=per.lam[mm],reference=lambda_exact(a,q),abs_error=abs(per.lam[mm]-lambda_exact(a,q))))
# FFT open op vs dense direct matrix at N=64
Nsmall=64; op=OpenOperator(a,Nsmall); rng=np.random.default_rng(3); v=rng.normal(size=Nsmall)+1j*rng.normal(size=Nsmall)
Lf=op.linear_open(v); J=kernel(a,Nsmall-1); Ld=np.zeros(Nsmall,dtype=complex)
for n in range(Nsmall):
 for k in range(Nsmall):
  if k!=n: Ld[n]+=J[abs(k-n)-1]*(v[k]-v[n])
for i,(x,y) in enumerate(zip(Lf,Ld)):
 rows.append(dict(test='fft_vs_direct',index=i,value=abs(x),reference=abs(y),abs_error=abs(x-y)))
pd.DataFrame(rows).to_csv(os.path.join(ROOT,'data/validation/operator_validation.csv'),index=False)

# Reference periodic RK4 convergence at Z=20.
def integrate_ref(h,Z=20):
 rr=PeriodicExactOperator(a,P); w=1+1e-6*np.cos(qc*np.arange(P)); steps=round(Z/h)
 for _ in range(steps): w=rk4_step(w,h,rr.rhs,rho)
 return w
hs=[0.08,0.04,0.02,0.01,0.005]; sol={h:integrate_ref(h) for h in hs}; ref=sol[0.005]
conv=[]
for h in hs[:-1]:
 e=float(np.sqrt(np.mean(abs(sol[h]-ref)**2)))
 conv.append(dict(h=h,error_vs_h005=e))
# observed order successive errors relative to finest proxy
for i in range(len(conv)-1):
 conv[i]['observed_order']=np.log(conv[i]['error_vs_h005']/conv[i+1]['error_vs_h005'])/np.log(2)
conv[-1]['observed_order']=np.nan
pd.DataFrame(conv).to_csv(os.path.join(ROOT,'data/validation/time_step_convergence.csv'),index=False)

# Open-domain key outputs under h refinement.
def pair_rhs(op):
 def rhs(x,rr):
  Lo=op.linear_open(x); L=Lo.copy(); L[1]+=op.kout*(1-x[1]); return 1j*(L+rr*(abs(x)**2-1)*x)
 return rhs
outs=[]
for h in [0.04,0.02,0.01]:
 per=PeriodicExactOperator(a,P); oo=OpenOperator(a,N); rhs=pair_rhs(oo)
 wr=1+1e-6*np.cos(qc*np.arange(P)); st=np.vstack([1+1e-6*np.cos(qc*np.arange(N))]*2).astype(complex); sl=central_slice(N,.25)
 steps=round(40/h); save=max(1,round(.1/h)); tlist=[]; eolist=[]
 for k in range(steps+1):
  if k%save==0:
   z=k*h; refn=periodic_repeat(wr,N); tlist.append(z); eolist.append(float(np.sqrt(np.mean(abs(st[0,sl]-refn[sl])**2))))
  if k<steps:
   wr=rk4_step(wr,h,per.rhs,rho); st=rk4_step(st,h,rhs,rho)
 t=np.array(tlist); e=np.array(eolist)
 for zq in [20,30,35,40]:
  j=int(np.argmin(abs(t-zq))); outs.append(dict(h=h,metric=f'field_error_Z{zq}',value=e[j]))
 for th in [1e-5,1e-4,1e-3]: outs.append(dict(h=h,metric=f'horizon_{th:g}',value=first_sustained_crossing(t,e,th,2)))
pd.DataFrame(outs).to_csv(os.path.join(ROOT,'data/validation/dynamic_step_sensitivity.csv'),index=False)
print('max periodic identity error',pd.DataFrame(rows).query("test=='periodization_identity'").abs_error.max())
print('max FFT-vs-direct error',pd.DataFrame(rows).query("test=='fft_vs_direct'").abs_error.max())
print(pd.DataFrame(conv))
print(pd.DataFrame(outs))
