import os,sys,time
import numpy as np, pandas as pd
from scipy.stats import wasserstein_distance
ROOT=os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0,os.path.join(ROOT,'src'))
from fidelity_model import *

a=.8; rho=.30; sigma=1e-3; R=64; seed=20260807; Nref=2048; Ns=[128,256,512]
h=.04; checkpoints=[15.,20.,25.,30.]; Zmax=max(checkpoints)
rng=np.random.default_rng(seed)
eta=(rng.standard_normal((R,Nref))+1j*rng.standard_normal((R,Nref)))/np.sqrt(2)
eta=eta-eta.mean(axis=1,keepdims=True)
master=1+sigma*eta
# Normalize each master realization to unit mean intensity, preserving paired restrictions.
master=master/np.sqrt(np.mean(np.abs(master)**2,axis=1,keepdims=True))
ref=master.copy(); per=PeriodicExactOperator(a,Nref)
ops={N:OpenOperator(a,N) for N in Ns}
states={}
for N in Ns:
 start=(Nref-N)//2; base=master[:,start:start+N]
 states[N]=np.stack([base.copy(),base.copy()],axis=0) # [closure,R,N]

def step_pair(x,op):
 def rhs(y,rr):
  # y shape 2,R,N
  Lo=op.linear_open(y)
  L=Lo.copy(); L[1]+=op.kout*(1-y[1])
  return 1j*(L+rr*(np.abs(y)**2-1)*y)
 return rk4_step(x,h,rhs,rho)

snapshots={}
steps=round(Zmax/h); cpsteps={round(z/h):z for z in checkpoints}
t0=time.time()
for k in range(steps+1):
 if k in cpsteps:
  z=cpsteps[k]; snapshots[z]={'ref':ref.copy(),'finite':{N:states[N].copy() for N in Ns}}
  print('checkpoint',z,'elapsed',time.time()-t0,flush=True)
 if k<steps:
  ref=rk4_step(ref,h,per.rhs,rho)
  for N in Ns: states[N]=step_pair(states[N],ops[N])


def paired_bootstrap_metrics(Ir, If, B, rng, chunk=100):
    """Paired realization bootstrap for W1 and Q99 bias.

    Ir and If have equal shape (R, L). For equal-size empirical samples,
    the one-dimensional W1 distance is the mean absolute difference of
    the sorted samples. Chunking vectorizes the bootstrap while keeping
    peak memory bounded.
    """
    Rloc=Ir.shape[0]
    bd=np.empty(B,float); bq=np.empty(B,float)
    pos=0
    while pos<B:
        nb=min(chunk,B-pos)
        ix=rng.integers(0,Rloc,size=(nb,Rloc))
        fr=Ir[ix].reshape(nb,-1)
        ff=If[ix].reshape(nb,-1)
        sr=np.sort(fr,axis=1); sf=np.sort(ff,axis=1)
        bd[pos:pos+nb]=np.mean(np.abs(sf-sr),axis=1)/np.mean(fr,axis=1)
        qr=np.quantile(fr,.99,axis=1); qf=np.quantile(ff,.99,axis=1)
        bq[pos:pos+nb]=(qf-qr)/qr
        pos+=nb
    return bd,bq

# Statistics and bootstrap by realization.
rows=[]; pdfrows=[]; paired=[]
B=500
for z,snap in snapshots.items():
 refw=snap['ref']
 for N in Ns:
  start=(Nref-N)//2; # compare central N/4 region in both systems
  sl=central_slice(N,.25); rsl=slice(start+sl.start,start+sl.stop)
  Ir=np.abs(refw[:,rsl])**2
  for ci,closure in enumerate(['open','corrected']):
   If=np.abs(snap['finite'][N][ci,:,sl])**2
   flat_r=Ir.ravel(); flat_f=If.ravel()
   dw=wasserstein_distance(flat_f,flat_r)/flat_r.mean()
   q99r=np.quantile(flat_r,.99); q99f=np.quantile(flat_f,.99); d99=(q99f-q99r)/q99r
   m4r=np.mean(flat_r**4)/(np.mean(flat_r)**4); m4f=np.mean(flat_f**4)/(np.mean(flat_f)**4)
   # per-realization maxima at checkpoint
   mr=Ir.max(axis=1); mf=If.max(axis=1); db=mf-mr
   # paired bootstrap over whole realizations (vectorized in memory-bounded chunks)
   bd,bq=paired_bootstrap_metrics(Ir,If,B,rng)
   rows.append(dict(alpha=a,rho=rho,sigma=sigma,R=R,Nref=Nref,N=N,Z=z,closure=closure,D_W=dw,D_W_ci_low=np.quantile(bd,.025),D_W_ci_high=np.quantile(bd,.975),Q99_ref=q99r,Q99_finite=q99f,Delta99=d99,Delta99_ci_low=np.quantile(bq,.025),Delta99_ci_high=np.quantile(bq,.975),M4_ref=m4r,M4_finite=m4f,paired_max_mean=float(db.mean()),paired_max_positive_fraction=float(np.mean(db>0))))
   for rr,(x,y) in enumerate(zip(mf,mr)):
    paired.append(dict(N=N,Z=z,closure=closure,realization=rr,max_finite=x,max_ref=y,delta_max=x-y))
   # common histogram edges reference + both finite range; density written for plot
   lo=min(flat_r.min(),flat_f.min()); hi=max(flat_r.max(),flat_f.max()); edges=np.linspace(lo,hi,81); centers=.5*(edges[:-1]+edges[1:])
   hr,_=np.histogram(flat_r,bins=edges,density=True); hf,_=np.histogram(flat_f,bins=edges,density=True)
   for c,pr,pf in zip(centers,hr,hf): pdfrows.append(dict(N=N,Z=z,closure=closure,intensity=c,pdf_ref=pr,pdf_finite=pf))

pd.DataFrame(rows).to_csv(os.path.join(ROOT,'data/stochastic/stochastic_pilot_summary.csv'),index=False)
pd.DataFrame(pdfrows).to_csv(os.path.join(ROOT,'data/stochastic/stochastic_pilot_pdfs.csv'),index=False)
pd.DataFrame(paired).to_csv(os.path.join(ROOT,'data/stochastic/stochastic_paired_maxima.csv'),index=False)
# Save master seed/metadata, not huge snapshots.
pd.DataFrame([dict(alpha=a,rho=rho,sigma=sigma,R=R,seed=seed,Nref=Nref,h=h,checkpoints=';'.join(map(str,checkpoints))) ]).to_csv(os.path.join(ROOT,'data/stochastic/stochastic_pilot_config.csv'),index=False)
print(pd.DataFrame(rows).to_string(index=False))
