import os,sys
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
ROOT=os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0,os.path.join(ROOT,'src'))
from fidelity_model import kernel
FOUT=os.path.join(ROOT,'figures'); DOUT=os.path.join(ROOT,'figure_data')
os.makedirs(FOUT,exist_ok=True); os.makedirs(DOUT,exist_ok=True)

# Journal-oriented defaults: title-free scientific panels, restrained grid, and
# line/marker redundancy so multi-curve plots remain readable in grayscale.
plt.rcParams.update({
    'font.size':9,'axes.labelsize':9,'legend.fontsize':7.5,'figure.dpi':160,
    'savefig.bbox':'tight','axes.linewidth':0.8,'lines.linewidth':1.25,
    'xtick.direction':'out','ytick.direction':'out'
})
COLORS=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#000000']
LINESTYLES=['-','--','-.',':',(0,(5,2)),(0,(3,1,1,1))]
MARKERS=['o','s','^','D','v','P']

def style(i):
    return dict(color=COLORS[i%len(COLORS)],ls=LINESTYLES[i%len(LINESTYLES)],marker=MARKERS[i%len(MARKERS)],
                mfc='white',mec=COLORS[i%len(COLORS)],ms=4.2,markevery=None)

def panel_labels(axs, labels=None):
    flat=np.asarray(axs).ravel()
    if labels is None: labels=[f'({chr(97+i)})' for i in range(len(flat))]
    for ax,lab in zip(flat,labels):
        ax.text(0.01,1.02,lab,transform=ax.transAxes,ha='left',va='bottom',fontweight='bold',fontsize=9,clip_on=False)

def clean_grid(ax, which='major'):
    ax.grid(True,which=which,alpha=.18,lw=.55)

def save(fig,num,df):
    df.to_csv(os.path.join(DOUT,f'figure{num:02d}.csv'),index=False)
    # Store the public release figures as vector PDF only.
    # Raster PNG duplicates are intentionally omitted to keep the repository compact.
    fig.savefig(os.path.join(FOUT,f'figure{num:02d}.pdf'))
    plt.close(fig)

# ---------- Figure 1: model and kernel ----------
rows=[]
fig,axs=plt.subplots(1,2,figsize=(7.1,2.7))
ax=axs[0]; ax.set_axis_off(); xs=np.linspace(.08,.92,9)
for x in xs: ax.add_patch(Circle((x,.58),.026,fill=False,lw=1.2))
for i in range(len(xs)-1): ax.add_patch(FancyArrowPatch((xs[i]+.03,.58),(xs[i+1]-.03,.58),arrowstyle='-',mutation_scale=8,lw=.9))
for j,rad in [(2,.16),(3,.24),(4,.32)]: ax.add_patch(FancyArrowPatch((xs[4],.61),(xs[4+j],.61),connectionstyle=f'arc3,rad={-rad}',arrowstyle='-',lw=.9,alpha=.8))
ax.add_patch(Rectangle((xs[2]-.04,.47),xs[6]-xs[2]+.08,.22,fill=False,ls='--',lw=1.0))
ax.text(.5,.38,'open truncation',ha='center')
ax.text(.5,.20,'tail correction approximates the omitted exterior\nby the homogeneous background',ha='center')
ax=axs[1]
r=np.arange(1,1001)
for i,a in enumerate([.8,1.2,1.6]):
    J=kernel(a,1000); kw=style(i); ax.loglog(r,J,label=fr'$\alpha={a}$',**kw)
    rows.extend([dict(panel='kernel',alpha=a,r=int(rr),J=float(jj)) for rr,jj in zip(r,J)])
ax.set_xlabel('separation $r$'); ax.set_ylabel('$J_r$'); ax.legend(frameon=False); clean_grid(ax,'both')
panel_labels(axs)
fig.tight_layout(); save(fig,1,pd.DataFrame(rows))

# ---------- Figure 2: dispersion and MI ----------
disp=pd.read_csv(os.path.join(ROOT,'data/theory/dispersion_mi.csv')); modes=pd.read_csv(os.path.join(ROOT,'data/theory/selected_modes.csv')); dynfit=pd.read_csv(os.path.join(ROOT,'data/dynamic/dynamic_fits.csv'))
fig,axs=plt.subplots(2,2,figsize=(7.1,5.1)); rows=[]
for i,a in enumerate([.8,1.2,1.6]):
    d=disp[(disp.alpha==a)&(disp.rho==.30)]
    kw=style(i); kw.pop('marker'); kw.pop('mfc'); kw.pop('mec'); kw.pop('ms'); kw.pop('markevery')
    axs[0,0].plot(d.q,d.lambda_inf,label=fr'$\alpha={a}$',**kw); rows.extend(d.assign(panel='dispersion').to_dict('records'))
axs[0,0].axhline(.30,ls='--',lw=.9,color='0.35'); axs[0,0].set(xlabel='$q$',ylabel=r'$\Lambda_\alpha(q)$'); axs[0,0].legend(frameon=False); clean_grid(axs[0,0])
for i,a in enumerate([.8,1.2,1.6]):
    d=disp[(disp.alpha==a)&(disp.rho==.30)]
    kw=style(i); kw.pop('marker'); kw.pop('mfc'); kw.pop('mec'); kw.pop('ms'); kw.pop('markevery')
    axs[0,1].plot(d.q,d.mi_gain,label=fr'$\alpha={a}$',**kw); rows.extend(d.assign(panel='mi_gain').to_dict('records'))
axs[0,1].set(xlabel='$q$',ylabel='$g(q)$'); axs[0,1].legend(frameon=False); clean_grid(axs[0,1])
mm=modes.copy(); axs[1,0].scatter(mm.q_star,mm.q_c,facecolors='white',edgecolors=COLORS[0],s=28)
lo=min(mm.q_star.min(),mm.q_c.min()); hi=max(mm.q_star.max(),mm.q_c.max()); axs[1,0].plot([lo,hi],[lo,hi],ls='--',lw=.9,color='0.35'); axs[1,0].set(xlabel='$q_*$',ylabel='$q_c$'); clean_grid(axs[1,0]); rows.extend(mm.assign(panel='selected_modes').to_dict('records'))
merged=dynfit.merge(modes[['alpha','rho','g_c']],on=['alpha','rho']); axs[1,1].scatter(merged.g_theory,merged.g_num,facecolors='white',edgecolors=COLORS[1],s=28)
lo=merged.g_theory.min()*.95; hi=merged.g_theory.max()*1.05; axs[1,1].plot([lo,hi],[lo,hi],ls='--',lw=.9,color='0.35'); axs[1,1].set(xlabel='theory $g(q_c)$',ylabel='numerical $g$'); clean_grid(axs[1,1]); rows.extend(merged.assign(panel='mi_validation').to_dict('records'))
panel_labels(axs)
fig.tight_layout(); save(fig,2,pd.DataFrame(rows))

# ---------- Figure 3: static closure scaling ----------
st=pd.read_csv(os.path.join(ROOT,'data/static/static_closure_scaling.csv')); loc=pd.read_csv(os.path.join(ROOT,'data/static/static_local_slopes.csv')); sf=pd.read_csv(os.path.join(ROOT,'data/static/static_fits.csv'))
fig,axs=plt.subplots(2,2,figsize=(7.1,5.1)); rows=[]
d=st[(st.alpha==.8)&(st.rho==.30)]
for i,(col,lab) in enumerate([('open_error','open'),('corrected_error','corrected'),('naive_periodic_error','naive periodic')]):
    axs[0,0].loglog(d.N,d[col],label=lab,**style(i))
axs[0,0].set(xlabel='$N$',ylabel='central operator defect'); axs[0,0].legend(frameon=False); clean_grid(axs[0,0],'both'); rows.extend(d.assign(panel='closure_scaling').to_dict('records'))
for i,a in enumerate([.8,1.2,1.6]):
    q=loc[(loc.alpha==a)&(loc.rho==.30)&(loc.closure=='open')]
    axs[0,1].semilogx(q.N,q.beta_local,label=fr'$\alpha={a}$',**style(i)); rows.extend(q.assign(panel='open_local_beta').to_dict('records'))
axs[0,1].set(xlabel='$N$',ylabel='local exponent'); axs[0,1].set_ylim(0.68,1.73); clean_grid(axs[0,1])
# Q1 layout fix: reserve a small in-axis band above the highest curve for the legend.
axs[0,1].legend(frameon=False,loc='upper center',ncol=3,fontsize=7,borderaxespad=.35,columnspacing=.8,handlelength=1.8)
openfit=sf[sf.closure=='open'].groupby('alpha',as_index=False).agg(beta=('beta','mean'),beta_min=('beta','min'),beta_max=('beta','max'))
open_yerr=np.vstack([openfit.beta-openfit.beta_min,openfit.beta_max-openfit.beta])
axs[1,0].errorbar(openfit.alpha,openfit.beta,yerr=open_yerr,fmt='o',mfc='white',mec=COLORS[0],ecolor=COLORS[0],capsize=2.5); axs[1,0].plot([.75,1.65],[.75,1.65],ls='--',lw=.9,color='0.35'); axs[1,0].set(xlabel=r'$\alpha$',ylabel=r'$\beta_{\rm open}$'); clean_grid(axs[1,0]); rows.extend(openfit.assign(panel='beta_vs_alpha').to_dict('records'))
last=loc[loc.closure=='corrected'].sort_values('N').groupby(['alpha','rho'],as_index=False).tail(1)
for i,(rho,g) in enumerate(last.groupby('rho')): axs[1,1].plot(g.alpha,g.beta_local,label=fr'$\rho={rho:.2f}$',**style(i))
a=np.linspace(.75,1.65,50); axs[1,1].plot(a,1+a,ls='--',lw=.9,color='0.35',label=r'$1+\alpha$')
axs[1,1].set(xlabel=r'$\alpha$',ylabel='largest-$N$ local exponent'); axs[1,1].legend(frameon=False); clean_grid(axs[1,1]); rows.extend(last.assign(panel='corrected_local_beta').to_dict('records'))
panel_labels(axs)
fig.tight_layout(); save(fig,3,pd.DataFrame(rows))

# ---------- Figure 4: dynamic amplification ----------
D=pd.read_csv(os.path.join(ROOT,'data/dynamic/dynamic_diagnostics.csv')); G=pd.read_csv(os.path.join(ROOT,'data/dynamic/error_growth_rates.csv'))
fig,axs=plt.subplots(2,2,figsize=(7.1,5.1)); rows=[]
d=D[(D.alpha==.8)&(D.rho==.30)]
for i,(N,g) in enumerate(d.groupby('N')):
    kw=style(i); kw['markevery']=max(1,len(g)//14)
    axs[0,0].semilogy(g.Z,g.open_field_error,label=f'N={N}',**kw); rows.extend(g.assign(panel='open_errors').to_dict('records'))
axs[0,0].set(xlabel='$Z$',ylabel=r'$E_{\rm field}$'); axs[0,0].legend(frameon=False,ncol=2); clean_grid(axs[0,0])
g=d[d.N==256]
for i,(col,lab) in enumerate([('open_field_error','open'),('corrected_field_error','corrected')]):
    kw=style(i); kw['markevery']=max(1,len(g)//14); axs[0,1].semilogy(g.Z,g[col],label=lab,**kw)
axs[0,1].set(xlabel='$Z$',ylabel=r'$E_{\rm field}$'); axs[0,1].legend(frameon=False); clean_grid(axs[0,1]); rows.extend(g.assign(panel='open_vs_corrected').to_dict('records'))
axs[1,0].scatter(G.g_theory,G.g_error,facecolors='white',edgecolors=COLORS[0],s=26)
lo=np.nanmin(G[['g_theory','g_error']].to_numpy()); hi=np.nanmax(G[['g_theory','g_error']].to_numpy()); axs[1,0].plot([lo,hi],[lo,hi],ls='--',lw=.9,color='0.35'); axs[1,0].set(xlabel=r'$g_{\rm MI}$',ylabel=r'$g_{\rm error}$'); clean_grid(axs[1,0]); rows.extend(G.assign(panel='g_error').to_dict('records'))
beta=float(sf[(sf.alpha==.8)&(sf.rho==.30)&(sf.closure=='open')].beta.iloc[0])
for i,z in enumerate([20,25,30]):
    q=d[np.isclose(d.Z,z)]
    axs[1,1].plot(q.N,q.open_field_error*q.N**beta,label=f'$Z={z:g}$',**style(i)); rows.extend(q.assign(panel='N_rescaling',scaled_error=q.open_field_error*q.N**beta).to_dict('records'))
axs[1,1].set_xscale('log',base=2); axs[1,1].set_xticks([128,256,512,1024],labels=['128','256','512','1024']); axs[1,1].minorticks_off(); axs[1,1].set(xlabel='$N$',ylabel=r'$E_{\rm field}N^{\beta}$'); axs[1,1].legend(frameon=False); clean_grid(axs[1,1])
panel_labels(axs)
fig.tight_layout(); save(fig,4,pd.DataFrame(rows))

# ---------- Figure 5: fidelity horizons ----------
H=pd.read_csv(os.path.join(ROOT,'data/dynamic/fidelity_horizons_all_thresholds.csv')); HF=pd.read_csv(os.path.join(ROOT,'data/dynamic/horizon_scaling_fits.csv'))
fig,axs=plt.subplots(2,2,figsize=(7.1,5.25)); rows=[]
h=H[(H.closure=='open')&(H.threshold==1e-5)]
for i,((a,rho),g) in enumerate(h.groupby(['alpha','rho'])):
    g=g.dropna(); kw=style(i); axs[0,0].plot(np.log(g.N),g.Z_f,label=fr'$\alpha={a},\rho={rho}$',**kw); rows.extend(g.assign(panel='horizons').to_dict('records'))
axs[0,0].set(xlabel=r'$\ln N$',ylabel='$Z_f$'); axs[0,0].set_ylim(15,90); clean_grid(axs[0,0])
# Q1 layout fix: expand the ordinate to create a clean legend band inside the panel.
axs[0,0].legend(frameon=False,loc='upper left',ncol=2,fontsize=6.4,borderaxespad=.45,columnspacing=.65,handlelength=1.6)
hf=HF[HF.threshold==1e-5]; axs[0,1].scatter(hf.predicted_slope_beta_over_g,hf.slope_Z_vs_lnN,facecolors='white',edgecolors=COLORS[0],s=28)
lo=min(hf.predicted_slope_beta_over_g.min(),hf.slope_Z_vs_lnN.min()); hi=max(hf.predicted_slope_beta_over_g.max(),hf.slope_Z_vs_lnN.max()); axs[0,1].plot([lo,hi],[lo,hi],ls='--',lw=.9,color='0.35'); axs[0,1].set(xlabel=r'predicted slope $\beta/g$',ylabel='measured slope'); clean_grid(axs[0,1]); rows.extend(hf.assign(panel='slope_test').to_dict('records'))
axs[1,0].scatter(np.arange(len(hf)),hf.slope_ratio,facecolors='white',edgecolors=COLORS[1],s=28); axs[1,0].axhline(1,ls='--',lw=.9,color='0.35'); labels=[f'{a:.1f}/{r:.2f}' for a,r in zip(hf.alpha,hf.rho)]; axs[1,0].set_xticks(np.arange(len(labels)),labels,rotation=45,ha='right'); axs[1,0].set(ylabel='measured / predicted slope'); clean_grid(axs[1,0]); rows.extend(hf.assign(panel='slope_ratio').to_dict('records'))
for i,(a,rho) in enumerate([(.8,.30),(1.2,.15)]):
    q=HF[(HF.alpha==a)&(HF.rho==rho)]
    axs[1,1].semilogx(q.threshold,q.slope_ratio,label=fr'$\alpha={a},\rho={rho}$',**style(i)); rows.extend(q.assign(panel='threshold_sensitivity').to_dict('records'))
axs[1,1].axhline(1,ls='--',lw=.9,color='0.35'); axs[1,1].set(xlabel='fidelity threshold',ylabel='slope ratio'); axs[1,1].legend(frameon=False); clean_grid(axs[1,1])
panel_labels(axs)
fig.tight_layout(); save(fig,5,pd.DataFrame(rows))

# ---------- Figure 6: propagation maps ----------
M=pd.read_csv(os.path.join(ROOT,'data/dynamic/propagation_map_a080_r030_N256.csv')); rows=M.assign(panel='propagation_map')
zs=np.sort(M.Z.unique()); ns=np.sort(M.n.unique())
def mat(col): return M.pivot(index='Z',columns='n',values=col).reindex(index=zs,columns=ns).to_numpy()
fig,axs=plt.subplots(2,2,figsize=(7.1,5.2),sharex=True,sharey=True)
for ax,col in zip(axs.flat[:3],['reference_intensity','open_intensity','corrected_intensity']):
    im=ax.imshow(mat(col),aspect='auto',origin='lower',extent=[ns.min(),ns.max(),zs.min(),zs.max()]); fig.colorbar(im,ax=ax,shrink=.8,label='$|w_n|^2$')
err=np.abs(mat('open_intensity')-mat('reference_intensity')); im=axs[1,1].imshow(err,aspect='auto',origin='lower',extent=[ns.min(),ns.max(),zs.min(),zs.max()]); fig.colorbar(im,ax=axs[1,1],shrink=.8,label=r'$|I_{\rm open}-I_{\rm ref}|$')
for ax in axs[:,0]: ax.set_ylabel('$Z$')
for ax in axs[1,:]: ax.set_xlabel('site $n$')
panel_labels(axs)
fig.tight_layout(); save(fig,6,rows)

# ---------- Figure 7: stochastic pilot ----------
S=pd.read_csv(os.path.join(ROOT,'data/stochastic/stochastic_pilot_summary.csv')); P=pd.read_csv(os.path.join(ROOT,'data/stochastic/stochastic_pilot_pdfs.csv')); PM=pd.read_csv(os.path.join(ROOT,'data/stochastic/stochastic_paired_maxima.csv'))
fig,axs=plt.subplots(2,2,figsize=(7.1,5.1)); rows=[]
for i,((N,cl),g) in enumerate(S.groupby(['N','closure'])):
    g=g.sort_values('Z'); kw=style(i); axs[0,0].plot(g.Z,g.D_W,label=f'{cl}, $N={N}$',**kw)
    if {'D_W_ci_low','D_W_ci_high'}.issubset(g.columns):
        axs[0,0].fill_between(g.Z,g.D_W_ci_low,g.D_W_ci_high,color=kw['color'],alpha=.08,linewidth=0)
    rows.extend(g.assign(panel='wasserstein').to_dict('records'))
axs[0,0].set(xlabel='$Z$',ylabel='$D_W$'); axs[0,0].legend(frameon=False,fontsize=5.9,ncol=2); clean_grid(axs[0,0])
for i,((N,cl),g) in enumerate(S.groupby(['N','closure'])):
    g=g.sort_values('Z'); kw=style(i); axs[0,1].plot(g.Z,g.Delta99,label=f'{cl}, $N={N}$',**kw)
    if {'Delta99_ci_low','Delta99_ci_high'}.issubset(g.columns):
        axs[0,1].fill_between(g.Z,g.Delta99_ci_low,g.Delta99_ci_high,color=kw['color'],alpha=.08,linewidth=0)
axs[0,1].axhline(0,ls='--',lw=.8,color='0.35'); axs[0,1].set(xlabel='$Z$',ylabel=r'$\Delta_{99}$'); clean_grid(axs[0,1]); rows.extend(S.assign(panel='delta99').to_dict('records'))
for i,cl in enumerate(['open','corrected']):
    q=P[(P.N==128)&(P.Z==30)&(P.closure==cl)]
    kw=style(i); kw.pop('marker'); kw.pop('mfc'); kw.pop('mec'); kw.pop('ms'); kw.pop('markevery'); axs[1,0].plot(q.intensity,q.pdf_finite,label=cl,**kw); rows.extend(q.assign(panel='pdf').to_dict('records'))
q=P[(P.N==128)&(P.Z==30)&(P.closure=='open')]; axs[1,0].plot(q.intensity,q.pdf_ref,ls='--',color='0.25',label='reference'); axs[1,0].set(xlabel='intensity $I$',ylabel='PDF'); axs[1,0].legend(frameon=False); clean_grid(axs[1,0])
q=PM[(PM.N==128)&(PM.Z==30)]
for i,(cl,g) in enumerate(q.groupby('closure')):
    vals=np.sort(g.delta_max); kw=style(i); kw.pop('marker'); kw.pop('mfc'); kw.pop('mec'); kw.pop('ms'); kw.pop('markevery'); axs[1,1].plot(vals,np.linspace(0,1,len(vals),endpoint=False),label=cl,**kw); rows.extend(g.assign(panel='paired_maxima').to_dict('records'))
axs[1,1].axvline(0,ls='--',lw=.8,color='0.35'); axs[1,1].set(xlabel=r'$B_N-B_{\rm ref}$',ylabel='empirical CDF'); axs[1,1].legend(frameon=False); clean_grid(axs[1,1])
panel_labels(axs)
fig.tight_layout(); save(fig,7,pd.DataFrame(rows))

# ---------- Figure 8: synthesis ----------
fig,axs=plt.subplots(2,2,figsize=(7.1,5.1)); rows=[]
axs[0,0].errorbar(openfit.alpha,openfit.beta,yerr=open_yerr,fmt='o',mfc='white',mec=COLORS[0],ecolor=COLORS[0],capsize=2.5); axs[0,0].plot([.75,1.65],[.75,1.65],ls='--',lw=.9,color='0.35'); axs[0,0].set(xlabel=r'$\alpha$',ylabel=r'$\beta_{\rm open}$'); clean_grid(axs[0,0]); rows.extend(openfit.assign(panel='spatial').to_dict('records'))
gsum=G.groupby(['alpha','rho'],as_index=False).agg(g_ratio=('g_ratio','median'))
axs[0,1].scatter(np.arange(len(gsum)),gsum.g_ratio,facecolors='white',edgecolors=COLORS[1],s=28); axs[0,1].axhline(1,ls='--',lw=.9,color='0.35'); labs=[f'{a:.1f}/{r:.2f}' for a,r in zip(gsum.alpha,gsum.rho)]; axs[0,1].set_xticks(np.arange(len(labs)),labs,rotation=45,ha='right'); axs[0,1].set(ylabel=r'$g_{\rm error}/g_{\rm MI}$'); clean_grid(axs[0,1]); rows.extend(gsum.assign(panel='temporal').to_dict('records'))
lo8=min(hf.predicted_slope_beta_over_g.min(),hf.slope_Z_vs_lnN.min()); hi8=max(hf.predicted_slope_beta_over_g.max(),hf.slope_Z_vs_lnN.max())
axs[1,0].scatter(hf.predicted_slope_beta_over_g,hf.slope_Z_vs_lnN,facecolors='white',edgecolors=COLORS[0],s=28); axs[1,0].plot([lo8,hi8],[lo8,hi8],ls='--',lw=.9,color='0.35'); axs[1,0].set(xlabel=r'$\beta/g$',ylabel=r'measured $\mathrm{d}Z_f/\mathrm{d}\ln N$'); clean_grid(axs[1,0]); rows.extend(hf.assign(panel='horizon').to_dict('records'))
q=S[S.Z==25].pivot_table(index='N',columns='closure',values='D_W').reset_index(); q['reduction']=1-q['corrected']/q['open']; axs[1,1].plot(q.N,q.reduction,**style(0)); axs[1,1].set_xscale('log',base=2); axs[1,1].set(xlabel='$N$',ylabel=r'$1-D_W^{\rm corr}/D_W^{\rm open}$'); clean_grid(axs[1,1]); rows.extend(q.assign(panel='correction_efficiency').to_dict('records'))
panel_labels(axs)
fig.tight_layout(); save(fig,8,pd.DataFrame(rows))
print('Created title-free, panel-labeled figures 1-8')
