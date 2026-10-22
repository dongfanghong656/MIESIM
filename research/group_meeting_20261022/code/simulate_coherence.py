#!/usr/bin/env python3
"""Reproduce coherent-amplitude counterexample and ensemble-intensity identity.

Run: python simulate_coherence.py --out-root /path/to/group_meeting_project
Dependencies: numpy, scipy, matplotlib. No input/experimental data required.
All coordinates are dimensionless illustrative sample coordinates. These are
mathematical forward-model checks, not a calibrated OCT simulation.
"""
from pathlib import Path
import argparse
import json
import os
import shutil
import tempfile
import numpy as np
from scipy.signal import fftconvolve
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'oct-coherence-mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SEED=20261022
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,
 'axes.spines.top':False,'axes.spines.right':False,'axes.grid':False,
 'svg.fonttype':'none','figure.facecolor':'white','axes.facecolor':'white',
 'savefig.facecolor':'white','axes.titlesize':15,'axes.labelsize':12})
COLORS=['#126B85','#D26942','#7657A2','#4E8447']

def c_psf(x):
    return np.exp(-x*x/(2*1.05**2))*np.exp(1j*.55*x*x)

def save(fig,path):
    fig.savefig(path.with_suffix('.svg'),bbox_inches='tight')
    fig.savefig(path.with_suffix('.png'),dpi=220,bbox_inches='tight')
    plt.close(fig)

def main(root):
    figs=root/'figures'; evidence=root/'evidence'
    figs.mkdir(parents=True,exist_ok=True);evidence.mkdir(parents=True,exist_ok=True)
    x=np.linspace(-5,5,1201);positions=[-.75,.75]
    h1=c_psf(x-positions[0]);h2=c_psf(x-positions[1])
    surrogate=np.abs(h1)+np.abs(h2)
    phases=[0,np.pi/2,np.pi]
    curves=[np.abs(h1+np.exp(1j*p)*h2) for p in phases]
    fig,axs=plt.subplots(1,2,figsize=(13.8,4.8),gridspec_kw={'width_ratios':[1,2.35]})
    ax=axs[0]
    ax.vlines(positions,0,[1,1],colors=COLORS[0],linewidth=3)
    ax.scatter(positions,[1,1],s=65,color=COLORS[0],zorder=3)
    ax.set(xlim=(-3,3),ylim=(0,1.35),xlabel='Object position (a.u.)',ylabel='Scatterer power',title='Same two scatterer powers')
    ax.text(0,.38,'Only relative phase changes',ha='center',fontsize=11)
    ax=axs[1]
    ax.plot(x,surrogate,'--',color='#333333',lw=2.3,label=r'Linear surrogate $|h|*|\eta|$')
    for curve,p,c,label in zip(curves,phases,COLORS,[r'$\Delta\theta=0$',r'$\Delta\theta=\pi/2$',r'$\Delta\theta=\pi$']):
        ax.plot(x,curve,color=c,lw=2,label=label)
    ax.set(xlim=(-4,4),ylim=(0,2.1),xlabel='Image coordinate (a.u.)',ylabel='OCT amplitude (a.u.)',title=r'Coherent result $|h*\eta|$ depends on phase')
    ax.legend(loc='upper right',fontsize=10,frameon=False)
    fig.subplots_adjust(wspace=.31,bottom=.2,top=.83)
    fig.text(.06,.02,'Illustrative noise-free calculation. Identical complex PSF and object powers in every curve.',fontsize=11,color='#444444')
    save(fig,figs/'coherent_cross_terms')
    counter=[]
    for p,c in zip(phases,curves):
        counter.append({'relative_phase_rad':float(p),'relative_L2_difference_from_amplitude_convolution':float(np.linalg.norm(c-surrogate)/np.linalg.norm(surrogate)),
          'center_amplitude':float(c[len(x)//2]),'linear_surrogate_center':float(surrogate[len(x)//2])})

    rng=np.random.default_rng(SEED)
    n=513;x2=np.linspace(-8,8,n);dx=float(x2[1]-x2[0]);kx=np.arange(-128,129)*dx
    h=c_psf(kx);h=h/np.sqrt(np.sum(np.abs(h)**2))
    variance=.12+.7*np.exp(-((x2+2)/.65)**2)+np.exp(-((x2-2.2)/.9)**2)
    truth=fftconvolve(variance,np.abs(h)**2,mode='same')
    amp_truth=np.sqrt(np.pi)/2*np.sqrt(truth)
    looks=[1,4,16,64,256,1024,4096,16384]
    power_sum=np.zeros(n);amp_sum=np.zeros(n);snapshots={};errs=[]
    last=0
    for target in looks:
        for start in range(last,target,128):
            size=min(128,target-start)
            eta=np.sqrt(variance/2)[None,:]*(rng.normal(size=(size,n))+1j*rng.normal(size=(size,n)))
            g=fftconvolve(eta,h[None,:],mode='same',axes=-1)
            power_sum+=(np.abs(g)**2).sum(axis=0)
            amp_sum+=np.abs(g).sum(axis=0)
        avg=power_sum/target;avg_a=amp_sum/target
        snapshots[target]=avg.copy()
        errs.append({'independent_looks':target,'relative_L2_intensity_error':float(np.linalg.norm(avg-truth)/np.linalg.norm(truth)),
                     'relative_L2_mean_amplitude_vs_Rayleigh_mean_error':float(np.linalg.norm(avg_a-amp_truth)/np.linalg.norm(amp_truth))})
        last=target
    fig,axs=plt.subplots(1,2,figsize=(13.8,4.8),gridspec_kw={'width_ratios':[2.1,1]})
    ax=axs[0]
    ax.plot(x2,snapshots[1],color='#CCCCCC',lw=1,label='1 independent look')
    ax.plot(x2,snapshots[64],color=COLORS[1],lw=1.6,label='64-look mean')
    ax.plot(x2,snapshots[16384],color=COLORS[0],lw=2,label='16,384-look mean')
    ax.plot(x2,truth,'--',color='#222222',lw=2,label=r'Theory $|h|^2*v$')
    ax.set(xlim=(-7,7),xlabel='Image coordinate (a.u.)',ylabel='Intensity (a.u.)',title='Ensemble intensity approaches a nonnegative linear model')
    ax.legend(fontsize=10,frameon=False,loc='upper right')
    ax=axs[1]
    e=np.array([v['relative_L2_intensity_error'] for v in errs]);ls=np.array(looks)
    ax.loglog(ls,e,'o-',color=COLORS[0],lw=2,ms=5)
    ax.loglog(ls,e[0]/np.sqrt(ls),'--',color='#777777',lw=1.5,label=r'$1/\sqrt{L}$ guide')
    ax.legend(fontsize=10,frameon=False)
    ax.set(xlabel='Independent looks L',ylabel='Relative L2 error',title='Finite-look sampling error')
    ax.grid(alpha=.18,which='major')
    fig.subplots_adjust(wspace=.32,bottom=.22,top=.83)
    fig.text(.06,.055,'Assumptions: zero-mean proper complex Gaussian scatterers, diagonal covariance, fixed PSF, independent looks.',fontsize=10.3,color='#444444')
    fig.text(.06,.015,'Repeated measurements of the same static speckle do not meet the independent-look assumption.',fontsize=10.3,color='#444444')
    save(fig,figs/'ensemble_intensity')
    for extension in ('.svg', '.png'):
        shutil.copyfile((figs/'ensemble_intensity').with_suffix(extension), (figs/'multilook_power').with_suffix(extension))
    amp_surrogate=fftconvolve(np.sqrt(np.pi)/2*np.sqrt(variance),np.abs(h),mode='same')
    metrics={'kind':'illustrative_forward_model_validation_not_experimental_OCT',
      'seed':SEED,'complex_psf':'exp(-x^2/(2*1.05^2))*exp(1j*0.55*x^2)',
      'coordinates':'dimensionless arbitrary units',
      'counterexample':counter,'ensemble':{'grid_n':n,'kernel_n':len(h),'boundary':'linear zero-padded convolution, cropped same',
      'normalization':'sum |h|^2 = 1','assumptions':['proper circular complex Gaussian scatterers','zero mean',
      'independent scatterer coefficients with known position-dependent variance','independent looks','same complex PSF in every look','no measurement noise'],
      'target_identity':'E|H eta|^2 = |H|^2 v','convergence':errs,
      'mean_amplitude_is_not_linear':{'correct_mean':'sqrt(pi)/2 * sqrt(|H|^2 v)',
      'wrong_mean':'|H| * (sqrt(pi)/2 * sqrt(v))',
      'relative_L2_wrong_mean_vs_correct_mean':float(np.linalg.norm(amp_surrogate-amp_truth)/np.linalg.norm(amp_truth))}},
      'warnings':['All figures are mathematical illustrations, not user measurements.',
      'Finite-look error is random and need not decrease monotonically for every random seed.',
      'Independent object microstates are explicitly generated; repeated static-speckle acquisitions do not provide them.',
      'The intensity identity targets scatterer covariance/power, not an arbitrary deterministic object amplitude.',
      'Spatially adjacent output pixels remain correlated even when input coefficients are independent.']}
    (evidence/'coherence_metrics.json').write_text(json.dumps(metrics,indent=2))
    np.savez_compressed(evidence/'coherence_plot_data.npz',x=x,amplitude_surrogate=surrogate,phase0=curves[0],phase90=curves[1],phase180=curves[2],ensemble_x=x2,scatterer_variance=variance,expected_intensity=truth,look1=snapshots[1],look64=snapshots[64],look16384=snapshots[16384])
    ix=np.rint(np.linspace(0,len(x)-1,161)).astype(int)
    ix2=np.rint(np.linspace(0,len(x2)-1,161)).astype(int)
    chart_data={
      'metadata':{'source':'simulate_coherence.py','kind':metrics['kind'],'seed':SEED,
       'sampling':'Exact subsampling of computed arrays, no smoothing or model refitting. 161 points per dense curve.',
       'coordinates':metrics['coordinates'],'warnings':metrics['warnings']},
      'coherent_cross_terms':{
       'scatterers':{'x':positions,'series':[{'name':'Identical scatterer powers','y':[1,1]}],
         'x_label':'Object position (a.u.)','y_label':'Scatterer power','render_note':'Two discrete delta-like scatterers; use stems/impulses, not a line connecting the two points.'},
       'amplitudes':{'x':x[ix].tolist(),'series':[
         {'name':'Linear surrogate |h|*|eta|','y':surrogate[ix].tolist(),'style':'dashed'},
         {'name':'Relative phase 0','y':curves[0][ix].tolist()},
         {'name':'Relative phase pi/2','y':curves[1][ix].tolist()},
         {'name':'Relative phase pi','y':curves[2][ix].tolist()}],
         'x_label':'Image coordinate (a.u.)','y_label':'OCT amplitude (a.u.)','display_xlim':[-4,4],'display_ylim':[0,2.1]}},
      'multilook_power':{
       'profiles':{'x':x2[ix2].tolist(),'series':[
         {'name':'1 independent look','y':snapshots[1][ix2].tolist()},
         {'name':'64-look mean','y':snapshots[64][ix2].tolist()},
         {'name':'16384-look mean','y':snapshots[16384][ix2].tolist()},
         {'name':'Theory |h|^2*v','y':truth[ix2].tolist(),'style':'dashed'}],
         'x_label':'Image coordinate (a.u.)','y_label':'Intensity (a.u.)','display_xlim':[-7,7]},
       'convergence':{'x':looks,'series':[
         {'name':'Relative L2 intensity error','y':e.tolist()},
         {'name':'1/sqrt(L) guide','y':(e[0]/np.sqrt(ls)).tolist(),'style':'dashed'}],
         'x_label':'Independent looks L','y_label':'Relative L2 error','x_scale':'log','y_scale':'log'},
       'assumptions':metrics['ensemble']['assumptions']}}
    (evidence/'coherence_chart_data.json').write_text(json.dumps(chart_data,indent=2))
    print(json.dumps({'counterexample':counter,'final_ensemble_error':errs[-1]},indent=2))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--out-root',type=Path,default=Path(__file__).resolve().parents[1]);args=parser.parse_args();main(args.out_root)
