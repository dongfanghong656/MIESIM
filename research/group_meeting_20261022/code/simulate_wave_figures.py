#!/usr/bin/env python3
"""Reproducible synthetic wave-optics figures for the 2026-10-22 group meeting.

Units are micrometres unless explicitly stated. No experimental data are used.
This is scalar, monochromatic, homogeneous-medium diffraction with a uniform
circular pupil. Exact angular-spectrum defocus is used; no paraxial expansion
is needed in the computation. Vector polarization and interfaces are omitted.
"""
from pathlib import Path
import json, os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/wave-fig-mpl')
import numpy as np
from scipy.special import j0
from scipy.interpolate import RegularGridInterpolator
from numpy.polynomial.legendre import leggauss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / 'figures'
EVID = ROOT / 'evidence'
FIG.mkdir(parents=True, exist_ok=True)
EVID.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, GREEN = '#176B91', '#D87438', '#4C8C6A'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,
    'axes.titlesize':16,'axes.labelsize':13,'figure.titlesize':21,
    'xtick.labelsize':11,'ytick.labelsize':11,'legend.fontsize':11,
    'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white',
    'svg.fonttype':'none','figure.facecolor':'white'})
METRICS = {'kind':'SYNTHETIC; no experimental data', 'units':'micrometres',
 'model':{'pupil':'uniform circular amplitude','propagation':'scalar exact angular-spectrum defocus',
          'medium_index':1.0,'omitted':['polarization/vector diffraction','interfaces','dispersion','noise except explicitly stated','sample scattering']},
 'figures':{}, 'verification':{}}

def save(fig, name, caption):
    fig.text(.01,.015,'SYNTHETIC  |  '+caption,ha='left',va='bottom',fontsize=10,color='#49525A')
    fig.savefig(FIG/(name+'.png'),dpi=160,bbox_inches='tight')
    fig.savefig(FIG/(name+'.svg'),bbox_inches='tight')
    plt.close(fig)
    print(name, flush=True)

def radial_psf(r, z, lam=.620, na=.4, n=1., order=400):
    """Hankel integral; peak at in-focus origin is one, no per-plane renormalization."""
    rho, w = leggauss(order); rho=(rho+1)/2; w=w/2
    J = j0((2*np.pi*na/lam)*np.asarray(r)[:,None]*rho[None,:])
    phase = 2*np.pi*n/lam * np.asarray(z)[None,:] * (np.sqrt(1-(na/n*rho[:,None])**2)-1)
    a = J @ ((2*rho*w)[:,None] * np.exp(1j*phase))
    return abs(a)**2

def fft_psf(lam=.620, na=.4, z=0., coeff=None, N=1024, pr=64, density=False):
    """Pupil discretization with same physical image grid if N/pr is fixed.
    Zernike coefficients are OPD micrometres, unit-RMS continuous-disk modes.
    """
    u=(np.arange(N)-N//2)/pr
    U,V=np.meshgrid(u,u,indexing='xy'); r2=U*U+V*V
    mask=(r2<=1).astype(float); theta=np.arctan2(V,U); r=np.sqrt(r2)
    # Mode polynomials only evaluated on disk for meaningful interpretation.
    zerns={'defocus':np.sqrt(3)*(2*r2-1),
           'astigmatism':np.sqrt(6)*(U*U-V*V),
           'coma':np.sqrt(8)*(3*r2-2)*U,
           'spherical':np.sqrt(5)*(6*r2*r2-6*r2+1)}
    W=np.zeros_like(U)
    for key, val in (coeff or {}).items(): W += val*zerns[key]
    phase = 2*np.pi/lam*W + 2*np.pi/lam*z*(np.sqrt(np.maximum(0,1-na*na*r2))-1)
    P=mask*np.exp(1j*phase)
    a=np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(P)))*(N*N)/mask.sum()
    h=abs(a)**2
    dx=lam*pr/(na*N); x=(np.arange(N)-N//2)*dx
    if density: h=h/(h.sum()*dx*dx)
    return x,h

def crop(x,h,lim):
    ix=np.flatnonzero(abs(x)<=lim)
    return x[ix],h[np.ix_(ix,ix)]

def fwhm_nonneg(x, h):
    h=h/h[0]; k=np.flatnonzero(h<=.5)[0]
    return float(2*np.interp(.5,[h[k],h[k-1]],[x[k],x[k-1]]))

def fig1():
    xx=np.linspace(-5,5,401); zz=np.linspace(-16,16,241)
    hz=radial_psf(abs(xx),zz)
    xy=np.linspace(-4,4,301); X,Y=np.meshgrid(xy,xy)
    rr=np.linspace(0,np.sqrt(32),801)
    hrad=radial_psf(rr,[-8,0,8])
    fig=plt.figure(figsize=(16,8.7)); gs=fig.add_gridspec(2,3,height_ratios=[1,1.1],hspace=.40,wspace=.25)
    for i,z in enumerate([-8,0,8]):
        ax=fig.add_subplot(gs[0,i]); v=np.interp(np.hypot(X,Y),rr,hrad[:,i])
        im=ax.imshow(v,extent=[-4,4,-4,4],origin='lower',cmap='magma',norm=LogNorm(1e-4,1))
        ax.set(title=f'z = {z:+d} µm',xlabel='x (µm)',ylabel='y (µm)')
    ax=fig.add_subplot(gs[1,:]); im=ax.imshow(hz.T,extent=[xx[0],xx[-1],zz[0],zz[-1]],origin='lower',aspect='auto',cmap='magma',norm=LogNorm(1e-4,1))
    ax.set(xlabel='x (µm)',ylabel='Defocus z (µm)',title='Through-focus section (y = 0); every plane uses the same intensity scale')
    fig.colorbar(im,ax=fig.axes,shrink=.75,pad=.02,label='Intensity / in-focus peak')
    fig.suptitle('A perfect finite-aperture microscope still defocuses',y=.99)
    save(fig,'wave_01_through_focus','λ = 620 nm, NA = 0.40, n = 1.0; logarithmic intensity; uniform circular pupil')
    # Same r/z values at twice the quadrature order.
    check=radial_psf(np.linspace(0,5,80),np.linspace(-16,16,30),order=800)
    base=radial_psf(np.linspace(0,5,80),np.linspace(-16,16,30),order=400)
    METRICS['verification']['hankel_400_vs_800_max_absolute_intensity_difference']=float(abs(check-base).max())
    METRICS['figures']['wave_01_through_focus']={'lambda_nm':620,'NA':.4,'z_um':[-8,0,8],'intensity_scale':'common, focus peak=1','quadrature_order':400}

def fig2():
    names=['defocus','astigmatism','coma','spherical']
    uv=np.linspace(-1,1,301); U,V=np.meshgrid(uv,uv); r2=U*U+V*V
    Z=[np.sqrt(3)*(2*r2-1),np.sqrt(6)*(U*U-V*V),np.sqrt(8)*(3*r2-2)*U,np.sqrt(5)*(6*r2*r2-6*r2+1)]
    fig,axs=plt.subplots(2,4,figsize=(16,8),layout='constrained')
    for i,name in enumerate(names):
        W=np.where(r2<=1,Z[i]*60,np.nan)
        im0=axs[0,i].imshow(W,extent=[-1,1,-1,1],origin='lower',cmap='RdBu_r',vmin=-175,vmax=175)
        axs[0,i].set(title=name.capitalize(),xlabel='Pupil u',ylabel='Pupil v' if i==0 else '')
        x,h=fft_psf(coeff={name:.060}); xc,hc=crop(x,h,2.7)
        im1=axs[1,i].imshow(hc,extent=[xc[0],xc[-1],xc[0],xc[-1]],origin='lower',cmap='magma',norm=LogNorm(1e-4,1))
        axs[1,i].set(xlabel='x (µm)',ylabel='y (µm)' if i==0 else '')
    fig.colorbar(im0,ax=axs[0,:],shrink=.75,label='OPD (nm)')
    fig.colorbar(im1,ax=axs[1,:],shrink=.75,label='Intensity / ideal peak')
    fig.suptitle('Equal wavefront RMS does not mean equal PSF shape',fontsize=21,y=.99)
    fig.get_layout_engine().set(rect=(0,.07,1,.85))
    save(fig,'wave_02_zernike_modes','60 nm OPD RMS in each mode; λ = 620 nm; NA = 0.40; no extra defocus')
    METRICS['figures']['wave_02_zernike_modes']={'OPD_RMS_nm':60,'modes':names,'normalization':'continuous unit disk RMS=1','intensity_scale':'all relative to ideal peak'}

def fig3():
    x,h0=fft_psf(); _,hd=fft_psf(z=4); _,ha=fft_psf(coeff={'coma':.06}); _,hc=fft_psf(z=4,coeff={'coma':.06})
    hsum=hd+ha-h0; residual=hc-hsum
    xc,ctrue=crop(x,hc,3.5); _,cadd=crop(x,hsum,3.5); _,cres=crop(x,residual,3.5)
    fig,axs=plt.subplots(1,3,figsize=(16,5.6),layout='constrained')
    im=axs[0].imshow(ctrue,origin='lower',extent=[xc[0],xc[-1],xc[0],xc[-1]],cmap='magma',vmin=0,vmax=.5)
    axs[0].set(title='Physical: h(defocus + coma)',xlabel='x (µm)',ylabel='y (µm)')
    lim=max(abs(cadd).max(),abs(cres).max())
    for ax,data,title in zip(axs[1:],[cadd,cres],['Proposed sum: h(d) + h(c) − h(0)','Physical response − proposed sum']):
        imd=ax.imshow(data,origin='lower',extent=[xc[0],xc[-1],xc[0],xc[-1]],cmap='RdBu_r',vmin=-lim,vmax=lim)
        ax.set(title=title,xlabel='x (µm)')
    fig.colorbar(im,ax=axs[0],shrink=.7,label='Intensity / ideal peak')
    fig.colorbar(imd,ax=axs[1:],shrink=.7,label='Signed intensity / ideal peak')
    fig.suptitle('Adding pupil phase does not add image intensity',fontsize=21,y=.99)
    fig.get_layout_engine().set(rect=(0,.16,1,.67))
    err=float(np.linalg.norm(residual)/np.linalg.norm(hc))
    fig.text(.5,.07,f'The proposed sum contains negative values (minimum {hsum.min():.3f}); relative L2 discrepancy = {err:.2f}',ha='center',fontsize=13)
    save(fig,'wave_03_intensity_coupling','z = 4 µm + 60 nm RMS coma, λ = 620 nm, NA = 0.40; common intensity scale')
    METRICS['figures']['wave_03_intensity_coupling']={'defocus_um':4,'coma_RMS_nm':60,'relative_L2_error':err,'minimum_unphysical_additive_prediction':float(hsum.min())}
    # Convergence of a non-symmetric case on exactly the same image sampling grid.
    xf,hf=fft_psf(z=4,coeff={'coma':.06},N=2048,pr=128)
    xc,c0=crop(x,hc,6); xcf,cf=crop(xf,hf,6)
    assert np.allclose(xc,xcf)
    METRICS['verification']['FFT_pupil_radius64_vs128_central12um_relativeL2']=float(np.linalg.norm(c0-cf)/np.linalg.norm(cf))
    METRICS['verification']['same_pupil_total_energy_relative_spread']=float(np.ptp([h0.sum(),hd.sum(),ha.sum(),hc.sum()])/h0.sum())

def fig4():
    l1,l2=.620,.930; s=l2/l1; d=4.
    fig,axs=plt.subplots(2,3,figsize=(16,9),layout='constrained')
    errors={}
    for j,(label,c) in enumerate([('Pure defocus',{}),('Defocus + fixed 60 nm RMS coma',{'coma':.06})]):
        x,h1=fft_psf(lam=l1,z=d,coeff=c,density=True)
        x2,h2=fft_psf(lam=l2,z=s*d,coeff=c,density=True)
        mapped=s*s*h2
        assert np.allclose(x2,s*x)
        xc,p1=crop(x,h1,3.5); _,p2=crop(x,mapped,3.5)
        vmax=max(p1.max(),p2.max())
        for ax,p,ti in zip(axs[j,:2],[p1,p2],['620 nm: h(x, z)','930 nm: s² h(sx, sz), s = 1.5']):
            im=ax.imshow(p,extent=[xc[0],xc[-1],xc[0],xc[-1]],origin='lower',vmin=0,vmax=vmax,cmap='magma')
            ax.set(title=ti,xlabel='Mapped x (µm)',ylabel=label+'\ny (µm)' if ax==axs[j,0] else '')
        mid=len(x)//2; ix=abs(x)<=3.5
        axs[j,2].plot(x[ix],h1[mid,ix],color=BLUE,lw=2,label='620 nm')
        axs[j,2].plot(x[ix],mapped[mid,ix],color=ORANGE,lw=2,ls='--',label='Scaled 930 nm')
        err=float(np.linalg.norm(h1-mapped)/np.linalg.norm(h1)); errors[label]=err
        axs[j,2].set(title=f'y = 0; full-plane relative L2 = {err:.3g}',xlabel='Mapped x (µm)',ylabel='Photon-density kernel (µm⁻²)')
        axs[j,2].legend(); axs[j,2].grid(alpha=.15)
    fig.suptitle('Wavelength similarity requires the pupil phase to stay the same',fontsize=21,y=.99)
    fig.get_layout_engine().set(rect=(0,.06,1,.85))
    save(fig,'wave_04_wavelength_scaling','λ: 620 → 930 nm; x and z both scale by 1.5; same n, NA and pupil amplitude; fixed physical OPD breaks similarity')
    METRICS['figures']['wave_04_wavelength_scaling']={'lambda_nm':[620,930],'scale':s,'defocus_um':[d,s*d],'relative_L2_errors':errors,'comparison':'s² h930(sx,sz) vs h620(x,z) after per-plane photon-density normalization'}

def fig5():
    f=np.linspace(0,.5,1000); sigma=2; H=np.exp(-2*np.pi*np.pi*sigma*sigma*f*f); alpha=1e-3
    fig,axs=plt.subplots(1,2,figsize=(15.8,6),layout='constrained')
    axs[0].semilogy(f,H,color=BLUE,lw=3); axs[0].axhline(np.sqrt(alpha),color=GREEN,ls='--',label='√α')
    axs[0].set(title='Optical transfer attenuates high frequencies',xlabel='Spatial frequency (cycles / sample)',ylabel='|H(f)|',ylim=(1e-9,1.2))
    axs[0].legend(); axs[0].grid(alpha=.2)
    axs[1].semilogy(f,1/H,color=ORANGE,lw=3,label='Direct inverse: 1 / |H|')
    axs[1].semilogy(f,H/(H*H+alpha),color=GREEN,lw=3,label='Regularized: |H| / (|H|² + α)')
    axs[1].set(title='The same inverse also amplifies noise',xlabel='Spatial frequency (cycles / sample)',ylabel='Gain applied to measurement noise',ylim=(1e-6,1e9))
    axs[1].legend(); axs[1].grid(alpha=.2)
    fig.suptitle('A known PSF does not make inversion stable',fontsize=21,y=.99)
    fig.get_layout_engine().set(rect=(0,.08,1,.78))
    save(fig,'wave_05_inverse_noise','1D Gaussian blur, σ = 2 samples; α = 10⁻³; illustrative transfer curves, not a microscope-specific MTF')
    METRICS['figures']['wave_05_inverse_noise']={'model':'1D continuous Gaussian sampled-frequency transfer','sigma_samples':2,'regularization_alpha':alpha,'max_direct_gain_in_plotted_band':float((1/H).max()),'max_regularized_gain':float((H/(H*H+alpha)).max())}
    chartdata={'x':f[::5].tolist(),'series':{'MTF':H[::5].tolist(),'inverse_noise_gain':(1/H)[::5].tolist(),'regularized_noise_gain':(H/(H*H+alpha))[::5].tolist()},'x_unit':'cycles/sample','y_axis':'logarithmic','synthetic':True,'caption':'1D Gaussian blur, sigma=2 samples, alpha=0.001; illustrative, not instrument MTF'}
    (EVID/'wave_05_chart_data.json').write_text(json.dumps(chartdata,indent=2),encoding='utf-8')

def sphere_convolved_profile(x,na,R=.1,lam=.620,order=8):
    """Integrate 3D point-source PSF over an isotropically emitting uniform solid sphere.
    Gaussian quadrature in radius and mu=cos(theta), uniform phi. h is symmetric
    in z for this ideal pupil. Interpolated radial/z table has independently
    controlled diffraction quadrature. Output radial profile.
    """
    rg=np.linspace(0,float(x.max()+2*R+.1),2201)
    zg=np.linspace(0,R,31)
    table=radial_psf(rg,zg,na=na,lam=lam,order=240)
    interp=RegularGridInterpolator((rg,zg),table,bounds_error=True)
    ar,aw=leggauss(order); rad=(ar+1)*R/2; rw=aw*R/2*(3*rad*rad/R**3)
    mu,mw=leggauss(order); mw=mw/2
    phi=np.arange(3*order)*(2*np.pi/(3*order))
    out=np.zeros_like(x)
    for r,w in zip(rad,rw):
        for u,v in zip(mu,mw):
            tx=r*np.sqrt(1-u*u)*np.cos(phi); ty=r*np.sqrt(1-u*u)*np.sin(phi)
            dist=np.sqrt((x[:,None]-tx[None,:])**2+ty[None,:]**2)
            z=np.full_like(dist,abs(r*u))
            vals=interp(np.stack([dist,z],axis=-1))
            out += w*v*vals.mean(axis=1)
    return out

def pixel_integrated_profile(x, radial_r, radial_h, pixel,order=16):
    u,w=leggauss(order); u=u*pixel/2; w=w/2; out=np.zeros_like(x)
    for a,wa in zip(u,w):
        for b,wb in zip(u,w):
            rr=np.sqrt((x-a)**2+b*b)
            out += wa*wb*np.interp(rr,radial_r,radial_h)
    return out

def fig6():
    fig,axs=plt.subplots(1,2,figsize=(16,6.9),layout='constrained')
    result=[]; chartsets=[]
    for i,(mag,na) in enumerate([(20,.4),(5,.14)]):
        p=3.76/mag; xmax=2 if mag==20 else 4.5
        xr=np.linspace(0,xmax+1,1801); xp=np.linspace(0,xmax,1201)
        hp=radial_psf(xr,[0],na=na)[:,0]
        hb=sphere_convolved_profile(xr,na=na,order=8)
        hi=pixel_integrated_profile(xp,xr,hb,pixel=p)
        pprof=np.interp(xp,xr,hp); bprof=np.interp(xp,xr,hb)
        f0=fwhm_nonneg(xp,pprof); fb=fwhm_nonneg(xp,bprof); fi=fwhm_nonneg(xp,hi)
        ax=axs[i]
        for yy,col,lbl,ls in [(pprof,BLUE,'Point emitter','-'),(bprof,ORANGE,'200 nm uniform sphere','--'),(hi,GREEN,'Sphere + square-pixel integration','-')]:
            ax.plot(np.r_[-xp[:0:-1],xp],np.r_[yy[:0:-1],yy]/yy[0],c=col,label=lbl,lw=2.5,ls=ls)
        # Centered sampling markers on the integrated response.
        xs=np.arange(-np.floor(xmax/p),np.floor(xmax/p)+1)*p
        hs=np.interp(abs(xs),xp,hi)/hi[0]
        ax.scatter(xs,hs,color=GREEN,s=24,zorder=3)
        ax.axhline(.5,color='#7D8388',lw=.7,ls=':')
        ax.set(xlabel='x (µm)',ylabel='Profile normalized to its own peak',title=f'{mag}×, NA = {na:.2f}, pixel = {p:.3f} µm',ylim=(-.02,1.07))
        ax.grid(alpha=.13); ax.legend(loc='upper right')
        ax.text(.5,-.20,f'FWHM point: {f0:.3f} µm   |   bead: {fb:.3f} µm\nFWHM bead + pixel: {fi:.3f} µm   |   change: {(fi/f0-1)*100:.2f}%',transform=ax.transAxes,ha='center',va='top',fontsize=12)
        result.append({'magnification':mag,'NA':na,'pixel_um':p,'fwhm_point_um':f0,'fwhm_bead_um':fb,'fwhm_bead_pixel_um':fi,'bead_relative_fwhm_change_percent':100*(fb/f0-1),'bead_pixel_relative_fwhm_change_percent':100*(fi/f0-1),'focal_peak_bead_over_point':float(hb[0]/hp[0]),'nyquist_max_pixel_um':.620/(4*na)})
        fullx=np.r_[-xp[:0:-1],xp]
        series={lbl: (np.r_[yy[:0:-1],yy]/yy[0])[::12].tolist() for lbl,yy in [('Point emitter',pprof),('200 nm sphere',bprof),('Sphere + pixel',hi)]}
        chartsets.append({'magnification':mag,'NA':na,'pixel_um':p,'x':fullx[::12].tolist(),'series':series,'x_unit':'micrometre','y_unit':'own-peak-normalized intensity'})
    fig.suptitle('Finite beads and camera pixels both broaden the measured core',fontsize=21,y=.99)
    fig.get_layout_engine().set(rect=(0,.20,1,.63))
    save(fig,'wave_06_bead_pixel','λ = 620 nm; 200 nm uniformly filled sphere; 3.76 µm sensor pixels; centered sampling; theoretical examples')
    # Increase spherical quadrature on a compact grid to check the integration.
    checkx=np.linspace(0,2,201)
    b8=sphere_convolved_profile(checkx,.4,order=8); b12=sphere_convolved_profile(checkx,.4,order=12)
    METRICS['verification']['sphere_order8_vs12_max_peak_relative_difference']=float(np.max(abs(b8-b12))/b12.max())
    METRICS['figures']['wave_06_bead_pixel']={'sphere_diameter_nm':200,'sensor_pixel_um':3.76,'detector_model':'square integration, no binning; actual magnification assumed equal to nominal','normalization':'own-peak used only for profile width display, not quantitative 3D kernel construction','results':result}
    (EVID/'wave_06_chart_data.json').write_text(json.dumps({'synthetic':True,'models':chartsets,'metrics':result},indent=2),encoding='utf-8')

if __name__=='__main__':
    for fun in [fig1,fig2,fig3,fig4,fig5,fig6]: fun()
    (EVID/'wave_figure_metrics.json').write_text(json.dumps(METRICS,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(METRICS['verification'],indent=2),flush=True)
