#!/usr/bin/env python
"""
Trivial example of aurora using Stan Solomon's GLOW Auroral model
code wrapping in Python by Michael Hirsch
"""
import numpy as np
import xarray
from datetime import datetime
from matplotlib.pyplot import figure
from pandas import DataFrame
from numpy import hstack,arange,append,array,rollaxis
from os import chdir
try:
    import seaborn
except ImportError:
    pass
#
from sciencedates import datetime2yeardoy as datetime2yd
#from pyiri90.runiri90 import runiri
#from msise00.runmsis import rungtd1d
import glowaurora
from glowaurora import glowfort
#
glowpath=glowaurora.__path__[0]

#def runglowaurora(eflux,e0,dt,glat,glon,f107a,f107,f107p,ap,mass):
def runglowaurora(params:dict, z_km:np.ndarray=None) -> xarray.Dataset:
    """ Runs Fortran GLOW program and collects results in friendly arrays with metadata. """
    #%% (-2) check/process user inputs
    assert isinstance(params['flux'],(float,int,np.ndarray))
    assert isinstance(params['E0'],   (float,'float32',int))
    assert isinstance(params['t0'],   (datetime,str))
    assert isinstance(params['glat'], (float,int))
    assert isinstance(params['glon'], (float,int))
#%% (-1) if no manual f10.7 and ap, autoload by date
    if not 'f107a' in params or params['f107a'] is None:
        if readmonthlyApF107 is None:
            raise ImportError(GRIDERR)
        f107Ap = readmonthlyApF107(params['t0'])
        params['f107a'] = params['f107p'] 
#= f107Ap['f107s']
        params['f107']  = f107Ap['f107o']
        params['Ap']    = (f107Ap['Apo'],)*7
    chdir(glowpath)
    yd,utsec = datetime2yd(params['t0'])[:2]

    glowfort.cglow.zz = z_km*1e5
    glowfort.cglow.znd[:]=0.

# Set other parameters and switches:
    glowfort.cglow.jlocal = 0
    glowfort.cglow.kchem = params['kchem'] #set in the main program
    glowfort.cglow.iscale = 0
    glowfort.cglow.xuvfac = 3.
    glowfort.cglow.hlybr = 0.
    glowfort.cglow.fexvir = 0.
    glowfort.cglow.hlya = 0.
    glowfort.cglow.heiew = 0.
#%% (1) setup flux at top of ionosphere
    ener,dE = glowfort.egrid()

    phitop = glowfort.maxt(params['flux'],params['E0'],ener, dE, itail=0, fmono=0, emono=0)

    phi = hstack((ener[:,None],dE[:,None],phitop[:,None]))

    glowfort.cglow.zz = z_km*1e5
    glowfort.cglow.znd[:]=0.
#%% (2) call MSIS
    #tselecopts = (1,)*25
    #dens,temp=rungtd1d(dt,z,glat,glon,f107a,f107,ap,mass,tselecopts) #no need for msis neutral densities get passed as parameters in from the main program
    glowfort.cglow.zo = params['O']/1e6
    glowfort.cglow.zn2 = params['N2']/1e6
    glowfort.cglow.zo2 = params['O2']/1e6
    glowfort.cglow.zrho = params['Total']/1e6
    glowfort.cglow.zns  = params['N']/1e6
    #glowfort.cglow.ztn  = temp['heretemp']
    glowfort.cglow.ztn  = params['NT'] #neutral temperature
#%% (3) call snoem
    """
    Call SNOEMINT to obtain NO profile from the Nitric Oxide Empirical
    Model (NOEM)
    """
    #glowfort.cglow.zno = glowfort.snoemint(dt.strftime('%Y%j'),
                               #glat,glon,f107,ap,z,temp['heretemp'])
    glowfort.cglow.zno=params['NO']/1e6
#%% (4a) call iri-90
    #outf,oarr = runiri(dt,z,glat,glon,f107,f107a,ap,mass=48)
    #chdir(glowpath) #need this since iri90 changes path
#%% (4b) store iri90 in COMMON blocks, after unit conversion
    glowfort.cglow.ze = params['ne']/1e6 # M-3 -> CM-3
    glowfort.cglow.ze[glowfort.cglow.ze<100.] = 100.

    glowfort.cglow.zti = params['Ti']
    i = glowfort.cglow.zti<glowfort.cglow.ztn
    glowfort.cglow.zti[i] = glowfort.cglow.ztn[i]

    glowfort.cglow.zte = params['Te']
    i = glowfort.cglow.zte<glowfort.cglow.ztn
    glowfort.cglow.zte[i] = glowfort.cglow.ztn[i]

    glowfort.cglow.zxden[2,:] = params['nO+']/1e6
    glowfort.cglow.zxden[5,:] = params['nO2+']/1e6
    glowfort.cglow.zxden[6,:] = params['nNO+']/1e6
#%% glow model

    ion,ecalc,photI,ImpI,isr = glowfort.auroran(z_km,yd,utsec,params['glat'],params['glon']%360,
                                             params['f107a'],params['f107'],phi)
#%% handle the outputs including common blocks
    zeta=glowfort.cglow.zeta.T #columns 11:20 are identically zero

    ver = DataFrame(index=z,
                    data=zeta[:,:11],
                    columns=[3371, 4278, 5200, 5577, 6300,7320,10400,3466,
                             7774, 8446,3726])
    photIon = DataFrame(index=z,
                   data=hstack((photI[:,None],ImpI[:,None],ecalc[:,None],ion)),
                    columns=['photoIoniz','eImpactIoniz','ne',
                    'nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+',
                    'nO','nO2','nN2','nNO'])

    isrparam = DataFrame(index=z,
                         data=isr,
                         columns=['ne','Te','Ti'])

    phitop = DataFrame(index=phi[:,0], #eV
                       data=phi[:,2],  #diffnumflux
                       columns=['diffnumflux'])
    zceta = glowfort.cglow.zceta.T

    return ver,photIon,isrparam,phitop,zceta
#%% plot
def plotaurora(phitop,ver,zceta,photIon,isr,dtime,glat,glon):
#%% incident flux at top of ionosphere
    ax = figure().gca()
    ax.plot(phitop.index,phitop['diffnumflux'])
    ax.set_title('Incident Flux',fontsize='x-large')
    ax.set_xlabel('Beam Energy [eV]',fontsize='large')
    ax.set_ylabel('Flux',fontsize='large')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(bottom=1e-4)
    ax.tick_params(axis='both',which='major',labelsize='medium')
#%% results of impacts
    fg = figure()
    axs = fg.subplots(1,4,sharey=True, figsize=(15,8))
    fg.suptitle('{} ({},{})'.format(dtime,glat,glon),fontsize='x-large')
    fg.tight_layout(pad=3.2, w_pad=0.3)

    ax = axs[0]
    ax.plot(ver.values,ver.index)
    ax.set_xlabel('VER',fontsize='large')
    ax.set_ylabel('altitude [km]',fontsize='large')
    ax.set_ylim(top=ver.index[-1],bottom=ver.index[0])
    ax.set_xscale('log')
    ax.set_xlim(left=1e-4)
    ax.legend(ver.columns)
    ax.set_title('Volume emission rate',fontsize='x-large')

    ax = axs[1]
    ax.plot(photIon[['photoIoniz','eImpactIoniz']],photIon.index)
    ax.set_xlabel('ionization',fontsize='large')
    ax.set_xscale('log')
    ax.set_xlim(left=1e-1)
    ax.legend(photIon.columns[:2])
    ax.set_title('Photo and e$^-$ impact ionization',fontsize='x-large')

    ax = axs[2]
    ax.semilogx(photIon[['ne','nO+(2P)','nO+(2D)','nO+(4S)','nN+','nN2+','nO2+','nNO+',
                    'nO','nO2','nN2','nNO']], photIon.index)
    ax.set_xlabel('Density',fontsize='large')
    ax.set_xscale('log')
    ax.set_xlim(left=1e-3)
    ax.legend(photIon.columns[2:])
    ax.set_title('Electron and Ion Densities',fontsize='x-large')

    ax = axs[3]
    ax.semilogx(isr[['Te','Ti']], isr.index)
    ax.set_xlabel('Temperature [K]',fontsize='large')
    ax.legend(isr.columns[1:])
    ax.set_title('Particle Temperature',fontsize='x-large')

    for a in axs:
        a.grid(True)
        a.tick_params(axis='both',which='major',labelsize='medium')

#%% total energy deposition vs. altitude
    fg = figure()
    axs = fg.subplots(1,2,sharey=True, figsize=(15,8))
    fg.suptitle('{} ({},{})'.format(dtime,glat,glon),fontsize='x-large')
    fg.tight_layout(pad=3.2, w_pad=0.3)

    ax = axs[0]
    tez = glowfort.cglow.tez
    ax.plot(tez,ver.index)
    ax.set_xscale('log')
    ax.set_xlim(left=1e-1)
    ax.set_ylim(top=ver.index[-1],bottom=ver.index[0])
    ax.set_xlabel('Energy Deposited',fontsize='large')
    ax.set_ylabel('Altitude [km]',fontsize='large')
    ax.set_title('Total Energy Depostiion',fontsize='x-large')
#%% e^- impact ionization rates from ETRANS
    ax = axs[1]
    sion = glowfort.cglow.sion
    sion = DataFrame(index=ver.index,data=sion.T,columns=['O','O2','N2'])
    ax.plot(sion,ver.index)
    ax.set_xscale('log')
    ax.set_xlim(left=1e-6)
    ax.set_xlabel('e$^-$ impact ioniz. rate',fontsize='large')
    ax.set_title('electron impact ioniz. rates',fontsize='x-large')
    #ax.legend(True)
#%% constituants of per-wavelength VER
#    zcsum = zceta.sum(axis=-1)

    ax = figure().gca()
    for zc in rollaxis(zceta,1):
        ax.plot(ver.index,zc)
    ax.set_xlabel('emission constituants',fontsize='large')
    #ax.legend(True)
