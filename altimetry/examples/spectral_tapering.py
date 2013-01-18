import numpy as np
import matplotlib.pyplot as plt
from altimetry.tools import spectrum as sp, detrend as detrend_fun

"""
    SPECTRAL_TAPERING
    @summary: Test spectral tapering on a white or red noise.
    @author: Renaud Dussurget, LER/PAC IFREMER
    @change: First created by RD, December 2012
"""

#Test parameters
################
N=2000 #Number of points in the time series
wsize=0.5*N #Size of windows used for tapering
noise='white' #Type of noise to analyse (white or red)
predetrend = False  #Remove any linear trend prior of the analysis 
detrend = True #Remove any linear trend from the input signal
normalise = False #Normalise the output spectrum by the overall energy content
add_trend = True #Adds a linear trend over the noise (overpasses the detrend option - first detrend, then add a linear trend)
trend_type = 2 #1 linear, 2 sine, 3 exp

#Generate white and red noise
#############################
wn=np.random.randn(N)
rn=np.cumsum(wn)
#if predetrend :
rn=detrend_fun(np.arange(N), rn)
wn=detrend_fun(np.arange(N), wn) 

if add_trend :
    cr = 100.0/N #These are the coefficients of the trends (should be based on random signal variability - greater for red noise than white). 
    cw = 20.0/N
    
    if trend_type == 1:
        rtrend=cr*np.arange(N)
        wtrend=cw*np.arange(N)
        trend_name='linear'
    elif trend_type == 2:
        rtrend=cr * N * np.cos((np.arange(N)/float(N)) * (np.pi/2)) #Change this factor to de/increase the trend wavelength
        wtrend=cw * N * np.cos((np.arange(N)/float(N)) * (np.pi/2))
        trend_name='sine (0->0.5pi)'
    elif trend_type == 3:
        rtrend=0
        wtrend=0
        trend_name='NONE'
        
    rn= rn + rtrend
    wn= wn + wtrend

if predetrend :
    rn=detrend_fun(np.arange(N), rn)
    wn=detrend_fun(np.arange(N), wn) 

#Setup plot parameters
######################
if noise == 'white':
    A=wn
    if normalise :
        eaxis=[1e-3,1e0,1e-2,1e2]
        paxis=[1e-3,1e0,1e-2,1e2]
        ezoom=np.copy(eaxis).tolist();ezoom[0]=1e-1;ezoom[1]=1e0
        pzoom=np.copy(paxis).tolist();pzoom[0]=1e-1;pzoom[1]=1e0
    else :
        eaxis=[1e-3,1e0,1e-1,1e3]
        paxis=[1e-3,1e0,1e-1,1e3]
        ezoom=np.copy(eaxis).tolist();ezoom[0]=1e-1;ezoom[3]=1e2
        pzoom=np.copy(paxis).tolist();pzoom[0]=1e-1;pzoom[3]=1e2
    noise='White'
elif noise == 'red' :
    A=rn
    if normalise :
        eaxis=[1e-3,1e0,1e-3,1e3]
        paxis=[1e-3,1e0,1e-3,1e3]
        ezoom=np.copy(eaxis).tolist();ezoom[0]=1e-1;ezoom[3]=1e-1
        pzoom=np.copy(paxis).tolist();pzoom[0]=1e-1;pzoom[3]=1e-1
    else :
        eaxis=[1e-3,1e0,1e-1,1e6]
        paxis=[1e-3,1e0,1e-1,1e6]
        ezoom=np.copy(eaxis).tolist();ezoom[0]=1e-1;ezoom[3]=1e2
        pzoom=np.copy(paxis).tolist();pzoom[0]=1e-1;pzoom[3]=1e2
    noise='Red'


#Run spectral analysis
######################
spec=sp.spectral_analysis(1.0, A, detrend=detrend,normalise=normalise)
F=spec['fq']
ESDref=spec['esd']
PSDref=spec['psd']
shamm=sp.spectral_analysis(1.0, A, tapering='hamming', overlap=0.5, wsize=wsize, detrend=False,normalise=normalise)
skais=sp.spectral_analysis(1.0, A, tapering='kaiser', overlap=0.75, wsize=wsize, detrend=False,normalise=normalise)
snone=sp.spectral_analysis(1.0, A, tapering='none',overlap=0.0, wsize=wsize, detrend=detrend,normalise=normalise)

#snone_at_all=sp.spectral_analysis(1.0, A, detrend=detrend,normalise=normalise)
#plt.loglog(snone_at_all['fq'],snone_at_all['psd'],'-c',linewidth=2);plt.loglog(snone_at_all['fq'],snone_at_all['psd'],'-b',linewidth=2);plt.show()

#Plot results
#############
if add_trend : plt.subplot(3,1,1);plt.plot(A);plt.plot(np.arange(N), rtrend if noise == 'Red' else wtrend,'-k');plt.title('Spectral analysis of {0} Noise  of length {1} pts + {2} trend'.format(noise,N,trend_name))
else : plt.subplot(3,1,1);plt.plot(A);plt.title('Spectral analysis of {0} Noise  of length {1} pts'.format(noise,N))
plt.subplot(3,2,3);plt.loglog(F,ESDref,'-r');plt.loglog(snone['fq'],snone['esd'],'-c',linewidth=2);plt.loglog(shamm['fq'],shamm['esd'],'-b');plt.loglog(skais['fq'],skais['esd'],'-g');plt.axis(eaxis);plt.ylabel('ESD (cm2)');leg=plt.legend(['{0} noise'.format(noise),'Segmented (N={0},WS={1},over={2})'.format(snone['params']['nwind'],snone['params']['wsize'],snone['params']['overlap']),'Hamming (N={0},WS={1},over={2})'.format(shamm['params']['nwind'],shamm['params']['wsize'],shamm['params']['overlap']),'Kaiser-Bessel(N={0},WS={1},over={2})'.format(skais['params']['nwind'],skais['params']['wsize'],skais['params']['overlap'])],prop={'size':8},labelspacing=0);plt.grid(True)
plt.subplot(3,2,5);plt.loglog(F,PSDref,'-r');plt.loglog(snone['fq'],snone['psd'],'-c',linewidth=2);plt.loglog(shamm['fq'],shamm['psd'],'-b');plt.loglog(skais['fq'],skais['psd'],'-g');plt.axis(paxis);plt.xlabel('Frequency (cph)');plt.ylabel('PSD (cm2.cph)');leg=plt.legend(['{0} noise'.format(noise),'Segmented (N={0},WS={1},over={2})'.format(snone['params']['nwind'],snone['params']['wsize'],snone['params']['overlap']),'Hamming (N={0},WS={1},over={2})'.format(shamm['params']['nwind'],shamm['params']['wsize'],shamm['params']['overlap']),'Kaiser-Bessel(N={0},WS={1},over={2})'.format(skais['params']['nwind'],skais['params']['wsize'],skais['params']['overlap'])],prop={'size':8},labelspacing=0);plt.grid(True)
plt.subplot(3,2,4);plt.loglog(F,ESDref,'-r');plt.loglog(snone['fq'],snone['esd'],'-c',linewidth=2);plt.loglog(shamm['fq'],shamm['esd'],'-b');plt.loglog(skais['fq'],skais['esd'],'-g');plt.axis(ezoom);plt.ylabel('ESD (cm2)');plt.grid(True);
plt.subplot(3,2,6);plt.loglog(F,PSDref,'-r');plt.loglog(snone['fq'],snone['psd'],'-c',linewidth=2);plt.loglog(shamm['fq'],shamm['psd'],'-b');plt.loglog(skais['fq'],skais['psd'],'-g');plt.axis(pzoom);plt.xlabel('Frequency (cph)');plt.ylabel('PSD (cm2.cph)');plt.grid(True)
plt.show()

