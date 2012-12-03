import numpy as np
import alti_tools as AT
import spectrum as sp
import matplotlib.pyplot as plt

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
detrend = True #Remove any linear trend from the input signal
normalise = True #Normalise the output spectrum by the overall energy content
add_trend = False #Adds a linear trend over the noise (overpasses the detrend option - first detrend, then add a linear trend)



#Generate white and red noise
#############################
wn=np.random.randn(N)
rn=np.cumsum(wn)
if detrend :
    rn=AT.detrend(np.arange(N), rn)
    wn=AT.detrend(np.arange(N), wn) 

if add_trend :
    rn= rn + (2./N)*np.arange(N)
    wn= wn + (2./N)*np.arange(N)

#Setup plot parameters
######################
if noise == 'white':
    A=wn
    if normalise :
        eaxis=[1e-3,1e1,1e-1,1e1]
        paxis=[1e-3,1e1,1e-0,1e2]
        ezoom=np.copy(eaxis).tolist();ezoom[0]=1e-1;ezoom[1]=1e0
        pzoom=np.copy(paxis).tolist();pzoom[0]=1e-1;pzoom[1]=1e0
    else :
        eaxis=[1e-3,1e1,1e-1,1e1]
        paxis=[1e-3,1e1,1e0,1e2]
        ezoom=np.copy(eaxis).tolist();ezoom[0]=1e-1;ezoom[1]=1e0
        pzoom=np.copy(paxis).tolist();pzoom[0]=1e-1;pzoom[1]=1e0
    noise='White'
elif noise == 'red' :
    A=rn
    if normalise :
        eaxis=[1e-3,1e1,1e-4,1e3]
        paxis=[1e-3,1e1,1e-3,1e4]
        ezoom=[1e-1,1e0,1e-3,1e-1]
        pzoom=[1e-1,1e0,1e-2,1e0]
    else :
        eaxis=[1e-3,1e1,1e-1,1e5]
        paxis=[1e-3,1e1,1e0,1e6]
        ezoom=np.copy(eaxis).tolist();ezoom[0]=1e-1;ezoom[1]=1e0
        pzoom=np.copy(paxis).tolist();pzoom[0]=1e-1;pzoom[1]=1e0
    noise='Red'


#Run spectral analysis
######################
spec=sp.spectral_analysis(1.0, A, detrend=detrend,normalise=normalise)
F=spec['fq']
ESDref=spec['esd']
PSDref=spec['psd']
shamm=sp.spectral_analysis(1.0, A, tapering='hamming', overlap=0.5, wsize=wsize, detrend=detrend,normalise=normalise)
skais=sp.spectral_analysis(1.0, A, tapering='kaiser', overlap=0.75, wsize=wsize, detrend=detrend,normalise=normalise)
snone=sp.spectral_analysis(1.0, A, tapering='none',overlap=0.0, wsize=wsize, detrend=detrend,normalise=normalise)

#Plot results
#############
if add_trend : plt.subplot(3,1,1);plt.plot(A);plt.title('Spectral analysis of {0} Noise  of length {1} pts + linear trend'.format(noise,N))
else : plt.subplot(3,1,1);plt.plot(A);plt.title('Spectral analysis of {0} Noise  of length {1} pts'.format(noise,N))
plt.subplot(3,2,3);plt.loglog(F,ESDref,'-r');plt.loglog(snone['fq'],snone['esd'],'-c',linewidth=2);plt.loglog(shamm['fq'],shamm['esd'],'-b');plt.loglog(skais['fq'],skais['esd'],'-g');plt.axis(eaxis);plt.ylabel('ESD (cm2)');leg=plt.legend(['{0} noise'.format(noise),'Segmented (N={0},WS={1},over={2})'.format(snone['params']['nwind'],snone['params']['wsize'],snone['params']['overlap']),'Hamming (N={0},WS={1},over={2})'.format(shamm['params']['nwind'],shamm['params']['wsize'],shamm['params']['overlap']),'Kaiser-Bessel(N={0},WS={1},over={2})'.format(skais['params']['nwind'],skais['params']['wsize'],skais['params']['overlap'])],prop={'size':8},labelspacing=0);plt.grid(True)
plt.subplot(3,2,5);plt.loglog(F,PSDref,'-r');plt.loglog(snone['fq'],snone['psd'],'-c',linewidth=2);plt.loglog(shamm['fq'],shamm['psd'],'-b');plt.loglog(skais['fq'],skais['psd'],'-g');plt.axis(paxis);plt.xlabel('Frequency (cph)');plt.ylabel('PSD (cm2.cph)');leg=plt.legend(['{0} noise'.format(noise),'Segmented (N={0},WS={1},over={2})'.format(snone['params']['nwind'],snone['params']['wsize'],snone['params']['overlap']),'Hamming (N={0},WS={1},over={2})'.format(shamm['params']['nwind'],shamm['params']['wsize'],shamm['params']['overlap']),'Kaiser-Bessel(N={0},WS={1},over={2})'.format(skais['params']['nwind'],skais['params']['wsize'],skais['params']['overlap'])],prop={'size':8},labelspacing=0);plt.grid(True)
plt.subplot(3,2,4);plt.loglog(F,ESDref,'-r');plt.loglog(snone['fq'],snone['esd'],'-c',linewidth=2);plt.loglog(shamm['fq'],shamm['esd'],'-b');plt.loglog(skais['fq'],skais['esd'],'-g');plt.axis(ezoom);plt.ylabel('ESD (cm2)');plt.grid(True);
plt.subplot(3,2,6);plt.loglog(F,PSDref,'-r');plt.loglog(snone['fq'],snone['psd'],'-c',linewidth=2);plt.loglog(shamm['fq'],shamm['psd'],'-b');plt.loglog(skais['fq'],skais['psd'],'-g');plt.axis(pzoom);plt.xlabel('Frequency (cph)');plt.ylabel('PSD (cm2.cph)');plt.grid(True)
plt.show()

