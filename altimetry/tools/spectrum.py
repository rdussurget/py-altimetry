# -*- coding: utf-8 -*-
import numpy as np
import scipy.fftpack as ft
#import alti_tools as AT

from scipy import stats

def get_kx(N,dx):
    """
    #+
    # GET_KX
    # @summary: Returns the frequencies to be used with FFT analysis
    #
    # @param N {type:numeric} : number of samples in data
    # @param dx {type:numeric} : sampling step
    #
    # @return:
    #    k: frequency
    #    L: length
    #    imx: index of maximum frequency (for separating positive and negative frequencies)
    #
    # @author: Renaud DUSSURGET (RD) - LER/PAC, Ifremer
    # @change: Created by RD, July 2012
    #-
    """
    
    # to work with dimensional axes
    L = N*dx
    
    #Odd or even
    odd = N&1 and True or False
    
    #Get frequency & wavenumber
    k=ft.fftfreq(N,d=dx)
#    k=f*N
    
    imx = (N-1)/2 if odd else N/2 #Maximum frequency
    
    #Convert into frequencies
#    k/=L
    
    return k, L, imx


def get_spec(dx,Vin,verbose=False,gain=1.0):
    """
    #+
    # GET_SPEC
    # @summary: Returns the spectrum of a regularly sampled dataset
    #
    # @param dq {type:numeric} : sampling interval (1D)
    # @param V {type:numeric} : data to analyse (1D).
    #
    # @note: NaN can not be used. 
    #
    # @return:
    #    psd: Power Spectral Density
    #    esd: Energy Spectral Density
    #    fq: frequency
    #    p: wavelength (period)
    #
    # @author: Renaud DUSSURGET (RD) - LER/PAC, Ifremer
    # @change: Created by RD, July 2012
    #     29/08/2012 : Changed the computation of frequencies and the spectral integration (spectrum is averaged at mid-width frequencies)
    #     30/11/2012 : Outstanding changes : corrected the spectral integration for computing psd and corrected the normalisation
    #-
    """
    V=Vin.copy()
    N=V.size

    # part 1 - estimate of the wavenumbers
    ######################################
    [k,L,imx] = get_kx(N,dx)
#    print "dx = %s" %(dx)
#    print "kx = ",k[0:3]
#    print "shape",k.shape

    # part 2 - estimate of the spectrum
    # fast fourier transform
    ##########################
    fft=ft.fft(V)/(gain)
    if verbose : print 'Check parseval theorem 1: SUM|Y(f)|²={0}, SUM|y(t)|²={1}'.format(((np.abs(fft)**2)/N).sum(),(V**2).sum()) 

#        Compute magnitude
#        CF.Wilks, 2006. Statistical Methods in the Atmopsheric Sciences. 2nd Edition. International Geophysics Series vol.96. Elsevier.
#        => Ck = [Ak^2 + Bk^2]^(1/2), where e^(iwt) = cos(wt) + i.sin(wt)
    a = fft.real
    b = fft.imag
    c = np.sqrt(a**2.0 + b**2.0) #This is the energy
    
    if verbose : print 'Check parseval theorem 2: SUM|Y(f)|²={0}, SUM|y(t)|²={1}'.format(((c**2)/N).sum(),(V**2).sum()) 
    
    #Normalise data
    c /= np.float32(N)
    
    if verbose : print 'Check parseval theorem 3: SUM|Y(f)|²={0}, SUM|y(t)|²={1}'.format(((c**2)*N).sum(),(V**2).sum()) 
    
    mean=fft.real[0]/N #Get average
    
    #Remove negative frequencies
    phase=np.angle(2*fft[1:imx]) #Get phase (not used yet)
    c = 2*c[1:imx-1] #Multiply by 2 (because of negative frequencies removal) - loses a bit of energy
    
    if verbose : print 'Check parseval theorem 4: SUM|Y(f)|²={0}, SUM|y(t)|²={1}'.format(((c**2)*(N/2.0)).sum(),((V-V.mean())**2).sum()) 
    
#    if not(N&1 and True or False) : print (np.sum(fft[1:N/2+1])-np.sum(fft[N/2+1:])), fft[N/2] #This is equal to fft[n/2] 
#    cbis=np.absolute(2*fft[0:imx-1]) #Slightly differing from previous one    
#    cbis /= np.float32(N)
   
    # integration of the spectrum
    # shift to have values centered on the intervals
    dk = k[1] - k[0]
    dk_half = dk/2  # half of the interval
        
    k = k[1:imx-1]  
    k_= k + dk_half  #Shift wavelengths by half the unit
    
    cesd = c**2 *(N/2) #Energy (ESD)
    csquared = c ** 2
    
    if verbose : print 'Check parseval theorem 5: SUM|Y(f)|²={0}, SUM|y(t)|²={1}'.format((csquared*(N/2)).sum(),((V-V.mean())**2).sum())
    
    #Spectral integration
    esd = k_*0.0
    psd = k_*0.0
    for i in np.arange(len(k_)):
        esd[i] = np.sum((csquared * (N/2.0))[(k > (k_[i]-dk)) & (k < (k_[i]+dk))])/2
        psd[i] = np.sum(csquared[(k > (k_[i]-dk)) & (k < (k_[i]+dk))]) #This is variance units integration
        
    psd = psd /dk 

    if verbose : print 'Check parseval theorem 6: SUM|Y(f)|²={0}, SUM|y(t)|²={1}'.format(esd.sum(),((V-V.mean())**2).sum())

    #Get frequencies and period  
    fq=k_
    p=1/fq
    
    return {'psd':psd,'esd':esd,'fq':fq,'p':p}
    
    # Normalisation (Danioux 2011)
#    specvar = specvar / (var.size)**2 / dk #No normalization!!

    #Power spectral density

    # Dimensionalization to be able to compare with the litterature
#    nbwave=k*1000.0/2.0 #/np.pi     # from rad/m to km/m
#    dk=dk*1000.0/2.0 #/np.pi
#    if fft: specvar=specvar*2.0*np.pi/1000.0

def spectral_analysis(dx,Ain,tapering=None,overlap=None,wsize=None,alpha=3.0,detrend=False,normalise=False):

    A=Ain.copy()

    #Check dimensions
    sh = A.shape
    ndims = len(sh)
    N = sh[ndims-1] #Time series are found along the last dimension
    
    #If vector, add one dimension
    if ndims == 1 :
        A = A.reshape((N,1))
        sh = A.shape
        ndims = len(sh)
    
    nr = sh[1] #Number of repeats  
    
#    gain=1.0 #Scaling gain... (used for tapering issues)
    
#    #Get the overall energy
#    spec=get_spec(dx, A[:,0])
#    F=spec['fq']   
#    Eref = ((A[:,0]-A[:,0].mean())**2).sum() #Get the reference energy level
#    ESDref=spec['esd']
#    SFactor=Eref/spec['esd'].sum()
#    ESDref*=SFactor
#    PSDref=spec['psd']*SFactor
#    print 'Check parseval theorem : SUM|Y(f)|²={0}, SUM|y(t)|²={1}'.format(spec['esd'].sum(),((A[:,0]-A[:,0].mean())**2).sum())
    
    
    #Apply tapering if asked
    ########################
    if tapering is not None:
        
        #Set tapering defaults
        overlap=0.50 if overlap is None else overlap
        wsize=0.5*N if wsize is None else wsize

        #Get time splitting (tapering) parameters
        #########################################
        a = np.float32(wsize)
        b = np.float32(overlap) 
        c = np.float32(N) 
        nn=(c - (a * b))/(a - (a * b)) #This is the number of segments
        ix = np.arange(nn) * ((1.0 - b) * a) #These are the starting points of each segments

        #Moving window
        ##############
        dum = np.zeros((wsize, nn),dtype=np.float64)
        for i in np.arange(nn): #looping through time to get splitted time series 
            dum[:,i] = AT.detrend(np.arange(wsize),A[ix[i] : ix[i] + wsize,0]) if detrend else A[ix[i] : ix[i] + wsize,0]
        
        #Set up tapering window
        #######################
        beta=np.pi*alpha
        hamm = np.hamming(wsize)
        hann = np.hanning(wsize)
        kbess = np.kaiser(wsize,beta)
        notaper = np.ones(wsize) #overpass tapering option
        gain=1.0
        
        if isinstance(tapering,bool) : which='hamm'
        elif isinstance(tapering,str) :
            if tapering.upper() == 'HAMMING' :
                which='hamm'
                gain=0.54
            elif tapering.upper() == 'HANNING' :
                which='hann'
                gain=0.50
            if tapering.upper() == 'KAISER' :
                which='kbess'
                gain=0.40
            if tapering.upper() == 'NONE' :
                which='notaper'
                gain=1.0
        elif isinstance(tapering,np.ndarray) : pass
        else :
            raise Exception('Bad value for tapering keyword')
        if not isinstance(tapering,np.ndarray) : exec('window='+which)
        else : window=tapering
        window = np.repeat(window,nn).reshape((wsize,nn))
    
        #Apply tapering on segmented data
        A=dum.copy()*window
        nr=nn
    
    #Run transform
    ###############
    for i in np.arange(nr):
        spec=get_spec(dx, A[:,i])
        if i == 0:
            esd = spec['esd']
            psd = spec['psd']
            fq = spec['fq']
        else : 
            esd = np.append(esd,spec['esd'])
            psd = np.append(psd,spec['psd'])
    
#    factor=((A[:,0]-A[:,0].mean())**2).sum()/spec['esd'].sum()
    
    #Average spectrum
    #################
    nf=len(fq)
    p=1./fq
    esd=esd.reshape(nr,nf)
    psd=psd.reshape(nr,nf)
    esd=(np.sum(esd,axis=0)/nr)#/gain
    psd=(np.sum(psd,axis=0)/nr)#/gain

    #Normalise by energy content    
    Scaling_Factor=len(fq)/esd.sum()
    esd*=Scaling_Factor
    psd*=Scaling_Factor
    
    if tapering is not None : return {'params':{'tapering':tapering is not None,'which':which,'wsize':int(wsize),'nwind':int(nn),'overlap':int(100.*overlap)},'psd':psd,'esd':esd,'fq':fq,'p':p}
    else : return {'params':{'tapering':tapering is not None},'psd':psd,'esd':esd,'fq':fq,'p':p}

#def get_cospec(dx,dy,var1,var2):
#
#    # part 1 - estimate of the wavenumbers
#    [kx,Lx] = get_kx(var1.size,dx)
#
#    # part 2 - estimate of the spectrum
#    # fast fourier transform
#    hat_phi1=np.fft.fft(var1)
#    hat_phi2=np.fft.fft(var2)
#    hat_phi1=hat_phi1*hat_phi1.conj()
#    hat_phi1=hat_phi1.real.copy()  # useless
#    hat_phi2=hat_phi2*hat_phi2.conj()
#    hat_phi2=hat_phi2.real.copy()  # useless
#    hat_cospec = hat_phi1 + hat_phi2
#    
#    # integration of the spectrum
#    # shift to have values centered on the intervals
#    dk = kx[1] - kx[0]
#    dk_half = dk/2  # half of the interval
#    k_= kx[0:min(var1.shape[1],var1.shape[0])/2] + dk
#    cospecvar = k_*0
#    for i in range(len(k_)):
#        # get indexes that satisfy two conditions
#        # integration over a spectral width
#        cospecvar[i] = sum( hat_cospec[ np.nonzero(np.logical_and(kx<=k_[i]+dk_half,kx>k_[i]-dk_half) ) ] )
#
#    # Normalisation (Danioux 2011)
#    cospecvar = cospecvar / (var1.shape[1]*var1.shape[0])**2 / dk
#
#    # Dimensionalization to be able to compare with the litterature
#    nbwave=k_*1000.0/2.0/np.pi     # from rad/m to km/m
#    dk=dk*1000.0/2.0/np.pi
#    cospecvar=cospecvar*2.0*np.pi/1000.0
##    if options.debug: print "\n Normalized co-spectrum : \n",cospecvar
#    return cospecvar,nbwave,dk


#def grid_track(dst,lat,lon,sla,remove_edges=True,backbone=None):
#    """
#    # GRID_TRACK
#    # @summary: This function allow detecting gaps in a set of altimetry data and rebin this data regularlyy, with informations on gaps.
#    # @param dst {type:numeric} : along-track distance.
#    # @param lat {type:numeric} : latitude
#    # @param lon {type:numeric} : longitude
#    # @param sla {type:numeric} : data
#    # @return:
#    #    outdst : resampled distance
#    #    outlon : resampled longitude
#    #    outlat : resampled latitude
#    #    outsla : resampled data
#    #    gaplen : length of the longest gap in data
#    #    ngaps : number of detected gaps in data
#    #    dx : average spatial sampling
#    #    interpolated : True when data was interpolated (empty bin)
#    #
#    # @author: Renaud DUSSURGET (RD) - LER/PAC, Ifremer
#    # @change: Created by RD, July 2012
#    #    29/08/2012 : Major change -> number of output variables changes (added INTERPOLATED), and rebinning modified
#    """
#    
#    #Find gaps in data
#    dx = dst[1:] - dst[:-1]
#    mn_dx = np.median(dx)
#    bins = np.ceil(dst.max() / mn_dx) + 1
#    range=(0/2.,mn_dx * bins) - mn_dx/2
#    hist,bin_edges=np.histogram(dst, bins=bins, range=range) #We have binned the data along a regular grid of size (bins) in the range (range)
#                                                             #Missing data is thus represented by no data in a given bin
#    
#    #Remove leading and trailing edges (and synchronise bin_edges)
#    if remove_edges == True :
#        while hist[0] == 0 : 
#            hist=np.delete(hist,[0])
#            bin_edges=np.delete(bin_edges,[0])
#        while hist[-1] == 0 :
#            hist=np.delete(hist,[len(hist)-1])
#            bin_edges=np.delete(bin_edges,[len(bin_edges)-1])
#    
#    #Get filled bins indices
#
#    ok = np.arange(len(hist)).compress(np.logical_and(hist,True or False))
#    empty = np.arange(len(hist)).compress(~np.logical_and(hist,True or False)) 
#    
#    outsla = np.repeat(np.NaN,len(hist))
#    outlon = np.repeat(np.NaN,len(hist))
#    outlat = np.repeat(np.NaN,len(hist))
#    outdst = bin_edges [:-1]+ mn_dx/2 #distances is taken at bins centers
#    outsla[ok] = sla
#    outlon[ok] = lon
#    outlat[ok] = lat
#    
#    #Fill the gaps if there are some
#    if len(empty) > 0 : 
#        #Interpolate lon,lat @ empty positions
#        outlon[empty] = AT.interp1d(ok, outlon[ok], empty, kind='cubic')
#        outlat[empty] = AT.interp1d(ok, outlat[ok], empty, kind='cubic')
#        outsla[empty] = AT.interp1d(ok, outsla[ok], empty, spline=True)
#        
#    
#    
#    #Get gap properties
#    ind=np.arange(len(hist))
#    dhist=(hist[1:] - hist[:-1])
#    st=ind.compress(dhist==-1)+1
#    en=ind.compress(dhist==1)
#    gaplen=(en-st) + 1
#    ngaps=len(st)
#    
#    #Get empty bin flag
#    interpolated=~hist.astype('bool')
#    
#    return outdst, outlon, outlat, outsla, gaplen, ngaps, dx, interpolated

#def grid_track_backbone(lat,lon,sla,backlat,backlon,fill=None):
#    """
#    # GRID_TRACK_BACKBONE
#    # @summary: This function allow detecting gaps in a set of altimetry data and rebin this data regularlyy, with informations on gaps.
#    # @param dst {type:numeric} : along-track distance.
#    # @param lat {type:numeric} : latitude
#    # @param lon {type:numeric} : longitude
#    # @param sla {type:numeric} : data
#    # @return:
#    #    outdst : resampled distance
#    #    outlon : resampled longitude
#    #    outlat : resampled latitude
#    #    outsla : resampled data
#    #    gaplen : length of the longest gap in data
#    #    ngaps : number of detected gaps in data
#    #    dx : average spatial sampling
#    #    interpolated : True when data was interpolated (empty bin)
#    #
#    # @author: Renaud DUSSURGET (RD) - LER/PAC, Ifremer
#    # @change: Created by RD, July 2012
#    #    29/08/2012 : Major change -> number of output variables changes (added INTERPOLATED), and rebinning modified
#    """
#    
#    #Find gaps in data
#    dst=AT.calcul_distance(backlat[0],backlon[0],lat,lon)
#    dstback=AT.calcul_distance(backlat,backlon)
#    dx = dstback[1:] - dstback[:-1]
#    mn_dx = np.median(dx)
#    
#    bins = np.ceil(dstback.max() / mn_dx) + 1
#    range=(0/2.,mn_dx * bins) - mn_dx/2
#    hist,bin_edges=np.histogram(dst, bins=bins, range=range) #We have binned the data along a regular grid of size (bins) in the range (range)
#                                                             #Missing data is thus represented by no data in a given bin
#    
#    #Get filled bins indices
#    ok = np.arange(len(hist)).compress(np.logical_and(hist,True or False))
#    empty = np.arange(len(hist)).compress(~np.logical_and(hist,True or False)) 
#    
#    outsla = np.repeat(np.NaN,len(hist))
#    outlon = np.repeat(np.NaN,len(hist))
#    outlat = np.repeat(np.NaN,len(hist))
#    outdst = bin_edges [:-1]+ mn_dx/2 #distances is taken at bins centers
#    outsla[ok] = sla
#    outlon[ok] = lon
#    outlat[ok] = lat
#    
#    #Fill the gaps if there are some
#    if (fill is not None) & (len(empty) > fill) : 
#        #Interpolate lon,lat @ empty positions
#        outlon[empty] = AT.interp1d(ok, outlon[ok], empty, kind='cubic')
#        outlat[empty] = AT.interp1d(ok, outlat[ok], empty, kind='cubic')
#        outsla[empty] = AT.interp1d(ok, outsla[ok], empty, spline=True)
#    
#    #Get empty bin flag
#    interpolated=(~hist.astype('bool'))
#    
#    return outdst, outlon, outlat, outsla, dx, interpolated

def get_slope(fq,spec,degree=1):
    """
    #+
    # GET_SLOPE
    # @summary: This function returns the spectral slope of a spectrum using a least-square regression 
    #
    # @param fq {type:numeric} : frequency
    # @param spec {type:numeric} : spectrum data
    #
    # @keyword degree {type:numeric}{default:1}: Degree of the least-square regression model 
    #
    # @return:
    #    slope : spectral slope (or model coefficients for a higher order model)
    #    intercept : Energy at unit frequency (1 cpkm)
    #
    # @author: Renaud DUSSURGET (RD) - LER/PAC, Ifremer
    # @change: Created by RD, August 2012
    #-
    """
    
    x = np.log10(fq).flatten()
    y = np.log10(spec).flatten()
    
    
    #1) Linear least-square regression
    if degree == 1 :
        (slope, intercept, rval, pval, err) = stats.linregress(x,y)

    #2) Least-square regression using a higher-order spectral model
    # -> Gives the same results as 
    #cf. http://pingswept.org/2009/01/24/least-squares-polynomial-fitting-in-python/                                                    
    else :
        A = np.vander(x, degree+1) # form the Vandermonde matrix 
        (coeffs, residuals, rank, sing_vals) = np.linalg.lstsq(A,y) # find the x that minimizes the norm of Ax-y
        (slope[:-1], intercept) = coeffs

    return slope,intercept

#+
# YULE_WALKER_REGRESSION : Estimation of an AR (autoregression) spectral model from data
#@todo: To be written
#-
def yule_walker_regression():
    pass