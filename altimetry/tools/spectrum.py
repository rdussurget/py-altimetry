# -*- coding: utf-8 -*-
import numpy as np
import scipy.fftpack as ft
from scipy import stats
from altimetry_tools import detrend as detrend_fun, grid_track
if __debug__ : import matplotlib.pyplot as plt

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


def get_spec(dx,Vin,verbose=False,gain=1.0,integration=True):
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
    
#    mean=fft.real[0]/N #Get average
    
    #Remove negative frequencies
#    phase=np.angle(2*fft[1:imx]) #Get phase (not used yet)
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
    k_= k[:-1] + dk_half  #Shift wavelengths by half the unit
    
    cesd = c**2 *(N/2) #Energy (ESD)
    csquared = c ** 2 
    
    if verbose : print 'Check parseval theorem 5: SUM|Y(f)|²={0}, SUM|y(t)|²={1}'.format((csquared*(N/2)).sum(),((V-V.mean())**2).sum())
    
    #Spectral integration
    if integration :
        esd = k_*0.0
        psd = k_*0.0
        for i in np.arange(len(k_)):
            esd[i] = np.sum((csquared * (N/2.0))[(k > (k_[i]-dk)) & (k < (k_[i]+dk))]) / 2.0
            psd[i] = np.sum(csquared[(k > (k_[i]-dk)) & (k < (k_[i]+dk))]) / (2.0**2) #This is variance units integration
        fq=k_
    else :
        esd = csquared
        psd = esd.copy()  /2.
        fq=k.copy()
        
    psd = psd / dk
#    psd = psd/ (N**2.0) / dk    # Normalisation (Danioux 2011) 

    if verbose : print 'Check parseval theorem 6: SUM|Y(f)|²={0}, SUM|y(t)|²={1}'.format(esd.sum(),((V-V.mean())**2).sum())

    #Get frequencies and period  
    p=1/fq
    
    return {'psd':psd,'esd':esd,'fq':fq,'p':p,'gain':gain}
    
    # Normalisation (Danioux 2011)
#    specvar = specvar / (var.size)**2 / dk #No normalization!!

    #Power spectral density

    # Dimensionalization to be able to compare with the litterature
#    nbwave=k*1000.0/2.0 #/np.pi     # from rad/m to km/m
#    dk=dk*1000.0/2.0 #/np.pi
#    if fft: specvar=specvar*2.0*np.pi/1000.0

def spectral_analysis(dx,Ain,tapering=None,overlap=None,wsize=None,alpha=3.0,detrend=False,normalise=False,integration=True):
    """
     Spectral_Analysis
     @summary: This function performs a spatial spectral analysis with different options on a time series of SLA profiles.
     @param dx {type:numeric} : sampling distance
     @param Ain {type:numeric} : 2D table of sla data with time along 2nd axis (NXxNT with NX the spatial length and NT the time length)
     @keyword tapering {type:string|bool|nd.array} : apply tapering to the data. <br \>
                    If this keyword is of type bool : apply hamming window. <br \>
                    If this keyword is a string : apply a hamming ('hamm'), hann ('hann'), kaiser-bessel ('kaiser'), kaiser-bessel ('blackman') or no ('none') tapering function. <br \>
                    If this keyword is an nd.array aobject : apply this array as taper.
     @keyword overlap {type:float} : overlap coefficient of the windows (0.75 means 75% overlap).
     @keyword wsize {type:numeric} : size of the sub-segments.
     @keyword normalise {type:bool,default:False} : If True, normalise the spectrum by its overall energy content.
     @keyword detrend {type:bool,default:False} : If True, removes a linear trend to the segmented signal (if tapered) or to the whole signal (if not tapered).
     @keyword integration {type:bool,default:False} : If True, integrate the spectrum between 2 frequencies. 
     @param sla {type:numeric} : data
     @return: a spectrum structrue with Energy Spectral Density ('esd'), Power Spectral Density ('PSD'), frequency ('fq'), wavelength ('p') and tapering parameters.
    
     @author: Renaud DUSSURGET (RD) - LER/PAC, Ifremer
     @change: Created by RD, December 2012
    """
    
    A=Ain.copy()

    #Check dimensions
    sh = A.shape
    ndims = len(sh)
    N = sh[0] #Time series are found along the last dimension
    
    #If vector, add one dimension
    if ndims == 1 :
        A = A.reshape((N,1))
        sh = A.shape
        ndims = len(sh)
    
    nr = sh[1] #Number of repeats  
    nt = nr
    
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
        nn=np.floor((c - (a * b))/(a - (a * b))) #This is the number of segments
        print 'Number of windows :{0}\nTotal windowed points : {1} ({2} missing)\nTotal points : {3}'.format(nn,nn*wsize,N - nn*wsize,N)
        
        ix = np.arange(nn) * ((1.0 - b) * a) #These are the starting points of each segments

        #Moving window
        ##############
        dum = np.zeros((wsize, nn, nr),dtype=np.float64)
        for j in np.arange(nr):
            for i in np.arange(nn): #looping through time to get splitted time series 
                dum[:,i,j] = detrend_fun(np.arange(wsize),A[ix[i] : ix[i] + wsize,j]) if detrend else A[ix[i] : ix[i] + wsize,j]
        
        #Set up tapering window
        #######################
        beta=np.pi*alpha
        hamm = np.hamming(wsize)
        hann = np.hanning(wsize)
        kbess = np.kaiser(wsize,beta)
        blackman = np.blackman(wsize)
        notaper = np.ones(wsize) #overpass tapering option
        gain=1.0
        
        if isinstance(tapering,bool) : which='hamm'
        elif isinstance(tapering,str) :
            if tapering.upper() == 'HAMMING' :
                which='hamm'
                gain=np.sum(hamm)/wsize #0.530416666667
            elif tapering.upper() == 'HANNING' :
                which='hann'
                gain=np.sum(hann)/wsize #0.489583333333
            elif tapering.upper() == 'KAISER' :
                which='kbess'
                gain=np.sum(kbess)/wsize #0.394170357504
            elif tapering.upper() == 'NONE' :
                which='notaper'
                gain=1.0
            elif tapering.upper() == 'BLACKMAN' :
                which='blackman'
                gain=np.sum(blackman)/wsize
            else : raise Exception('Unknown taper {0}'.format(tapering))
        elif isinstance(tapering,np.ndarray) : pass
        else :
            raise Exception('Bad value for tapering keyword')
        if not isinstance(tapering,np.ndarray) : exec('window='+which)
        else : window=tapering
        window = np.repeat(window,nn*nr).reshape((wsize,nn,nr))
    
        #Apply tapering on segmented data
        A=dum.copy()*window
        A=A.reshape(wsize,nr*nn) #Reshapa matrix
        nr=nn*nr
    else :
        if detrend :
            for i in np.arange(nr): A[:,i] = detrend_fun(np.arange(N),A[:,i]) if detrend else A[:,i]
        gain=1.0
    
    #Run transform
    ###############
    for i in np.arange(nr):
        spec=get_spec(dx, A[:,i],integration=integration,gain=gain)
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
    
    psd = psd * (gain**0.5)
#    print gain, np.sqrt(gain), gain **2, gain*0.5, gain/2.
#    esd=(np.sum(esd,axis=0))#/gain
#    psd=(np.sum(psd,axis=0))#/gain


    #Normalise by energy content    
    Scaling_Factor=len(fq)/esd.sum()
    if normalise :
        esd*=Scaling_Factor
        psd*=Scaling_Factor
    
    if tapering is not None : return {'params':{'tapering':tapering is not None,'which':which,'wsize':int(wsize),'nwind':int(nn),'overlap':int(100.*overlap),'gain':gain},'psd':psd,'esd':esd,'fq':fq,'p':p}
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


def preprocess(lat,lon,sla,N_min=None,per_min=15.0,max_gap=None,remove_edges=True,interp_over_continents=False,truncate_if_continents=True):
    sh=sla.shape
    nt=sh[0]
    nx=sh[1]
    dumsla=sla.copy()
    
    continent=False #This is a keyword to know if the pass went over a land mass
    
    #Remove profiles with less than 3 points
    ok=np.where(sla.mask.sum(axis=1) < (nx -3))[0]
    dumsla=dumsla[ok,:]
    nt=len(ok)
    
    for i in np.arange(nt):
        fg=~dumsla.mask[i,:]
        dst, dumlon, dumlat, dsla, lgaps, n, edges, int = grid_track(lat[fg], lon[fg], dumsla[i,:][fg],remove_edges=False,backbone=[lon,lat],interp_over_continents=interp_over_continents)
        if dsla.mask.sum() > 0:
            pass
        if (len(dumlon) > len(lon)) & (i == 0) :
            continent=True
            lendiff = len(dumlon) - len(lon)
            print '[WARNING] : Pass goes over a land mass, changing the track size from {0} to {1}'.format(nx,nx+lendiff)
            nx+=lendiff
#        dumsla[i,:]=dsla
        dumslaout=dsla.reshape((1,len(dsla))) if i == 0 else np.ma.concatenate([dumslaout,dsla.reshape((1,len(dsla)))],axis=0)
        if i == 0 : gaplen = [lgaps]
        else : gaplen.append(lgaps)
        ngaps = n if i == 0 else np.append(ngaps,n) 
    dumsla=dumslaout.copy()
        
    if max_gap is not None:
        
        #Remove profiles with long gaps
        gapmax=np.array([np.max(g) if len(g) > 0 else 0 for g in gaplen])
        id1 = np.where(gapmax <= max_gap)[0]
        dumsla=dumsla[id1,:]
        
        #Remove profiles with not enough coverage
#        per=100 * dumsla.mask.sum(axis=0) / np.float(nt)
        per = 100. * dumsla.mask.sum(axis=1)/np.float(nx)
        if N_min is None :
            N_min = np.round((0.01* (100- per_min)) * nx)
#        if per_min is None : per_min = 100* (1 - N_min/np.float(nx))
        
        id2 = np.where( per <= per_min)[0]
        dumsla=dumsla[id2,:]
        
        #At this point track edges are removed
        dumsla,id3=get_segment(dumsla,N_min,remove_edges=remove_edges,truncate_if_continents=truncate_if_continents)
        
        return dumsla, id1[id2[id3]]
        
    else :
        return dumsla, ngaps, gaplen

def get_segment(sla,N,last=True,mid=None,first=None,remove_edges=True,truncate_if_continents=True):
    
    #Set defaults
    if first is not None :
        last=None
        mid=None
    elif mid is not None :
        last=None
        first=None
    
    dumsla=sla.copy()
    nx=sla.shape[1]
    nt=sla.shape[0]
    
    #Save input mask
    dumsla.data[dumsla.mask]=dumsla.fill_value
    mask = np.ma.array(dumsla.mask.copy(),mask=np.zeros(sla.shape,dtype=bool))
    dumsla.mask[:]=False
    
    #Get edges
    if remove_edges : xid=np.ma.array(np.repeat(np.arange(nx),nt).reshape(nx,nt).transpose(),mask=mask.data)
    else : xid=np.ma.array(np.repeat(np.arange(nx),nt).reshape(nx,nt).transpose(),mask=np.zeros(sla.shape,dtype=bool))
    left=xid.min(axis=1)
    right=xid.max(axis=1)
    
    #Shift towards end
    if last :
        st=(right-N).astype(int) +1
        en=(right).astype(int) + 1
    elif mid :
        midpt=nx/2
        rlag=right-midpt
        llag=midpt-left
        odd = np.int(N)&1 and True or False
        if odd : nr=nl=np.int(N)/2
        else :
            nr=np.int(N)/2  + 1
            nl=np.int(N)/2
        for i,jk in enumerate(zip(*(llag,rlag))):
            j,k=jk  
        print 'THIS MUST BE COMPLETED!!!'
    elif first :
        st=(left).astype(int)
        en=(left+N).astype(int)
    
    if not remove_edges :
        st[st < 0] = 0
        en[en > nx] = nx
    
    for i in np.arange(nt) :
        dumsla.mask[i,:st[i]]=True
        dumsla.mask[i,en[i]:]=True
        mask.mask[i,:st[i]]=True
        mask.mask[i,en[i]:]=True
    
    #Update nt
    cycempty=dumsla.mask.sum(axis=1) == N
    ind=np.arange(nt)[~cycempty]
    nt=(~cycempty).sum()
    
    #Reform stuff
    dumsla=dumsla.compressed().reshape(nt,N)
    mask=mask.compressed().reshape(nt,N)
    
    if truncate_if_continents :
        empty=mask.sum(axis=0) == nt
        if empty.sum() > 0 :
            dumsla=dumsla[:,~empty]
            mask=mask[:,~empty]
            print '[WARNING] Points over land mass - removed {} pts'.format(empty.sum())
    
    return np.ma.array(dumsla,mask=mask), ind

def get_slope(fq,spec,degree=1,frange=None,threshold=0.):
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
    
    sh=spec.shape
    ndims=len(sh)
    
    if ndims == 1 :
        
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
    else :
        
        x = np.log10(fq[(fq > np.min(frange)) & (fq <= np.max(frange))]).flatten()
        y = np.log10(spec[(fq > np.min(frange)) & (fq <= np.max(frange)),:])
        nx = len(x)
        nt = sh[1]
        degree = 1
        
        out = []
        for i in np.arange(nt):
            (slope, intercept, rval, pval, err) = stats.linregress(x,y[:,i])
            flag = (y.mask[:,i].sum(dtype=float)/nx) <= threshold
            if i == 0 : out.append((slope,intercept,flag))
            else : out.append((slope,intercept,flag))

        return out

#+
# YULE_WALKER_REGRESSION : Estimation of an AR (autoregression) spectral model from data
#@todo: To be written
#-
def yule_walker_regression():
    pass