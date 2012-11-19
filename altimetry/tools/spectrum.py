import numpy as np
import alti_tools as atools

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
    
    #Get wavenumber
    k=np.fft.fftfreq(N)*N
    
    imx = (N-1)/2 if odd else N/2 #Maximum frequency
    
    #Convert into frequencies
    k/=L
    
    return k, L,imx


def get_spec(dx,var):
    """
    #+
    # GET_SPEC
    # @summary: Returns the spectrum of a regularly sampled dataset
    #
    # @param dq {type:numeric} : sampling interval (1D)
    # @param var {type:numeric} : data to analyse (1D).
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
    #-
    """
    
    N=var.size

    # part 1 - estimate of the wavenumbers
    ######################################
    [k,L,imx] = get_kx(N,dx)
#    print "dx = %s" %(dx)
#    print "kx = ",k[0:3]
#    print "shape",k.shape

    # part 2 - estimate of the spectrum
    # fast fourier transform
    ##########################
    fft=np.fft.fft(var)
    
#        Compute magnitude
#        CF.Wilks, 2006. Statistical Methods in the Atmopsheric Sciences. 2nd Edition. International Geophysics Series vol.96. Elsevier.
#        => Ck = [Ak^2 + Bk^2]^(1/2), where e^(iwt) = cos(wt) + i.sin(wt)
    a = fft.real
    b = fft.imag
    c = a**2 + b**2
    
    #Checking parseval's Theorem
    #-->                     (c/N).sum() == (var**2).sum()
    #OR : np.linalg.norm(fft)/np.sqrt(N) == np.linalg.norm(var)
    
    #Get phase (not used yet)
    phase=np.angle(fft)

    # integration of the spectrum
    # shift to have values centered on the intervals
    dk = k[1] - k[0]
    dk_half = dk/2  # half of the interval
    
    ist=1
    
    k = k[ist:imx]  #Remove negative frequencies and zero frequency (mean)
    k_= k #k + dk  #Shift wavelengths by half the unit
    
    esd = c[ist:imx]
    
    
    #Spectral integration
    var = k_*0.0
    for i in np.arange(len(k_)):
#        mesd[i] = np.mean(c[(k > (k_[i]-dk)) & (k < (k_[i]+dk))])
        var[i] = 2*dk*np.sum(c[(k > (k_[i]-dk)) & (k < (k_[i]+dk))]) #This is variance units integration

    #Get frequencies and period  
    fq=k_
    p=1/fq
    
    psd = esd / fq

#from http://www.mathworks.com/matlabcentral/newsreader/view_thread/283145

#N = 1024 % 1024
#dt = 0.001 % 1e-3
#Fs = 1/dt % 1000 fft period
#T = N*dt % 1.024 ifft period
#
#t = 0:dt:T-dt;
#
#f0 = 30 % 30
#T0 = 1/f0 % 0.033333 signal period
#Nc = T/T0 % 30.72 noninteger number of cycles
#y = sin(2*pi*f0*t);
#
#% Power and Energy in the Time Domain
#
#p = y.^2 ; % Instantaneous Power
#P = sum(p) % 511.04 Total Power
#Pav = mean(p) % 0.49906 Average Power
#e = dt*p; % Instantaneous Energy
#Et = sum(e) % 0.51104 Total Energy
#Etav = mean(e) % 0.00049906 Average Energy
#
#http://groups.google.com/group/comp.soft-sys.matlab/
#msg/009339c581e63467?hl=en
#
#% Power and Energy in the Freqency Domain
#
#df = Fs/N % 0.97656 = 1/T
#f = 0:df:Fs-df; % Unipolar freq interval
#fb = f - df*ceil((N-1)/2); % Bipolar freq interval
#Y = fft(y); % Unipolar freq spectrum
#Yb = fftshift(Y); % Bipolar freq spectrum
#
#% Note Y and Yb are interchangeable in most of the following
#
#absYb = abs(Yb);
#PSDb = absYb.^2/N; % Power Spectral Density
#
#figure
#plot(fb,PSDb)
#axis([-Fs/2 Fs/2 0 1.1*max(PSDb)])
#title('Power Spectrum')
#
#% Check Parseval's Theorem (See Wikipedia)
#
#check = sum(p)-sum(PSDb) % -4.5475e-013
#
#You can determine the appropriate expressions
#for the frequency domain calculations of
#P, Pav, Ef and Efav
#
#http://groups.google.com/group/comp.soft-sys.matlab/
#msg/2222327db2ea7f51
#
#
#Hope this helps.
#
#Greg 
    # Normalisation (Danioux 2011)
#    specvar = specvar / (var.size)**2 / dk #No normalization!!

    #Power spectral density

    # Dimensionalization to be able to compare with the litterature
#    nbwave=k*1000.0/2.0 #/np.pi     # from rad/m to km/m
#    dk=dk*1000.0/2.0 #/np.pi
#    if fft: specvar=specvar*2.0*np.pi/1000.0

    return {'psd':psd,'esd':esd,'var':var,'fq':fq,'p':p}


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


def grid_track(dst,lat,lon,sla,remove_edges=True,backbone=None):
    """
    # GRID_TRACK
    # @summary: This function allow detecting gaps in a set of altimetry data and rebin this data regularlyy, with informations on gaps.
    # @param dst {type:numeric} : along-track distance.
    # @param lat {type:numeric} : latitude
    # @param lon {type:numeric} : longitude
    # @param sla {type:numeric} : data
    # @return:
    #    outdst : resampled distance
    #    outlon : resampled longitude
    #    outlat : resampled latitude
    #    outsla : resampled data
    #    gaplen : length of the longest gap in data
    #    ngaps : number of detected gaps in data
    #    dx : average spatial sampling
    #    interpolated : True when data was interpolated (empty bin)
    #
    # @author: Renaud DUSSURGET (RD) - LER/PAC, Ifremer
    # @change: Created by RD, July 2012
    #    29/08/2012 : Major change -> number of output variables changes (added INTERPOLATED), and rebinning modified
    """
    
    #Find gaps in data
    dx = dst[1:] - dst[:-1]
    mn_dx = np.median(dx)
    bins = np.ceil(dst.max() / mn_dx) + 1
    range=(0/2.,mn_dx * bins) - mn_dx/2
    hist,bin_edges=np.histogram(dst, bins=bins, range=range) #We have binned the data along a regular grid of size (bins) in the range (range)
                                                             #Missing data is thus represented by no data in a given bin
    
    #Remove leading and trailing edges (and synchronise bin_edges)
    if remove_edges == True :
        while hist[0] == 0 : 
            hist=np.delete(hist,[0])
            bin_edges=np.delete(bin_edges,[0])
        while hist[-1] == 0 :
            hist=np.delete(hist,[len(hist)-1])
            bin_edges=np.delete(bin_edges,[len(bin_edges)-1])
    
    #Get filled bins indices

    ok = np.arange(len(hist)).compress(np.logical_and(hist,True or False))
    empty = np.arange(len(hist)).compress(~np.logical_and(hist,True or False)) 
    
    outsla = np.repeat(np.NaN,len(hist))
    outlon = np.repeat(np.NaN,len(hist))
    outlat = np.repeat(np.NaN,len(hist))
    outdst = bin_edges [:-1]+ mn_dx/2 #distances is taken at bins centers
    outsla[ok] = sla
    outlon[ok] = lon
    outlat[ok] = lat
    
    #Fill the gaps if there are some
    if len(empty) > 0 : 
        #Interpolate lon,lat @ empty positions
        outlon[empty] = atools.interp1d(ok, outlon[ok], empty, kind='cubic')
        outlat[empty] = atools.interp1d(ok, outlat[ok], empty, kind='cubic')
        outsla[empty] = atools.interp1d(ok, outsla[ok], empty, spline=True)
        
    
    
    #Get gap properties
    ind=np.arange(len(hist))
    dhist=(hist[1:] - hist[:-1])
    st=ind.compress(dhist==-1)+1
    en=ind.compress(dhist==1)
    gaplen=(en-st) + 1
    ngaps=len(st)
    
    #Get empty bin flag
    interpolated=~hist.astype('bool')
    
    return outdst, outlon, outlat, outsla, gaplen, ngaps, dx, interpolated

def grid_track_backbone(lat,lon,sla,backlat,backlon,fill=None):
    """
    # GRID_TRACK_BACKBONE
    # @summary: This function allow detecting gaps in a set of altimetry data and rebin this data regularlyy, with informations on gaps.
    # @param dst {type:numeric} : along-track distance.
    # @param lat {type:numeric} : latitude
    # @param lon {type:numeric} : longitude
    # @param sla {type:numeric} : data
    # @return:
    #    outdst : resampled distance
    #    outlon : resampled longitude
    #    outlat : resampled latitude
    #    outsla : resampled data
    #    gaplen : length of the longest gap in data
    #    ngaps : number of detected gaps in data
    #    dx : average spatial sampling
    #    interpolated : True when data was interpolated (empty bin)
    #
    # @author: Renaud DUSSURGET (RD) - LER/PAC, Ifremer
    # @change: Created by RD, July 2012
    #    29/08/2012 : Major change -> number of output variables changes (added INTERPOLATED), and rebinning modified
    """
    
    #Find gaps in data
    dst=atools.calcul_distance(backlat[0],backlon[0],lat,lon)
    dstback=atools.calcul_distance(backlat,backlon)
    dx = dstback[1:] - dstback[:-1]
    mn_dx = np.median(dx)
    
    bins = np.ceil(dstback.max() / mn_dx) + 1
    range=(0/2.,mn_dx * bins) - mn_dx/2
    hist,bin_edges=np.histogram(dst, bins=bins, range=range) #We have binned the data along a regular grid of size (bins) in the range (range)
                                                             #Missing data is thus represented by no data in a given bin
    
    #Get filled bins indices
    ok = np.arange(len(hist)).compress(np.logical_and(hist,True or False))
    empty = np.arange(len(hist)).compress(~np.logical_and(hist,True or False)) 
    
    outsla = np.repeat(np.NaN,len(hist))
    outlon = np.repeat(np.NaN,len(hist))
    outlat = np.repeat(np.NaN,len(hist))
    outdst = bin_edges [:-1]+ mn_dx/2 #distances is taken at bins centers
    outsla[ok] = sla
    outlon[ok] = lon
    outlat[ok] = lat
    
    #Fill the gaps if there are some
    if (fill is not None) & (len(empty) > fill) : 
        #Interpolate lon,lat @ empty positions
        outlon[empty] = atools.interp1d(ok, outlon[ok], empty, kind='cubic')
        outlat[empty] = atools.interp1d(ok, outlat[ok], empty, kind='cubic')
        outsla[empty] = atools.interp1d(ok, outsla[ok], empty, spline=True)
    
    #Get empty bin flag
    interpolated=(~hist.astype('bool'))
    
    return outdst, outlon, outlat, outsla, dx, interpolated

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