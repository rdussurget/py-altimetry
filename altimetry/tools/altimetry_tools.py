# -*- coding: utf-8 -*-
import numpy as np

try:
 import seawater.gibbs
except ImportError, e:
 pass # module doesn't exist, deal with it.
 
from altimetry.tools import deriv, calcul_distance, interp1d, loess
if __debug__ : import matplotlib.pyplot as plt

def track_orient(x,y,orient=False):
#    ;Calculate track orientation (on 3 points lagrangian interpolation - see DERIV help page)
    
    getorient=orient
    
    dxy=deriv(x,y)
    
    dx=np.median(deriv(x))
    dy=np.median(deriv(y))

    orient = np.arctan(dxy) if dx > 0 else np.arctan(dxy) + np.pi #Rotate backward (180�)
                                                                  #orientation when track goes
                                                                  #westward
  
    orient = np.remainder(orient,2.0*np.pi)
#    if (medorient < 0) : orient += 2.0*np.pi
#    if (medorient > 2.0*np.pi) THEN orient -= 2D*!dpi
    #Sense is positive for ascending tracks
    sense = True if np.median(orient) < np.pi else False
  
    if getorient : return sense, orient
    else : return sense

def genweights(p,q,dt,error=None , n=False):
    """
    ;+
    ;   GENWEIGTHS : return optimal weigthing coefficients from Powell and Leben 2004<br />
    ;   translated from Matlab genweigths.m program to IDL<br /><br />
    ;  
    ;   Reference : Powell, B. S., et R. R. Leben (2004), An Optimal Filter for <br />
    ;   Geostrophic Mesoscale Currents from Along-Track Satellite Altimetry, <br />
    ;   Journal of Atmospheric and Oceanic Technology, 21(10), 1633-1642.
    ;   
    ; @author Renaud DUSSURGET, LEGOS/CTOH
    ; @history Created Sep. 2009 from genweights.m (Brian Powell (c) 2004, <br />
    ;   University of Colorado, Boulder)<br />
    ;-
    """
    
    p = np.abs(p)
    q = np.abs(q)
  
    #check inputs
    if (-p > q) : raise "genweights : P must be lesser than q"
  
    #Build matrices
    N = p + q
    T = N + 1
    A = np.matrix(np.zeros((T,T)))
    A[T-1,:] = np.append(np.repeat(1.0,N),0)
    sn = np.arange(T) - p
    sn = sn.compress(sn != 0)
    for i in np.arange(len(sn)) :
        A[i,:] = np.append(((1./sn)*(-sn[i]/2.)),sn[i]**2.*dt**2./4.) #Eq.11 (PL) 
        A[i,i] = -1.
        
    B = np.zeros(T)
    B[N] = 1.0
  
    #Compute the coefficients
    cn=np.dot(A.I,B)##B
#    cn=cn.transpose()
    cn=np.array([i for i in cn.flat])
    cn = cn[0:N] #Check the indices
  
    #Compute the error
    error = np.sqrt(np.sum(cn.transpose()/(sn*dt))**2. + np.sum( (cn.transpose()/(sn*dt))**2. ) );
  
    return cn, sn if n else cn

def optimal_slope(cn,n,dt,z,i):
    """
    ;+
    ;
    ;   OPTIMAL_SLOPE : Function to compute the slope using a Powell et Leben (2004) <br />
    ;   filter with given characteristics.<br /><br />
    ;   
    ;   Reference : Powell, B. S., et R. R. Leben (2004), An Optimal Filter for <br />
    ;   Geostrophic Mesoscale Currents from Along-Track Satellite Altimetry, <br />
    ;   Journal of Atmospheric and Oceanic Technology, 21(10), 1633-1642.
    ;
    ; @param cn {in}{required}{type:NUMERIC} Filter coefficients
    ; @param n {in}{required}{type:NUMERIC} Filter indices reffered to central point<br />
    ;          (from -p to q - 0 is the center)
    ; @param dt {in}{required}{type:NUMERIC} Mean along-track grid step
    ; @param z {in}{required}{type:NUMERIC}{default:12} Along-track data of height<br />
    ;          anomalies (all dataset for current track)
    ; @param i {in}{required}{type:NUMERIC}{default:12} Index of points to use with<br />
    ;          current computation wihtin track.
    ;          
    ; @returns Optimal slope. Unit depends on distance and height anomalies units.
    ;
    ;
    ;
    ; @author Renaud DUSSURGET, LEGOS/CTOH
    ; @history Created Sep. 2009 from genweights.m (Brian Powell (c) 2004, <br />
    ;   University of Colorado, Boulder)
    ;
    ;
    ; @example dh = optimal_slope(cn,n,dt,z,id) : Compute the Optimal Slope using <br />
    ;   points in z[id], with a dt step, coefficients cn (of indices n).
    ;
    ;-
    """

#Return optimal slope
    #FIX : This version uses a regular grid (uses n * dt)
    #      We should rather use an irregular grid (real array of along-track distances)
  
    #; 1st member EQ 0 if filter centered
    dh = np.nansum((-cn)/(n * dt)) * z[i] \
        + np.nansum(cn / (n * dt) * z[i+n])
    return dh

#FUNCTION powell_leben_filter_KM, lon, lat, z, p=p, q=q, verbose=verbose
def powell_leben_filter_km(*args,**kwargs):
    """
    ;+
    ;
    ;   POWEL_LEBEN_FILTER_KM : Compute geostrophic speeds from a sea level dataset <br />
    ;   using a Powell et Leben (2004) along-track filtering to maintain measurement<br />
    ;    noise constant.<br /><br />
    ;
    ;   Reference : Powell, B. S., et R. R. Leben (2004), An Optimal Filter for <br />
    ;   Geostrophic Mesoscale Currents from Along-Track Satellite Altimetry, <br />
    ;   Journal of Atmospheric and Oceanic Technology, 21(10), 1633-1642.
    ;
    ; @param lon {in}{optional}{type:NUMERIC} longitude in degrees
    ; @param lat {in}{optional}{type:NUMERIC} latitude in degrees
    ; @param z {in}{required}{type:NUMERIC} sea level surface. Can either be absolute<br />
    ;          values (SSH) or relative values (SLA). This MUST be given in METERS.
    ; @param p {in}{optional}{type:NUMERIC}{default:12} Filter half-width BEFORE<br />
    ;          center of filtering window in KILOMETERS.
    ; @param q {in}{optional}{type:NUMERIC}{default:12} Filter half-width AFTER <br />
    ;          center of filtering window in KILOMETERS.
    ;          
    ; @returns Geostrophic velocity component, positive eastward
    ;
    ;
    ;
    ; @author Renaud DUSSURGET, LEGOS/CTOH
    ; @history Created Sep. 2009 from genweights.m (Brian Powell (c) 2004, <br />
    ;   University of Colorado, Boulder)<br />
    ;   Modified May 2010 to be compliant with 20Hz datasets (p & n can vary).<br />
    ;     Warining may also be issued for data with holes within the width of the<br />
    ;     window.<br />
    ;   Modified June 2010 to include filtering window width in KM instead of nb. of<br />
    ;     points (Equivalent btw. 1Hz and 20Hz data).<br />
    ;
    ; @uses CALCUL_DISTANCE, EXIST, GENWEIGTHS, SETINTERSECTION, SETUNION, <br />
    ;   OPTIMAL_SLOPE, GRAVITY, CORIOLIS, TRACK_ORIENT
    ;
    ; @example dummy1=powell_leben_filter_KM(lon,lat,sla,p=11,q=11) :<br />
    ;   Return along-track velocity anomalies using a 11km by 11km filter window.
    ;
    ;-
    """
    
    
    #From Powell and Leben, 2004. JAOT.
    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    lon = args[0]
    lat = args[1]
    nu = args[2]
    dst = calcul_distance(lat,lon) * 1e3 #distance in meters
    
    n=nu.size

    #We use a p+q filter (default 2 1Hz points on each side ~ 12km)
    #-> Set p & q distances to defaults if they don't exist    
    p = kwargs.pop('p',12.)
    q = kwargs.pop('q',12.)
    verbose = kwargs.pop('verbose',False)
        
    dt = np.median(dst[1:] - dst[:-1])   #Median point-to-point distance in meters
      
    #Get the corresponding number of points rounded to the nearest integer
    pn=np.round((p*1e3)/dt).astype(int)
    qn=np.round((q*1e3)/dt).astype(int)
    
    if pn + qn > n : raise 'Filtering window is too large wrt array length'
     
    #Compute weights
    #;;;;;;;;;;;;;;;
    
    #Initialise all possible configurations of filters
    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    cn = []#=np.array([])   #array of weigth arrays
    nout= [] #np.array([]) #array of point indices in filter arrays
    pout=np.ndarray(pn+qn+1) #array containg number of points before any filter
    qout=np.ndarray(pn+qn+1) #array containg number of points after any filter
      
      
    #Compute each filter configuration
    #1)Start with increasing number of points before point
    #2)Compute full resolution filter (p & q points before and after filter)
    #3)Finish with a deacreasing number of points after filter
    
    #1)Weights before
    for i in np.arange(pn):
        w,n=genweights(i,qn,dt,n=True)
        cn.append(w)
        nout.append(n)
        pout[i]=i
        qout[i]=qn
        
    #2)Centered weights
    w,n=genweights(pn,qn,dt,n=True)
    cn.append(w)
    nout.append(n)
    pout[pn]=pn
    qout[pn]=qn
        
    #3)Weights after
    for i in np.arange(pn+1, pn+1+qn):
        w,n=genweights(pn,pn+qn-i,dt,n=True)
        cn.append(w)
        nout.append(n)
        pout[i]=pn
        qout[i]=pn+qn-i
        
    #find track edges and gaps
    #;;;;;;;;;;;;;;;;;;;;;;;;;
    empty = np.where(np.isnan(nu))[0]
    ok = np.where(~np.isnan(nu))[0]
    emptycnt = empty.size
    okcnt = ok.size
    
    st = np.min(ok) #start point
    en = np.max(ok) #end of track
        
    #Initiate dh for slope calculation, Eq. (5)
    dh = np.repeat(np.nan,nu.size)
      
    #Loop over valid points in track
    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    for i in np.arange(okcnt):
        
        #find valid points within filter width around filter center
        #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        cur = np.arange(pn+qn+1)+ok[i]-pn               #Points around filter center +-p/q points
        cur = np.sort(list(set(cur).intersection(set(ok))))  #intersection of filter window with valid points
        curcnt = cur.size
        
        #check points currently selected
        #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        if (curcnt > 1) :
            #Be careful -> This is intended to select points following directly central point regularly (without gaps within window)
            #We'll see that for a hi-rate case, we'll be able to accept gaps within filter window
            #FIX:this might be buggy...
            a = np.where(np.append(0,cur[1:] - cur[:-1]) == 1)[0]  #select following valid points
            b = np.where(np.append(cur[:-1] - cur[1:],0) ==  -1)[0] #select preceeding valid points
            aftcnt = a.size
            befcnt = b.size
            
#          b = WHERE([cur[:-1] - cur[1:],0] EQ -1,befcnt) ; select preceeding valid points
            #If there are points at least before OR after central point, keep going, otherwise switch to next point
            if (aftcnt != 0) & (befcnt != 0):
                cur = cur[np.sort(list(set(a).union(b)))] #Get the union of valid points before and after central point
                curcnt = cur.size
            else : curcnt =0
        else : curcnt =0
        
        #Adapt filter length to current case
        #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        #(We know how many points are following directly before and after)
        
        if (curcnt != 0) :
            
            #Get number of points before and after central point
            #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            nbef = ok[i] - cur[0]
            naft = cur[curcnt - 1] - ok[i]
          
            #If we have less than p or q on each side, we'll set full filter width to
            #the greater size, and reduce the smaller size
            #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            if (naft > nbef) & (naft+nbef > qn) : naft=pn
            if (naft < nbef) & (naft+nbef > pn) : nbef=qn
          
            #Select adapted filter in filter list
            #(given by pout and qout)
            #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            which=(np.where((nbef == pout) & (naft == qout)))[0]
            nwhich = which.size
          
            #If we have found a filter, we use it, otherwise we issue a warning message and flag the point
            if (nwhich > 0) : dh[ok[i]] = optimal_slope(cn[which],nout[which],dt,nu,ok[i])
            else :
                if (verbose) : print "[WARNING] No available filter for point no.{0}".format(ok[i])
                dh[ok[i]] = np.NaN #If no suitable filter, flag the data
       
    #Compute currents with given slope
    #;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    #Compute gravitational & coriolis forces
    g = gravity(lat)
    f = coriolis(lat)

    #Get geostrophic velocities
    try :
        ug = - (g*dh) / (f)
    except ValueError :
        print 'error'
  
    #Inverse sign of ug for descending tracks as Coriolis is oriented to the right
    #northward
    if (not track_orient(lon,lat)) : #descending tracks
        ug *=-1
    
    return ug


def coriolis(lat):

    sideral_day = 86164.1 #sideral day in seconds
    omega = (2. * np.pi) / sideral_day #earth rotation rate
  
    #Compure Coriolis frequency and returns it
    return 2. * omega * np.sin(np.deg2rad(lat))

def gravity(*args):
    try :
        return seawater.gibbs.grav(*args)
    except :
        return np.repeat(6.67384,len(args[0]))

def geost_1d(*args,**kwargs) : #(lon,lat,nu): OR (dst,nu)
    """
    ;+
    ;
    ;   GEOST_1D : Compute geostrophic speeds from a sea level dataset <br />
    ;
    ;   Reference : Powell, B. S., et R. R. Leben (2004), An Optimal Filter for <br />
    ;   Geostrophic Mesoscale Currents from Along-Track Satellite Altimetry, <br />
    ;   Journal of Atmospheric and Oceanic Technology, 21(10), 1633-1642.
    ;
    ; @param lon {in}{optional}{type:NUMERIC} longitude in degrees
    ; @param lat {in}{optional}{type:NUMERIC} latitude in degrees
    ; @param dst {in}{optional}{type:NUMERIC} along-track distance.
    ; @param z {in}{required}{type:NUMERIC} sea level surface. Can either be absolute<br />
    ;          values (SSH) or relative values (SLA). This MUST be given in METERS.
    ; @keyword strict {in}{optional}{type:BOOLEAN} If True, compute gradient at mid-distance.
    ; @keyword pl04 {in}{optional}{type:BOOLEAN} If True, use the Powell & Leben 2004 method.
    ;          
    ; @returns Geostrophic velocity component, positive eastward
    ;
    ;
    ;
    ; @author Renaud DUSSURGET, LEGOS/CTOH
    ; @history Created Sep. 2009 from genweights.m (Brian Powell (c) 2004, <br />
    ;   University of Colorado, Boulder)<br />
    ;   Modified May 2010 to be compliant with 20Hz datasets (p & n can vary).<br />
    ;     Warining may also be issued for data with holes within the width of the<br />
    ;     window.<br />
    ;   Modified June 2010 to include filtering window width in KM instead of nb. of<br />
    ;     points (Equivalent btw. 1Hz and 20Hz data).<br />
    ;
    ; @uses CALCUL_DISTANCE, EXIST, GENWEIGTHS, SETINTERSECTION, SETUNION, <br />
    ;   OPTIMAL_SLOPE, GRAVITY, CORIOLIS, TRACK_ORIENT
    ;
    ; @example dummy1=geost_1D(lon,lat,sla,pl04=True,p=11,q=11) :<br />
    ;   Return along-track velocity anomalies using a 11km by 11km Powell & Leben 2004 filter window <br />
    ;          dummy2=geost_1D(dst,sla,strict=True) :<br />
    ;   Return along-track velocity anomalies computed at mid-distance <br />
    ;
    ;-
    """
    
    lon = args[0]
    lat = args[1]
    
    dst = args[2] if len(args) == 4 else calcul_distance(lat,lon) * 1e3 #distance in meters
    nu = args [3] if len(args) == 4 else args[2]
    
    isVector = len(np.shape(nu)) == 1 
    
    #Reshape nu if vector
    if isVector : nu=np.reshape(nu,(len(nu),1))
    
    nt = np.shape(nu)[1] if not isVector else 1
    sh = nu.shape
    nufilt=np.ma.array(np.empty(sh),mask=True,dtype=nu.dtype)
    
    pl04 = kwargs.pop('pl04',False)
    filter = kwargs.pop('filter', None)
    strict = kwargs.pop('strict',False)
    verbose = kwargs.pop('verbose',False)
    
    if filter is not None :
        for t in np.arange(nt) : 
            nufilt[:,t] =loess(nu[:,t],dst,filter*1e3)
        nu=nufilt
    
    if pl04 :
        ug = np.ma.array(np.empty(sh),mask=True,dtype=nu.dtype)
        for t in np.arange(nt) :
            ug[:,t] = powell_leben_filter_km(lon,lat,nu[:,t],verbose=verbose,**kwargs)
        if isVector : ug=ug.flatten()
        return ug
    
    #If strict option is set to True, compute gradients at mid-distance between points
    if strict :
        lon = (lon[1:] - lon[:-1])/2. + lon[0:-1]
        lat = (lat[1:] - lat[:-1])/2. + lat[0:-1]
    
    
    #Compute gravitational & coriolis forces
    if strict : sh = (sh[0]-1,sh[1])
    g = np.repeat(gravity(lat),nt).reshape(sh)
    f = np.repeat(coriolis(lat),nt).reshape(sh)
    
    #Compute SSH 1st derivative
#    dh = deriv(dst,nu) #(deriv is very bad...)
    
    
    dh = np.ma.array(np.empty(sh),mask=True,dtype=nu.dtype)
    for t in np.arange(nt) :
        dh[:,t] = (nu[1:,t] - nu[:-1,t])/(dst[1:] - dst[:-1]) if strict else deriv(dst,nu[:,t])
      
    #Compute geostrophy
#    print f
#    print g
#    print dh
    ug = - (g*dh) / (f)
  
    #Inverse sign of ug for descending tracks as Coriolis is oriented to the right
    #northward
    if (not track_orient(lon,lat)) : #descending tracks
        ug *=-1
    
    if isVector : ug=ug.flatten()
    
    return (lon,lat,ug) if strict else ug


    
def grid_track(lat,lon,sla,remove_edges=None,backbone=None,interp_over_continents=True):
    """
    # GRID_TRACK
    # @summary: This function allow detecting gaps in a set of altimetry data and rebin this data regularly, with informations on gaps.
    # @param lat {type:numeric} : latitude
    # @param lon {type:numeric} : longitude
    # @param sla {type:numeric} : data
    # @return:
    #    outdst : resampled distance
    #    outlon : resampled longitude
    #    outlat : resampled latitude
    #    outsla : resampled data
    #    gaplen : length of the longest gap in data
    #    dx : average spatial sampling
    #    interpolated : True when data was interpolated (empty bin)
    #
    # @author: Renaud DUSSURGET (RD) - LER/PAC, Ifremer
    # @change: Created by RD, July 2012
    #    29/08/2012 : Major change -> number of output variables changes (added INTERPOLATED), and rebinning modified
    #    06/11/2012 : Included in alti_tools lib
    #    19/12/2012 : Added backbone option (reproject along the backbone grid)
    """
    
    
    
    #Find gaps in data
    if backbone is not None :
        backlon=backbone[0]
        backlat=backbone[1]
        ascending=track_orient(lon,lat)
        dst=calcul_distance(backlat[0],backlon[0],lat,lon)
        if ascending : dst[lat < backlat[0]]*=-1
        if not ascending : dst[lat > backlat[0]]*=-1
        dstback=calcul_distance(backlat,backlon)
        dx = dstback[1:] - dstback[:-1]
        mn_dx = np.median(dx)
        bins = np.round(dstback.max() / mn_dx)+1
        range=(0/2.,mn_dx * bins) - mn_dx/2
#        dfback=list(set(dstback).difference(set(dst)))
        bhist,bbin_edges=np.histogram(dstback, bins=bins, range=range)
        continent=np.where(bhist==0)[0]
        if remove_edges is None : remove_edges=False
    else :
        dst=calcul_distance(lat,lon)
        #Find gaps in data
        dx = dst[1:] - dst[:-1]
        mn_dx = np.median(dx)
        bins = np.ceil(dst.max() / mn_dx) + 1
        range=(0/2.,mn_dx * bins) - mn_dx/2
        
        if remove_edges is None : remove_edges=True
    
    
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
    
    nH =len(hist)
    
    #Get filled bins indices
    ok = np.arange(len(hist)).compress(np.logical_and(hist,True or False))
    empty = np.arange(len(hist)).compress(~np.logical_and(hist,True or False)) 
    
    if isinstance(sla,np.ma.masked_array) : outsla = np.ma.masked_array(np.repeat(sla.fill_value,nH),mask=np.ones(nH,dtype=bool),dtype=sla.dtype)
    else : outsla = np.ma.masked_array(np.repeat(np.ma.default_fill_value(1.0),nH),mask=np.ones(nH,dtype=bool),dtype=np.float32)
    if isinstance(sla,np.ma.masked_array) : outlon = np.ma.masked_array(np.repeat(lon.fill_value,nH),mask=np.ones(nH,dtype=bool),dtype=lon.dtype)
    else : outlon = np.ma.masked_array(np.repeat(np.ma.default_fill_value(1.0),nH),mask=np.ones(nH,dtype=bool),dtype=np.float32)
    if isinstance(sla,np.ma.masked_array) : outlat = np.ma.masked_array(np.repeat(lat.fill_value,nH),mask=np.ones(nH,dtype=bool),dtype=lat.dtype)
    else : outlat = np.ma.masked_array(np.repeat(np.ma.default_fill_value(1.0),nH),mask=np.ones(nH,dtype=bool),dtype=np.float32)
    
    outdst = bin_edges [:-1]+ mn_dx/2 #distances is taken at bins centers
    outsla[ok] = sla
    outlon[ok] = lon
    outlat[ok] = lat
    
    #Remove land mass point if asked
    if not interp_over_continents :
        sempty=np.sort(np.array(list(set(empty).difference(set(continent)))))
    else : sempty=empty.copy()
    
    #Fill the gaps if there are some
    if len(empty) > 0 : 
        #Interpolate lon,lat @ empty positions
        outlon[empty] = interp1d(ok, outlon[ok], empty, kind='cubic', fill_value=lon.fill_value)
        outlat[empty] = interp1d(ok, outlat[ok], empty, kind='cubic', fill_value=lat.fill_value)
    if len(sempty) > 0 :outsla[sempty] = interp1d(ok, outsla[ok], empty, kind=0, fill_value=sla.fill_value) #0-th order Spline interpolation
    
    outlon.mask[outlon.data == outlon.fill_value] = outlon.fill_value
    outlat.mask[outlat.data == outlat.fill_value] = outlat.fill_value
    outsla.mask[outsla.data == outsla.fill_value] = outsla.fill_value
    
    
    #Get gap properties
    ind=np.arange(len(hist))
    dhist=(hist[1:] - hist[:-1])
    
    #Remove edge gaps
    if (dhist!=0).sum()> 0 :
        if (dhist[dhist!=0])[0] == 1 : dhist[(np.arange(nH)[dhist!=0])[0]]=0 #This do not count start and end of track as gaps
    if (dhist!=0).sum()> 0 :
        if (dhist[dhist!=0])[-1] == -1 : dhist[(np.arange(nH)[dhist!=0])[-1]]=0
    st=ind[dhist==-1]+1
    en=ind[dhist==1]
    gaplen=(en-st) + 1
    ngaps=len(st)
    gapedges=np.array([st,en])
    
    #Get empty bin flag
    interpolated=~hist.astype('bool')
    
    return outdst, outlon, outlat, outsla, gaplen, ngaps, gapedges, interpolated

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
    dst=calcul_distance(backlat[0],backlon[0],lat,lon)
    dstback=calcul_distance(backlat,backlon)
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
        outlon[empty] = interp1d(ok, outlon[ok], empty, kind='cubic')
        outlat[empty] = interp1d(ok, outlat[ok], empty, kind='cubic')
        outsla[empty] = interp1d(ok, outsla[ok], empty, spline=True)
    
    #Get empty bin flag
    interpolated=(~hist.astype('bool'))
    
    return outdst, outlon, outlat, outsla, dx, interpolated


def grid_time(time,remove_edges=False):
    """
    # GRID_TIME
    # @summary: This function allow to regularly resample time as an array in an altimetry dataset
    # @param time {type:numeric} : time
    # @return:
    #    outtime : resampled time
    #    gaplen : length of the longest gap in data
    #    ngaps : number of detected gaps in data
    #    dx : average spatial sampling
    #    interpolated : True when data was interpolated (empty bin)
    #
    # @author: Renaud DUSSURGET (RD) - LER/PAC, Ifremer
    # @change: Created by RD, July 2012
    #    29/08/2012 : Major change -> number of output variables changes (added INTERPOLATED), and rebinning modified
    #    06/11/2012 : Included in alti_tools lib
    """
    
    #Find gaps in data
    nt=time.shape[0]
    nx=time.shape[1]

    #Get series with the longest record
    xnt=np.isfinite(time).sum(axis=0)
    xid = np.argmax(xnt)
    
    #Make a linear regression of time over each cycle
    w = np.linalg.lstsq(np.array([np.arange(nt)[~time.mask[:,xid]],np.ones(nt)[~time.mask[:,xid]]]).T,time[:,xid][~time.mask[:,xid]])[0]
    t = w[0]*np.arange(nt)+w[1]
    
#    dft = t[1:] - t[:-1]
#    mn_dt = np.median(dft)
#    bins = np.ceil(t.max() / mn_dt) + 1
#    range=(0/2.,mn_dt * bins) - mn_dt/2
#    hist,bin_edges=np.histogram(time, bins=bins, range=range) #We have binned the data along a regular grid of size (bins) in the range (range)
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
#    outtime = np.repeat(np.NaN,len(hist))
##    outtime = bin_edges [:-1]+ mn_dx/2 #distances is taken at bins centers
#    outtime[ok] = time
#    
#    #Fill the gaps if there are some
#    if len(empty) > 0 : 
#        #Interpolate lon,lat @ empty positions
#        outtime[empty] = interp1d(ok, outtime[ok], empty, kind='cubic')
#        
#    #Get gap properties
#    ind=np.arange(len(hist))
#    dhist=(hist[1:] - hist[:-1])
#    st=ind.compress(dhist==-1)+1
#    en=ind.compress(dhist==1)
#    gaplen=(en-st) + 1
#    ngaps=len(st)
#    gapedges=np.array([st,en])
#    
#    #Get empty bin flag
#    interpolated=~hist.astype('bool')
#    
#    return outtime, gaplen, ngaps, gapedges, interpolated
    return t

def fill_gaps(lat,lon,sla,mask,remove_edges=False):
    """
    # FILL_GAPS
    # @summary: This function allow interpolating data in gaps, depending on gap size. Data must be regularly gridded
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
    #    06/11/2012 : Included in alti_tools lib
    """
    
    dst=calcul_distance(lat,lon)
    
    #Find gaps in data
    dx = dst[1:] - dst[:-1]
    mn_dx = np.median(dx)
    nx=len(sla)

    flag=~mask
    
    #Get filled bins indices
    outsla = sla.copy()
    outlon = lon.copy()
    outlat = lat.copy()
    outind = np.arange(nx)
    
    #Replace missing data on edges by the latest valid point
    first=np.where((flag))[0].min()
    last=np.where((flag))[0].max()
    if remove_edges :
        outsla=outsla[first:last+1]
        outlon=outlon[first:last+1]
        outlat=outlat[first:last+1]
        outind=outind[first:last+1]
        mask=mask[first:last+1]
        flag=flag[first:last+1]
    else:
        outsla[0:first] = outsla[first]
        outsla[last:] = outsla[last]
    
    #Get gap properties
    hist=np.ones(nx,dtype=int)
    hist[outsla.mask]=0
    while hist[0] == 0 : 
        hist=np.delete(hist,[0])
    while hist[-1] == 0 :
        hist=np.delete(hist,[len(hist)-1])
    
    ind=np.arange(len(hist))
    dhist=(hist[1:] - hist[:-1])
    st=ind.compress(dhist==-1)+1
    en=ind.compress(dhist==1)
    gaplen=(en-st) + 1
    ngaps=len(st)
    gapedges=np.array([st,en])
    
    
    ok = np.where(flag)[0]
    empty = np.where(mask)[0] 
    
    #Fill the gaps if there are some
    if len(empty) > 0 : 
        #Interpolate lon,lat @ empty positions
        outsla[empty] = interp1d(ok, outsla[ok], empty)
    
    #Get empty bin flag
    interpolated=~hist.astype('bool')
    
    return outsla, outlon, outlat, outind, ngaps, gapedges, gaplen, interpolated

def detrend(X,Z,deg=1):
    ndims=len(Z.shape)
    isVector=False
    if ndims == 1:
        Z=np.reshape(Z,(1,Z.size))
        isVector=True
    notMa=False
    if not isinstance(Z,np.ma.masked_array):
        notMa=True
        Z=np.ma.array(Z,mask=np.zeros(Z.shape))
    nt=Z.shape[0]
    valid=(~Z.mask).sum(axis=1)
    a=np.arange(deg+1)
    for t in np.arange(nt)[valid > 0]:
        fit=np.polyfit(X[~Z[t,:].mask],Z[t,:][~Z[t,:].mask], deg)
        for d in a : Z[t,:]-=np.power(X,a[::-1][d])*fit[d]
    if isVector : Z=Z.reshape(Z.size) 
    return Z.data if notMa else Z
