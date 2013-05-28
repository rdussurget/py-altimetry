import numpy as np
import matplotlib.pyplot as plt

import altimetry.tools as atools
from altimetry.tools.spectrum import get_spec, get_slope
from altimetry.data import alti_data

from scipy import stats


if __name__ == "__main__" :
    
    limit=np.array([36.5,2.,44.5,9.]) # NW Med
#    limit=np.array([41.0,5.,44.5,9.]) # Ligure
    limit=np.array([43.,3.,44.,9.]) # Track 9 coastal...

#    track_list=np.array([9])
    
    track_list_in = None
#    track_list_in = [146]
#    track_list_in = [9]
    
    sampling_rate=5. #5Hz sampling for PISTACH Products
#    sampling_rate=1. #Standard data
    
    sst_dev=2.5
    verbose=4
    vminsla=-0.15
    vmaxsla=0.15
    figsize=[23,15]
    
    trange_str = ['17/01/2012','07/05/2012']
    trange=[atools.cnes_convert(trange_str[0])[0],atools.cnes_convert(trange_str[1])[0]]
    
#    length_threshold = 600
#    N_min = 120 #(J2)
#    N_min = 90  #(C2)
#    N_min = 50  #LPC zone
    N_min = 15.

    N_min *= sampling_rate

    gaplen_max = 3 * sampling_rate
    
    
    
#    alti_pattern = "C:\\VMShared/data/alti/regional/europe/j2_cf/nrt_europe*.nc"
    alti_pattern = "C:\\VMShared\\data\\alti\\regional\\mersea-dt\\j2_cf\\dt_mersea*.nc"
#    alti_pattern = "C:\\VMShared\\data\\alti\\regional\\mersea\\j2_cf\\nrt_mersea*.nc"

#    alti_pattern = "C:\\VMShared\\data\\alti\\regional\\mersea\\c2_cf\\nrt_mersea*.nc"
#    alti_pattern = "C:\\VMShared\\data\\alti\\regional\\mersea-dt\\c2_cf\\dt_mersea*.nc"
#    alti_pattern = "C:\\VMShared/data/alti/regional/europe/c2_cf/nrt_europe*.nc"

#    alti_pattern='C:\\VMShared\\data\\alti\\regional\\PISTACH\\PISTACH_rtkRED3_NWMED\\PISTACH_L3_Product_NWMED_RED3_*.nc' #tr9_5hz.nc'
#    alti_pattern='C:\\VMShared\\data\\alti\\regional\\PISTACH\\PISTACH_rtkRED3_NWMED\\PISTACH_L3_Product_NWMED_RED3_tr9_5hz.nc'
#    alti_pattern='C:\\VMShared\\data\\alti\\regional\\PISTACH\\PISTACH_rtkRED3_NWMED\\PISTACH_L3_Product_NWMED_RED3_tr146_5hz.nc'

    alti_pattern='C:\\VMShared\\data\\alti\\regional\\PISTACH\\PISTACH_rtkMLE4_NWMED\\PISTACH_L3_Product_NWMED_MLE4_tr9_5hz.nc'
#    alti_pattern='C:\\VMShared\\data\\alti\\regional\\PISTACH\\PISTACH_rtkMLE4_NWMED\\PISTACH_L3_Product_NWMED_MLE4_tr146_5hz.nc'

#    alti=atools.alti_data(alti_pattern,limit=limit,verbose=verbose,time_range=trange)
    alti=alti_data(alti_pattern,limit=limit,verbose=verbose)
    
    alti.create_Variable('sla', alti.sla_filtered_21pts, {'time':alti._dimensions['time']})
    
    #update time range with real value
    trange_str = alti.time_range()[0]
    
    cycle_list = alti.cycle_list()
    
    slamat=np.array([],dtype=np.float64)
    dxout=np.array([],dtype=np.float64)
    
    empty=np.array([],dtype=np.bool)
    
    npts=np.array([],dtype=np.int)
    trnb=np.array([],dtype=np.int)
    cycnb=np.array([],dtype=np.int)
    
#    pmap=atools.plot_map(0,0,0,limit=limit)
#    cmap=plt.get_cmap('jet') 
    
    Npts_all = np.array([])
    
    #for cryosat, find long args (crossing med)
    
    for cycle in cycle_list :
        if track_list_in is None :
            print 'Warning : track_list is not yet set . Being defined from 1st cycle data'
            track_list=alti.track_list(alti.cycle == cycle)
        else : track_list=track_list_in
            
        if cycle == np.min(cycle_list) :
            orient = np.repeat(-1,len(track_list))
        
        for track in track_list:
            
            sla=alti.sla.compress((alti.track == track) & (alti.cycle == cycle))
            flag=abs(sla) <= 1.
            sla=sla.compress(flag)
            
            lon = alti.lon.compress((alti.track == track) & (alti.cycle == cycle)).compress(flag)
            lat = alti.lat.compress((alti.track == track) & (alti.cycle == cycle)).compress(flag)
            dst = atools.calcul_distance(lat,lon)
            
            Npts_all=np.append(Npts_all,len(sla))
            
            if orient[np.where(track == track_list)[0]] == -1:
                orient[np.where(track == track_list)[0]] = atools.track_orient(lon,lat)
            
#            print len(sla)
            
            #Check if current track is long enough
#            if dst.max() >= 500 :
            if len(sla) >= N_min :
                
                #Regrid SLA into the regular grid defined by hist
                dumdst, dumlon, dumlat, dumsla, gaplen, ngaps, dx, interpolated = atools.grid_track(dst,lat,lon,sla)
                dxout=np.append(dxout,dx)
                
                #Get some stats on raw tracks
                trnb=np.append(trnb,track)
                cycnb=np.append(cycnb,cycle)
                npts=np.append(npts,len(sla))
                
                
                #Extract a segment of length N (Centered on window)
                st=(len(dumdst) - N_min) /2
                en = st+N_min
                
                if len(dumsla[st:en]) != N_min :print "problem"
                
                #If gaps are detected, use interpolated values if only gaps are shorter than 3 points
                if ngaps > 0 : 
                    if gaplen.max() <= gaplen_max :
                        slamat=np.append(slamat,dumsla[st:en])
                        empty=np.append(empty,interpolated[st:en])
                else : 
                    slamat=np.append(slamat,dumsla[st:en])
                    empty=np.append(empty,interpolated[st:en])
                        
#                    psd,esd,fq,p = get_spec(mn_dx, sla,fft=True)
#                    plt.loglog(p,psd)
#                    plt.axis([1000.,10.,1e-3,1e4])
#                    plt.show()

                print 'cycle {0} - track {1}: {2} gap'.format(cycle,track,ngaps)

    #output matrix
    shape = (slamat.size/N_min,N_min)
    slamat=np.reshape(slamat,shape)
    empty=np.reshape(empty,shape)
    
    #Get SLA hovmoller
    rawsla = slamat.flatten()
    rawsla[np.arange(rawsla.size).compress(empty.flatten())]=np.NaN
    rawsla = np.reshape(rawsla,shape)
    
    plt.figure(figsize=[8,8])
#    plt.axis([10000.,1.,1e0,1e8])
    plt.axis([100.,1.,1e-1,1e5])
    plt.xlabel('Spatial wavelength (km)')
    plt.ylabel('Power spectral density (cm2.cpkm-1)')
    plt.title('{0} \n [{1}], {2} \n {3} passes of {4} elements '.format(alti.filelist[0],','.join('{0}'.format(x) for x in limit),' - '.join('{0}'.format(x) for x in trange_str),shape[0],N_min))
    
    mn_dx = np.median(dxout)
    psdmat = np.array([],dtype=np.float64)
    
    for i in np.arange(shape[0]) :
        spec = get_spec(mn_dx*1e3, slamat[i,:])
        psdmat = np.append(psdmat,spec['psd'])
        if i == 0 :
            fq=spec['fq']
            nfq=len(fq)
            
    #Unit conversions
    psdmat /=10. #convert from m2.cpm-1 to cm2.cpkm-1 (cpm = 1e3cpkm; m2 = 1e4cm2)
    psdmat=np.reshape(psdmat,(shape[0],nfq)) 
    
    #Mean spectrum
    mn_psd = np.mean(psdmat,axis=0)
#    mn_psd = np.median(psdmat,axis=0)
    
    fq*=1e3      #convert to m-1
    p = 1/fq
    
    #Compute spectral slope using linear regression (70-250km - Xu & Fu 2011)
    #########################################################################
    
    #1) Linear regression
    fqsub = fq[(p >= 70.) & (p <= 250.)]
    mnspecsub = mn_psd[(p >= 70.) & (p <= 250.)]
    f = np.repeat(fqsub,shape[0]).reshape((fqsub.size,shape[0])).transpose()
    specsub = psdmat[:,(p >= 70.) & (p <= 250.)]
    
    (slope,intercept) = get_slope(fqsub, mnspecsub)
    (slope2,intercept2) = get_slope(f, specsub)
    
    poly=np.poly1d([slope,intercept])
    poly2=np.poly1d([slope2,intercept2])
    
    #Get noise slope
    (nslope,nintercept) = get_slope(fq[(p >= 15.) & (p <= 30.)], mn_psd[(p >= 15.) & (p <= 30.)])
    npoly=np.poly1d([nslope,nintercept])
    
    #Compute signal/noise limit:
    # Es . k^(-sigma) = En.k^(-nu) --> Es : LS regression intercept, En : same for noise, sigma : signal slope, nu: noise slope, k:frequency
    #==> k = exp( (np.log(En) - np.log(Es)) / (nu-sigma) )
    klim = np.exp( (np.log(nintercept) - np.log(intercept)) / (nslope - slope) )  
    
    
#    #white noise?
#    wnoise=0.024*np.random.rand(2*len(fq))
#    pn,en,ff,pp = get_spec(mn_dx*1e3, wnoise)
    
    for i in np.arange(shape[0]) :
        plt.loglog(p,psdmat[i,:],'.k')
    
    
    plt.loglog(p,mn_psd,'-r',linewidth=4)
    
    if ~(np.isnan(slope) | np.isnan(slope2) | np.isnan(nslope)) : 
        xtxt = np.round(np.log10(fq.max()))-1
        ytxt = 0.1*(np.log10(plt.get(plt.gca(),'ylim')).max() - np.log10(plt.get(plt.gca(),'ylim')).min()) + poly(xtxt)
        nxtxt = np.round(np.log10(1/20.))
        nytxt = 0.1*(np.log10(plt.get(plt.gca(),'ylim')).max() - np.log10(plt.get(plt.gca(),'ylim')).min()) + npoly(nxtxt)
        
        
        plt.loglog(p, 10.0**poly(np.log10(fq)),':m',linewidth=2)
        plt.text(10.0**(-xtxt), 10.0**ytxt, '{0:4.1f}'.format(slope), color='m')
        plt.loglog(p[p < 50.], 10.0**npoly(np.log10(fq[p < 50.])),':b',linewidth=2)
        plt.text(10.0**(-nxtxt), 10.0**nytxt, '{0:4.1f}'.format(nslope), color='b')
    
    
    plt.gca().xaxis.grid(True,linestyle=':')
    plt.gca().yaxis.grid(True,linestyle=':')
    plt.show()
    
    import scipy.io as io
#    outfile='Q:\\Documents\\Post_Doc\\rapports\\rapport_ete_012\\J2_spec.{0}.sav'.format(track_list_in[0])
    outfile='Q:\\Documents\\Post_Doc\\rapports\\rapport_ete_012\\spectras\\mat{0}.HF.{1:.0f}-{2:.0f}'.format(alti.filelist[0],limit[0],limit[2])
    print outfile
    io.savemat(outfile,{'fq':fq,'spec':psdmat,'mn_spec':mn_psd})
    
    
    print 'done'

    