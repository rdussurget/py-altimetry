import numpy as np
import matplotlib.pyplot as plt

import alti_tools as atools

from spectrum import get_kx, get_spec, grid_track

if __name__ == "__main__" :
    
    limit=np.array([36.5,4.,44.5,9.]) # NW Med
    limit=np.array([36.5,2.,44.5,9.]) # NW Med
#    track_list=np.array([9])
    
#    track_list_in = None
#    track_list_in = [9,146]
    track_list_in = [9]
    
    sst_dev=2.5
    verbose=0
    vminsla=-0.15
    vmaxsla=0.15
    figsize=[23,15]
    
    trange_str = ['17/01/2012','07/05/2012']
    trange=[atools.cnes_convert(trange_str[0])[0],atools.cnes_convert(trange_str[1])[0]]
    
#    length_threshold = 600
    N_min = 120 #(J2)
    N_min = 90  #(C2)
    
#    alti_pattern = "C:\\VMShared/data/alti/regional/europe/j2_cf/nrt_europe*.nc"
    alti_pattern = "C:\\VMShared\\data\\alti\\regional\\mersea-dt\\j2_cf\\dt_mersea*.nc"
#    alti_pattern = "C:\\VMShared\\data\\alti\\regional\\mersea\\j2_cf\\nrt_mersea*.nc"

#    alti_pattern = "C:\\VMShared\\data\\alti\\regional\\mersea\\c2_cf\\nrt_mersea*.nc"
#    alti_pattern = "C:\\VMShared\\data\\alti\\regional\\mersea-dt\\c2_cf\\dt_mersea*.nc"
#    alti_pattern = "C:\\VMShared/data/alti/regional/europe/c2_cf/nrt_europe*.nc"

#    alti=atools.alti_data(alti_pattern,limit=limit,verbose=verbose,time_range=trange)
    alti=atools.alti_data(alti_pattern,limit=limit,verbose=verbose)
    
    #update time range with real value
    trange_str = alti.time_range()[0]
    
    cycle_list = alti.cycle_list()
    
    slamat=np.array([],dtype=np.float64)
    dxout=np.array([],dtype=np.float64)
    
#    pmap=atools.plot_map(0,0,0,limit=limit)
#    cmap=plt.get_cmap('jet') 
    
    #for cryosat, find long args (crossing med)
    for cycle in cycle_list :
        if track_list_in is None :
            print 'Warning : track_list is not yet set . Being defined from 1st cycle data'
            track_list=alti.track_list(alti.cycle == cycle)
        else : track_list=track_list_in
        
        mndst=0
        for track in track_list:
            
            sla=alti.sla.compress((alti.track == track) & (alti.cycle == cycle))
            flag=abs(sla) <= 1.
            sla=sla.compress(flag)
            
            lon = alti.lon.compress((alti.track == track) & (alti.cycle == cycle)).compress(flag)
            lat = alti.lat.compress((alti.track == track) & (alti.cycle == cycle)).compress(flag)
            dst = atools.calcul_distance(lat,lon)
            
#            print len(sla)
            
            #Check if current track is long enough
#            if dst.max() >= 500 :
            if len(sla) >= N_min :
                
                #Regrid SLA into the regular grid defined by hist
                dumdst, dumlon, dumlat, dumsla, gaplen, ngaps, dx = grid_track(dst,lat,lon,sla)
                dxout=np.append(dxout,dx)
                
                #Extract a segment of length N
                st=(len(dumdst) - N_min) /2
                en = st+N_min
                
#                if ngaps == 0 :
                if ngaps > 0 : 
                    if gaplen.max() <= 3 :
                        slamat=np.append(slamat,dumsla[st:en])
                else : 
                    slamat=np.append(slamat,dumsla[st:en])
                        
#                    psd,esd,fq,p = get_spec(mn_dx, sla,fft=True)
#                    plt.loglog(p,psd)
#                    plt.axis([1000.,10.,1e-3,1e4])
#                    plt.show()

                print 'cycle {0} - track {1}: {2} gap'.format(cycle,track,ngaps)

    #output matrix
    shape = (slamat.size/N_min,N_min)
    slamat=np.reshape(slamat,shape)
    
    plt.figure(figsize=[8,8])
    plt.axis([1000.,10.,1e0,1e6])
    plt.title('J2 spectrum \n [{0}], {1} \n {2} passes of {3} elements '.format(','.join('{0}'.format(x) for x in limit),' - '.join('{0}'.format(x) for x in trange_str),shape[0],N_min))
    
    mn_dx = np.median(dxout)
    psdmat = np.array([],dtype=np.float64)
    
    for i in np.arange(shape[0]) :
        psd,esd,fq,p = get_spec(mn_dx*1e3, slamat[i,:])
        nfq=len(fq)
        psdmat = np.append(psdmat,psd)

    
    dst = mn_dx * np.arange(N_min)
    fq = np.fft.fftfreq(N_min, d=mn_dx)[0:nfq]
    fq = fq + (fq[1] - fq[0])
    p = 1/fq
    psdmat=np.reshape(psdmat,(shape[0],nfq))
    
    psdmat /=10. #convert from m2.cpm-1 to cm2.cpkm-1 (cpm = 1e3cpkm; m2 = 1e4cm2) 
    
    mn_psd = np.mean(psdmat,axis=0)
    
    for i in np.arange(shape[0]) :
        plt.loglog(p,psdmat[i,:],'.k')
    
    plt.loglog(p,mn_psd,'-r',linewidth=4)
    plt.gca().xaxis.grid(True,linestyle=':')
    plt.gca().yaxis.grid(True,linestyle=':')
    plt.show()
    
    print 'done'

    