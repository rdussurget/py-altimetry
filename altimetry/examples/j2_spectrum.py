'''
Created on 30 avr. 2013

@author: rdussurget
'''
import numpy as np
import matplotlib.pyplot as plt

from altimetry.data import alti_data
from altimetry.tools.spectrum import spectral_analysis, preprocess
import altimetry.tools as AT
from altimetry.tools.altimetry_tools import geost_1d
from altimetry.tools.others import deriv
from altimetry.tools.map_tools import plot_map

verbose=1
integrate=True
detrend=False
max_gap=3
per_min=15. #in %
t=222 #track number
#limit=[40.5,6,45,9]
limit=[40.5,6,43.5,10]

trange_str = ['12/07/2008','08/07/2012']


if __name__ == '__main__':
    
    trange,tdatetime=AT.cnes_convert(trange_str)
    
    #Load SLA data
    alti_pattern = '/archive/MSA/rdussurget/sauvegarde_13032013/VMShared/data/alti/regional/mfstep-dt/j2_cf/dt_upd_med_j2_sla_vxxc_*.nc'
    alti=alti_data(alti_pattern,limit=limit,verbose=verbose,time_range=trange.tolist(),track=t)
    alti.update_with_slice(alti.slice('track', t))
    alti.reorder()  
    
    cycle_list=alti.cycle_list()
    tracks=alti.track_list()
    trange_str = alti.time_range()[0]
    
    nt=len(cycle_list)
    
    sla=alti.sla*100.0
    lon=alti.lon
    lat=alti.lat
    dst=AT.calcul_distance(lat,lon)
    dx=np.median(AT.deriv(dst))
    N=len(lon)
    
    slaft=np.ma.empty(sla.shape)
    ugpl=slaft.copy()
    ugpllo=slaft.copy()
    ugfd=slaft.copy()
    ugfds=slaft.copy()
    
    #Filter small scales
    for i in np.arange(nt):
        slaft[i,:]=AT.loess(sla[i,:], dst, 40.)
        ugpl[i,:] = geost_1d(lon,lat,sla[i,:]/100.,pl04=True,p=20.,q=20.)
        ugpllo[i,:] = geost_1d(lon,lat,slaft[i,:]/100.,pl04=True,p=20.,q=20.)
        ugfd[i,:] = geost_1d(lon,lat,slaft[i,:]/100.,pl04=False)
        dumlon,dumlat,ugfds[i,:-1] = geost_1d(lon,lat,slaft[i,:]/100.,pl04=False,strict=True)

    ugpl.mask=np.isnan(ugpl)
    ugpllo.mask=np.isnan(ugpllo)
    
    plt.subplot(2,2,1)
    pmap=plot_map(0,0,0,limit=limit,resolution='i',)
    pmap.plot(lon,lat,np.mean(sla,0))
    pmap.across_track_arrow(lon, lat, (np.nansum(ugpl,0)/np.nansum(np.isfinite(ugpl),0)),scale=5)
    pmap.across_track_arrow(lon, lat, (np.nansum(ugpl,0)/np.nansum(np.isfinite(ugpl),0)),scale=5,ref=0.1,color='r')
    pmap.title('Powell & Leben')
    
    plt.subplot(2,2,2)
    pmap=plot_map(0,0,0,limit=limit,resolution='i',title='Loess 40km + Powell & Leben')
    pmap.plot(lon,lat,np.mean(sla,0))
    pmap.across_track_arrow(lon, lat, (np.nansum(ugpllo,0)/np.nansum(np.isfinite(ugpllo),0)),scale=5)
    pmap.across_track_arrow(lon, lat, (np.nansum(ugpllo,0)/np.nansum(np.isfinite(ugpllo),0)),scale=5,ref=0.1,color='r')
    
    plt.subplot(2,2,3)
    pmap=plot_map(0,0,0,limit=limit,resolution='i',title='Loess 40km + 3pts Fin. Diff.')
    pmap.plot(lon,lat,np.mean(sla,0))
    pmap.across_track_arrow(lon, lat, (np.nansum(ugfd,0)/np.nansum(np.isfinite(ugfd),0)),scale=5)
    pmap.across_track_arrow(lon, lat, (np.nansum(ugfd,0)/np.nansum(np.isfinite(ugfd),0)),scale=5,ref=0.1,color='r')

    plt.subplot(2,2,4)
    pmap=plot_map(0,0,0,limit=limit,resolution='i',title='Loess 40km + 2pts Fin. Diff.')
    pmap.plot(lon,lat,np.mean(sla,0))
    pmap.across_track_arrow(lon, lat, (np.nansum(ugfds,0)/np.nansum(np.isfinite(ugfds),0)),scale=5)
    pmap.across_track_arrow(lon, lat, (np.nansum(ugfds,0)/np.nansum(np.isfinite(ugfds),0)),scale=5,ref=0.1,color='r')


    pmap.show()

    plt.suptitle('MFSTEP J2 #222')
    plt.subplot(1,4,1)
    plt.pcolormesh(lat,alti.pass_time()/365.25 + 1950,ugpl*100,vmin=-10,vmax=20);plt.xlabel('Latitude (N)');plt.ylabel('Date');plt.title('Powell & Leben')
    plt.subplot(1,4,2)
    plt.pcolormesh(lat,alti.pass_time()/365.25 + 1950,ugpllo*100,vmin=-10,vmax=20);plt.xlabel('Latitude (N)');plt.title('Loess 40km + Powell & Leben')
    plt.subplot(1,4,3)
    plt.pcolormesh(lat,alti.pass_time()/365.25 + 1950,ugfd*100,vmin=-10,vmax=20);plt.xlabel('Latitude (N)');plt.title('Loess 40km + 3pts Fin. Diff.')
    plt.subplot(1,4,4)
    plt.pcolormesh(lat,alti.pass_time()/365.25 + 1950,ugfds*100,vmin=-10,vmax=20);plt.xlabel('Latitude (N)');plt.title('Loess 40km + 2pts Fin. Diff.')
    plt.show()
    
    plt.suptitle('MFSTEP J2 #222')
    plt.subplot(2,1,1)
    plt.plot(lat,np.mean(sla,0),'-r')
    plt.plot(lat,np.mean(slaft,0),'-b')
    plt.legend(['Raw SLA','Loess 40km'],loc=3)
    plt.ylabel('SLA (m)')
    plt.subplot(2,1,2)
    plt.plot(lat,(np.nansum(ugpl,0)/np.nansum(np.isfinite(ugpl),0))*100,'-r')
    plt.plot(lat,(np.nansum(ugpllo,0)/np.nansum(np.isfinite(ugpllo),0))*100,':r')
    plt.plot(lat,(np.nansum(ugfd,0)/np.nansum(np.isfinite(ugfd),0))*100,'-g')
    plt.plot(dumlat,(np.nansum(ugfds,0)/np.nansum(np.isfinite(ugfds),0))[:-1]*100,':g')
    plt.legend(['Powell & Leben','Loess 40km + Powell & Leben','Finite Diff (3 pts)','Finite Diff (2 pts)'],loc=4)
    plt.xlabel('Latitude (N)')
    plt.ylabel('Ugeo (m/s)')
    plt.show()

    #Preprocess SLA dataset
    #----------------------
#    slapp,ind=preprocess(lat, lon, sla, per_min=per_min, max_gap=max_gap)
    slapp,ind=preprocess(lat, lon, sla, per_min=per_min, max_gap=max_gap, remove_edges=True)
    npp=slapp.shape[1]
    slapp2=sla.copy()
    slapp2.mask[:]=True
    slapp2[ind,0:npp]=slapp
    
    #Run spectral analysis
    #---------------------
    spec = spectral_analysis(dx, slapp.transpose(), integration=integrate)
    AR5 = spectral_analysis(dx, slapp.transpose(), integration=integrate, ARspec=5)
    
    #Plot preprocessed datasets
    #---------------------------
    plt.subplot(1,2,1);plt.pcolormesh(np.arange(N),cycle_list,sla);plt.title('MFSTEP DT Upd vxxc')
    plt.subplot(1,2,2);plt.pcolormesh(np.arange(npp),cycle_list,slapp2);plt.title('PyValid preprocessing')
    bar=plt.colorbar();bar.set_label('SLA (cm)')
    plt.show()
    
    
    
    #Plot SLA spectra
    plt.axis([1000.,10.,1e1,1e4])
    plt.xlabel('Spatial wavelength (km)')
    plt.ylabel('Power spectral density (cm2.cpkm-1)')
    plt.title('Spatial spectrum of AVISO Residuals\nJason-2 #{0} (1Hz)'.format(t))

#    gain=np.mean(spec['psd']/AR5['psd'])

    plt.loglog(spec['p'],spec['psd'], '-oc', linewidth=4, markersize=4)
    plt.loglog(AR5['p'],AR5['psd']/2.0, '-ob', linewidth=2, markersize=2)
    
    plt.gca().xaxis.grid(True, linestyle=':')
    plt.gca().yaxis.grid(True, linestyle=':')
    plt.legend(['FFT spectrum','AR(5) spectrim'])
    plt.show()
    
        
    print 'done'
   
