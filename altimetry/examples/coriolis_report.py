import altimetry.data.hydro as htools
import matplotlib.pyplot as plt
#from mayavi import mlab
#import chaco.shell as plt 
import numpy as np
import altimetry.tools as atools


#limit=np.array([42.,6.,44.5,11.])
limit=np.array([41.,4.,44.5,12.])
verbose=0

#cor_extract='C:\\VMShared\\data\\hydro\\IMEDIA\\extraction_Coriolis\\DataSelection_20120510_104737_39401\\'
cor_extract='C:\\VMShared\\data\\hydro\\IMEDIA\\extraction_Coriolis\\DataSelection_20120921_112034_50901\\'

#argo=htools.glider_data(cor_extract+'\\argo\\argo-profiles-6900504.nc',limit=limit,verbose=4)
#argotraj=htools.buoy_data(cor_extract+'\\argo\\argo-trajectory-6900504.nc',limit=limit,verbose=4)

argo=htools.glider_data(cor_extract+'argo\\argo-profiles*.nc',limit=limit,verbose=verbose)
##argotraj=htools.buoy_data(cor_extract+'\\argo\\argo-trajectory*.nc',limit=limit,verbose=verbose)
#dft=htools.buoy_data('C:\\VMShared\\data\\hydro\\IMEDIA\\extraction_Coriolis\\DataSelection_20120510_104737_39401\\drifters\\drifting-buoys-*.nc',limit=limit,verbose=verbose)
gld=htools.glider_data(cor_extract+'glider*\*.nc',limit=limit,verbose=verbose)
tsg=htools.buoy_data('C:\\VMShared\\data\\hydro\\IMEDIA\\extraction_Coriolis\\DataSelection_20120510_104737_39401\\tsg\\tsg-*.nc',limit=limit,verbose=verbose)
xbt=htools.glider_data(cor_extract+'xbt\\*.nc',limit=limit,verbose=verbose)

#dft.summary(all=True,fig=cor_extract+'\\drifters_',legend='Drifters',col='.k', endpoint='*y',ms=4,fontsize=12,textcolor='m',)
argo.summary(all=True,fig=cor_extract+'\\argo_',col='o-c',legend='Argo', endpoint='', ms=4,linewidth=1,fontsize=12,textcolor='m')
gld.summary(all=True,fig=cor_extract+'\\glider_', col='.r',endpoint='*-r', ms=4,fontsize=12,textcolor='m',legend='Glider')
tsg.summary(all=True,fig=cor_extract+'\\TSG_',legend='TSG',col='.-b', ms=4,fontsize=12,linewidth=1,endpoint='*g',textcolor='m')
xbt.summary(all=True,fig=cor_extract+'\\xbt_',legend='xbt',col='og', endpoint='*g', ms=4,linewidth=1,fontsize=12,textcolor='m')

##Get dates
#start=min([argo.date_surf.min(),dft.date.min(),gld.date_surf.min(),tsg.date.min(),xbt.date_surf.min()])
#end=max([argo.date_surf.max(),dft.date.max(),gld.date_surf.max(),tsg.date.max(),xbt.date_surf.max()])
start=min([argo.date_surf.min(),gld.date_surf.min(),tsg.date.min(),xbt.date_surf.min()])
end=max([argo.date_surf.max(),gld.date_surf.max(),tsg.date.max(),xbt.date_surf.max()])

scr_size=atools.get_screen_size('i')

pmap=atools.plot_map(0,0,0,limit=limit)
pmap.figure(figsize=scr_size)
p1,pmap=argo.map(np.ones(len(argo.id_surf),dtype='bool'), pmap=pmap, col='o-c', endpoint='', ms=8,linewidth=3,fontsize=12,textcolor='m',show=False)
#p2,pmap=dft.map(np.ones(dft.count,dtype='bool'), pmap=pmap, col='.k', endpoint='*y',ms=8,fontsize=12,textcolor='m',show=False)
p3,pmap=gld.map(np.ones(len(gld.id_surf),dtype='bool'), pmap=pmap, col='.r',endpoint='*-r', ms=8,fontsize=12,textcolor='m',show=False)
p4,pmap=tsg.map(np.ones(tsg.count,dtype='bool'), pmap=pmap, col='.-b', ms=8,fontsize=12,linewidth=3,endpoint='*g',textcolor='m',show=False)
p5,pmap=xbt.map(np.ones(len(xbt.id_surf),dtype='bool'), pmap=pmap, col='og', endpoint='*g', ms=8,linewidth=3,fontsize=12,textcolor='m',show=False)

#pmap.plot(dft.lon,dft.lat,'.k',ms=10)
#pmap.plot(argo.lon_surf,argo.lat_surf,'oc',ms=10)
#pmap.plot(gld.lon_surf, gld.lat_surf, '.r', ms=10)
#pmap.plot(tsg.lon, tsg.lat, '.b',ms=10)
#pmap.plot(xbt.lon_surf, xbt.lat_surf, 'og', ms=10)
pmap.legend([p1,p3,p4,p5],['drifters','argo floats','glider profiles','TSG transect','XBT profiles'])
pmap.title('Coriolis DB extraction from {0} to {1}'.format(atools.cnes_convert(start)[0][0],atools.cnes_convert(end)[0][0]))

#Scatter plot
#############
#pmap.scatter(dft.lon,dft.lat,dft.temp,vmin=12,vmax=15,s=5,edgecolor='none')
#pmap.scatter(gld.lonsurf, gld.latsurf, gld.tempsurf,vmin=12,vmax=15,s=5,edgecolor='none')
#pmap.scatter(tsg.lon, tsg.lat, tsg.temp,vmin=12,vmax=15,s=5,edgecolor='none')
#pmap.scatter(xbt.lonsurf, xbt.latsurf, xbt.tempsurf,vmin=12,vmax=15,s=5,edgecolor='none')
#pmap.colorbar()

pmap.setup_map()
#pmap.show()
pmap.savefig(cor_extract+'\\extraction_map.png')


##Old stuff
###########
##dft=htools.buoy_data('C:\\VMShared\\data\\hydro\\IMEDIA\\buoys\\positime_buoy_new_100512.dat')
##dft=htools.buoy_data('C:\\VMShared\\data\\hydro\\IMEDIA\\extraction_Coriolis\\DataSelection_20120510_104737_39401\\drifters\\drifting-buoys-61850.nc',verbose=4)
##pmap=atools.plot_map(0,0,0,limit=limit)
##pmap.scatter(dft.lon,dft.lat,dft.temp,vmin=12,vmax=15,s=5,edgecolor='none')
##pmap.colorbar()
##pmap.show()
#
##tsg=htools.TSG_data('C:\\VMShared\\data\\hydro\\IMEDIA\\TSG\\positime_boat.dat')
##
##pmap=atools.plot_map(0,0,0,limit=limit)
###pmap.scatter(tsg.lon,tsg.lat,tsg.temp,vmin=5,vmax=10,s=5,edgecolor='none')
##pmap.scatter(tsg.lon,tsg.lat,tsg.psal,vmin=38,vmax=38.5,s=5,edgecolor='none')
##pmap.colorbar()
##pmap.show()


#gld=htools.glider_data(cor_extract+'\glider\glider-profiles-58970.nc',limit=limit,verbose=4)
#flag=gld.get_file('glider-profiles-58970.nc')
#plt.figure(0)
#gld.contour_transect(gld.dist.compress(flag), gld.deph.compress(flag), gld.temp.compress(flag),vmin=13,vmax=14,marker=None)
#plt.figure(1)
#gld.plot_transect(gld.dist.compress(flag), gld.deph.compress(flag), gld.temp.compress(flag),vmin=13,vmax=14,marker=None)
##gld.plot_transect(gld.dist.compress(flag), gld.deph.compress(flag), gld.rho.compress(flag)-1000.,vmax=32,vmin=28)
##plt.colorbar()
#plt.show()
