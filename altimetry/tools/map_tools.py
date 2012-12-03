# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from mpl_toolkits.basemap import Basemap
import numpy as np

class plot_map(Basemap):
    def __init__(self,lon,lat,z,
                 resolution='i',
                 xoffset=.05,
                 yoffset=.05,
                 length=50.,
                 scale_lon=None,
                 scale_lat=None,
                 s=25,
                 c='k',
                 limit=None,
                 edgecolor='none'):
    
#        super(plot_map,self).__init__()#lon,lat,s=25,
#                 c=None,
#                 limit=None,
#                 resolution='i',
#                 xoffset=0.05,
#                 yoffset=0.05,
#                 length=50.,
#                 scale_lon=None,
#                 scale_lat=None)

        #Startup sequence
        #################
        
        #Force fields
        lon=np.array(lon)
        lat=np.array(lat)
        
        #Load default values (requested by Basemap.__init__()
        if limit is None : limit=np.array([lat.min() - 0.5, lon.min() - 0.5, lat.max() + 0.5, lon.max() + 0.5])
        else : limit=np.array(limit)
        
        clon = limit[[1, 3]].mean()
        clat = limit[[0, 2]].mean()
        
        #Init sequence
        ##############
        
        
        # Init Basemap class
        Basemap.__init__(self,projection='tmerc', resolution=resolution, llcrnrlat=limit[0], llcrnrlon=limit[1], urcrnrlat=limit[2], urcrnrlon=limit[3], lat_0=clat, lon_0=clon)
        
        #Set map characteristics fields
        self.limit=limit
        self.clon=clon
        self.clat=clat
        self.mcenter = self(self.clon, self.clat) #central position
        cparallel = np.array(self(self.limit[[1, 3]], [self.clat, self.clat])) #Get central paralell and meridian coordinates
        cmeridian = np.array(self(self.limit[[0, 2]], [self.clon, self.clon]))
        self.mllcorner = self(self.limit[1], self.limit[0])
        self.width = np.sqrt(np.diff(cparallel[:, 0]) ** 2 + np.diff(cparallel[:, 1]))
        self.heigth = np.sqrt(np.diff(cmeridian[:, 0]) ** 2 + np.diff(cmeridian[:, 1]))
        self.mxoffset = xoffset * self.width
        self.myoffset = yoffset * self.heigth
        
        if scale_lon is None :
            self.scale_lon = 0.9 * self.width + self.mllcorner[0]
            self.scale_lat = 0.1 * self.heigth + self.mllcorner[1]
        else :
            self.scale_lon = self(scale_lon, self.clat)
            self.scale_lat = self(self.clon, scale_lat)
        
        #Set plot characterisitics fields
        self.resolution=resolution
        self.length=length

        #Draw map 
        if np.size(z) > 1 : return self.setup_map(lon, lat, z, s=s, edgecolor=edgecolor)
        else :
            if z != 0 :  return self.setup_map(lon, lat, z, s=s, edgecolor=edgecolor)
        
#        self.show() #RQ : Show erase everything!!       
#        self.m.drawmapscale(self.scale_lon, self.scale_lat, self.mcenter[0], self.mcenter[1], length=self.length)

    
    #Redefined methods
    ##################
    def scatter(self,*args,**kwargs):

        mlon, mlat = self(np.array(args[0]),np.array(args[1]))
        
        if len(args) == 3 :
            c = np.array(args[2])
            if 'c' in kwargs :del kwargs['c']
        
        if 'edgecolor' in kwargs :
            edgecolor=kwargs['edgecolor']
            del kwargs['edgecolor']
        else : edgecolor = 'none'

        return Basemap.scatter(self,mlon,mlat,c=c,edgecolor=edgecolor,**kwargs)
#        self.colorbar(orientation='horizontal')
    
    def pcolormesh(self,*args,**kwargs):
        
        if len(args) == 1 : raise Exception('Lon/Lat arguments must be provided')
        
        lon=np.array(args[0]) 
        lat=np.array(args[1])
        
        if lon.ndim != 2 :
            lon,lat=np.meshgrid(lon,lat)
            lon=lon.transpose()
            lat=lat.transpose()
        
        mlon, mlat = self(lon,lat)
        image=np.ma.masked_array(args[2]) if not isinstance(args[2], np.ma.masked_array) else args[2]
#        image=np.ma.masked_array(args[2])
        return Basemap.pcolormesh(self,mlon,mlat,image,**kwargs)
    
    def contour(self,*args,**kwargs):
        
        lon=np.array(args[0]) 
        lat=np.array(args[1])
        
        if lon.ndim != 2 :
            lon,lat=np.meshgrid(lon,lat)
            lon=lon.transpose()
            lat=lat.transpose()
        
        mlon, mlat = self(lon,lat)
        
        if len(args) == 3 : return Basemap.contour(self,mlon,mlat,args[2],**kwargs)
        else : return Basemap.contour(self,mlon,mlat,args[2],args[3],**kwargs)
        
    def contourf(self,*args,**kwargs):
        
        lon=np.array(args[0]) 
        lat=np.array(args[1])
        
        if lon.ndim != 2 :
            lon,lat=np.meshgrid(lon,lat)
            lon=lon.transpose()
            lat=lat.transpose()
        
        mlon, mlat = self(lon,lat)
        
        if len(args) == 3 : return Basemap.contourf(self,mlon,mlat,args[2],**kwargs)
        else : return Basemap.contourf(self,mlon,mlat,args[2],args[3],**kwargs)
    
#    def clabel(self,*args,**kwargs):
#        Basemap.cl
    
    def drawmapscale(self, lon, lat, lon0, lat0, length, barstyle='simple', 
        units='km', fontsize=9, yoffset=None, labelstyle='simple', 
        fontcolor='k', fillcolor1='w', fillcolor2='k', ax=None, 
        format='%d', zorder=None):
        
#        x1,y1 = 0.3*m.xmax, 0.25*m.ymax
#        mlon,mlat = self(lon,lat,inverse=True)
#        mlon0,mlat0 = self(lon0,lat0,inverse=True)
        return Basemap.drawmapscale(self, lon, lat, lon0, lat0, length, barstyle=barstyle, units=units, fontsize=fontsize, yoffset=yoffset, labelstyle=labelstyle, fontcolor=fontcolor, fillcolor1=fillcolor1, fillcolor2=fillcolor2, ax=ax, format=format, zorder=zorder)
    
    def arrow(self,*args,**kwargs):
        lon,lat=self(args[0],args[1])
        dx=args[2]
        dy=args[3]
        ax=plt.gca()
        rg=np.arange(np.size(lon))
        if rg.max() > 0:
            for i in rg : ax.arrow(lon[i], lat[i], dx[i], dy[i], **kwargs)
        else : ax.arrow(lon, lat, dx, dy, **kwargs)
    
    def across_track_arrow(self,lon,lat,z, scale=None, ref=None, angle=0.0, **kwargs):
        
        if scale is None : scale = calcul_distance(self.clat,self.clon,self.clat,self.clon+1)*1e3
        else : scale *= calcul_distance(self.clat,self.clon,self.clat,self.clon+1)*1e3
        
        if ref is None : sense,orient = track_orient(lon, lat, orient=True)
        else :
            sense = True
            orient = (np.pi/2.)
            
        across_track = orient + angle - (np.pi/2.) if sense else orient + angle + (np.pi/2.)

        if isinstance(z, np.ma.masked_array) : z = z.data

        dx = scale * z * np.cos(across_track)
        dy = scale * z * np.sin(across_track)
        
        if ref is None : self.arrow(lon,lat,dx,dy,**kwargs)
        else :
            reflon = (self.lonmax - self.lonmin)*0.80 + self.lonmin
            reflat = (self.latmax - self.latmin)*0.15 + self.latmin
            reflatarr = (self.latmax - self.latmin)*0.2 + self.latmin
            self.arrow(reflon,reflatarr,scale*ref,0.0,**kwargs)
            if type(ref) == bool : self.text(reflon, reflat, '{0:.1f} m/s'.format(z[0]))
            else : self.text(reflon, reflat, '{0:.1f} m/s'.format(ref))
        
    
    def quiver(self,*args,**kwargs):
        lon=np.array(args[0]) 
        lat=np.array(args[1])
        
        if lon.ndim != 2 :
            lon,lat=np.meshgrid(lon,lat)
            lon=lon.transpose()
            lat=lat.transpose()
        
        mlon, mlat = self(lon,lat)
        
        u=np.ma.masked_array(args[2]) if not isinstance(args[2], np.ma.masked_array) else args[2]
        v=np.ma.masked_array(args[3]) if not isinstance(args[3], np.ma.masked_array) else args[3]
        
        return Basemap.quiver(self,mlon,mlat,u,v,**kwargs)  
    
    #Following is not working
    def plot(self,*args,**kwargs):
        
        mlon, mlat = self(args[0],args[1])
        
        #Get marker color
       
            
        #override if passed directly through 3rd argument
        if len(args) == 3 :
            c = args[2]
        else :
            if 'c' in kwargs :
                c = kwargs['c']
                del kwargs['c']
            else : c = '.k'
            
        #Override plot ms keyword - becomes the same syntax as scatter
        if 's' in kwargs :
            kwargs['ms'] = kwargs['s']
            del kwargs['s']
        
        if isinstance(c,str) : return Basemap.plot(self,mlon,mlat,c,**kwargs)
#        self.colorbar(orientation='horizontal')
    
#    def plot(self,lon,lat,
#             s=25,
#             c='k',
#             edgecolors=None):
##             *args):
#
##        if len(args) == 1 : 
##        z=args[0]
#                
#        mlon, mlat = self(lon, lat)
##        if c is None : 
##        if len(args) == 1 :self.m.scatter(mlon, mlat, s=s, c=z, edgecolors=edgecolors)
##        else : self.m.scatter(mlon, mlat, s=s, c=c, edgecolors=edgecolors)
#        self.scatter(mlon, mlat, c=c, s=s, edgecolors=edgecolors)
        
#        if not hasattr(self, 'cbar') and len(args) == 1 : self.cbar = plt.colorbar(orientation='horizontal')
    
    def get_cursor(self,*args):
        c= zip(*(pylab.ginput(n=0,mouse_add=1, mouse_pop=2, mouse_stop=3)))
        lonpt,latpt = self(c[0],c[1],inverse=True)
        for enum in enumerate(zip(*(lonpt,latpt))):
            dst=calcul_distance(enum[1][1],enum[1][0],args[0],args[1])
            i=np.argmin(dst)
            xdum = args[0][i]
            ydum = args[1][i]
            if len(args) == 3 : zdum = args[2][i]
            if enum[0] == 0 :
                xout=xdum
                yout=ydum
                if len(args) == 3 : zout = zdum
            else :
                xout=np.append(xout,xdum)
                yout=np.append(yout,ydum)
                if len(args) == 3 : zout = np.append(zout,zdum)
        if len(args) == 3 : return xout, yout, zout
        else : return xout, yout
                
    def get_cursor_id(self,x,y):
        c= zip(*(pylab.ginput(n=0,mouse_add=1, mouse_pop=2, mouse_stop=3)))
        lonpt,latpt = self(c[0],c[1],inverse=True)
        for enum in enumerate(zip(*(lonpt,latpt))):
            dst=calcul_distance(enum[1][1],enum[1][0],x,y)
            if enum[0] == 0 : i=[np.argmin(dst).tolist()]
            else : i.append(np.argmin(dst))
        return i
    
    def figure(self,*args,**kwargs):
        return plt.figure(*args,**kwargs)
        
    def legend(self,*args,**kwargs):
        plt.legend(*args,**kwargs)
        
    def text(self,lon, lat, string,**kwargs):
        mlon,mlat=self(lon,lat)
        return plt.text(mlon,mlat,string,**kwargs)
    
    def title(self,s,**kwargs):
        plt.title(s,**kwargs)
    
    def colorbar(self,*args,**kwargs):
        return plt.colorbar(*args,**kwargs)
        
    def setup_map(self,londel=1.,latdel=1.,*args,**kwargs):
        """
        SETUP_MAP : Draw the map and coastlines.
        
        @note: To be compliant with self.legend() -> setup_map() must be called after legend
        """
        
        s = kwargs['s'] if kwargs.has_key('s') else 25
        edgecolor = kwargs['edgecolor'] if kwargs.has_key('edgecolor') else 'none'
        
        #Plot Sequence
        ##############
        self.drawcoastlines()
        cs=self.drawbathy(**kwargs)
        
        self.drawparallels(np.arange(np.floor(self.limit[0]), np.ceil(self.limit[2]), latdel), labels=[1, 0, 0, 0])
        self.drawmeridians(np.arange(np.floor(self.limit[1]), np.ceil(self.limit[3]), londel), labels=[0, 0, 0, 1])
                
        if len(args) == 3 :
            lon=args[0]
            lat=args[1]
            z=args[2]
            if np.size(z) > 1 : self.scatter(lon,lat,z,s=s,edgecolor=edgecolor)
            else : 
                if type(z) != type(str()) : z = '.k'
                self.plot(lon,lat,z,ms=s)
        return cs
    
    def drawbathy(self,fname=None,V=[-250,-2000],Nlevels=None,step=None,colors='k',linestyles='dotted', linewidths=2,**kwargs):

        #If 1 point within, go for menor        
        if np.sum(in_limits(self.limit[[1,3]],self.limit[[0,2]],[35,0,45,25])[1]) == 1 :
            bathy='MENOR'
        else : bathy = 'ETOPO'
        
        if fname is None :
            bat=load_bathy(bathy=bathy,limit=self.limit)
        else :
            bat=load_bathy(fname=fname,limit=self.limit)
        
        glon,glat=np.meshgrid(np.squeeze(bat['lon']), np.squeeze(bat['lat']))
        
        mn=0.01
        mx=np.nanmin(bat['Z'].data)
        if Nlevels is not None :
            if step is None : step = (mx - mn) / Nlevels
            V = np.arange(mn,mx,step,dtype=np.float)
        
        return self.contour(glon,glat,bat['Z'],V,colors=colors,linestyles=linestyles,linewidths=linewidths,**kwargs)
        
        
        
    def show(self):
        plt.show()
        
    def savefig(self,*args,**kwargs):
        plt.savefig(*args,**kwargs)
    
#    def set_local(self,VarName,VarValue):
#        if VarName not in locals() : exec('self.'+VarName+' = VarValue')

def load_bathy(bathy='MENOR',**kwargs):
    
    if bathy == 'ETOPO' :
        fname='C:\\VMShared\\data\\spare_products\\bathy\\ETOPO2v2g_f4.nc'
        ncf=nc(fname)
        if kwargs.has_key('lon') : kwargs['x']=recale(kwargs.pop('lon'),degrees=True)
        if kwargs.has_key('lat') : kwargs['y']=kwargs.pop('lat')
        if kwargs.has_key('limit') :
            limit=kwargs.pop('limit')
            kwargs['x']=tuple(recale(limit[[1,3]],degrees=True))
            kwargs['y']=tuple(recale(limit[[0,2]],degrees=True))
        str=ncf.read(**kwargs)
        str['lon']=np.ma.array(recale(str.pop('x'),degrees=True))
        str['lat']=np.ma.array(str.pop('y'))
        str['Z']=np.ma.array(str.pop('z'))
        
        return str
        
    elif bathy == 'MENOR' :
        fname='C:\\VMShared\\data\\spare_products\\bathy\\bathy_menor.mat'
        m_dict=OrderedDict( (('lon_menor',[]), ('lat_menor',[]), ('H0',[])) )
        str=io.loadmat(fname,m_dict)
        
        str['lon']=np.ma.array(np.squeeze(str['lon_menor']))
        str['lat']=np.ma.array(np.squeeze(str['lat_menor']))
        str['Z']=np.ma.array(str['H0'])
        
        del str['lon_menor']
        del str['lat_menor']
        del str['H0']
    
        return str
    else : raise Exception('Unknown bathymetry') 