# -*- coding: utf-8 -*-
'''
;@author: rdussurg
;
; CNES_convert: converts a date from a specified time reference frame to another
;
;@param in {in}{required}{type=STRING|NUMERIC} input date (either a scalar or a vector)
;@keyword julian {in}{optional}{type=BOOLEAN} True if output date is in julian format
;@keyword calendar {in}{optional}{type=BOOLEAN} True if output date is in calendar format
;@keyword nasa {in}{optional}{type=BOOLEAN} True if reference frame should be NASA time (days since 01/01/1958) instead of CNES time (days since 01/01/1950)
;@keyword string {in}{optional}{type=BOOLEAN} True if output date should be formatted as a string instead of a 3 column matrix (days, months, years)
;@keyword out {out}{optional}{type=STRING|NUMERIC} Used to return output variable
;@keyword quiet {in}{optional}{type=BOOLEAN} Turn off comments
;@keyword calformat {out}{optional}{type=STRING} Calendar format used for conversion. dd : day, mm : month, yyyy:year; with any non-alphanumeric separator
;
;@examples 1) Converts an array of calendar dates into julian dates<br />
;   IDL> CNES_convert, ['01/01/2007','02/01/2007'], /julian, out=out<br />
;   IDL> print, out<br />
;   20819<br />
;   20820<br />
;   example 2) Converts an array of julian dates into string formatted calendar dates<br />
;   IDL> CNES_convert, [20819,20820], /calendar, out=out, /string<br />
;   IDL> print, out<br />
;   1/1/2007 2/1/2007<br />
;   example 3) Converts a NASA julian date into a calendar data<br />
;   IDL> CNES_convert, 17897, /calendar, /nasa<br />
;   Calendar date : 1/1/2007<br />
;   example 4) Converts a calendar date into a julian date<br />
;   IDL> CNES_convert, '01/01/2007', /julian<br />
;   Julian day : 20819 (CNES time)
;
;@todo Fix bug with NaN values : CNES_convert, !VALUES.D_NAN returns "19/12/****"
;
;@author : Renaud DUSSURGET (RD), LEGOS/CTOH, Toulouse
;          Renaud DUSSURGET (RD@IFR), IFREMER LER/PAC, La Seyne/Mer
;@history
;   - Jan. 2009 : First release<br />
;   - 9 Feb. 2009 : changed default options and fixed logical errors<br />
;   - March 2009 : Working on arrays and added OUT keyword to get output data<br /> 
;   - April 2009 : quiet mode<br />
;   - May 2009 : Automatic recognition of output format (calendar or julian)<br />
;   - May 2010 : Added CALFORMAT option for date conversion in any format
;   - Apr. 2012 : Ported to Python by RD@IFR (Created on 24 avr. 2012)
;
'''

#Loading libraries
#import sys
import datetime
import numpy as np
import scipy as sc
from scipy import interpolate, io as interpolate, io
from matplotlib import pyplot, pylab as plt, pylab
from mpl_toolkits.basemap import Basemap
from bisect import bisect_left
from Tkinter import Tk
from netCDF4 import Dataset as ncfile
import bisect
import seawater.gibbs
from collections import OrderedDict

import glob
import os

import hydro_tools as htools
import esutils_stat as es
from nctools import nc, load_ncVar



def in_limits(lon, lat, limit):
    
    limit = recale_limits(limit)
    
    lats = limit[0]
    lons = limit[1]
    late = limit[2]
    lone = limit[3]

    ind=[]
    
    lon = recale(lon,degrees=True)

    # find out lon between S and E (start and end)
    if lons < lone:
        lon_flag = (lon >= lons) & (lon <= lone)
        
    # find out lon outside S and E (between E and S)
    else:
        lon_flag = (lon >= lons) | (lon <= lone)

    # nothing remains
    if not np.sometrue(lon_flag):
        print "* WARNING : No lon [", lons, ",", lone, "]"
        return ind, lon_flag

    # find lat constraints
    lat_flag = (lat >= lats) & (lat <= late)
 
    # nothing remains
    if not np.sometrue(lat_flag):
        print "* WARNING : No Lat [", lats, ",", lats, "]  with lon  [", lons, ",", lone, "]"
        return ind, lat_flag
    
    #Construct flag and index arrays
    if (len(lon_flag) == len(lat_flag)) :
        flag = lon_flag & lat_flag
        ind = np.arange(len(lon)).compress(flag)
    else :
        flag = (lon_flag,lat_flag)
        ind = (np.arange(len(lon)).compress(flag[0]),np.arange(len(lat)).compress(flag[1]))
    

    # compress indexes
#    ind = ind.compress(flag)
    
    return ind, flag

def recale_limits(limitin, zero_2pi=True, minpi_pi=None): 
    
    #defaults
    if (minpi_pi is not None) : zero_2pi != minpi_pi
    
    limitout=limitin
    if isinstance(limitin,np.ndarray) : limitout[[1,3]]=recale(limitin[[1,3]], zero_2pi=zero_2pi, minpi_pi=minpi_pi,degrees=True) #Always in degrees!
    else :
        limitout[1]=recale(limitin[1], zero_2pi=zero_2pi, minpi_pi=minpi_pi,degrees=True)
        limitout[3]=recale(limitin[3]-1e2*np.MachAr().resolution, zero_2pi=zero_2pi, minpi_pi=minpi_pi,degrees=True)

    return limitout

def recale(angle, degrees=False, zero_2pi=True, minpi_pi=None) :

    #defaults
    if (minpi_pi is not None) : zero_2pi = not minpi_pi
    radians = not degrees
  
    modulus= 360.0 if degrees else 2.0 * np.pi
  
    limits=np.array([0.0, 2.0*np.pi])
    if minpi_pi : limits-=np.pi
    if degrees : limits = np.rad2deg(limits)
    
    angle = np.mod(angle + limits[0],limits[0]+limits[1]) - limits[0]
    
#    over=angle >= limits[1]
#    under=angle < limits[0]
#    angle+=under*modulus
#    angle-=over*modulus
  
    return angle


def rgb_int2tuple(rgbint):
    return (rgbint // 256 // 256 % 256, rgbint // 256 % 256, rgbint % 256)

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
        
def cnes_convert(argin,
                 julian=True,
                 calendar=False,
                 matlab=False,
                 #nasa=False,
                 #string=False,
                 #calformat=False,
                 epoch=None,
                 fromReference=None,
                 verbose=False):
    
    #Scalar to vector conversion
    if np.size(argin) == 1 :
        if type(argin) is not list : argin = [argin]
    
    if type(argin[0]) == str : julian = True
    else : calendar = True
    
#    print type(argin[0])
    
    if calendar is True : julian = False
    if julian is True : calendar = False

    if epoch is None : epoch = datetime.datetime(1950, 01, 01)
    if matlab is True : epoch = datetime.datetime(1, 1, 1)
    
    if julian is True :
        if verbose is True : print("julian is true")
#        srlist=[np.array(x.split("/",2) for x in argin,dtype=int)]
        narg = np.size(argin)
        
        strlist  =[x.split("-") for x in argin]
        datelist = []
        for x in strlist :
            if len(x) == 1 : datelist.append(datetime.datetime.strptime(x[0].strip(),"%d/%m/%Y") )
            else : datelist.append(datetime.datetime.strptime('{0}-{1}'.format(x[0].strip(),x[1].strip()),"%d/%m/%Y-%H:%M"))  
#        datelist = [datetime.datetime.strptime(x,"")""]
#        strlist = np.array([x.split("/", 2) for x in argin], dtype=int)
#        datelist = [datetime.datetime(strlist[x, 2], strlist[x, 1], strlist[x, 0]) for x in np.arange(narg)]
        return np.array([(datelist[x] - epoch).days for x in np.arange(narg)], dtype=float)
        
    if calendar is True :
        if verbose is True : print("caldendar is true")
        datelist = [epoch.toordinal() + x for x in  argin]
        if verbose is True : print(datelist)
#        days = [np.int(x) for x in datelist]
#        decim = np.array(datelist) - np.array(days)
#        return [datetime.date.fromordinal(y) for y in datelist]
        a=datetime.datetime.fromordinal(int(datelist[0]))
        dateObjList=[datetime.datetime.fromordinal(int(y)) + datetime.timedelta(seconds=(y - np.floor(y))*86400.0) for y in datelist]
        dateStr=[Obj.strftime('%Y/%m/%d') for Obj in dateObjList]
        return dateStr,dateObjList
#        dd,mm,yyyy=np.array(argin.split("/",2),dtype=int)
#        varout = np.array((datetime.datetime(yyyy,mm,dd) - epoch).days,dtype=float)
        
        

def calcul_distance(*args):
    """
    ;+
    ; CALCUL_DISTANCE : fonction permettant de calculer la distance entre deux points sur la sphere
    ;
    ; @Author : Mathilde FAILLOT (LEGOS/CTOH), 2005
    ; @History :
    ;    - 2005 : First release, source:<br />
    ;             Weisstein, Eric W. "Great Circle." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/GreatCircle.html<br /> 
    ;    - 2008 : Renaud DUSSURGET, optimization + vectors + Mercurial<br />
    ;    - Feb. 2009 : RD, debug on repeated positions removal<br />
    ;    - March 2010 : RD, force coordinates to FLOATS if Floating-Point operand errors<br />
    ;    - Nov. 2010 : 2 params option for calculating distance vector
    ;    - Apr. 2012 : Converted to Python
    ;-
    """
#def calcul_distance(norepeat=False, *args):
    
    lon_a_in = args[1]
    lat_a_in = args[0]
    
    if np.size(lon_a_in) != np.size(lat_a_in) : raise 'Error : arguments are not the same size'
    
    if len(args) == 2 :
        lon_b_in = lon_a_in
        lat_b_in = lat_a_in
        lon_a_in = get_zero_element(lon_a_in)
        lat_a_in = get_zero_element(lat_a_in)
    elif len(args) == 4 :
        lat_b_in = args[2]
        lon_b_in = args[3]
        if np.size(lon_b_in) != np.size(lat_b_in) : raise 'Error : arguments are not the same size'
    else :
        print "ERROR"
        return -1

    #   Define constants
    
    #Degree to radians conversion
    #deg_rad = np.pi/360.0
    
    #Earth radius (km)
    rt = 6378.137
    
    #;Remove repeated positions
    #IF (KEYWORD_SET(NOREPEAT)) THEN id = WHERE(~((lon_b_in EQ lon_a_in[0]) * (lat_b_in EQ lat_a_in[0])),okCnt) ELSE okcnt = 1 
    #IF (okCnt EQ 0) THEN RETURN lon_a_in, DOUBLE(0.)

    #Degree to radians conversion
    lon_a = np.deg2rad(lon_a_in)
    lat_a = np.deg2rad(lat_a_in)
    lon_b = np.deg2rad(lon_b_in)
    lat_b = np.deg2rad(lat_b_in)


    #;   calcul de la distance entre les deux points consideres
#    interm = (np.cos(lon_b) * np.cos(lat_b)) * (np.cos(lon_a) * np.cos(lat_a)) + \
#      (np.sin(lon_a) * np.cos(lat_a)) * (np.sin(lon_b) * np.cos(lat_b)) + \
#      np.sin(lat_a) * np.sin(lat_b)
    
#    interm = (np.cos(lon_b) * np.cos(lat_b)) * (np.cos(lon_a) * np.cos(lat_a)) + \
#      (np.sin(lon_a) * np.cos(lat_a)) * (np.sin(lon_b) * np.cos(lat_b)) + \
#      np.sin(lat_a) * np.sin(lat_b)

    interm = np.cos(lat_a) * np.cos(lat_b) * np.cos(lon_a - lon_b) + np.sin(lat_a) * np.sin(lat_b) #Wolfram Math World definition
    
    dist = rt * np.arccos(interm)
    
    #Remove computation errors
    fgerr=(lon_b == lon_a) & (lat_b == lat_a)
    if fgerr.sum() > 0 : dist[fgerr] = 0.
    
    return dist

def get_zero_element(array):
    try : array = array.flat.next()
    except StopIteration : pass
    return array

def cumulative_distance(lat, lon,dist_int=None):
    '''
    ;+
    ; CUMULATIVE_DISTANCE : permet de calculer la distance le long d'une ligne.
    ;
    ; @Author : Renaud DUSSURGET, LEGOS/CTOH
    ; @History :
    ;    - Feb. 2009 : First release (adapted from calcul_distance)
    ;-
    '''

    #   Define constants
    
    #Degree to radians conversion
    #deg_rad = np.pi/360.0
    
    #Earth radius (km)
    rt = 6378.137
    
    #;Remove repeated positions
    #IF (KEYWORD_SET(NOREPEAT)) THEN id = WHERE(~((lon_b_in EQ lon_a_in[0]) * (lat_b_in EQ lat_a_in[0])),okCnt) ELSE okcnt = 1 
    #IF (okCnt EQ 0) THEN RETURN lon_a_in, DOUBLE(0.)

    nelts=lon.size

    #Degree to radians conversion
    lon_a = np.deg2rad(lon[0:nelts - 1])
    lon_b = np.deg2rad(lon[1:nelts])
    lat_a = np.deg2rad(lat[0:nelts - 1])
    lat_b = np.deg2rad(lat[1:nelts])
    
    interm = np.cos(lat_a) * np.cos(lat_b) * np.cos(lon_a - lon_b) + np.sin(lat_a) * np.sin(lat_b) #Wolfram Math World definition
    
    dist_int=np.append(0,rt*np.arccos(interm))
    
    return dist_int.cumsum()


def distance_matrix(lat_a, lon_a, lat_b, lon_b):
    
    na = np.size(lat_a)
    nb = np.size(lat_b)
    
    #Earth radius (km)
    rt = 6378.137
    
#    lon_a_tile = np.tile(np.deg2rad(lon_a), (nb, 1))
#    lat_a_tile = np.tile(np.deg2rad(lat_a), (nb, 1))
#    lon_b_tile = np.tile(np.deg2rad(lon_b).transpose(), (na, 1)).transpose()
#    lat_b_tile = np.tile(np.deg2rad(lat_b).transpose(), (na, 1)).transpose()

#    plt.figure()
#    plt.imshow(lon_a_tile,interpolation='NEAREST',aspect='auto')
#    plt.show()
#    
#    plt.figure()
#    plt.imshow(lon_b_tile,interpolation='NEAREST',aspect='auto')
#    plt.show()
    
#    interm = np.cos(lat_a_tile) * np.cos(lat_b_tile) * np.cos(lon_a_tile - lon_b_tile) + np.sin(lat_a_tile) * np.sin(lat_b_tile)
    
    #lower memory usage for huge matrices
    dum = np.tile(np.deg2rad(lat_a), (nb, 1))
    dist_matrix = np.cos(dum)
    interm = np.sin(dum)
    
    del dum
    dum =  np.tile(np.deg2rad(lat_b).transpose(), (na, 1)).transpose()
    dist_matrix *= np.cos(dum)
    interm *= np.sin(dum)
    
    del dum
    dum = np.tile(np.deg2rad(lon_a), (nb, 1))
    dum -= np.tile(np.deg2rad(lon_b).transpose(), (na, 1)).transpose()
    dist_matrix *= np.cos(dum)
    
    del dum
    
    dist_matrix += interm
    
    del interm
        
    dist_matrix = rt * np.arccos(dist_matrix)
    
    return dist_matrix


def nearest(t, x):
    adiff=np.abs(t-np.float(x))
    i=np.argmin(adiff)
#    i = bisect_left(t, x)
#    if i == len(t) : i-=1
#    if t[i] - x > 0.5:
#        i-=1
    return i



def bytscl(array, max=None , min=None , nan=0, top=255 ):
    """
    see http://star.pst.qub.ac.uk/idl/BYTSCL.html
    note that IDL uses slightly different formulae for bytscaling floats and ints. 
    here we apply only the FLOAT formula...
    """
    if max is None: max = np.nanmax(array)
    if min is None: min = np.nanmin(array)
    return np.maximum(np.minimum(((top+1.0)*(array-min)/(max-min)).astype(np.int16), top),0)

def scale_vector(x,bottom,top):
    return (x - bottom) / (top-bottom)

#Screen size widgets
def get_screen_size(units='p'):
    units='1'+units
    width = Tk().winfo_fpixels(str(Tk().winfo_screenwidth())+'p')/Tk().winfo_fpixels(units)
    height = Tk().winfo_fpixels(str(Tk().winfo_screenheight())+'p')/Tk().winfo_fpixels(units)
    return np.array([width,height])
def get_screen_dpi(units='1i'):
    return Tk().winfo_fpixels(units)

def where_list(list1,list2):
    index=[]
    for i in list1 :
        try:
            index.append(list2.index(i))
        except ValueError:
            index.append(-1)

    return index
#    return [list2.index(i) for i in list1] 


#Load alti data
###############
class alti_data(htools.hydro_data) :
    def __init__(self,file_pattern,time_range=None,**kwargs):
        
        if time_range is not None :
            ls=glob.glob(file_pattern)
            filelist=[os.path.basename(p) for p in ls]
            st=[f.split('_')[-3] for f in filelist]
            en=[f.split('_')[-2] for f in filelist]
            if (~pylab.is_numlike(st) & ~pylab.is_numlike(en)) :
                self.warning('Date range cannot be found - extract all files from file_pattern')
            else :
                jst = [cnes_convert('{0}/{1}/{2}'.format(s[-2:],s[-4:-2],s[0:4]))[0] for s in st]
                jen = [cnes_convert('{0}/{1}/{2}'.format(s[-2:],s[-4:-2],s[0:4]))[0] for s in en]
                dt=np.fix(np.median(np.array(jst[1:]) - np.array(jst[:-1])))
    
                hist,R= es.histogram(np.array(jst),binsize=dt,use_weave=False,rev=True,min=time_range[0] - dt/2.,max=time_range[1] + dt/2.)
                dumind = histogram_indices(hist, R)
                ind=np.array([])
                for i in dumind : ind=np.append(ind,i)
                ind=ind.tolist()
                file_pattern=np.array(ls)[ind]
            
        htools.hydro_data.__init__(self,file_pattern,**kwargs)
        self.set_sats()
    
    def set_sats(self):
        self.sat=[]
        for enum in enumerate(self.filelist):
            self.sat=np.append(self.sat,[enum[1].split('_')[2]]*self.filelist_count[enum[0]])
    
    def read(self,filename,datatype=None,**kwargs):
        
        fname,extension = os.path.splitext(filename)
        
        if os.path.basename(filename).count('.') > os.path.basename(filename).count('_'): delim='.'
        else : delim = '_'
        
        #Get data type
        if datatype is None :
            if os.path.basename(filename).split(delim)[0] == 'ctoh' : datatype='CTOH'
            if os.path.basename(filename).split(delim)[0] == 'PISTACH' : datatype='PISTACH'
            if os.path.basename(filename).split(delim)[0] == 'nrt' : datatype='NRT'
            if os.path.basename(filename).split(delim)[0] == 'dt' : datatype='DT'
        else :
            datatype=datatype.upper() #Setup default as AVISO dt
        
        if (datatype == 'DT') | (datatype == 'NRT') | (datatype == 'PISTACH') :
            outStr=self.read_SLACF(filename,datatype=datatype,**kwargs)
            self.update_fid_list(os.path.basename(filename),outStr['_dimensions']['time'])
        elif (datatype == 'CTOH') :
            outStr=self.read_CTOH(filename,**kwargs)
            self.update_fid_list(os.path.basename(filename),outStr['_dimensions']['time'])
        if (datatype == 'RES') :
            outStr=self.read_SLARES(filename,**kwargs)
            self.update_fid_list(os.path.basename(filename),outStr['_dimensions']['Data'])
       
        
        return outStr
    
    def read_SLACF(self,filename,params=None,force=False,timerange=None,datatype=None,**kwargs):
        
        """
        READ_SLACF : Read AVISO Along-Track SLA CF-compliant products
        
        @return: outStr {type:dict} Output data structure containing all recorded parameters as specificied by NetCDF file PARAMETER list.
        @author: Renaud Dussurget
        """
        
        #Open file
        self._filename = filename
        self._ncfile = ncfile(self._filename, "r")
        
        #Get delimiter
        if os.path.basename(filename).count('.') > os.path.basename(filename).count('_'): delim='.'
        else : delim = '_'
        
        #Gat sat name
        splitted=os.path.basename(filename).split(delim)
        if (datatype == 'DT') | (datatype == 'NRT') : sat_name = splitted[2] if splitted[0] == 'nrt' else splitted[3]
        if datatype == 'PISTACH' : sat_name = 'J2'
        
     
        #Get list of recorded parameters:
        par_list=[i.encode() for i in self._ncfile.variables.keys()]
        for i in ['time','longitude','latitude'] : par_list.pop(par_list.index(i))
        nparam=len(par_list)
        
        self.message(2,'Recorded parameters : '+str(nparam)+' -> '+str(par_list))
     
        lon = load_ncVar('longitude',nc=self._ncfile,**kwargs)
        lon['data'] = recale(lon['data'], degrees=True, zero_2pi=True) #shift longitudes
        lat = load_ncVar('latitude',nc=self._ncfile,**kwargs)

        
        #Extract within limits
        ind, flag = in_limits(lon['data'],lat['data'],limit=self.limit)
        dim_lon = lon['_dimensions']
        lat = lat['data'].compress(flag)
        lon = lon['data'].compress(flag)
        dist=cumulative_distance(lat, lon)
        
        sz=np.shape(lon)
        ndims=np.size(sz)
            
        date = load_ncVar('time',nc=self._ncfile,time=ind,**kwargs)
        dimStr = date['_dimensions']
        date=date['data']
        
        outStr={'_dimensions':dimStr,'lon':lon,'lat':lat,'date':date}
        
        for param in par_list :
            dumVar = load_ncVar(param,nc=self._ncfile,time=ind,**kwargs) #Load variables
            dimStr=dumVar['_dimensions']
            
            #update dimensions
            curDim = [str(dimname) for dimname in dimStr.keys()[1:]] #[str(dimname) for dimname in self._ncfile.variables['LONGITUDE'].dimensions]
            curDimval = [dimStr[dim] for dim in curDim] #[len(self._ncfile.dimensions[dimname]) for dimname in curDim]
            flag = [(np.array(dimname) == outStr['_dimensions'].keys()).sum() == 0 for dimname in curDim] #find dimensions to update
            dimUpdate = np.array(curDim).compress(flag)
            for enum in enumerate(dimUpdate) : 
                self.message(3, 'Appending dimensions {0}:{1} to dataStructure'.format(enum[1],np.array(curDimval).compress(flag)[enum[0]]))
                outStr['_dimensions'].update({enum[1]:np.array(curDimval).compress(flag)[enum[0]]}) #Append new dimension
                outStr['_dimensions']['_ndims']+=1 #update dimension counts
            
            cmd = 'dumStr = {\''+param.lower()+'\':dumVar[\'data\']}'
            self.message(4, 'exec : '+cmd)
            exec(cmd)
            outStr.update(dumStr)
        
        id=np.repeat(sat_name,sz)
        
        
        outStr.update({'id':id})
        self._ncfile.close()
        
        return outStr

    def read_SLARES(self,filename,params=None,force=False,timerange=None,datatype=None,**kwargs):
        
        """
        READ_SLARES : Read AVISO Along-Track SLA older format products (RES files)
        
        @return: outStr {type:dict} Output data structure containing all recorded parameters as specificied by NetCDF file PARAMETER list.
        @author: Renaud Dussurget
        """
        
        #Open file
        self._filename = filename
        self._ncfile = ncfile(self._filename, "r")
        
        #Get delimiter
        if os.path.basename(filename).count('.') > os.path.basename(filename).count('_'): delim='.'
        else : delim = '_'
        
#        #Gat sat name
#        splitted=os.path.basename(filename).split(delim)
#        if (datatype == 'DT') | (datatype == 'NRT') : sat_name = splitted[2] if splitted[0] == 'nrt' else splitted[3]
#        if datatype == 'PISTACH' : sat_name = 'J2'
        
     
        #Get list of recorded parameters:
        par_list=[i.encode() for i in self._ncfile.variables.keys()]
        for i in ['BeginDates','Longitudes','Latitudes'] : par_list.pop(par_list.index(i))
        nparam=len(par_list)
        
        self.message(2,'Recorded parameters : '+str(nparam)+' -> '+str(par_list))
     
        lon = load_ncVar('Longitudes',nc=self._ncfile,**kwargs)
        lon['data'] = recale(lon['data'], degrees=True, zero_2pi=True) #shift longitudes
        lat = load_ncVar('Latitudes',nc=self._ncfile,**kwargs)

        
        #Extract within limits
        ind, flag = in_limits(lon['data'],lat['data'],limit=self.limit) #This ve
        dim_lon = lon['_dimensions']
        lat = lat['data'].compress(flag)
        lon = lon['data'].compress(flag)
        dist=cumulative_distance(lat, lon)
        
        sz=np.shape(lon)
        ndims=np.size(sz)
            
#        cyc = load_ncVar('Cycles',nc=self._ncfile,**kwargs)
        date = load_ncVar('BeginDates',nc=self._ncfile,**kwargs) #This vector is no longer in phase with indexed data.
        dimStr = date['_dimensions']
        dimStr.update({dim_lon.keys()[1]:len(lon)})
        dimStr['_ndims']+=1
        date=date['data']
        
        outStr={'_dimensions':dimStr,'lon':lon,'lat':lat,'date':date}
        
        for param in par_list :
            dumVar = load_ncVar(param,nc=self._ncfile,Data=ind,**kwargs) #Load variables
            dimStr=dumVar['_dimensions']
            
            #update dimensions
            curDim = [str(dimname) for dimname in dimStr.keys()[1:]] #[str(dimname) for dimname in self._ncfile.variables['LONGITUDE'].dimensions]
            curDimval = [dimStr[dim] for dim in curDim] #[len(self._ncfile.dimensions[dimname]) for dimname in curDim]
            flag = [(np.array(dimname) == outStr['_dimensions'].keys()).sum() == 0 for dimname in curDim] #find dimensions to update
            dimUpdate = np.array(curDim).compress(flag)
            for enum in enumerate(dimUpdate) : 
                self.message(3, 'Appending dimensions {0}:{1} to dataStructure'.format(enum[1],np.array(curDimval).compress(flag)[enum[0]]))
                outStr['_dimensions'].update({enum[1]:np.array(curDimval).compress(flag)[enum[0]]}) #Append new dimension
                outStr['_dimensions']['_ndims']+=1 #update dimension counts
            
            cmd = 'dumStr = {\''+param.lower()+'\':dumVar[\'data\']}'
            self.message(4, 'exec : '+cmd)
            exec(cmd)
            outStr.update(dumStr)
        
        id=np.repeat(self._ncfile.__dict__['Mission'].lower(),outStr['_dimensions']['Data'])
        
        
        outStr.update({'id':id})
        self._ncfile.close()
        
        return outStr


    def read_CTOH(self,filename,params=None,force=False,timerange=None,datatype=None,**kwargs):
        
        """
        READ_CTOH : Read CTOH data products
        
        @return: outStr {type:dict} Output data structure containing all recorded parameters as specificied by NetCDF file PARAMETER list.
        @author: Renaud Dussurget
        """
        
        #Open file
        self._filename = filename
        self._ncfile = ncfile(self._filename, "r")
        
        #Get delimiter
        delim = '.'
        
        #Gat sat name
        splitted=os.path.basename(filename).split(delim)
        sat_name = splitted[4]
     
        #Get list of recorded parameters:
        par_list=[i.encode() for i in self._ncfile.variables.keys()]
        for i in ['cycle','lon','lat'] : par_list.pop(par_list.index(i))
        nparam=len(par_list)
        
        self.message(2,'Recorded parameters : '+str(nparam)+' -> '+str(par_list))
     
        lon = load_ncVar('lon',nc=self._ncfile,**kwargs)
        lon['data'] = recale(lon['data'], degrees=True, zero_2pi=True) #shift longitudes
        lat = load_ncVar('lat',nc=self._ncfile,**kwargs)

        
        #Extract within limits
        ind, flag = in_limits(lon['data'],lat['data'],limit=self.limit)
        dim_lon = lon['_dimensions']
        lat = lat['data'].compress(flag)
        lon = lon['data'].compress(flag)
        dist=cumulative_distance(lat, lon)
        
        sz=np.shape(lon)
        ndims=np.size(sz)
        
        #Get date table and convert it to dimension
        date = load_ncVar('time',nc=self._ncfile,nbpoints=ind,**kwargs)
        
        dimStr = date['_dimensions']
        date=date['data']
        
        outStr={'_dimensions':dimStr,'lon':lon,'lat':lat,'date':date}
        
        for param in par_list :
            dumVar = load_ncVar(param,nc=self._ncfile,time=ind,**kwargs) #Load variables
            dimStr=dumVar['_dimensions']
            
            #update dimensions
            curDim = [str(dimname) for dimname in dimStr.keys()[1:]] #[str(dimname) for dimname in self._ncfile.variables['LONGITUDE'].dimensions]
            curDimval = [dimStr[dim] for dim in curDim] #[len(self._ncfile.dimensions[dimname]) for dimname in curDim]
            flag = [(np.array(dimname) == outStr['_dimensions'].keys()).sum() == 0 for dimname in curDim] #find dimensions to update
            dimUpdate = np.array(curDim).compress(flag)
            for enum in enumerate(dimUpdate) : 
                self.message(3, 'Appending dimensions {0}:{1} to dataStructure'.format(enum[1],np.array(curDimval).compress(flag)[enum[0]]))
                outStr['_dimensions'].update({enum[1]:np.array(curDimval).compress(flag)[enum[0]]}) #Append new dimension
                outStr['_dimensions']['_ndims']+=1 #update dimension counts
            
            cmd = 'dumStr = {\''+param.lower()+'\':dumVar[\'data\']}'
            self.message(4, 'exec : '+cmd)
            exec(cmd)
            outStr.update(dumStr)
        
        id=np.repeat(sat_name,sz)
        
        
        outStr.update({'id':id})
        self._ncfile.close()
        
        return outStr


    def track_list(self,*args):
        noargs = len(args) == 0
        return np.unique(self.track) if noargs else np.unique(self.track.compress(args[0]))
    
    def cycle_list(self,*args):
        noargs = len(args) == 0
        return np.unique(self.cycle) if noargs else np.unique(self.cycle.compress(args[0]))

    def reorder(self):
        #update time range with real value
        trange_str = self.time_range()[0]
        cycle_list = self.cycle_list()
        track_list=self.track_list()
            
        #Detect recurrent  lon/lat/track triplets
        triplets = np.unique(zip(*(self.lon, self.lat, self.track)))
        
        #Sort by track (we do not sort using other columns to preserve descending/ascending order)
        triplets=np.array(sorted(triplets, key=operator.itemgetter(2)))
        
        lon = triplets[:,0]
        lat = triplets[:,1]
        track = triplets[:,2]
        cycle=cycle_list
        
        N=len(lon)
        ntracks=len(track_list)
        ncycles=len(cycle_list)
        
        ind = np.arange(N)
        
    #    dst = atools.calcul_distance(lat, lon)
        
        #Get local index on nominal tracks
        xid = np.array([np.where((ln == lon) & (self.lat[i] == lat))[0][0] for i,ln in enumerate(self.lon) ])
        tid = np.array([np.where(c == cycle_list)[0][0] for c in self.cycle])
        
        #Get object attributes to reform
        varSize = np.array([np.size(self.__dict__[k]) for k in self.__dict__.keys()])
        varSize == self._dimensions['time']
        par_list=np.array(self.__dict__.keys())[varSize == self._dimensions['time']].tolist()
        
        #Refine param list
        for par in ['lon','lat','track','cycle'] :
            par_list.pop(par_list.index(par))
            exec('self.{0}={0}'.format(par))
            
        #Set new dimensions array
        self._dimensions={'_ndims':2,'npoints':N,'ncycles':ncycles}
        
        #Init output matrices using object fields
        for par in par_list :
            cmd = '{0} = np.ma.array(np.zeros((ncycles,N)),mask=np.ones((ncycles,N),dtype=bool),dtype=self.{0}.dtype)'.format(par)
            exec(cmd)
            cmd = '{0}[tid,xid]=self.{0}'.format(par)
            exec(cmd)
            exec('self.{0}={0}'.format(par))
            
#        #Rework date array
#        self.time=self.date
#        
#        self.time.m
        
        
#        ntime=len(self.time)
#
#        #Regrid time
#        rtime=time - time[0]
#        mn_dt = np.median(atools.deriv(rtime)) #1st time derivative approximation
#        mn_dt = np.median((rtime[1:] - rtime[:-1])[np.floor(rtime[1:] - rtime[:-1]) == np.floor(mn_dt)]) #Refined value
#        bins = np.ceil(rtime.max() / mn_dt) + 1
#        range=(0/2.,mn_dt * bins) - mn_dt/2
#        time = np.arange(0,rtime.max(),mn_dt)+time[0] #There is a long term trend in this due to error on mn_dt
#        hist,bin_edges=np.histogram(rtime, bins=bins, range=range)
#        
#        ntime = len(hist)
        
def bilinear_interpolation(x, y, points):
    '''Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

        >>> bilinear_interpolation(12, 5.5,
        ...                        [(10, 4, 100),
        ...                         (20, 4, 200),
        ...                         (10, 6, 150),
        ...                         (20, 6, 300)])
        165.0
        
        @note : source -> http://stackoverflow.com/questions/8661537/how-to-perform-bilinear-interpolation-in-python
    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    points = sorted(points)               # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

#    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
#        raise ValueError('points do not form a rectangle')
#    if not x1 <= x <= x2 or not y1 <= y <= y2:
#        raise ValueError('(x, y) not within the rectangle')

    try :
        if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
            raise ValueError('points do not form a rectangle')
        if not x1 <= x <= x2 or not y1 <= y <= y2:
            raise ValueError('(x, y) not within the rectangle')
        return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
           ) / ((x2 - x1) * (y2 - y1) + 0.0)
    except ValueError : return np.nan

def interp1d(x,Z,xout,spline=False,kind='linear',**kwargs):
    """    
    INTERP1D : Interpolate values from a 1D vector at given positions
    
    @param x: 1st dimension vector of size NX
    
    @author: Renaud DUSSURGET, LER/PAC, Ifremer La Seyne
    """

    linear = not spline
    
    
    nx=len(x)
    
    if linear :
        
        try :
            f = scipy.interpolate.interp1d(x, Z, kind=kind,bounds_error=False)
            Zout = f(xout)
        except RuntimeError : Zout = np.repeat(np.NaN,nx)

    else :
        tck = scipy.interpolate.splrep(x,Z,s=0)
        try : Zout = scipy.interpolate.splev(xout,tck,der=0)
        except RuntimeError : Zout = np.repeat(np.NaN,nx)

    return Zout
    

def interp2d(x,y,Z,xout,yout,**kwargs):
    """    
    INTERP2D : Interpolate values from a 2D matrix along a 1D trajectory
    
    @param x: 1st dimension vector of size NX
    @param y: 2nd dimension vector of size NY
    
    @author: Renaud DUSSURGET, LER/PAC, Ifremer La Seyne
    """
#    gx, gy = np.meshgrid(x,y)
#    gz = np.reshape(Z.transpose(),Z.size)
#     
#    points = zip(*(np.reshape(gx,gx.size),np.reshape(gy,gy.size)))
#    xi = zip(*(xout,yout))
#    
#    Zout = sc.interpolate.griddata(points, gz, xi, **kwargs)
#    return Zout

    gx = np.reshape(np.repeat(x,y.size),(x.size,y.size))
    gy = np.reshape(np.repeat(y,x.size),(y.size,x.size)).transpose((1,0))

    gz = Z.flatten()
     
    points = zip(*(gx.flatten(),gy.flatten())) 
    xi = zip(*(xout,yout))
    
    try : Zout = scipy.interpolate.griddata(points, gz, xi, **kwargs)
    except RuntimeError : Zout = np.NaN
#    Zout = sc.interpolate.griddata(points, gz, xi, **kwargs)
    return Zout
    



def interp3d(x,y,t,Z,xout,yout,tout,**kwargs):
    
    nx=x.size
    ny=y.size
    nt=t.size
    
    #Turn around matrix to (nx,ny,nt)
    order=(Z.shape.index(nx),Z.shape.index(ny),Z.shape.index(nt))
    Z=Z.transpose(order)
    
    N=xout.size
    
    for enum in enumerate(zip(*(xout,yout,tout))):
        time=enum[1][2]
        lon=enum[1][0]
        lat=enum[1][1]
        idx = np.arange(bisect.bisect_left(x,lon)-1,bisect.bisect_left(x,lon)+1)
        idy = np.arange(bisect.bisect_left(y,lat)-1,bisect.bisect_left(y,lat)+1)
        idt = np.arange(bisect.bisect_left(t,time)-1,bisect.bisect_left(t,time)+1)
        
        idx=idx.compress((idx >= 0) & (idx < nx))
        idy=idy.compress((idy >= 0) & (idy < ny))
        idt=idt.compress((idt >= 0) & (idt < nt))
        nxid=idx.size
        nyid=idy.size
        ntid=idt.size
        
#        print enum[0]
        
        if ntid == 0 :
            res = np.NaN
        elif ((ntid == 1) & (nxid == 2) & (nyid == 2)):
            zi = np.squeeze(Z[idx,:,:][:,idy,:][:,:,idt])
#            res = interp2d(x[idx], y[idy], zi, [lon], [lat])
            res = bilinear_interpolation(lon, lat, zip(*(x[idx[[0,0,1,1]]],y[idy[[0,1,0,1]]], np.reshape(zi,4))))
        elif ((ntid == 2) & (nxid == 2) & (nyid == 2)) :
            zi = Z[idx,:,:][:,idy,:][:,:,idt]
#            a=interp2d(x[idx], y[idy], zi[:,:,0], [lon], [lat])
#            b=interp2d(x[idx], y[idy], zi[:,:,1], [lon], [lat])
            a = bilinear_interpolation(lon, lat, zip(*(x[idx[[0,0,1,1]]],y[idy[[0,1,0,1]]], np.reshape(zi[:,:,0],4))))
            b = bilinear_interpolation(lon, lat, zip(*(x[idx[[0,0,1,1]]],y[idy[[0,1,0,1]]], np.reshape(zi[:,:,1],4))))
            res = ((b-a)/(t[idt[1]] - t[idt[0]]))*(time - t[idt[0]]) + a #linear interpolation
#            res = interp3d_core(x, y, t[idt], zi, [enum[1][0]], [enum[1][1]], [time])
        else : res = np.NaN
        if enum[0] == 0 : Zout=res
        else : Zout=np.append(Zout,res)  
        
    return Zout

def interp3d_core(x,y,t,Z,xout,yout,tout,**kwargs):
    """    
    INTERP3D : Interpolate values from a 3D matrix along a 1D trajectory
    
    @param x: 1st dimension vector of size NX
    @param y: 2nd dimension vector of size NY
    @param t: 3rd dimension vector of size NT
    
    @author: Renaud DUSSURGET, LER/PAC, Ifremer La Seyne
    """
    
    #this below can take a very LOOOOOONG time
    gx = np.reshape(np.repeat(x,y.size*t.size),(x.size,y.size,t.size))
    gy = np.reshape(np.repeat(y,x.size*t.size),(y.size,x.size,t.size)).transpose((1,0,2))
    gt = np.reshape(np.repeat(t,x.size*y.size),(t.size,x.size,y.size)).transpose((1,2,0))

    gz = Z.flatten()
     
    points = zip(*(gx.flatten(),gy.flatten(),gt.flatten())) 
    xi = zip(*(xout,yout,tout))
    
    Zout = sc.interpolate.griddata(points, gz, xi, **kwargs)
    return Zout    

def cart2polar(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return r, theta

def polar2cart(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

import operator
def lagrange(x, x_values, y_values):
    def _basis(j):
        p = [(x - x_values[m])/(x_values[j] - x_values[m]) for m in xrange(k) if m != j]
        return reduce(operator.mul, p)
    assert len(x_values) != 0 and (len(x_values) == len(y_values)), 'x and y cannot be empty and must have the same length'
    k = len(x_values)
    return sum(_basis(j)*y_values[j] for j in xrange(k))

#
#def lagrange(x):
#    tmp = scipy.poly1d([0])
#    result=scipy.poly1d([0])
#    
#    for i in x.keys():
#        numerator=scipy.poly1d([1])
#        denom = 1.0
#        for j in x.keys():
#            if (i != j):
#                tmp = scipy.poly1d([1,-j])
#                numerator = numerator * tmp
#                denom = denom * (i - j)
#        tmp = (numerator/denom) * x.get(i)
#        result = result + tmp
#
#return result

#from scipy.interpolate import lagrange, BarycentricInterpolator, barycentric_interpolate
def deprecated_deriv(*args):
    if len(args) == 1 :
        y = args[0]
        x = np.arange(len(y))
    if len(args) == 2 :
        x = args[0]
        y = args[1]
    
    dx = x[1:] - x[:-1]
    dy = y[1:] - y[:-1]
    
#    d2 = lagrange(np.arange(len(x)),np.arange(len(dx))+0.5,dy/dx)

    dx = np.append(np.array(2*dx[1] - dx[0]),dx)
    dx = np.append(dx,np.array(2*dx[-1] - dx[-2]))
    
    dy = np.append(np.array(2*dy[1] - dy[0]),dy)
    dy = np.append(dy,np.array(2*dy[-1] - dy[-2]))
    
#    dxint=lagrange(np.arange(nx),np.arange(nx-1)+0.5,dx)
#    dyint = lagrange(np.arange(nx),np.arange(nx-1)+0.5,dy)
    
#    d = BarycentricInterpolator(dx)
#    dxint = lagrange(np.arange(len(dx)-1) + 0.5,np.arange(len(dx)),dx)
#    dyint = lagrange(np.arange(len(dx)-1) + 0.5,np.arange(len(dx)),dy)
    d = lagrange(np.arange(len(x)) + 0.5,np.arange(len(dx)),dy/dx)
    
#    p = barycentric_interpolate(np.arange(len(dx))+0.5,dy/dx,np.arange(len(dx)))
#    p(np.arange(len(dx)))
#    p(np.arange(len(dx)))
    
    return d


def deriv(*args):
    """
    ; Copyright (c) 1984-2009, ITT Visual Information Solutions. All
    ;       rights reserved. Unauthorized reproduction is prohibited.
    ;
    
    ;+
    ; NAME:
    ;    DERIV
    ;
    ; PURPOSE:
    ;    Perform numerical differentiation using 3-point, Lagrangian 
    ;    interpolation.
    ;
    ; CATEGORY:
    ;    Numerical analysis.
    ;
    ; CALLING SEQUENCE:
    ;    Dy = Deriv(Y)         ;Dy(i)/di, point spacing = 1.
    ;    Dy = Deriv(X, Y)    ;Dy/Dx, unequal point spacing.
    ;
    ; INPUTS:
    ;    Y:  Variable to be differentiated.
    ;    X:  Variable to differentiate with respect to.  If omitted, unit 
    ;        spacing for Y (i.e., X(i) = i) is assumed.
    ;
    ; OPTIONAL INPUT PARAMETERS:
    ;    As above.
    ;
    ; OUTPUTS:
    ;    Returns the derivative.
    ;
    ; COMMON BLOCKS:
    ;    None.
    ;
    ; SIDE EFFECTS:
    ;    None.
    ;
    ; RESTRICTIONS:
    ;    None.
    ;
    ; PROCEDURE:
    ;    See Hildebrand, Introduction to Numerical Analysis, Mc Graw
    ;    Hill, 1956.  Page 82.
    ;
    ; MODIFICATION HISTORY:
    ;    Written, DMS, Aug, 1984.
    ;    Corrected formula for points with unequal spacing.  DMS, Nov, 1999.
    ;-
    ;
    ; on_error,2              ;Return to caller if an error occurs
    """
    x = args[0]
    n = x.size
    if n < 3 : raise 'Parameters must have at least 3 points'


    if (len(args) == 2) :
        y=args[1]
        if n != y.size : raise 'Vectors must have same size'
    
        #;df/dx = y0*(2x-x1-x2)/(x01*x02)+y1*(2x-x0-x2)/(x10*x12)+y2*(2x-x0-x1)/(x20*x21)
        #; Where: x01 = x0-x1, x02 = x0-x2, x12 = x1-x2, etc.
    
        if isinstance(x,np.ma.masked_array) :  x = x.data                # Convert masked arrays to classic arrays
        if not isinstance(x.data[0],np.float) : x.astype(np.float) #;If not floating type, ensure floating... 
        
        x12 = x - np.roll(x,-1)                                      #;x1 - x2
        x01 = np.roll(x,1) - x                                       #;x0 - x1
        x02 = np.roll(x,1) - np.roll(x,-1)                           #;x0 - x2
    
        d = np.roll(y,1) * (x12 / (x01*x02)) \
            + y * (1./x12 - 1./x01) \
            - np.roll(y,-1) * (x01 / (x02 * x12))                    #Middle points
        
        
        #Formulae for the first and last points:
        d[0] = y[0] * (x01[1]+x02[1])/(x01[1]*x02[1]) \
            - y[1] * x02[1]/(x01[1]*x12[1]) \
            + y[2] * x01[1]/(x02[1]*x12[1])                          #;First point
        n2 = n-2
        d[n-1] = -y[n-3] * x12[n2]/(x01[n2]*x02[n2]) \
            + y[n-2] * x02[n2]/(x01[n2]*x12[n2]) \
            - y[n-1] * (x02[n2]+x12[n2]) / (x02[n2]*x12[n2])             #;Last point

    #Equally spaced point case
    else :
        d = (np.roll(x,-1) - np.roll(x,1))/2.
        d[0] = (-3.0*x[0] + 4.0*x[1] - x[2])/2.
        d[n-1] = (3.*x[n-1] - 4.*x[n-2] + x[n-3])/2.

    return d

def track_orient(x,y,orient=False):
#    ;Calculate track orientation (on 3 points lagrangian interpolation - see DERIV help page)
    
    getorient=orient
    
    dxy=deriv(x,y)
    
    dx=np.median(deriv(x))
    dy=np.median(deriv(y))

    orient = np.arctan(dxy) if dx > 0 else np.arctan(dxy) + np.pi #Rotate backward (180°)
                                                                  #orientation when track goes
                                                                  #westward
  
    orient = np.remainder(orient,2.0*np.pi)
#    if (medorient < 0) : orient += 2.0*np.pi
#    if (medorient > 2.0*np.pi) THEN orient -= 2D*!dpi
    #Sense is positive for ascending tracks
    sense = True if np.median(orient) < np.pi else False
  
    if getorient : return sense, orient
    else : return sense

def coriolis(lat):

    sideral_day = 86164.1 #sideral day in seconds
    omega = (2. * np.pi) / sideral_day #earth rotation rate
  
    #Compure Coriolis frequency and returns it
    return 2. * omega * np.sin(np.deg2rad(lat))  

def gravity(*args):
    return seawater.gibbs.grav(*args)
    
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
    ; @returns Geostrophic velocity component, positive to the right of the track
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
    
    pl04 = kwargs.pop('pl04') if kwargs.has_key('pl04') else False
    strict = kwargs.pop('strict') if kwargs.has_key('strict') else False
    
    if pl04 : return powell_leben_filter_km(lon,lat,nu,verbose=True,**kwargs)
    
    #If strict option is set to True, compute gradients at mid-distance between points
    if strict :
        lon = (lon[1:] - lon[:-1])/2. + lon[0:-1]
        lat = (lat[1:] - lat[:-1])/2. + lat[0:-1]
    
    
    #Compute gravitational & coriolis forces
    g = gravity(lat)
    f = coriolis(lat)
    
    #Compute SSH 1st derivative
#    dh = deriv(dst,nu) #(deriv is very bad...)
    
    dh = (nu[1:] - nu[:-1])/(dst[1:] - dst[:-1]) if strict else deriv(dst,nu)
      
    #Compute geostrophy
#    print f
#    print g
#    print dh
    ug = - (g*dh) / (f)
  
    if (not track_orient(lon,lat)) : ug*=-1
    
    return (lon,lat,ug) if strict else ug

def loess(h,x,cut):
    """
    #===========================================================================
    # ;+
    # ; LOESS : Filter a signal H along X with a cutoff frequency CUT<br /><br />
    # ; 
    # ; Reference : Schlax, M. G., et D. B. Chelton (1992), Frequency Domain Diagnostics<br />
    # ;   for Linear Smoothers, Journal of the American Statistical Association, <br />
    # ;   87(420), 1070-1081. 
    # ; 
    # ; @param h {in}{required}{type:NUMERIC} Array of height anomalies
    # ; @param x {in}{required}{type:NUMERIC} Array of corresponding distances
    # ; @param cut {in}{required}{type:NUMERIC} Cutoff distance in units of X.
    # ; @param nan {in}{optional}{type:LOGICAL}{default:FALSE} Use NaNs<br />
    # ;          
    # ; @returns Hfilt (H filtered at cutoff CUT)
    # ; 
    # ; 
    # ; @uses exist
    # ; 
    # ; @author Renaud DUSSURGET, LEGOS/CTOH
    # ; @history Created Jan. 2010 by R.DUSSURGET from loess.c develloped by Mog2d team<br />
    # ;          May 2010 : Take NaNs into account for computation<br />
    # ;                     Checked the filters formulation (Schlax and Chelton, 1992) <br /><br />
    # ;
    # ;  LEGACY :<br />
    # ;  void dloess1d_irregular(int n, double l_c, double *h, double *x, double mask,double *lf,int *rstatus)<br /><br /> 
    # ;  
    # ;  ./varReduction.cpp:          dloess1d_irregular(n,cut,h,t,(double) satData.mask,lf,&status);<br />
    # ;  ./functions.extraction.cpp:  dloess1d_irregular(n,cut,h2,x,mask,lf,&status);<br /><br />
    # ;  
    # ;  rq: n-> nb points dans série<br />
    # ;  l_c -> cut<br />
    # ;  *h -> sla<br />
    # ;  *x -> dimensional array<br />
    # ;  mask -> valeur du masque<br />
    # ;  *lf -> partie basse freq<br />
    # ;  *status -> retour<br />
    # ;-
    #===========================================================================
    """
          
    n=np.size(h)

    #If masked_array --> Get data and replace masked values by NaNs
    if isinstance(h,np.ma.masked_array): h=mask2NaN(h)
    else : h=np.ma.masked_array(h,mask=np.zeros(n,dtype=bool))

    #Start filtering if there are more than 1 point
    if n == 1 : lf = h[0]
    else : pass
    
    #Get the half-width of filter window
    l_c=cut/2.0          
    
    #Initiate weights and outputs
    w=np.zeros(n)
    lf = np.repeat(np.NaN,n)
    
    #Get valid points in serie
    flag = ~np.isnan(h)
    fcnt = flag.sum()
    fnt = np.arange(n).compress(flag)
    
    x = np.ma.masked_array(x,mask=(np.isnan(h) + h.mask) )
#        WHERE(FINITE(h),fcnt)

    #Vectorized form
    q = np.reshape(np.repeat(x,n),(n,n))
    q = (np.abs(q - q.transpose())/l_c).reshape(n*n)
    
    s = (1.0 - q*q*q)
    
#    outOfFilter_flag = q > 1.0
    outOfFilter_flag = np.greater(q,1)
    outCnt = np.nansum(outOfFilter_flag)
    outOfFilter = np.arange(n*n).compress(outOfFilter_flag)
    if (outCnt > 0) :
        s[outOfFilter]=0.0
        s.mask[outOfFilter]=True
    
    w=(s*s*s).reshape((n,n))  
    sumvar=np.nansum(w,1)

    #Compute current height (heights times weigths, divided by sum of weights)
    try :
        lf=np.nansum(w*np.repeat(h,n).reshape((n,n)),0)/sumvar
    except (ValueError):
        lf=np.nansum(w*np.repeat(h,n).reshape((n,n)),0)/sumvar #Retry

    return lf


def loess_inline(h,x,cut,nan=True):
    """
    #This is a raw, inline version of the loess filtering function. Runs much more slowlier.
    """
    
    n=np.size(h)

    #Start filtering if there are more than 1 point
    if n == 1 : lf = h[0]
    else : pass
    
    #Get the half-width of filter window
    l_c=cut/2.0          
    
    #Initiate weights and outputs
    w=np.zeros(n)
    lf = np.repeat(np.NaN,n)
    
    #Get valid points in serie
    flag = ~np.isnan(h)
    fcnt = flag.sum()
    fnt = np.arange(n).compress(flag)
#        WHERE(FINITE(h),fcnt)

    #Loop from 1st to last valid point
    for i in fnt :

        icur=i
      
        #Get Q (distance from central point divided by half-window width)
        q=(np.abs((x-x[icur])/l_c))
            
        #Compute S from Q³
        s = (1.0 - q*q*q)
        
     
        #Set all values out of filter window to 0
        outOfFilter_flag = q > 1.0
        outCnt = outOfFilter_flag.sum() 
        outOfFilter = np.arange(n).compress(outOfFilter_flag)
        if (outCnt > 0) : s[outOfFilter]=0.0
      
        #Compute weights (S³ - Tricube)
        w=s*s*s  
        sumvar=np.nansum(w)

        #Compute current height (heights times weigths, divided by sum of weights)
        lf[icur]=np.nansum(w*h)/sumvar

    return lf

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
    ; @returns Geostrophic velocity component, positive to the right of the track
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
    p = kwargs.pop('p') if kwargs.has_key('p') else 12.
    q = kwargs.pop('q') if kwargs.has_key('q') else 12.
    
    verbose = kwargs.pop('verbose') if kwargs.has_key('verbose') else False
        
    dt = np.median(dst[1:] - dst[:-1])   #Median point-to-point distance in meters
      
    #Get the corresponding number of points rounded to the nearest integer
    pn=np.fix(p*1e3/dt).astype(int)
    qn=np.fix(q*1e3/dt).astype(int)
    
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
    if (not track_orient(lon,lat)) : ug*=-1
    
    return ug

def rms(array):
    """
    ;+
    ; RMS : Returns root-mean-squared deviation of an array
    ; 
    ; @param array {in}{required}{type:NUMERIC} 2-dimensionnal array. Time dimension<br />
    ;   should be the last dimension.
    ;   
    ; @returns RMS array (1 value per time serie)
    ; 
    ; @author Renaud DUSSURGET, LEGOS/CTOH
    ;-
    """
    #;time should be the last dimension
  
    n=array.size
    
#    IF sz[ndims-1] EQ 1 THEN RETURN, 0.
    
    nans = np.where(~np.isnan(array))
    cnt=len(nans)
    
    if (cnt == 0) : return np.NaN
    
    nval=np.nansum(~np.isnan(array))
    
    mn = np.nansum(array) / nval

    return np.sqrt(np.nansum((array - mn)**2.0)/nval)
#    RETURN, SQRT((TOTAL((array - REBIN(mn,sz))^2D,ndims,/DOUBLE,/NAN) / nval))
#;  RETURN, SQRT((TOTAL((array - REBIN(mn,sz))^2D,ndims,/DOUBLE,/NAN) / sz[ndims -1]))


def histogram_indices(hist,R):
    ind = []
    for k in np.arange(len(hist)) : ind.append(R[R[k] : R[k+1]]) if hist[k] > 0 else ind.append([])
    return ind
#    for k in notempty : ind.append(R[R[k] : R[k+1]])

def mask2NaN(array):
    n=array.size
    if array.mask.size != n :
        array.mask = np.zeros(n,dtype=bool)
#        raise np.ma.MaskError("[mask2NaN]Error : mask length is not consistent with data")
#        array.mask = np.ones(n,dtype=bool)
    array.data[np.arange(n).compress(array.mask)] = np.NaN
    return array

def nanargmin(array,axis=None):
    shape=array.shape
    nx=shape[0]
    ny=shape[1]
    if axis is None : return np.nanargmin()
    
def grid_track(lat,lon,sla,remove_edges=True,backbone=None):
    """
    # GRID_TRACK
    # @summary: This function allow detecting gaps in a set of altimetry data and rebin this data regularlyy, with informations on gaps.
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
        outlon[empty] = interp1d(ok, outlon[ok], empty, kind='cubic')
        outlat[empty] = interp1d(ok, outlat[ok], empty, kind='cubic')
        outsla[empty] = interp1d(ok, outsla[ok], empty, spline=True)
        
    
    
    #Get gap properties
    ind=np.arange(len(hist))
    dhist=(hist[1:] - hist[:-1])
    st=ind.compress(dhist==-1)+1
    en=ind.compress(dhist==1)
    gaplen=(en-st) + 1
    ngaps=len(st)
    gapedges=np.array([st,en])
    
    #Get empty bin flag
    interpolated=~hist.astype('bool')
    
    return outdst, outlon, outlat, outsla, gaplen, ngaps, gapedges, interpolated

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


def cnes2modis(cnes_date,YYYYDDDHHMM=None,YYYYDDD=True):
    """
    #+
    # CNES2MODIS : Convert CNES julian days to MODIS data format (YYYYDDD)
    # 
    # @author: Renaud DUSSURGET (LEGOS/CTOH)
    # @history: Created by RD on 2/12/2011
    #           Adapted to Python on 29/10/2012 by RD (now at LER PAC/IFREMER)
    #
    #-
    """

    #Setup defaults
    if YYYYDDDHHMM is None : YYYYDDDHHMM = False
    if YYYYDDD : YYYYDDDHHMM = False
    if YYYYDDDHHMM : YYYYDDD = False

    str,obj=cnes_convert(cnes_date)
    dat=[o.strftime("%Y%j%H%M") for o in obj]
    
    return dat

def modis2cnes(modis_date):
    """
    #+
    # MODIS2CNES : Convert MODIS date format to CNES JULIAN DAYS (YYYYDDD or YYYYDDDHHMM)
    # 
    # @author: Renaud DUSSURGET (LER PAC/IFREMER)
    # @history: Created by RD on 29/10/2012
    #
    #-
    """

    if not isinstance(modis_date,list) : modis_date=[modis_date]

    if len(modis_date[0]) ==  7 :
        obj=[datetime.datetime.strptime(d,"%Y%j") for d in modis_date] 
    else :
        obj=[datetime.datetime.strptime(d,"%Y%j%H%M") for d in modis_date]
    
#    str,obj=cnes_convert(modis_date)
    dat=[o.strftime("%d/%m/%Y-%H:%M") for o in obj]
    
    return cnes_convert(dat)

def modis_filename2modisdate(modis_fname):
    """
    #+
    # MODIS_FILENAME2DATE : Convert MODIS file name to MODIS date
    # 
    # @author: Renaud DUSSURGET (LER PAC/IFREMER)
    # @history: Created by RD on 29/10/2012
    #
    #-
    """
    
    if not isinstance(modis_fname,list) : modis_fname=[modis_fname]
    
    return [os.path.splitext(os.path.basename(m))[0][1:12] for m in modis_fname]

def modis_filename2cnes(modis_fname):

    modis_date=modis_filename2modisdate(modis_fname)
    return modis2cnes(modis_date)
    
def detrend(X,Z,deg=1):
    if Z.shape == 1 :
        nt = 1
        valid=np.array([True])
    else :
        nt=Z.shape[0]
        valid=(~Z.mask).sum(axis=1)
    a=np.arange(deg+1)
    for t in np.arange(nt)[valid > 0]:
        fit=np.polyfit(X[~Z[t,:].mask],Z[t,:][~Z[t,:].mask], deg)
        for d in a : Z[t,:]-=np.power(X,a[::-1][d])*fit[d]
    return Z
    