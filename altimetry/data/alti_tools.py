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
import scipy.interpolate
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from mpl_toolkits.basemap import Basemap
from bisect import bisect_left
from Tkinter import Tk
from netCDF4 import Dataset as ncfile
import bisect
import seawater.gibbs

import glob
import os

import hydro_tools as htools
import esutils_stat as es


def in_limits(lon, lat, limit):
    
    lats = limit[0]
    lons = limit[1]
    late = limit[2]
    lone = limit[3]

    ind=[]

    # find out lon between S and E
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
    
    flag = lon_flag & lat_flag
    
    # compress lat list according to lon constraint 
    ind = np.arange(len(lon))
    ind = ind.compress(flag)

    # compress indexes
#    ind = ind.compress(flag)
    
    return ind, flag

def recale_limits(limitin, zero_2pi=True, minpi_pi=None): 
    
    #defaults
    if (minpi_pi is not None) : zero_2pi != minpi_pi
    
    limitout=limitin
    if isinstance(limitin,np.ndarray) : limitout[[1,3]]=recale(limitin[[1,3]], zero_2pi=zero_2pi, minpi_pi=minpi_pi,degrees=True) #Always in degrees!
    else :
        for i in [1,3] : limitout[i]=recale(limitin[i], zero_2pi=zero_2pi, minpi_pi=minpi_pi,degrees=True)

    return limitout

def recale(angle, degrees=False, zero_2pi=True, minpi_pi=None) :

    #defaults
    if (minpi_pi is not None) : zero_2pi = not minpi_pi
    radians = not degrees
  
    modulus= 360.0 if degrees else 2.0 * np.pi
  
    limits=np.array([0.0, 2.0*np.pi])
    if minpi_pi : limits-=np.pi
    if degrees : limits = np.rad2deg(limits)
    
    over=angle > limits[1]
    under=angle <= limits[0]
    angle+=under*modulus
    angle-=over*modulus
  
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
        if np.size(z) > 1 : self.setup_map(lon, lat, z, s=s, edgecolor=edgecolor)
        else :
            if z != 0 :  self.setup_map(lon, lat, z, s=s, edgecolor=edgecolor)
        
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
        
    def show(self):
        plt.show()
        
    def savefig(self,*args,**kwargs):
        plt.savefig(*args,**kwargs)
    
#    def set_local(self,VarName,VarValue):
#        if VarName not in locals() : exec('self.'+VarName+' = VarValue')
        
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
        strlist = np.array([x.split("/", 2) for x in argin], dtype=int)
        datelist = [datetime.datetime(strlist[x, 2], strlist[x, 1], strlist[x, 0]) for x in np.arange(narg)]
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
        
        
    

#    args={lambda:julian:None,
#          lambda:calendar:None,
#          lambda:nasa:None,
#          lambda:string:None,
#          lambda:quiet:None,
#          lambda:calformat:None}

#  ON_ERROR, 2
#
#  ;Set up default
#  IF ~KEYWORD_SET(julian)  AND KEYWORD_SET(calendar) THEN BEGIN
#    julian = 0 & calendar = 1
#  ENDIF
#  IF ~KEYWORD_SET(calendar) AND KEYWORD_SET(julian) THEN BEGIN
#    julian = 1 & calendar = 0
#  ENDIF
#  IF ~KEYWORD_SET(calendar) AND ~KEYWORD_SET(julian) THEN BEGIN
#    type = SIZE(in,/TYPE)
#    IF ((type GE 1 AND type LE 5) OR (type GE 12) ) THEN BEGIN
#      julian = 0 & calendar = 1    
#    ENDIF ELSE IF (type EQ 7) THEN BEGIN
#      julian = 1 & calendar = 0
#    ENDIF
#  ENDIF
#  IF ~KEYWORD_SET(nasa) THEN nasa = 0
#  IF ~(KEYWORD_SET(string)) THEN string = 0
#  IF ~KEYWORD_SET(quiet) THEN quiet = 0
#  
#  IF (~exist(calformat)) THEN calformat="dd/mm/yyyy"
#
#  IF (N_ELEMENTS(in) EQ 1) $
#    THEN vec = 0 $
#    ELSE vec = 1
#
#  ;Detect output type
#  IF (julian) AND ~(calendar) $
#    THEN jul = 1 $
#    ELSE IF ~(julian) AND (calendar) $
#      THEN jul = 0 $datetime.date.fromordinal(argin)
#      ELSE jul = 1
#
#  ;Get calendar format and ordering
#  calsep=STREGEX(CALFORMAT,"[^a-zA-Z0-9_]",/EXTRACT)
#  
#  calorder=[STREGEX(CALFORMAT,'dd'),STREGEX(CALFORMAT,'mm'),STREGEX(CALFORMAT,'yyyy')]
#  IF (MIN(CALORDER) LT 0) THEN MESSAGE, STRING(CALFORMAT, FORMAT='(%"Invalid format %s")')
#  calorder=SORT(calorder)
#  juldayorder=calorder[[1,0,2]]
#  
#
#  ;Get reference time, 
#  IF ~(nasa) $
#    THEN ref = JULDAY(01,01,1950) $
#    ELSE ref = JULDAY(01,01,1958)
#    
#  SWITCH (vec * 2^1 + jul * 2^0) OF
#    ;Vector + julian output
#    3: BEGIN
#      date = MAKE_ARRAY(3,N_ELEMENTS(in))
#      FOR i = 0L, N_ELEMENTS(in) - 1 DO date[*,i] =  STRSPLIT(in[i],calsep,/EXTRACT)
#;      out = JULDAY(date[1,*],date[0,*],date[2,*]) - ref
#      out = JULDAY(date[juldayorder[0],*], date[juldayorder[1],*], date[juldayorder[2],*]) - ref
#      BREAK
#    END
#    ;Vector + calendar output
#    2:BEGIN
#      CALDAT, ref + in, mm, dd, yyyy
#      out=TRANSPOSE(([[dd],[mm],[yyyy]])[*,[calorder]]) ;WARNING TRANSPOSITION AND CONCATENATION
#                                                        ;SUBJECT TO INPUT VECTOR DIMENSIONS!!!get arguments
#      IF (string) $
#        THEN out = STRTRIM(out[calorder[0],*],2)+calsep+STRTRIM(out[calorder[1],*],2)+calsep+STRTRIM(out[calorder[2],*],2)
#      BREAK
#    END
#    ;Scalar + julian output
#    1: BEGIN
#      date = STRSPLIT(in,calsep,/EXTRACT)
#;      out = JULDAY(date[1],date[0],date[2]) - ref
#      out = JULDAY(date[juldayorder[0]], date[juldayorder[1]], date[juldayorder[2]]) - ref
#      IF (~quiet) THEN PRINT, out, FORMAT='(%"Julian day : %d ",$)'
#      IF ~(nasa) $
#        THEN IF (~quiet) THEN PRINT, "(CNES time)" $
#        ELSE IF (~quiet) THEN PRINT, "(NASA time)"
#      BREAK
#    END
#    ;Scalar + calendar output
#    0: BEGIN
#      CALDAT, ref + in, mm, dd, yyyy
#      out=TRANSPOSE(([[dd],[mmdef main():
#                                ;strout=STRTRIM(out[calorder[0]],2)+calsep+STRTRIM(out[calorder[1]],2)+calsep+STRTRIM(out[calorder[2]],2)
#      IF (string) $
#        THEN out = STRTRIM(out[calorder[0]],2)+calsep+STRTRIM(out[calorder[1]],2)+calsep+STRTRIM(out[calorder[2]],2)
#      IF (~quiet) THEN PRINT, dd,mm,yyyy, FORMAT='(%"Calendar date : '+STRJOIN((['%d','%d','%4d'])[calorder],calsep)+'")'
#      BREAK
#    END
#  ENDSWITCH


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
    i = bisect_left(t, x)
    if i == len(t) : i-=1
    if t[i] - x > 0.5:
        i-=1
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
    
    def read(self,filename,**kwargs):
        
        fname,extension = os.path.splitext(filename)
        outStr=self.read_slaext(filename,**kwargs)
        self.update_fid_list(os.path.basename(filename),outStr['_dimensions']['time'])
        
        return outStr
    
#    def read_CorrSSH(self,filename,limit=None):
#        
#        #Open file
#        self._filename = filename
#        self._ncfile = ncfile(self._filename, "r")
##        lon = np.array(self._ncfile.variables['Longitudes'][:])
##        lat = np.array(self._ncfile.variables['Latitudes'][:])
#        lon = self._ncfile.variables['Longitudes'][:]
#        lat = self._ncfile.variables['Latitudes'][:]
#        
#        sat=os.path.basename(filename).split('_')[2]
#        
#        #Extract within limits
#        id, flag = in_limits(lon,lat,limit)
#        lat = lat.compress(flag)
#        lon = lon.compress(flag)
#        
#        #reconstruct track vector
#        tracklist=(np.ma.masked_array(self._ncfile.variables['Tracks'][:])).compressed()
#        ntracks=tracklist.size
#        nbpoints=(np.ma.masked_array(self._ncfile.variables['NbPoints'][:])).compressed()
#        tracks=[]
#        for i in np.arange(ntracks) :
#            tracks=np.append(tracks,np.repeat(tracklist[i],nbpoints[i]))
#        
#        tracks=tracks.compress(flag)
#            
##        dum=np.reshape([np.arange(i) for i in nbpoints],nbpoints.cumsum().max())
##        nbpoints.cumsum().max()
#              
#        nid=len(id)
#        sla = (self._ncfile.variables['CorSSH'][id] - self._ncfile.variables['mean_sea_surface'][id]) * 100. #Turn everything in cm
#        date = np.float64(self._ncfile.variables['TimeDay'][id,0])+(np.float64(self._ncfile.variables['TimeSec'][id,0])/86400.0) +np.float64(self._ncfile.variables['TimeMicroSec'][id,0])/(86400*1e6)
#
#        self._ncfile.close()
#        
#        return {'_sat':sat,'_nbpoints':nid,'lon':lon,'lat':lat,'date':date,'tracks':tracks,'sla':sla}

    def read_slaext(self,filename,params=None,force=False,timerange=None,**kwargs):
        
        """
        READ_SLAEXT : Read AVISO Along-Track SLA regional products
        
        @return: outStr {type:dict} Output data structure containing all recorded parameters as specificied by NetCDF file PARAMETER list.
        @author: Renaud Dussurget
        """
        
        #Open file
        self._filename = filename
        self._ncfile = ncfile(self._filename, "r")
        
        #Gat sat name
        splitted=os.path.basename(filename).split('_')
        sat_name = splitted[2] if splitted[0] == 'nrt' else splitted[3]
        
     
        #Get list of recorded parameters:
        par_list=[i.encode() for i in self._ncfile.variables.keys()]
        for i in ['time','longitude','latitude'] : par_list.pop(par_list.index(i))
        nparam=len(par_list)
        
        self.message(1,'Recorded parameters : '+str(nparam)+' -> '+str(par_list))
     
        lon = self.load_ncVar('longitude',**kwargs)
        lon['data'] = recale(lon['data'], degrees=True, zero_2pi=True) #shift longitudes
        lat = self.load_ncVar('latitude',**kwargs)

        
        #Extract within limits
        ind, flag = in_limits(lon['data'],lat['data'],limit=self.limit)
        dim_lon = lon['_dimensions']
        lat = lat['data'].compress(flag)
        lon = lon['data'].compress(flag)
        dist=cumulative_distance(lat, lon)
        
        sz=np.shape(lon)
        ndims=np.size(sz)
            
        date = self.load_ncVar('time',time=ind,**kwargs)
        dimStr = date['_dimensions']
        date=date['data']
        
        outStr={'_dimensions':dimStr,'lon':lon,'lat':lat,'date':date}
        
        for param in par_list :
            dumVar = self.load_ncVar(param,time=ind,**kwargs) #Load variables
            dimStr=dumVar['_dimensions']
            
            #update dimensions
            curDim = [str(dimname) for dimname in dimStr.keys()[1:]] #[str(dimname) for dimname in self._ncfile.variables['LONGITUDE'].dimensions]
            curDimval = [dimStr[dim] for dim in curDim] #[len(self._ncfile.dimensions[dimname]) for dimname in curDim]
            flag = [(np.array(dimname) == outStr['_dimensions'].keys()).sum() == 0 for dimname in curDim] #find dimensions to update
            dimUpdate = np.array(curDim).compress(flag)
            for enum in enumerate(dimUpdate) : 
                self.message(2, 'Appending dimensions {0}:{1} to dataStructure'.format(enum[1],np.array(curDimval).compress(flag)[enum[0]]))
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
            f = scipy.interpolate.interp1d(x, Z, kind=kind)
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
    
    try : Zout = sc.interpolate.griddata(points, gz, xi, **kwargs)
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
    
def geost_1d(*args,**kwargs) : #(lon,lat,dst,nu):
    
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
    lf=np.nansum(w*np.repeat(h,n).reshape((n,n)),0)/sumvar

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
    ug = - (g*dh) / (f)
  
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
        raise np.ma.MaskError("[mask2NaN]Error : mask length is not consistent with data")
#        array.mask = np.ones(n,dtype=bool)
    array.data[np.arange(n).compress(array.mask)] = np.NaN
    return array

def nanargmin(array,axis=None):
    shape=array.shape
    nx=shape[0]
    ny=shape[1]
    if axis is None : return np.nanargmin()
    