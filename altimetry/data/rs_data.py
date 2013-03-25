import os
import glob

import datetime

import numpy as np
#import matplotlib.pyplot as plt
from pyhdf.SD import SD, SDC
#from altimetry.tools.spatial_tools import in_limits
#from altimetry.tools.interp_tools import interp1d,interp2d2d
#from altimetry.tools.dates import cnes_convert, modis2cnes
#from altimetry.tools.others import nearest
from altimetry.tools import in_limits,interp1d,interp2d2d,cnes_convert, modis2cnes, nearest

#import alti_tools as atools


#Read HDF data
##############
def modis_sst(file_name,
              limit=None,
              flagLevel=None,
              param='sst'):
    
    #Setup flag tables for L2 data
    flagTable=np.zeros(32,dtype=bool)
    if flagLevel == 0 : flagTable[[1]]=True
    elif flagLevel == 1 :flagTable[[  1,  3,  5,            15,16,      19,         25,26]]=True
    elif flagLevel == 2 : flagTable[[0,1,  3,4,5,8,9,10,  14,15,16,      19,21,22,23,25,26]]=True
    elif flagLevel >= 3 : flagTable[[0,1,  3,4,5,8,9,10,  14,15,16,      19,21,22,23,25,26]]=True
    flags=np.where(flagTable)[0]

    #Read MODIS HDF4 data
    f = SD(file_name, SDC.READ)
    fattr=f.attributes()
    
    #Get dimensions
    nScans=fattr['Number of Scan Lines']
    sCtl=fattr['Number of Scan Control Points']
    pCtl=fattr['Number of Pixel Control Points']
    nPix=fattr['Pixels per Scan Line']

    #Load coordinates 
    #UPDATE THIS SECTION TO ALLOW CLEAN LOAD OF THE IMAGE...
    lonsel=f.select('longitude')
    latsel=f.select('latitude')
#    
#    pCtl_ind=np.arange(pCtl,dtype=np.float32)
#    p_ind=(pCtl-1)*np.arange(nPix,dtype=np.float32)/(nPix-1)
#    
#    sCtl_ind=np.arange(sCtl,dtype=np.float32)
#    s_ind=(sCtl-1)*np.arange(nScans,dtype=np.float32)/(nScans-1)
    
#    dum=interp2d2d(sCtl_ind, pCtl_ind, lonsel.get(), s_ind, p_ind)
#    dumlon=interp1d(p_ind, lonsel.get(), pCtl_ind, spline=True)
#    dumlat=interp1d(p_ind, lonsel.get(), pCtl_ind, spline=True)

#    shlon=shlat=lon.dimensions().values()

    #Shrink image
    #############
    
    
    #scene dimensions
    
    info_sst=f.datasets()[param]
    dnames=info_sst[0]
    d=info_sst[1]
    
    lonvec=lonsel.get().reshape(d[0]*d[1])
    latvec=latsel.get().reshape(d[0]*d[1])
    
    
    #Get points within domain
    if limit is not None :
       indvec,flagvec=in_limits(lonvec, latvec, limit)
    
    flagmat=flagvec.reshape(d[0],d[1])
    rowsum=np.sum(flagmat, 0)
    colsum=np.sum(flagmat, 1)
    yflag=rowsum >= 1
    xflag=colsum >= 1
    xid=np.arange(d[0])
    xid=xid.compress(xflag)
    xcnt=np.int(xid.size)
    yid=np.arange(d[1])
    yid=yid.compress(yflag)
    ycnt=np.int(yid.size)
    

    if xcnt == 0: raise Exception('Error : no longitude within limits')
    if ycnt == 0: raise Exception('Error : no latitude within limits')
    
    #Get start points
    xst=np.int(xid.min())
    yst=np.int(yid.min())
    
    

    #Shrink lon & lat
    lon_var=lonsel.get(start=[xst,yst], count=[xcnt,ycnt])
    lat_var=latsel.get(start=[xst,yst], count=[xcnt,ycnt])
    
    #Load SST image
    ###############
    sst=f.select(param)
    attr=sst.attributes()
    slope=attr['slope']
    intercept=attr['intercept']
    flagValue=attr['bad_value_scaled']
    sst_var=sst.get(start=[xst,yst], count=[xcnt,ycnt]) #Shrink sst image
    
    #Compute mask
    if (param == 'sst') or (param == 'sst4') :
        fg=f.select('qual_'+param)
        fg_var=fg.get(start=[xst,yst], count=[xcnt,ycnt])
    
        if flagLevel is None :
            mask=sst_var == flagValue
        else :
            mask= (sst_var == flagValue) | (fg_var >= flagLevel)
    elif param == 'chlor_a'  :
        fg=f.select('l2_flags')
        fg_var=fg.get(start=[xst,yst], count=[xcnt,ycnt]).flatten()
#        dumvar=[False]*(xcnt*ycnt)#np.zeros((xcnt,ycnt,32),dtype=str).reshape((xcnt*ycnt*32))
        dumfg=[[np.int(b) for b in np.binary_repr(f,32)[::-1]] for f in fg_var] #Rq : bits should be read from end to start
        dumvar=np.sum(np.array(dumfg)[:,flags],1) >= 1
#        for i,f in enumerate(fg_var) :
#            dumvar[i] =  (np.array([b for b in np.binary_repr(f,32)])[flags] == '1').any()
        mask=np.reshape(dumvar,(xcnt,ycnt))
    
    #Check flags
#    plt.bar(np.arange(1,33),np.sum(dumfg,0)[::-1]/np.float64(xcnt*ycnt)); plt.show()
        
    sst_var=np.ma.masked_array(sst_var*slope + intercept, mask=mask,type=float)
    
    
    #Image reprojection to avoid bow tie effect
#    import pyresample as pr
#    swath_def=pr.geometry.SwathDefinition(lons=lon_var, lats=lat_var)
#    lons,lats=np.meshgrid(np.arange(lon_var.min(),lon_var.max(),0.005), np.arange( lat_var.min(),lat_var.max(),0.005))
#    grid_def = pr.geometry.GridDefinition(lons=lons, lats=lats)
#    
#    area_id = 'NWMEd'
#    area_name = 'NWMed'
#    proj_id = 'cyl'
#    proj4_args = '+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +units=m'
#    x_size = 601
#    y_size = 351
#    area_extent = (5., 41., 11., 44.5)
#    area_def = pr.utils.get_area_def(area_id, area_name, proj_id, proj4_args, x_size, y_size, area_extent )
#    res = pr.kd_tree.resample_gauss(swath_def, sst_var.data, grid_def, radius_of_influence=10000.,sigmas=500.,fill_value=None)
    

    
    #Quality control
#    sstfg=f.select('qual_sst')
#    l2fg=f.select('l2_flags')
    

    #Return output dictionary
    return {'lon':lon_var,'lat':lat_var,'sst':sst_var}

def modis_oc(file_name,
             limit=None,
             flagLevel=None,
            param='sst'):
#  ;Example of setting up flag level for non SST images
#  ;--> cf http://oceancolor.gsfc.nasa.gov/DOCS/Ocean_Level-2_Data_Products.pdf
#  flagtable=REPLICATE(0B,32) ;32 bits flag
#  ;flaglevel[[0,1,2,3,4,5,8,9,10,12,14,15,16,17,18,19,21,22,23,25,26,28,29]]=1B ;RQ: remove 1 for table correspondance in documentation (table 3)
#  
#  IF (~exist(sfglevel)) THEN sfglevel=1
#  
#  case (flaglevel) of
#    1: flagTable[[  1,  3,  5,            15,16,      19,         25,26]]=1B
#    2: flagTable[[0,1,  3,4,5,8,9,10,  14,15,16,      19,21,22,23,25,26]]=1B
#    3: flagTable[[0,1,  3,4,5,8,9,10,  14,15,16,      19,21,22,23,25,26]]=1B
#    else: flagTable[[0,1,  3,4,5,8,9,10,  14,15,16,      19,21,22,23,25,26]]=1B
#  endcase
#
#;   flaglevel[[  1,  3,  5,       12,   15,16,      19,         25,26]]=1B ;low value flag
#;  ;flaglevel[[0,1,  3,4,5,8,9,10,12,14,15,16,      19,21,22,23,25,26]]=1B ;operational flag (strong removal)
    #Read MODIS HDF4 data
    f = SD(file_name, SDC.READ)
    fattr=f.attributes()

    #Load coordinates
    lon=f.select('longitude')
    lat=f.select('latitude')

    #Shrink image
    #############
    
    
    #scene dimensions
    
    info_sst=f.datasets()[param]
    dnames=info_sst[0]
    d=info_sst[1]
    
    lonvec=lon.get().reshape(d[0]*d[1])
    latvec=lat.get().reshape(d[0]*d[1])
    
    
    #Get points within domain
    if limit is not None :
       indvec,flagvec=in_limits(lonvec, latvec, limit)
    
    flagmat=flagvec.reshape(d[0],d[1])
    rowsum=np.sum(flagmat, 0)
    colsum=np.sum(flagmat, 1)
    yflag=rowsum >= 1
    xflag=colsum >= 1
    xid=np.arange(d[0])
    xid=xid.compress(xflag)
    xcnt=int(xid.size)
    xst=int(xid.min())
    yid=np.arange(d[1])
    yid=yid.compress(yflag)
    ycnt=int(yid.size)
    yst=int(yid.min())

    #Shrink lon & lat
    lon_var=lon.get(start=[xst,yst], count=[xcnt,ycnt])
    lat_var=lat.get(start=[xst,yst], count=[xcnt,ycnt])
    
    
    #Load SST image
    ###############
    sst=f.select(param)
    attr=sst.attributes()
    slope=attr['slope']
    intercept=attr['intercept']
    flagValue=attr['bad_value_scaled']
    sst_var=sst.get(start=[xst,yst], count=[xcnt,ycnt]) #Shrink sst image
    
    #Compute mask
    fg=f.select('qual_'+param)
    fg_var=fg.get(start=[xst,yst], count=[xcnt,ycnt])
    
    if flagLevel is None :
        mask=sst_var == flagValue
    else :
        mask= (sst_var == flagValue) | (fg_var >= flagLevel)
    
    sst_var=np.ma.masked_array(sst_var*slope + intercept, mask=mask,type=float)
      


#Date conversion tools
######################

def modis_date_convert(argin,
                       julian=True,
                       calendar=False,
                       matlab=False,
                       #nasa=False,
                       #string=False,
                       #calformat=False,
                       epoch=None,
                       verbose=False):
    
    #Scalar to vector conversion
    if np.size(argin) == 1 :
        if type(argin) is not list : argin = [argin]
        
    #Check conversion direction
    #MODIS to CNES (or calendar) -> forward
    #CNES to MODIS               -> backward
    forward = type(argin[0]) is str
        
    
    if calendar is True : julian = False
    if julian is True : calendar = False
    
    if epoch is None : epoch = datetime.date(1950, 01, 01)
    
    if forward :
        if len(argin[0]) == 7 : strptime_fmt='%Y%j'
        if len(argin[0]) == 11 : strptime_fmt='%Y%j%H%M'
        if len(argin[0]) == 13 : strptime_fmt='%Y%j%H%M%S'
        moddate=[datetime.datetime.strptime(x,strptime_fmt) for x in argin]
        if julian is True :
            return [x.toordinal() +                                                                     #Number of days since reference
                    (datetime.timedelta(seconds=x.second,minutes=x.minute,hours=x.hour)).seconds/86400. #Time lag in seconds within day
                    - epoch.toordinal() + 1 for x in moddate]                                           #Change reference to epoch
        if calendar is True :
            return [x.strftime('%d/%m/%Y') for x in moddate]
    else :
        _,dateObj=cnes_convert(argin)
        if type(argin[0]) is float : strftime_fmt="%Y%j%H%M%S"
        else : strftime_fmt="%Y%j"
        
        return [O.strftime(strftime_fmt) for O in dateObj]

def modis_date_from_filename(argin,
                             julian=True,
                             calendar=False,
                             matlab=False,
                             epoch=None):
    
    if np.size(argin) == 1 :
        if type(argin) is not list : argin = [argin]
    
    if calendar is True : julian = False
    if julian is True : calendar = False
    
    if epoch is None : epoch = datetime.date(1950, 01, 01)
    
    return modis_date_convert([x[1:12] for x in argin],julian=julian,calendar=calendar)

def modis_filename2cnes(modis_fname):

    modis_date=modis_filename2modisdate(modis_fname)
    return modis2cnes(modis_date)

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

#Data manipoulation object
##########################

class image_data:
    
    def __init__(self,file_pattern,limit=None,current=0,param='sst'):
        
        if limit is None : self.limit = [-90,-180,90,180]
        else : self.limit=limit
        
        self.param=param
        
        ls=glob.glob(file_pattern)
        dirname=[os.path.dirname(p) for p in ls]
        filename=[os.path.basename(i) for i in ls]
        self.datelist=np.array(modis_date_from_filename(filename,julian=True))
        self.filelist=[dirname[i] + "/" + filename[i] for i in range(len(filename))]
        self.current=current

        self.sst_dev=1.5        
        self.cache_size=9
        
        #Sort file list by date
        sort_order=np.argsort(self.datelist)
        self.datelist=self.datelist[sort_order]
        filelist=[self.filelist[i] for i in sort_order]
        self.filelist=filelist
        
        #Setup cache
        self.cache_index=np.ma.masked_array(np.zeros(self.cache_size),mask=np.repeat(True,self.cache_size))
        self.cached=self.cache_index.compressed().size
        self.cached_lon=[]
        self.cached_lat=[]
        self.cached_sst=[]

    def show_image(self,pmap,logscale=False,vmin=None,vmax=None, colorbar=True,title=None,cbar_label=None,cmap=None):
        
        mnsst=self.sst.mean()
        var_sst=np.std(self.sst)
        
        if vmin is None :
            vmin = mnsst - self.sst_dev * var_sst
            if np.abs(vmin - mnsst) >= 2 : vmin=mnsst - 2
        
        if vmax is None :
            vmax = mnsst + self.sst_dev * var_sst
            if np.abs(vmax - mnsst) >= 2 : vmax=mnsst + 2
            
        if title is None :
            title=''
        
        title+='\nSST chosen->{0:%Y/%m/%d}\nshown->{1:%Y/%m/%d - %H:%M} (offset : {2:d} days)\nfile->{3}'.format(self.chosen_date,self.datetime,-self.offset.days,os.path.basename(self.filename))
        
#        if logscale: pmap.pcolormesh(self.lon,self.lat,np.log10(self.sst),vmin=vmin,vmax=vmax)
        pmap.pcolormesh(self.lon,self.lat,self.sst,vmin=vmin,vmax=vmax,cmap=cmap)
        if colorbar is True : pmap.colorbar(label=cbar_label)
        
        pmap.title(title)

    def set_current(self,date,**kwargs):
        try :
            if len(date) > 0: date = date[0]
        except TypeError : date = date
        
        self.current=nearest(self.datelist, date) #return nearest value index
        self.offset=-datetime.timedelta(days=self.datelist[self.current] - date)
        self.filename=self.filelist[self.current]
        self.date=self.datelist[self.current] # date in days
        days=divmod(self.date,1)
        self.datetime=datetime.datetime.fromordinal(int(days[0])) + datetime.timedelta(datetime.datetime(1950,1,1).toordinal() + days[1] - 1)
        hours=divmod(datetime.timedelta(days=divmod(self.date,1)[1]).seconds,3600)
        minutes=divmod(hours[1],60)
        self.hours=hours[0]
        self.minutes=minutes[0]
        
        self.chosen_date=self.datetime+self.offset
        self.chosen_date_modis=modis_date_convert(date)[0]
        self.chosen_date_cal=cnes_convert(date)[0][0]
        self.load(**kwargs)
        
    def load(self,**kwargs):
        
        cache_id=np.arange(self.cache_size)
        
        #Check if current date is already cached
        if self.cache_index.compress(self.cache_index == self.current).size == 1 :
            cache_id=cache_id.compress(self.cache_index == self.current)

            
        #If no data found in cache, load it
        ###################################
        else :
            
            #Check cache first :
            #Cache is full -> remove first element
            if self.cached == self.cache_size :
                print "clearing oldest cache element"
                self.cached_lon.pop(0)
                self.cached_lat.pop(0)
                self.cached_sst.pop(0)
                
                #update cache index
                ###################
                cached=self.cache_index[1:]
                toappend=np.ma.array(0,mask=True)
                self.cache_index=np.append(cached,toappend)
                self.cache_index.mask=np.append(cached.mask,toappend.mask)
                self.cached=self.cache_index.compressed().size
                
#            #Check if there is a cache space left empty
#            empty=cache_id.compress(self.cache_index.mask == True)
            
            #Load data and fill cache
            print "Loading "+os.path.basename(self.filename)+" - modis"
            modis=modis_sst(self.filename,limit=self.limit,param=self.param,**kwargs)
            self.cached_lon+=[modis['lon'][:]]
            self.cached_lat+=[modis['lat'][:]]
            self.cached_sst+=[modis['sst'][:]]
            
            #Update cache index
#            if self.cached == self.cache_size :
#                print 'problem!'
            self.cache_index.mask[self.cached]=False
            self.cache_index[self.cached]=self.current
            self.cached=self.cache_index.compressed().size

            #Set current position within cache
            cache_id=cache_id.compress((self.cache_index == self.current).compressed())
            
            
        #Update cached variables to be returned
        #######################################
        self.lon=self.cached_lon[cache_id]
        self.lat=self.cached_lat[cache_id]
        self.sst=self.cached_sst[cache_id]
        self.date=self.datelist[cache_id]
        
        print 'Cache size : {0} (shape : {2}) - loading {1}'.format(self.cached,cache_id,np.shape(self.cached_lon))
        
