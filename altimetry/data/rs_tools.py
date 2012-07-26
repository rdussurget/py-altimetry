import sys
import os
import glob

import datetime

import numpy as np
import matplotlib.pyplot as plt
from pyhdf.SD import SD, SDC

import alti_tools as atools


#Read HDF data
##############
def modis_sst(file_name,
              limit=None,
              flagLevel=None,
              param='sst'):

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
       indvec,flagvec=atools.in_limits(lonvec, latvec, limit)
    
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
    if param == 'sst' :
        fg=f.select('qual_'+param)
        fg_var=fg.get(start=[xst,yst], count=[xcnt,ycnt])
    
        if flagLevel is None :
            mask=sst_var == flagValue
        else :
            mask= (sst_var == flagValue) | (fg_var >= flagLevel)
    else : mask = np.zeros((xcnt,ycnt),dtype='bool')
        
    sst_var=np.ma.masked_array(sst_var*slope + intercept, mask=mask,type=float)
      

    
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
       indvec,flagvec=atools.in_limits(lonvec, latvec, limit)
    
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
        _,dateObj=atools.cnes_convert(argin)
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

    def show_image(self,pmap,vmin=None,vmax=None, colorbar=True,title=None):
        
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
        
        pmap.pcolormesh(self.lon,self.lat,self.sst,vmin=vmin,vmax=vmax)
        if colorbar is True : pmap.colorbar()
        
        pmap.title(title)

    def set_current(self,date):
        self.current=atools.nearest(self.datelist, date) #return nearest value index
        self.offset=-datetime.timedelta(days=self.datelist[self.current] - date)
        self.filename=self.filelist[self.current]
        self.date=self.datelist[self.current] # date in days
        days=divmod(self.date,1)
        self.datetime=datetime.datetime.fromordinal(int(days[0])) + datetime.timedelta(datetime.datetime(1950,1,1).toordinal() + days[1])
        hours=divmod(datetime.timedelta(days=divmod(self.date,1)[1]).seconds,3600)
        minutes=divmod(hours[1],60)
        self.hours=hours[0]
        self.minutes=minutes[0]
        
        self.chosen_date=self.datetime+self.offset
        self.chosen_date_modis=modis_date_convert(date)[0]
        self.chosen_date_cal=atools.cnes_convert(date)[0][0]
        self.load()
        
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
            modis=modis_sst(self.filename,limit=self.limit,flagLevel=1,param=self.param)
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
            cache_id=cache_id.compress(self.cache_index == self.current)
            
            
        #Update cached variables to be returned
        #######################################
        self.lon=self.cached_lon[cache_id]
        self.lat=self.cached_lat[cache_id]
        self.sst=self.cached_sst[cache_id]
        self.date=self.datelist[cache_id]
        
        print 'Cache size : {0} (shape : {2}) - loading {1}'.format(self.cached,cache_id,np.shape(self.cached_lon))
        
