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

import numpy as np
from netCDF4 import Dataset as ncfile
import operator
import glob
import os


import hydro as htools #This object is based on the hydro_data object
from altimetry.externals import esutils_stat as es
from altimetry.tools import cnes_convert, histogram_indices, recale, in_limits, cumulative_distance , nctools
from collections import OrderedDict

if __debug__ : import matplotlib.pyplot as plt


#Load alti data
###############
class alti_data(htools.hydro_data) :
    def __init__(self,file_pattern,time_range=None,output_is_dict=True,**kwargs):
        
        if time_range is not None :
            ls=np.array(glob.glob(file_pattern))
            if len(ls) == 0: raise Exception('File not found : {0}'.format(file_pattern))
            filelist=[os.path.basename(p) for p in ls]
            st=[f.split('_')[-3] for f in filelist]
            
#            print st

            jst = np.array([cnes_convert('{0}/{1}/{2}'.format(s[-2:],s[-4:-2],s[0:4]))[0][0] for s in st])            
            
            #Time sort
            order=np.argsort(jst)
            jst=jst[order]
            ls=ls[order]
            
            dft=jst[1:] - jst[:-1]
            dt=np.fix(np.median(dft[dft != 0])) #If dft == 0 : then duplicates!

            hist,R= es.histogram(jst,binsize=dt,use_weave=False,rev=True,min=time_range[0] - dt/2.,max=time_range[1] + dt/2.)
            dumind = histogram_indices(hist, R)
            ind=np.array([])
            for i in dumind : ind=np.append(ind,i)
            ind=ind.tolist()
            file_pattern=ls[ind]
            
        htools.hydro_data.__init__(self,file_pattern,output_is_dict=output_is_dict,**kwargs)
        self.set_sats()
    
    def set_sats(self):
        self.sat=[]
        for enum in enumerate(self.filelist):
            self.sat=np.append(self.sat,[enum[1].split('_')[2]]*self.filelist_count[enum[0]])
        self.sat=np.ma.array(self.sat,mask=False) 
    
    def read(self,filename,datatype=None,slaext=False,**kwargs):
        
        fname,extension = os.path.splitext(filename)
        
        if os.path.basename(filename).count('.') > os.path.basename(filename).count('_'): delim='.'
        else : delim = '_'
        
        #Get data type
        if datatype is None :
            if os.path.basename(filename).split(delim)[0] == 'ctoh' : datatype='CTOH'
            if os.path.basename(filename).split(delim)[0] == 'PISTACH' : datatype='PISTACH'
            if os.path.basename(filename).split(delim)[0] == 'nrt' : datatype='NRT'
            if os.path.basename(filename).split(delim)[0] == 'dt' : datatype='DT'
#        else :
#            datatype='RAW' #Setup default as raw NetCDF file
        
        if (datatype == 'DT') | (datatype == 'NRT') | (datatype == 'PISTACH') :
            if slaext : outStr=self.read_slaext(filename,datatype=datatype,**kwargs)
            else : outStr=self.read_sla(filename,datatype=datatype,**kwargs)
            self.update_fid_list(os.path.basename(filename),outStr['_dimensions']['time'])
        elif (datatype == 'CTOH') :
            outStr=self.read_CTOH(filename,**kwargs)
            self.update_fid_list(os.path.basename(filename),outStr['_dimensions']['time'])
        else: #Setup default as raw NetCDF file
            outStr=self.read_nc(filename,**kwargs)
            self.update_fid_list(os.path.basename(filename),outStr['_dimensions'][outStr['_dimensions'].keys()[1]])
        
        
        return outStr
    
    def read_sla(self,filename,params=None,force=False,timerange=None,datatype=None,**kwargs):
        
        """
        READ_SLA : Read AVISO Along-Track products
        
        @return: outStr {type:dict} Output data structure containing all recorded parameters as specificied by NetCDF file PARAMETER list.
        @author: Renaud Dussurget
        """
        
        from time import time
        from datetime import timedelta
#        a = time()
        
        self.message(2,'Reading AVISO DT data ({0})'.format(datatype))
        
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
        for i in ['BeginDates','Longitudes','Latitudes'] : par_list.pop(par_list.index(i))
        nparam=len(par_list)
        
        self.message(2,'Recorded parameters : '+str(nparam)+' -> '+str(par_list))
     
        lon = self.load_ncVar('Longitudes',**kwargs)
        lon['data'] = recale(lon['data'], degrees=True, zero_2pi=True) #shift longitudes
        lat = self.load_ncVar('Latitudes',**kwargs)
        
#        lon = self.load_ncVar('Longitudes',Longitudes=(self.limit[1],self.limit[3]),**kwargs)
#        lon['data'] = recale(lon['data'], degrees=True, zero_2pi=True) #shift longitudes
#        lat = self.load_ncVar('Latitudes',Latitudes=(self.limit[0],self.limit[1]),**kwargs)

        
        #Extract within limits
        ind, flag = in_limits(lon['data'],lat['data'],limit=self.limit)
        dim_lon = lon['_dimensions']
        lat = lat['data'].compress(flag)
        lon = lon['data'].compress(flag)
        dist=cumulative_distance(lat, lon)
        
        sz=np.shape(lon)
        ndims=np.size(sz)
        
        #Get dates
        stDate = self.load_ncVar('BeginDates',**kwargs)['data']
        dumVar = self.load_ncVar('Cycles',**kwargs)
        nbCyc = dumVar['data']
        Ncycs = dumVar['_dimensions']['Cycles']
        Ntra = dumVar['_dimensions']['Tracks']
        nbTra = self.load_ncVar('Tracks',**kwargs)['data']
        
#        if np.size(stDate) == 1 : stDate.
        
        DeltaT = self._ncfile.variables['DeltaT'][0]  / 86400. #* self._ncfile.variables['DeltaT'].scale_factor
        npts = self.load_ncVar('NbPoints',**kwargs)['data']
        dumind=np.cumsum(npts)
        
        #Loop 1
#        date=np.ma.array([],mask=[])
#        cycles=np.ma.array([],mask=[])
#        tracks=np.ma.array([],mask=[])
#        for i in xrange(Ncycs) :
#            np.ma.concatenate((nbTra,nbTra))
#        
#        for i,nc in enumerate(nbCyc.data.flatten()):
#            N=npts[i]
#            curInd=np.array(list(set(xrange(dumind[i]-N,dumind[i]) if N > 0 else []).intersection(ind)))
#            ncur=len(curInd)
#            date=np.ma.concatenate((date,(curInd - dumind[0])*DeltaT+stDate.flatten()[i]))
#            cycles=np.ma.concatenate((cycles,np.ma.array((nbCyc.data.flatten()[i],)*ncur)))#,mask=np.repeat(nbCyc.mask[i][j],ncur))))
#            tracks=np.ma.concatenate((tracks,np.ma.array((nbTra.data.flatten()[i],)*ncur)))
        
        #Loop 2
        
        import operator
        
        date = ()
        cycles = ()
        tracks = ()
        
#        rowind = (0,)*Ntra
        
#        nind=0
#        for i in xrange(Ncycs): nind+=(npts*(~nbCyc.mask.T[i])).sum()
        
        indcopy=ind.copy()
#        npts_copy=npts.copy()
        npts[npts.mask]=0
        dumind[dumind.mask]=0
        
        nbTra_copy=nbTra.copy()
        
        for i in np.arange(1,Ncycs,1.0,dtype=int) :
            nbTra=np.ma.concatenate((nbTra,nbTra_copy))
            npts=np.ma.concatenate((npts,tuple((~nbCyc.T[i].mask)*1*npts)))
#            rowind+=(i,)*Ntra      
        
        npts=npts.reshape(nbCyc.shape[::-1]).T
        nbTra=nbTra.reshape(nbCyc.shape[::-1]).T
#        rowind=np.reshape(rowind,nbCyc.shape[::-1]).T
        
#        npts.mask=nbCyc.mask
        nbTra.mask=nbCyc.mask
        
        npts=npts.flatten()
        nbTra=nbTra.flatten()
#        rowind=rowind.flatten()
        
        nbCyc_flatten=nbCyc.flatten()
        nbTra_flatten=nbTra.flatten()
        stDate_flatten=stDate.flatten()
        
#        nind=0
        
        outInd=[]
        
        for i,nc in enumerate(nbCyc.data.flatten()):
            N=npts[i]
            Nprev=npts[i-Ncycs] if i >= (Ncycs) and np.remainder(float(i),Ncycs) == 0 else 0
            indcopy-=Nprev #if rowind[i] == 0 else 0
            curInd=tuple(sorted(set(xrange(N) if N > 0 else []).intersection(indcopy)))
            ncur=len(curInd)
#            nind+=ncur
            outInd+=map(operator.sub, curInd,(( (curInd[0] if len(curInd) > 0 else 0) - (outInd[-1] +1 if len(outInd) > 0 else 0) - len(ind)*(np.remainder(float(i),Ncycs)),)*ncur))
            curInd=tuple(map(operator.mul, curInd, (DeltaT,)*ncur))     
            date+=tuple(map(operator.add, curInd, (stDate_flatten[i],)*ncur))
            cycles+=(nbCyc_flatten[i],)*ncur
            tracks+=(nbTra_flatten[i],)*ncur
        
        date=np.ma.masked_array(date,mask=False)
        cycles=np.ma.masked_array(cycles,mask=False)
        tracks=np.ma.masked_array(tracks,mask=False)
                
        #Loop 3
#        date=np.ma.array([],mask=[])
#        cycles=np.ma.array([],mask=[])
#        tracks=np.ma.array([],mask=[])
#        for j in xrange(Ncycs) :
#            for i,N in enumerate(npts.data) :
##                curFg=(ind >= dumind[i]-N) & (ind <= dumind[i])
#                curInd=np.array(list(set(xrange(dumind[i]-N,dumind[i]) if N > 0 else []).intersection(ind)))
#                ncur=len(curInd)
#                date=np.ma.concatenate((date,(curInd - dumind[0])*DeltaT+stDate[i][j]))
#                cycles=np.ma.concatenate((cycles,np.ma.array((nbCyc.data[i][j],)*ncur)))#,mask=np.repeat(nbCyc.mask[i][j],ncur))))
#                tracks=np.ma.concatenate((tracks,np.ma.array((nbTra.data[i],)*ncur)))#,mask=np.repeat(nbCyc.mask[i][j],ncur))))
        
        outInd=np.array(outInd,dtype=int)
        
        #Check output index
#         print outInd.shape[0],npts.cumsum().max()
        
        nt=len(date)
        date.mask=(False,)*nt
        cycles.mask=date.mask
        tracks.mask=date.mask
#        date=date.reshape((Ncycs,)+(npts.sum(),)).T
#        mask=date.mask
#        date=date.compressed()
#        cycles=cycles.reshape((Ncycs,)+(npts.sum(),)).T.compressed()
#        tracks=tracks.reshape((Ncycs,)+(npts.sum(),)).T.compressed()
        
        
#        lon=np.repeat(lon,Ncycs)
#        lat=np.repeat(lat,Ncycs)
#        mask=~lon.mask
        
        dimStr = dim_lon
        dimStr.pop('Data')
        nrec=len(date)
        dimStr.update({'time':nrec})
        
        for i in ['DeltaT','NbPoints','Cycles','Tracks','DataIndexes'] : par_list.pop(par_list.index(i))
        
        outStr={'_dimensions':dimStr,'lon':lon,'lat':lat,'date':date,'cycle':cycles,'track':tracks}
        
        
        
        for param in par_list :
            a = time()
            dumVar = self.load_ncVar(param,Data=ind,**kwargs) #Load variables            
            runtime = time() - a
#            print 'runtime:', timedelta(seconds=runtime)
            
            dimStr=dumVar['_dimensions']
            dimStr.pop('Cycles')
            dimStr.pop('Data')
            dimStr['time']=nrec
            dimStr['_ndims']=len(dimStr.keys())-1
            
            #update dimensions
            curDim = [str(dimname) for dimname in dimStr.keys()[1:]] #[str(dimname) for dimname in self._ncfile.variables['LONGITUDE'].dimensions]
            curDimval = [dimStr[dim] for dim in curDim] #[len(self._ncfile.dimensions[dimname]) for dimname in curDim]
            flag = [(np.array(dimname) == outStr['_dimensions'].keys()).sum() == 0 for dimname in curDim] #find dimensions to update
            dimUpdate = np.array(curDim).compress(flag)
            for enum in enumerate(dimUpdate) : 
                self.message(3, 'Appending dimensions {0}:{1} to dataStructure'.format(enum[1],np.array(curDimval).compress(flag)[enum[0]]))
                outStr['_dimensions'].update({enum[1]:np.array(curDimval).compress(flag)[enum[0]]}) #Append new dimension
                outStr['_dimensions']['_ndims']+=1 #update dimension counts
            
#            dumStr = {param.lower() : dumVar['data']}
            dumStr = {param.lower() : dumVar['data'].flatten()[outInd]}
#            cmd = 'dumStr = {\''+param.lower()+'\':dumVar[\'data\']}'
#            self.message(4, 'exec : '+cmd)
#            exec(cmd)
            outStr.update(dumStr)
        
        id=np.repeat(sat_name,outStr['_dimensions']['time'])
        
        
        outStr.update({'id':id})
        self._ncfile.close()
        
        #Checkit [len(outStr[k]) for k in outStr.keys()]
        return outStr
    
    def read_slaext(self,filename,params=None,force=False,timerange=None,datatype=None,**kwargs):
        
        """
        READ_SLAEXT : Read AVISO Along-Track SLA regional products
        
        @return: outStr {type:dict} Output data structure containing all recorded parameters as specificied by NetCDF file PARAMETER list.
        @author: Renaud Dussurget
        """
        
        self.message(2,'Reading SLAext data ({0})'.format(datatype))
        
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

    def read_CTOH(self,filename,params=None,force=False,timerange=None,datatype=None,**kwargs):
        
        """
        READ_SLAEXT : Read AVISO Along-Track SLA regional products
        
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
     
        lon = self.load_ncVar('lon',**kwargs)
        lon['data'] = recale(lon['data'], degrees=True, zero_2pi=True) #shift longitudes
        lat = self.load_ncVar('lat',**kwargs)

        
        #Extract within limits
        ind, flag = in_limits(lon['data'],lat['data'],limit=self.limit)
        dim_lon = lon['_dimensions']
        lat = lat['data'].compress(flag)
        lon = lon['data'].compress(flag)
        dist=cumulative_distance(lat, lon)
        
        sz=np.shape(lon)
        ndims=np.size(sz)
        
        #Get date table and convert it to dimension
        date = self.load_ncVar('time',nbpoints=ind,**kwargs)
        
#        mn = np.min(date['data'],1)
#        mx = np.max(date['data'],1)
#        trange=mx - mn
#        mxid = np.argmax(trange)
#        mxtrange = trange[mxid] #Get maximum time range
#        deltat = (date['data'][mxid]))[0,1:data[mxid].mes-1] - (*(data[mxid].data))[0,0:data[mxid].mes-2]
        
#        ;Get theoretical time range
#          FOR i=0, N_ELEMENTS(data) - 1 DO mn = (i EQ 0)? (*(data[i].data))[0,0] : [mn,(*(data[i].data))[0,0]]; Get minimum time at each time series
#          FOR i=0, N_ELEMENTS(data) - 1 DO mx = (i EQ 0)? (*(data[i].data))[0,data[i].mes - 1] : [mx,(*(data[i].data))[0,data[i].mes - 1]] ;same for max time
#          trange=[mx - mn] ;Get time range for each time series
#          mxtrange=max(trange,mxid) ;Get maximum time range
#          deltat=(*(data[mxid].data))[0,1:data[mxid].mes-1] - (*(data[mxid].data))[0,0:data[mxid].mes-2]
#          hist=HISTOGRAM(deltat,LOCATIONS=loc,BINSIZE=1);Put deltat into 1 day boxes
#          dum=MAX(hist,histmx) ;Get modal class
#          tstep=MAX(deltat[WHERE(deltat LT loc[histmx] + 1.AND deltat GT loc[histmx])]) ;delta is the time spelled btw.
#                                                                                        ;max time from previous pass and
#                                                                                        ;min time of following pass.
#                                                                                        ;Thus the greatest time is the closest
#                                                                                        ;to the complete cycle
#          
#          ;tstep=MIN(DERIV((*(data[mxid].data))[0,*]))
#        ;  nt=ROUND(( (*(data[mxid].data))[0,data[mxid].mes - 1] - (*(data[mxid].data))[0,0] ) / tstep) + 1
#          nt=LONG(ROUND(( (*(data[mxid].data))[0,data[mxid].mes - 1] - (*(data[mxid].data))[0,0] ) / tstep) + 1)
#        ;  theoretical_time=Scale_Vector(DINDGEN(nt), MEDIAN(mn),MEDIAN(mx))
#          
#          theoretical_time=DINDGEN(nt)*tstep + MIN(mn) ;This is the theoretical time from reference date MIN(min)
        
        
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
    
    def read_nc(self,filename,**kwargs):
        #Set filename
        self._filename = filename
        
        #Read data from NetCDF
        obj=nctools.nc(verbose=self.verbose,limit=self.limit,use_local_dims=True)
        outStr=obj.read(filename,**kwargs)
        
        #Remove attributes already existing in data object
        for a in self.__dict__.keys():
            if outStr.has_key(a) and not a.startswith('_') :
                outStr.pop(a)
                self.message(4, 'Attribute {0} already exists'.format(a))
        
        return outStr


    def track_list(self,*args):
        noargs = len(args) == 0
        return np.unique(self.track) if noargs else np.unique(self.track.compress(args[0]))
    
    def cycle_list(self,*args):
        noargs = len(args) == 0
        return np.unique(self.cycle) if noargs else np.unique(self.cycle.compress(args[0]))

    def reorder(self,*args,**kwargs):
        #update time range with real value
        trange_str = self.time_range()[0]
        cycle_list = self.cycle_list()
        track_list=self.track_list()
        
        #Solve precision problem (for unique)!
        precision=np.finfo(self.lon.dtype).precision-2 #9
        self.lon=self.lon.round(decimals=precision)
        self.lat=self.lat.round(decimals=precision)
            
        #Detect recurrent  lon/lat/track triplets
        self.message(2, 'Detect recurrent lon/lat/track triplets')
        triplets = np.unique(zip(*(self.lon, self.lat, self.track)))
        
        
        #Sort by track (we do not sort using other columns to preserve descending/ascending order)
        triplets=np.ma.array(sorted(triplets, key=operator.itemgetter(2)))
        
        #Get main variables
        lon = triplets[:,0]
        lat = triplets[:,1]
        tracknb = triplets[:,2]
        
        #!!!!
        #PROBLEM WITH UNIQUE!!!! DUPLICATES ARE STILL PRESENT! 
        
        track=track_list
        cycle=cycle_list
        
        N=len(lon)
        ntracks=len(track_list)
        ncycles=len(cycle_list)
        
        ind = np.arange(N)
        
    #    dst = atools.calcul_distance(lat, lon)
        
        #Get local index on nominal tracks (THIS IS long)
        self.message(2, 'Computing space and time indices')
        xid = np.array([np.where((ln == lon) & (self.lat[i] == lat))[0][0] for i,ln in enumerate(self.lon) ])
        tid = np.array([np.where(c == cycle_list)[0][0] for c in self.cycle])
        
        #PROBLE
        
        #Get object attributes to reform
        varSize = np.array([np.size(self.__dict__[k]) for k in self.__dict__.keys()])
        par_list=np.array(self.__dict__.keys())[varSize == self._dimensions['time']].tolist()
        
        #Refine param list
        for par in ['lon','lat','track','cycle'] :
            par_list.pop(par_list.index(par))
            self.__setattr__(par,locals()[par])
            
        #Set new dimensions array
        self._dimensions=OrderedDict({'_ndims':3,'cycle':ncycles,'record':N,'track':ntracks})
        record=np.ma.array(ind,mask=np.zeros(N,dtype=bool))
        
        #Update lon, lat, track and cycle arrays
        lon.mask=np.zeros(N,dtype=bool)
        lat.mask=np.zeros(N,dtype=bool)
        tracknb.mask=np.zeros(N,dtype=bool)
        track.mask=np.zeros(ntracks,dtype=bool)
        cycle.mask=np.zeros(ncycles,dtype=bool)
        
        lon.__setattr__('_dimensions',{'_ndims':1,'record':N}); lon.__setattr__('long_name','longitude')
        lat.__setattr__('_dimensions',{'_ndims':1,'record':N}); lat.__setattr__('long_name','latitude')
        tracknb.__setattr__('_dimensions',{'_ndims':1,'record':ntracks}); track.__setattr__('long_name','track_number')
        track.__setattr__('_dimensions',{'_ndims':1,'track':ntracks}); track.__setattr__('long_name','track_list')
        cycle.__setattr__('_dimensions',{'_ndims':1,'cycle':ncycles}); cycle.__setattr__('long_name','cycle_number')
        record.__setattr__('_dimensions',{'_ndims':1,'record':N}); record.__setattr__('long_name','record_index')
        
        self.lon=lon
        self.lat=lat
        self.tracknb=tracknb
        self.track=track
        self.cycle=cycle
        self.record=record
        
        #Init output matrices using object fields
        for par in par_list :
            self.message(2, 'Reforming {0}'.format(par))
            locals()[par]=np.ma.array(np.zeros((ncycles,N)),mask=np.ones((ncycles,N),dtype=bool),dtype=self.__getattribute__(par).dtype)
            locals()[par][tid,xid]=self.__getattribute__(par)
            locals()[par].__setattr__('_dimensions',{'_ndims':2,'cycle':ncycles,'record':N})
            if hasattr(self.__dict__[par],'__dict__') :
                #Remove attributes already in output variable
                attrStr=OrderedDict()
                null=[attrStr.update({a:self.__dict__[par].__dict__[a]}) for a in self.__dict__[par].__dict__.keys() if not locals()[par].__dict__.has_key(a)]
                locals()[par].__dict__.update(attrStr)
            self.__setattr__(par,locals()[par])
        
        par_list=np.append(par_list,['lon','lat','tracknb','track','cycle','record'])
        self.par_list=par_list #Add record dimension
    
    def pass_time(self):
        date=self.date
        nt=self._dimensions['cycle']
        N=self._dimensions['record']
        for t in np.arange(nt):
            poly=np.polyfit(np.arange(N)[~date.mask[t,:]], date[t,:][~date.mask[t,:]], 1)
            date[t,:][date.mask[t,:]]=poly[0]*np.arange(N)[date.mask[t,:]] + poly[1]
        date.mask=False
        return date[:,N/2]
    
#    def update_with_slice(self,flag):
#        
#        N=flag.sum()
#        
#        #Get object attributes to update
#        varSize = np.array([np.size(self.__dict__[k]) for k in self.__dict__.keys()])
#        par_list=np.array(self.__dict__.keys())[varSize == self._dimensions['time']].tolist()
#        
#        self._dimensions['time']=N
#        
#        for par in par_list :
#            self.__setattr__(par,self.__dict__[par][flag])
#            if (hasattr(self.__dict__[par], '_dimensions')) : self.__dict__[par]._dimensions['time']=self._dimensions['time']
        
    def ncstruct(self):
        par_list = self.par_list.tolist()
        dimStr=self._dimensions
        dimlist = dimStr.keys()
        outStr = OrderedDict({'_dimensions':dimStr})
        
        if (np.array(par_list) == 'sat').sum() : par_list.pop(par_list.index('sat')) #Remove satellite info
        varlist=np.append(np.array(dimlist[1:]),np.array(par_list))
        
        for d in varlist :
            self.message(2, 'Updating output structure with {0}'.format(d))
            curDim=getattr(self.__getattribute__(d),'_dimensions',None)
            attributes = [a for a in self.__getattribute__(d).__dict__.keys() if not a.startswith('_')]
#            attributes = np.append(attributes,'_dimensions')
            outStr[d]={'_dimensions':self.__getattribute__(d)._dimensions}
            outStr[d].update({'data':self.__getattribute__(d)})
            for a in attributes : outStr[d].update({a:self.__getattribute__(d).__getattribute__(a)})
        
        return outStr
    
    def write_nc(self,filename,clobber=False,**kwargs):
        obj=nctools.nc(verbose=self.verbose,limit=self.limit,use_local_dims=True)
        ncalti=self.ncstruct() #Get an netcdf structure from data
        res=obj.write(ncalti,filename,clobber=clobber,**kwargs) #Save processed datase
    
    def push_nc(self,*args,**kwargs):
        obj=nctools.nc(verbose=self.verbose,limit=self.limit,use_local_dims=True)
        res=obj.push(*args,**kwargs)
        









    




