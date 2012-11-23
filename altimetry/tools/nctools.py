'''
Created on 7 sept. 2012

@author: rdussurg
'''

import datetime
import numpy as np
import scipy.io as io
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

#from grid.regridding import griddata, xy2xy

from netCDF4 import Dataset as ncfile
import glob
import os

import alti_tools as atools
from collections import OrderedDict

class nc :
    
    def __init__(self, file_pattern, limit=[-90., 0., 90., 360.], verbose=0, zero_2pi=False, use_local_dims=False, **kwargs):
        
        #Init system variables
#        if limit is None : limit=[-90.,0.,90.,360.]
        self.zero_2pi = zero_2pi
        self.limit = np.array(atools.recale_limits(limit, zero_2pi=self.zero_2pi))
        self.verbose = verbose
        self.fileid = np.array([])

        self.count = 0
        self.size = 0

        self.use_local_dims=use_local_dims

        #Setup file list
        if isinstance(file_pattern, str) : ls = glob.glob(file_pattern)
        else :
            ls = file_pattern.tolist()
            file_pattern = file_pattern[0]
        
        if len(ls) == 0 :
            self.Error('File pattern not matched : ' + file_pattern)
                
        self.filelist = [os.path.basename(i) for i in ls]
        self.filelist_count = [0] * len(self.filelist)
        enum = list(enumerate(ls))
        enum = zip(*enum)
        
        self.fid_list = np.array(enum[0])
        self.dirname = os.path.dirname(file_pattern)
        
        self.par_list = np.array([])
        self.dim_list = np.array([])
        self._dimensions = {'_ndims':0}
    
    def read(self,**kwargs):  
        #Loop over data files
        #####################
        for i in np.arange(len(self.fid_list)) :
            
            #Read data file
            ###############
            filename = self.dirname+os.sep+self.filelist[i]
            self.message(0, "Loading " + os.path.basename(filename))
            
            res = self.load(filename, **kwargs) #read() function is specific of each class
#            self.update_dataset(res) #update class with loaded data
#            
#            self.check_variables()
        
#        self.update()

        return(res)

    def load(self, filename, params=None, force=False, depthrange=None, timerange=None, **kwargs):
        """
        READ_ARGONC : Argo NetCDF reader
        
        @return: outStr {type:dict} Output data structure containing all recorded parameters as specificied by NetCDF file PARAMETER list.
        @author: Renaud Dussurget
        """
        
        if (params is not None) & isinstance(params,str): params=[params]
        
        #Open file
        self._filename = filename
        ncf = ncfile(self._filename, "r")
     
        #Get list of recorded parameters:
        dum = ncf.variables.keys()
        nparam = np.shape(dum)[0]
    #        par_list=np.array([''.join(ncf.variables.keys()[0,i,:].compressed()) for i in np.arange(nparam)])
        par_list = np.array(['{0}'.format(v) for v in ncf.variables.keys()])
        
        #remove empty items and update nparam
        par_list = par_list.compress([len(par) != 0 for par in par_list])
        nparam = par_list.size
        
        #Get dimensions
        ncdimlist = np.array(['{0}'.format(d) for d in ncf.dimensions.keys()])
        ndims = len(ncdimlist)
        dimStr = OrderedDict()
        dimStr.update({'_ndims':ndims})
    
        
        #Check for the presence of strategic dimensions
        checkedDims = np.array(['lon', 'lat', 'time', 'depth'])
        existDim = -np.ones(4)
        if not self.use_local_dims :
            for i,d in enumerate(ncdimlist) :
                if (d.lower().startswith('lon')) | (d.lower().find('longitude') != -1) : existDim[0]=i
                if (d.lower().startswith('lat')) | (d.lower().find('latitude') != -1) : existDim[1]=i
                if (d.lower().startswith('time')) | (d.lower().startswith('date')) : existDim[2]=i
                if (d.lower().startswith('lev')) | (d.lower().startswith('dep')) : existDim[3]=i
                    
#            existDim[0] = np.where([d.lower().startswith('lon') | (d.lower().find('longitude') != -1) for d in ncdimlist])[0]
#            existDim[1] = np.where([d.lower().startswith('lat') | (d.lower().find('latitude') != -1) for d in ncdimlist])[0]
#            existDim[2] = np.where([(d.lower().startswith('time')) | (d.lower().startswith('date')) for d in ncdimlist])[0]
#            existDim[3] = np.where([(d.lower().startswith('lev')) | (d.lower().startswith('dep')) for d in ncdimlist])[0]
        
        identified = existDim > -1
        
#        checkedDims[identified]=checkedDims[identified][existDim.compress(identified).astype(int)]
#        checkedDims=checkedDims[identified]
        
#        for cn, vn in enumerate(checkedDims) : dimStr.update({cn:len(ncf.dimensions[vn])})
        
        #Update dimension structure with identified dimensions
        #Load dimensional variables
        #TODO : Add scaling here in case...
        for i,d in enumerate(existDim) : 
            if identified[i] :
                dimStr.update({checkedDims[i]:len(ncf.dimensions[ncdimlist[d]])})
                cmd = checkedDims[i] + '=load_ncVar(\'' + ncdimlist[d] + '\',nc=ncf)'
                self.message(4, 'exec : ' + cmd)
                exec(cmd)

        for i,d in enumerate(ncdimlist) :
            if not identified[i] :
                dimStr.update({d:len(ncf.dimensions[d])})
                if ncf.variables.has_key(d) :
                    cmd = d + '=load_ncVar(\'' + d + '\',nc=ncf)'
                    self.message(4, 'exec : ' + cmd)
                    exec(cmd)
                #If the variable associated to the dimension do not exist, generate it
                else :
                    self.message(1, '[WARING] Netcdf file not standard - creating data for {0} dimnsion'.format(d))
                    ndim=len(ncf.dimensions[d])
                    var = {'_dimensions':{'_ndims':1,'record':ndim}, 'data':np.arange(ndim)}
                    cmd = d + '=var'
                    self.message(4, 'exec : ' + cmd)
                    exec(cmd)
        
#        #Update dimension structure with identified dimensions
#        for cn, vn in zip(*(checkedDims[identified], ncdimlist[identified])) : dimStr.update({cn:len(ncf.dimensions[vn])})
#        for vn in ncdimlist[~identified] : dimStr.update({vn:len(ncf.dimensions[vn])})
        
        
        #Load dimensional variables
        #TODO : Add scaling here in case...
#        for cn, vn in zip(*(checkedDims[identified], ncdimlist[identified])) :
##            cmd = 'self.'+cn+'=ncf.variables[\''+vn+'\'][:]'
#            cmd = cn + '=load_ncVar(\'' + vn + '\',nc=ncf)'
#            self.message(4, 'exec : ' + cmd)
#            exec(cmd)
#            
#        for vn in ncdimlist[~identified] :
##            cmd = 'self.'+vn+'=ncf.variables[\''+vn+'\'][:]'
#            cmd = vn + '=load_ncVar(\'' + vn + '\',nc=ncf)'
#            self.message(4, 'exec : ' + cmd)
#            exec(cmd)
            
        #Update dimlist with dimensions present in the object
#        dimlist = np.append(checkedDims[identified], ncdimlist[~identified])
        dimlist=[]
        for d in dimStr.keys() :
            if not d.startswith('_') : dimlist = np.append(dimlist,d)
#        dimlist=[(d if not d.startswith('_') else None) for d in dimStr.keys()]
        
        if params is not None :
            if force : par_list = [i.upper() for i in params]
            else :par_list = list(set(params).intersection(par_list))
        else : par_list = par_list.tolist()
        
        self.message(1, 'Recorded parameters : ' + str(nparam) + ' -> ' + str(par_list))
     
        
        #Extract within limits
        if (existDim[0] > -1) & (existDim[1] > -1):
            llind, flag = atools.in_limits(lon['data'],lat['data'], limit=self.limit)
            lon = atools.recale(lon['data'].compress(flag),degrees=True)
            lat = lat['data'].compress(flag)
            dimStr['lon']=len(lon)
            dimStr['lat']=len(lat)
#            self.message(4, 'self.lon & self.lat updated')
        
        if (existDim[2] > -1):
            if (timerange is not None) : timeflag = (time['data'] >= np.min(timerange)) & (time['data'] <= np.max(timerange))
            else : timeflag = np.ones(len(time['data']), dtype=bool)
            if timeflag.sum() == 0 : self.Error('No data within specified depth range (min/max = {0}/{1})'.format(np.min(time), np.max(time)))
            time = time['data'].compress(timeflag)
            dimStr['time']=len(time)
#            self.message(4, 'self.lon & self.lat updated')
        
        #Extract within depth range
        if (existDim[3] > -1):
            if (depthrange is not None) : depthflag = (depth['data'] >= np.min(depthrange)) & (depth['data'] <= np.max(depthrange))
            else : depthflag = np.ones(len(depth['data']), dtype=bool)
            if depthflag.sum() == 0 : self.Error('No data within specified depth range (min/max = {0}/{1})'.format(np.min(depth), np.max(depth)))
            depth = depth['data'].compress(depthflag)
            dimStr['depth']=len(depth)
        
        outStr = OrderedDict()
        outStr.update({'_dimensions':dimStr})
        if (existDim[0] > -1) : outStr.update({'lon':lon})
        if (existDim[1] > -1) : outStr.update({'lat':lat})
        if (existDim[2] > -1) : outStr.update({'time':time})
        if (existDim[3] > -1) : outStr.update({'depth':depth})
        
        
        #Update object with remaining dimensions
        for d in dimlist.compress([not outStr.has_key(f) for f in dimlist]) :
            cmd = 'outStr.update({\''+d+'\':'+d+'[\'data\']})'
            self.message(4, 'exec : '+cmd)
            exec(cmd)
        
        ncdimStr=outStr.copy()
        #Get dimension lengths
        shape=()
        for d in dimlist:
            exec('shape+=(len('+d+'),)')
        
        ndims = np.size(shape)
     
#        #Create dimension structure
#        curDim = [str(dimname) for dimname in dimStr.keys()[1:]] #[str(dimname) for dimname in ncf.variables['LONGITUDE'].dimensions]
#        curDimval = [dimStr[dim] for dim in curDim] #[len(ncf.dimensions[dimname]) for dimname in curDim]
#   
#        outStr={'_dimensions':{'_ndims':ndims,'nbpoints':sz[0]},'lon':lon,'lat':lat,'date':date}
        
#        for d in dimlist : outStr.update({d:self.__dict__[d]})

        #Sort NCDIMLIST to match DIMLIST
#        ncdimlist[np.arange(len(identified))[np.arange(len(identified))]]=ncdimlist[existDim[identified].tolist()]
        ncdimlist[np.sort(existDim.astype(np.int)[identified])]=ncdimlist[existDim[identified].tolist()]


        #Setup kwargs with current dimensionnal properties
        for d, ncd in zip(*(dimlist,ncdimlist)):
            if not kwargs.has_key(ncd) :
                if kwargs.has_key(d) :
                    
                    kwargs.update({ncd:kwargs[d]})
                    del kwargs[d]
                else :
                    kwargs.update({ncd:(ncdimStr[d].min(),ncdimStr[d].max())})

        
        for param in par_list :
#            dumVar = load_ncVar(param,  nc=ncf, lon=llind[0], lat=llind[1], time=np.arange(len(time['data'])).compress(timeflag),**kwargs) #Load variables
#            dumVar = load_ncVar(param,  nc=ncf, longitude=(self.limit[1],self.limit[3]), latitude=(self.limit[0],self.limit[2]), time=(self.time.min(),self.time.max()),**kwargs) #Load variables
            
            
            dumVar = load_ncVar(param,  nc=ncf, **kwargs) #Load variables
            dimStr = dumVar['_dimensions']
            
            #update dimensions
            curDim = [str(dimname) for dimname in dimStr.keys()[1:]] #[str(dimname) for dimname in ncf.variables['LONGITUDE'].dimensions]
            curDimval = [dimStr[dim] for dim in curDim] #[len(ncf.dimensions[dimname]) for dimname in curDim]
            
            curDim = dimlist[atools.where_list(curDim, ncdimlist.tolist())] #Convert to object dimension names
            
##            curDim = [str(dimname) for dimname in ncf.variables[param].dimensions]
##            curDimval = [len(ncf.dimensions[dimname]) for dimname in curDim]
#            flag = [(np.array(dimname) == outStr['_dimensions'].keys()).sum() == 0 for dimname in curDim] #find dimensions to update
#            dimUpdate = np.array(curDim).compress(flag)
#            for enum in enumerate(dimUpdate) : 
#                self.message(2, 'Appending dimensions {0}:{1} to dataStructure'.format(enum[1], np.array(curDimval).compress(flag)[enum[0]]))
#                outStr['_dimensions'].update({enum[1]:np.array(curDimval).compress(flag)[enum[0]]}) #Append new dimension
#                outStr['_dimensions']['_ndims'] += 1 #update dimension counts
            
            cmd = 'dumStr = {\'' + param + '\':dumVar[\'data\']}'
            self.message(4, 'exec : ' + cmd)
            exec(cmd)
            outStr.update(dumStr)
            
            cmd = 'self.'+param+'='
        
        ncf.close()
        return outStr
        
        
    def message(self, MSG_LEVEL, str):
        """
         MESSAGE : print function wrapper. Print a message depending on the verbose level
         
         @param MSG_LEVEL {in}{required}{type=int} level of the message to be compared with self.verbose
         
         @example self.log(0,'This message will be shown for any verbose level') 
        """
        
        if MSG_LEVEL <= self.verbose : print(str)
        
    def Error(self, ErrorMsg):    
        raise Exception(ErrorMsg)


def load_ncVar(varName, nc=None, **kwargs):
        
        if (nc is None) : raise Exception('No Netcdf file passed')
        
        #Load variable
        var = nc.variables[varName]
        
        var.set_auto_maskandscale(False)
        
        #Load dimensions
        varDim = [str(dim) for dim in var.dimensions]
        varDimval = [len(nc.dimensions[dimname]) for dimname in varDim]
        
        ind_list = [] #Init index list
        dims = {'_ndims':0} #Init dimensions
        
        dstr=[]
        shape=()

        #Construct index list
        #looping on variable dimension list
        for vid,vn in enumerate(varDim) :
            
            #No indexation on current dimension
            if not kwargs.has_key(vn) :
                dstr=np.append(dstr,':')
                
                sz=np.long(varDimval[vid])
#                ind_list.append(range(varDimval[vid])) # if no restriction in kargs then equivalent to [:]
#                dims.update({vn:varDimval[vid]})
            
            #Data is indexed along current dimension
            else :
                drange=kwargs[vn]
                if len(drange) == 2 : drange = drange + (1,)
                if nc.variables.has_key(vn) : #Check if current dimension exists
                    dumvar = nc.variables[vn][:]
                else :
                    dumvar = np.arange(len(nc.dimensions[vn]))
                if vn.startswith('lon') : dumvar=atools.recale(dumvar,degrees=True)
                fg=(dumvar >= drange[0]) & (dumvar <= drange[1])
                if fg.sum() == 0 :
                    #retry switrhcing lon/lat
                    dumvar=atools.recale(dumvar,degrees=True)
                    drange=tuple(atools.recale(drange,degrees=True).astype(np.long))
                    fg=(dumvar >= drange[0]) & (dumvar <= drange[1])
                if fg.sum() == 0 : raise IndexError('{0} {1} is not matching given dimensions {2}'.format(vn,(np.nanmin(nc.variables[vn][:]),np.nanmax(nc.variables[vn][:])),drange))
                if len(fg) == 1 :
                    dstr=np.append(dstr,':')
                    sz=1L
                elif len(fg) == 0:
                    sz=0L
                else :
                    dumind=np.arange(varDimval[vid]).compress(fg)
                    bg=dumind[0]
                    en=dumind[-1]+1
                    st=drange[2]
                    dstr=np.append(dstr,'{0}:{1}:{2}'.format(bg,en,st))
                    sz = np.long(np.mod(np.float(en-bg-1)/st,np.float(en-bg)) + 1.)
            
            dims.update({vn:sz})
            shape = shape + (sz,)
#                if isinstance(dumind, np.ndarray) : dumind = dumind.tolist() #Rq: tolist() can take a very long time to run on large arrays
#                if type(dumind) is not list : dumind = [dumind] 
#                ind_list.append(dumind)
#                dims.update({vn:len(dumind)})
        
#        #check index list
#        sz = [np.size(i) for i in ind_list]
        
        dstr=','.join(dstr)
#        print len(shape)
        
        cmd = 'varOut = var[{0}]'.format(dstr)
        exec(cmd)
        
        #find empty variables
#        if not (atools.where_list([0], shape)[0] == -1) : varOut = var[[0]][[]]
#        else : varOut = var[ind_list]
        
        #Mask it!
        if var.__dict__.has_key('_FillValue') : mask = varOut == var._FillValue
        elif var.__dict__.has_key('missing_value') : mask = varOut == var.missing_value
        else : mask = np.zeros(varOut.shape, dtype='bool')
        
        #Scale it
        #note : we do not use the *= or += operators to force casting to scaling attribute types
        if var.__dict__.has_key('scale') : varOut = varOut * var.scale
        elif var.__dict__.has_key('scale_factor') : varOut = varOut * var.scale_factor
        if var.__dict__.has_key('add_offset') : varOut = varOut + var.add_offset
        
        #Set masks properly
        if isinstance(varOut, np.ndarray) : varOut = np.ma.masked_array(varOut, mask=mask)
        elif isinstance(varOut, np.ma.masked_array) : var.mask = mask
        else : raise 'This data type {} has not been defined - code it!'.format(type(varOut))
            
        #Build up output structure
        outStr = {'_dimensions':dims, 'data':varOut}
        dims.update({'_ndims':len(dims.keys()[1:])})
        
        return outStr
#            ind_list=[[]] 
    
