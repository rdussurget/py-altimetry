'''
NCTOOLS
@summary: Netcdf data object, to help loading and writing data
@change Created on 7 sept. 2012
@author: rdussurg
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

#from grid.regridding import griddata, xy2xy

from netCDF4 import Dataset as ncfile
import glob
import os
from altimetry.tools import recale, in_limits, where_list, recale_limits, get_caller
#import altimetry.data.alti_tools as atools
from collections import OrderedDict
from warnings import warn

class nc :
    
    def __init__(self, limit=[-90., 0., 90., 360.], verbose=0, zero_2pi=False, use_local_dims=False, **kwargs):
        
        #Init system variables
#        if limit is None : limit=[-90.,0.,90.,360.]
        self.zero_2pi = zero_2pi
        self.limit = np.array(recale_limits(limit, zero_2pi=self.zero_2pi))
        self.verbose = verbose
        self.fileid = np.array([])

        self.count = 0
        self.size = 0

        self.use_local_dims=use_local_dims
    
    def write(self, data, outfile, **kwargs):
        
        #Open file
        if self.verbose == 1 : self.message(1, 'Writing data file {}'.format(os.path.basename(outfile)))
        elif self.verbose > 1 : self.message(2, 'Writing data file {}'.format(outfile))
        root_grp=ncfile(outfile, 'w', format='NETCDF4', clobber=True)
#        root_grp.description = 'nctools.write() file'
        
        #Get attributes
        if data.has_key('_attributes'):
            self.message(2, 'Adding attributes data')
            attrStr=data.pop('_attributes')
            self.message(4, 'Attribute list (N={0}) :[{1}]'.format(len(attrStr.keys()), ','.join(attrStr.keys())) )
            root_grp.setncatts(attrStr)
        
        #Get dimensions
        dimStr=data.pop('_dimensions')
        ndims=dimStr.pop('_ndims')
        dimlist = dimStr.keys()
        dimVal = dimStr.values()
        
        #Get variables
        parlist = data.keys()
        
        
        # Set up dimensions
        
        #Put dimensions
        for d in dimlist :
            self.message(2, 'Adding D {0}={1}'.format(d,dimStr[d]))
            root_grp.createDimension(d, dimStr[d])
        
        #Loop over variables
        for p in parlist :
            
            #Get dimensions for current variable
            if not data[p].has_key('_dimensions') : self.Error('_dimension attribute is not set for variaple'+p)
            pardim=data[p].pop('_dimensions')
            if isinstance(pardim,dict) :pardim=tuple(pardim.keys()[1:]) if pardim.has_key("_ndims") else tuple(pardim.keys())
            elif isinstance(pardim,list) : pardim = tuple(pardim)
            elif isinstance(pardim,tuple) : pass
            else : self.Error('_dimensions must be dict, list or tuple - not {0}'.type(pardim)) 
            
            #Convert to numpy array if scalar or non numpy
            if not hasattr(data[p]['data'],'__iter__') or not hasattr(data[p]['data'], 'dtype'):  data[p]['data']=np.array(data[p]['data'])

            if not (data[p]['data'].dtype == '|S6') and not (data[p]['data'].dtype == '|S2') :
                self.message(2, 'Adding V {0} (dims={{{1}}},attr={{{2}}})'.
                             format(p,
                                    ', '.join(['\'{0}\':{1}'.format(d,dimStr[d]) for d in pardim]),
                                    ', '.join(['\'{0}\':{1}'.format(d,data[p][d]) for d in data[p].keys() if (d != '_dimensions') and (d != 'data')]) )
                             )
                if hasattr(data[p]['data'],'fill_value') :
                    locals()[p] = root_grp.createVariable(p,
                                                   data[p]['data'].dtype,
                                                   pardim,
                                                   fill_value=data[p]['data'].fill_value)
                else :
                    locals()[p] = root_grp.createVariable(p,
                                                   data[p]['data'].dtype,
                                                   pardim)
                locals()[p][:]=data[p].pop('data').transpose(tuple(range(len(pardim))[::-1])) #Transpose data before writing it into file

                
                #Update with attribute list
                locals()[p].setncatts(data[p])
#                for a in data[p].keys() :
#                    if not hasattr(locals()[p],a) :
#                        locals()[p].setncattr(a,data[p][a])

        
        self.message(2, 'Closing file')
        root_grp.close()

        return True
        
        
    
    def read(self,file_pattern,**kwargs): 
        
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
        self._dimensions = OrderedDict({'_ndims':0})
         
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
                    var = {'_dimensions':{'_ndims':1,d:ndim}, 'data':np.arange(ndim)}
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
            llind, flag = in_limits(lon['data'],lat['data'], limit=self.limit)
            lon = recale(lon['data'].compress(flag),degrees=True)
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
        
        #Create output data structure
        outStr = OrderedDict()
        outStr.update({'_dimensions':dimStr})
        
        if (existDim[0] > -1) : outStr.update({'lon':lon})
        if (existDim[1] > -1) : outStr.update({'lat':lat})
        if (existDim[2] > -1) : outStr.update({'time':time})
        if (existDim[3] > -1) : outStr.update({'depth':depth})
        
        
        #Update object with remaining variables
        for d in dimlist.compress([not outStr.has_key(f) for f in dimlist]) :
#            cmd = 'outStr.update({\''+d+'\':'+d+'[\'data\']})'
            cmd = 'outStr.update({\''+d+'\':'+d+'})'
            self.message(4, 'exec : '+cmd)
            exec(cmd)
        
        ncdimStr=outStr.copy()
        #Get dimension lengths
        shape=()
        for d in dimlist:
            exec('shape+=(np.shape('+d+'[\'data\']),)')
        
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
                    kwargs.update({ncd:(ncdimStr[d]['data'].min(),ncdimStr[d]['data'].max())})
#            else :
#                outStr['NbLatitudes']['data']        
        for param in par_list :
#            dumVar = load_ncVar(param,  nc=ncf, lon=llind[0], lat=llind[1], time=np.arange(len(time['data'])).compress(timeflag),**kwargs) #Load variables
#            dumVar = load_ncVar(param,  nc=ncf, longitude=(self.limit[1],self.limit[3]), latitude=(self.limit[0],self.limit[2]), time=(self.time.min(),self.time.max()),**kwargs) #Load variables
            
            
            dumVar = load_ncVar(param,  nc=ncf, **kwargs) #Load variables
            dimStr = dumVar['_dimensions']
            
            #update dimensions
            curDim = [str(dimname) for dimname in dimStr.keys()[1:]] #[str(dimname) for dimname in ncf.variables['LONGITUDE'].dimensions]
            curDimval = [dimStr[dim] for dim in curDim] #[len(ncf.dimensions[dimname]) for dimname in curDim]
            
            curDim = dimlist[where_list(curDim, ncdimlist.tolist())] #Convert to object dimension names
            
##            curDim = [str(dimname) for dimname in ncf.variables[param].dimensions]
##            curDimval = [len(ncf.dimensions[dimname]) for dimname in curDim]
#            flag = [(np.array(dimname) == outStr['_dimensions'].keys()).sum() == 0 for dimname in curDim] #find dimensions to update
#            dimUpdate = np.array(curDim).compress(flag)
#            for enum in enumerate(dimUpdate) : 
#                self.message(2, 'Appending dimensions {0}:{1} to dataStructure'.format(enum[1], np.array(curDimval).compress(flag)[enum[0]]))
#                outStr['_dimensions'].update({enum[1]:np.array(curDimval).compress(flag)[enum[0]]}) #Append new dimension
#                outStr['_dimensions']['_ndims'] += 1 #update dimension counts
            
#            cmd = 'dumStr = {\'' + param + '\':dumVar[\'data\']}'
            cmd = 'dumStr = {\'' + param + '\':dumVar}'
            self.message(4, 'exec : ' + cmd)
            exec(cmd)
            outStr.update(dumStr)
            
            #Update output dimensions with extracted dimensions
            for ddum in dumStr[param]['_dimensions'].keys()[1:] :
                if outStr['_dimensions'].get(ddum) != dumStr[param]['_dimensions'][ddum] : outStr['_dimensions'][ddum]=dumStr[param]['_dimensions'][ddum]
            
            cmd = 'self.'+param+'='
        
        ncf.close()
        return outStr
    
    def message(self, MSG_LEVEL, str):
        """
         MESSAGE : print function wrapper. Print a message depending on the verbose level
         
         @param MSG_LEVEL {in}{required}{type=int} level of the message to be compared with self.verbose
         
         @example self.log(0,'This message will be shown for any verbose level')
         @author: Renaud DUSSURGET (RD), LER PAC/IFREMER
         @change: Added a case for variables with missing dimensions
         
        """
        caller=get_caller()
        if MSG_LEVEL <= self.verbose : print('[{0}.{1}()] {2}'.format(__name__,caller.co_name,str))
        
    def Error(self, ErrorMsg):    
        raise Exception(ErrorMsg)
    
    def attributes(self, filename, **kwargs):
        """
        ATTRIBUTRES: Get attributes of a NetCDF file
        
        @return: outStr {type:dict} Attribute structure.
        @author: Renaud Dussurget
        """
        
        #Open file
        self._filename = filename
        ncf = ncfile(self._filename, "r")
     
        #Get list of recorded parameters:
        keys = ncf.__dict__.keys()
        outStr = OrderedDict()
        for a in keys: outStr.update({a:ncf.__getattr__(a)})

        return outStr
        

def load_ncVar(varName, nc=None, **kwargs):
        
        if (nc is None) : raise Exception('No Netcdf file passed')
        
        var = nc.variables[varName]
        
        var.set_auto_maskandscale(False)
        
        #Load dimensions
        varDim = [str(dim) for dim in var.dimensions] #Revert the dimensions indices for numpy
#        varDim = [str(dim) for dim in var.dimensions][::-1] #Revert the dimensions indices for numpy
        missDim=len(varDim) == 0
        if (missDim): warn('No dimension found')
        else : varDimval = [len(nc.dimensions[dimname]) for dimname in varDim]
        
        ind_list = [] #Init index list
        dims = OrderedDict({'_ndims':0}) #Init dimensions
        
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
                if vn.startswith('lon') : dumvar=recale(dumvar,degrees=True)
                fg=(dumvar >= drange[0]) & (dumvar <= drange[1])
                if fg.sum() == 0 :
                    #retry switrhcing lon/lat
                    dumvar=recale(dumvar,degrees=True)
                    drange=tuple(recale(drange,degrees=True).astype(np.long))
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
        
        dstr=','.join(dstr) #invert dimension list for numpy
#        dstr=','.join(dstr[::-1]) #invert dimension list for numpy
        if missDim : cmd = 'varOut = var[:]'
        else : cmd = 'varOut = var[{0}]'.format(dstr)
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
        
        #Switch dimensions
        varOut=np.transpose(varOut,tuple(range(len(dims.keys()[1:]))[::-1]))
        
        #Build up output structure
        dims.update({'_ndims':len(dims.keys()[1:])})
        outStr = {'_dimensions':dims, 'data':varOut}
        
        #Add variable attributes
        for A in var.__dict__.keys():
            outStr[A]=var.getncattr(A)
        
        return outStr
#            ind_list=[[]] 
    
