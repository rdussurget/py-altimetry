import os
import glob
import inspect
import fnmatch

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as ncfile
import copy

from altimetry.tools.nctools import load_ncVar, load_ncVar_v2
try:
    import seawater.gibbs as gsw
    import seawater.csiro as csw
except ImportError:
    pass # module doesn't exist, deal with it.

#import alti_tools as atools
from scipy import interpolate
from warnings import warn
from altimetry.tools import recale_limits, in_limits, cumulative_distance, calcul_distance, \
    where_list, \
    cnes_convert, \
    plot_map, \
    get_caller
from collections import OrderedDict


class hydro_data(object):
    
    def __init__(self,file_pattern,limit=None,verbose=1,round=True,zero_2pi=True,output_is_dict=True,**kwargs):
        
        #Init system variables
#        if limit is None : limit=[-90.,0.,90.,360.]
        self.zero_2pi=zero_2pi
        
        self.limit_set=False
        if limit is None :
            limit=[-90.,0.,90.,360.]
        else : self.limit_set = True
        
        self.limit = np.array(recale_limits(limit, zero_2pi=self.zero_2pi)) 
            
        self.verbose = verbose
        self.fileid = np.array([])

        self.count=0
        self.size=0


        #Setup file list
        if isinstance(file_pattern, str) : ls=glob.glob(file_pattern)
        else :
            ls = file_pattern.tolist()
            file_pattern=file_pattern[0]
        
        if len(ls) == 0 :
            self.Error('File pattern not matched : '+file_pattern)
                
        self.filelist=[os.path.basename(j) for i,j in enumerate(ls)]
        self.filelist_count = [0]*len(self.filelist)
        enum=list(enumerate(ls))
        enum = zip(*enum)
        
        self.fid_list=np.array(enum[0])
        self.dirname=os.path.dirname(file_pattern)
        
        self.par_list=np.array([])
        self.dim_list=np.array([])
        self._dimensions={'_ndims':0}
        
        #Loop over data files
        #####################
        for i in np.arange(len(self.fid_list)) :
            
            #Read data file
            ###############
            filename = enum[1][i]
            self.message(1,"Loading "+os.path.basename(filename))
            
            self.current_file=enum[0][i]
            
            res=self.read(filename,output_is_dict=output_is_dict,**kwargs) #read() function is specific of each class
            self.update_dataset(res) #update class with loaded data
            
            self.check_variables()
        
        if not self.limit_set : self.limit=self.extension(round=round)
#        self.update()

    def update_dataset(self,dataStr,flatten=False):

        #Load keys and dimensions
        #########################
        dataDim = dataStr.pop('_dimensions')
        attrStr = dataStr.pop('_attributes',{})
        ndims = dataDim.pop('_ndims')
        dimensions = [dataDim.keys(),dataDim.values()]
        
        keys = dataStr.keys()
        self.message(2, 'Loaded variables : '+str(keys))
        
        #Check what is the current variable type
        isStructure = True if isinstance(dataStr[keys[0]],dict) else False
        
#        datalen = [np.size(dataStr[key]) for key in keys]
        datalen = [list(np.shape(dataStr[key]['data'])) for key in keys] if isStructure else [list(np.shape(dataStr[key])) for key in keys] #####!!!!!! WARNING!!!! NEW FEATURE TO TEST
        if isStructure :
            varDim = [list(dataStr[key]['_dimensions'])[1:][::-1] for key in keys]
            ind = [where_list(vDim,dimensions[0]) for vDim in varDim] #Dimensions indices from actual variables' dimensions
            #Check dimension lengths
            dimOk = np.array([enum[1][0] == dimensions[1][ind[enum[0]][0]] for enum in enumerate(datalen)])
            if (~dimOk).sum() > 0 : 
                notOk = np.where(~dimOk)[0]
                self.Error('Problem with {0} variables : {1}'.format(len(notOk),','.join(np.array(dataStr.keys())[notOk])))
        else :
            ind = [where_list(dlen,dimensions[1]) for dlen in datalen] #Dimensions indices from variable length
            if (np.array(ind).sum() == -1)!= 0 : self.Error('At least one variable have not been properly defined')
            
        dimname = [np.array(dimensions[0])[i].tolist() for i in ind]  #Get correspondance between data structure dimensions and variables
                
        curDim, nself=self.get_currentDim()
        createDim=np.array([np.array([w == -1 for w in where_list(j, curDim[0])]) for i,j in enumerate(dimname) ])
        createDim=np.squeeze(createDim)
#        curInd = atools.where_list(dimname_reduced,curDim[0]) #Get correspondance between data structure dimensions and object dimensions
        
#        createDim = (np.array(curInd) == -1) #Get dimensions to be created   
        toCreate = np.array([not self.__dict__.has_key(key) for key in keys])
        
        updateDim=[]  
        
        self.message(2, 'Updating object with '+str(['{0}({1}:{2})'.format(i[0],i[1],i[2]) for i in zip(*(keys,dimname,datalen))]))
        
        for enum in enumerate(keys) :
            
            ind=enum[0]
            key=enum[1]
            
            #Load variable
            ##############
#            var=dataStr.get(key)
            dum=dataStr.get(key).get('data') if isStructure else dataStr.get(key) 
            if flatten :
                if isStructure :
                    dum['data']=dum['data'].flatten()
                else : dum=dum.flatten
#            if isinstance(dum,np.ma.masked_array) : 
            
#            #Get associated dimensions
#            ##################################
#            datalen = datalen[ind]#[len(dataStr[key]) for key in keys]
#            ind = atools.where_list([datalen],dimensions[1])[0]
#            if (ind == -1) : self.Error('Dimensions of current variable ('+key+') have not been properly defined')
#            dimname = dimensions[
            
            #Initialize variable if required
#            if toCreate :
#            updateDim.append(self.create_Variable(key, dum, dimensions={dimname[ind]:datalen[ind]},toCreate=toCreate[ind],createDim=createDim[ind]))
            updateDim.append(self.create_Variable(key, dum, dimensions=dict(zip(dimname[ind],datalen[ind])),toCreate=toCreate[ind],createDim=createDim[ind]))

        
        
        #Final sequence
        
        zipped_upd=zip(*(np.hstack(dimname)[~np.hstack(createDim)],np.hstack(datalen)[~np.hstack(createDim)]))
        updateDim_List = np.array(list(set(tuple(i) for i in np.array(zipped_upd,dtype='|S16').tolist()))) #2D unique
#        updateDim_List = np.unique(np.array(zipped_upd,dtype='|S16')) #[str(i) for i in datalen]
#        if updateDim_List.size > 0 : updateDim_List.resize((2,updateDim_List.size/2))
#        updateDim_List = np.unique(zip(*(np.array(dimname)[~createDim],np.array(datalen)[~createDim]))) #[str(i) for i in datalen]

        zipped_dims=zip(*(np.hstack(dimname)[np.hstack(createDim)],np.hstack(datalen)[np.hstack(createDim)]))
        createDim_list = np.array(list(set(tuple(i) for i in np.array(zipped_dims,dtype='|S16').tolist()))) #2D unique
#        clist, inv = np.unique(np.array(zipped_dims,dtype='|S16'),return_inverse=True) #RQ : THIS WILL FAIL IF NUMBERS HAVE MORE THAN 16 DIGITS #[str(i) for i in datalen]
#        if createDim_list.size > 0 : createDim_list.resize((2,createDim_list.size/2))
        
#        createDim_list = np.unique(zip(*(np.array(dimname)[createDim],np.array(datalen)[createDim]))) #[str(i) for i in datalen]
        
        for dname,dim in createDim_list :
            self.create_Dim(dname, np.int(dim))
        
        for dname,dim in updateDim_List:
            self.update_Dim(dname, np.int(dim))
        
            
    #                curInd = atools.where_list([l_datalen],curDim[1])[0]
    #                createDim = (curInd == -1) #Get dimensions to be created 
            
#            else :
#                pass
            
        
        
#        #No data present in the object -> initialisation phase
#        if self._dimensions['_ndims'] == 0 :
#            self._dimensions=szstr
#            
#        #Append or create data sets.
#        else :
#        
#            pass
#        #Get current dimensions
#        self._dimensions['_ndims']
#        
#        sz = np.array([szstr.get(szstr.keys()[j]) for j in np.arange(ndims)]) #reforms a vector containing dimensions
#       for j in ndims : sz[i] = szstr.get(szstr.keys()[j])
#            
#        #Fill class looping over keys
#        #############################
#        keys=dataStr.keys()
#        self.message(4, 'Loaded variables : '+str(keys))
#        for key in keys :
#            
#            #Load variable and flatten 2D matrix into vector (for concatenation)
#            cmd='dum=dataStr[\''+key+'\'][:].flatten()'
#            self.message(4,'exec : '+cmd)
#            exec(cmd)
#            self.update_variable(key,dum)
#                
#                #Init variables
#                if i == 0 : cmd='self.'+key+'=dum'
#                
#                #Append data to existing variables
#                else :
#                    #Use np.ma.concatenate if masked array
#                    #and np.append for other cases
#                    if isinstance(dum,np.ma.core.MaskedArray) : cmd = 'self.'+key+'=np.ma.concatenate((self.'+key+',dum))'
#                    else : cmd='self.'+key+'=np.append(self.'+key+',dum)'
#                self.message(4,'exec : '+cmd)
##                exec(cmd)
#
#    def update_variable(self,name,data):
#        
#        nelts=len(data)
#        
#        
#        
#        #Test wether the variable exists or not
#        if self.__dict__.has_key(name) :
#        
#            
#        
#        #variable do not exist
#        else :
#        
#        pass
#
#    #Update stuff        
    
    def check_variables(self):
        """
        CHECK_VARIABLES : Forces variables to respect dimensions
        """
        
        self.count = len(self.fileid)
        self.size = np.size([np.size(self.__dict__.get(par)) for par in self.par_list])
        
        infos = zip(*(self.par_list,self.dim_list))
#        curDim, nself = self.get_currentDim()
        
        for enum in enumerate(infos):
            varSize = np.size(self.__dict__.get(enum[1][0]))
            varShape = np.shape(self.__dict__.get(enum[1][0]))[::-1] #Data and Netcdf dimensions are inverted (not allways?)
            dimSize = tuple([(self._dimensions).get(d) for d in enum[1][1]])
            masked = isinstance(self.__dict__.get(enum[1][0]), np.ma.masked_array)
            
            #Check mask consistency (mask length should be the same as data)
            if masked :
                if (self.__dict__.get(enum[1][0]).mask.size != self.__dict__.get(enum[1][0]).data.size): raise np.ma.core.MaskError("Mask length is not consistent with data")
            
            #Check dimensions
            self.message(4, 'checking variables -> {0}(N={1}) - {2}:{3}'.format(enum[1][0],varSize,enum[1][1],dimSize))
            for n,sh in enumerate(varShape) :
                if (sh > dimSize[n]) :
                    self.Error('Object variable {0} greater than corresponding dimension ({1})'.format(enum[1][0],enum[1][1]))
                elif (sh < dimSize[n]):
                    self.message(3, 'Variable {0}(N={1}) being extended to match dimension {2}:{3}'.format(enum[1][0],varSize,enum[1][1],dimSize))  
    #                self.__dict__[enum[1][0]] = np.ma.concatenate((self.__dict__[enum[1][0]], np.ma.masked_array(np.repeat(np.nan,dimSize - varSize),mask=np.zeros(dimSize - varSize,dtype='bool'))))
                    self.__dict__[enum[1][0]] = np.ma.masked_array( np.append(self.__dict__[enum[1][0]].data,np.repeat(self.__dict__[enum[1][0]].fill_value if hasattr(self.__dict__[enum[1][0]],'fill_value') else np.NaN,dimSize[n] - varSize)),
                                                                    mask=np.append(self.__dict__[enum[1][0]].mask,np.ones(dimSize[n] - varSize,dtype='bool')) )
    
    def create_Dim(self, name,value):
        if not self._dimensions.has_key(name) :
            self.message(3, 'Create dimension {0}:{1}'.format(name,value))
            self._dimensions[name]=value
            self._dimensions['_ndims']=len(self._dimensions) - 1
        else :
            self.message(3, 'Dimension {0} already exists'.format(name))
        
    def update_Dim(self,name,value):
        oldVal=self._dimensions[name]
        self._dimensions[name] += value
        self.message(2, 'Updating dimension {0} (from {1} to {2})'.format(name,oldVal,self._dimensions[name]))
    
    def update_fid_list(self,filename,N):
        self.filelist_count[self.filelist.index(filename)] = N
        fid=self.fid_list.compress([enum[1][0] == os.path.basename(filename) for enum in enumerate(zip(*(self.filelist,self.fid_list)))])
        self.__dict__.update({'fileid' :np.append(self.fileid,np.repeat(fid,N))})
    
    def delete_Variable(self,name):
        self.message(1,'Deleting variable {0}'.format(name))
        self.__dict__.pop(name)
        self.par_list=self.par_list[self.par_list != name]
        
    
    def create_Variable(self,name,value,dimensions,toCreate=None,createDim=None,extend=True):
        """
        create_Variable : This function adds a variable to the current data object
        """
        
        #Check variable name
        ####################
        #This allows to remove impossible variable names
        #!!!! This is not a good solution
        name=name.replace('.','_')
        
        #Check if data is structured or not
        isStructure = True if isinstance(value,dict) else False
        
        #Get dimensions
        dimName = np.array(dimensions.keys())
        dimVal = np.array(dimensions.values())
        keys=np.array(self._dimensions.keys())
        
#        if createDim is None : createDim = self._dimensions.has_key(dimName[0])
        createDim = np.array([not self._dimensions.has_key(dim) for dim in dimName]) if createDim is None else np.array(createDim)
        if toCreate is None : toCreate = np.sum(self.par_list == name) == 0
        
        self.message(3,'Loading {0} ({1}:{2}) from {3}'.format(name,dimName,dimVal,os.path.basename(self._filename)))

        #Cast variable into masked array first
        ######################################
        if (not isinstance(value['data'],np.ma.core.MaskedArray) if isStructure else not isinstance(value,np.ma.core.MaskedArray)) :
            value = np.ma.masked_array(value['data'],mask=np.zeros(tuple(dimVal),dtype='bool')) if isStructure else np.ma.masked_array(value,mask=np.zeros(tuple(dimVal),dtype='bool'))
            self.message(4,'Casting variable to np.ma.MaskedArray')
        
        #Restructure dataset if structure
        if isStructure :
            dumvalue=value.pop('data')
            for a in value.keys() :
                dumvalue.__setattr__(a,value[a])
            value=dumvalue
        
        
        curDim, nself=self.get_currentDim()
        
        
        curInd=np.array(where_list(dimName,curDim[0]))
        curDimVal=np.array(where_list(dimVal,curDim[1]))
        
        existDims= (curInd != -1)
        createDim = (curInd == -1)
        createInd = np.where(createDim)[0]
        
        appendDim=existDims & (curDimVal == -1)
        appendInd=curInd[appendDim]
        
        
#        curInd = set(atools.where_list(dimVal,curDim[1])).intersection(set(atools.where_list(dimName,curDim[0]))) 

        #Get dims to be created
        #######################
        
        

        #Choose case between all different solutions :
        ##############################################
        # 1: create a new variable with at least 1 new dimension
        # 2: extend -> create a new variable using existing dimensions
        # 3: append exisiting variable with data
        # 4: impossible case ?
        if createDim.any() & toCreate :           
            #Create Variable
            self.message(4,'Create variable '+name)
#            self.__setattr__(name,value)
#            cmd='self.'+name+'=value'

            #Append variable infos to object
            self.par_list=np.append(self.par_list,name)
            dimlist_cp=self.dim_list.tolist()
            dimlist_cp.append(dimName.tolist())
            self.dim_list=np.array(dimlist_cp) #np.append(self.dim_list,dimName.tolist())
        
            updateDim=False
        
        elif (not createDim.any()) & toCreate :
            
            #extend variable
            if extend : value = np.ma.masked_array(np.append(np.zeros(curDim[1][curInd]),value.data),mask=np.append(np.ones(curDim[1][curInd],dtype='bool'),value.mask))
            
            self.message(4,'Extend variable '+name)
#            self.__setattr__(name,value)
#            cmd='self.'+name+'=value'
#            self.message(4,'exec : '+cmd)

            #Append variable infos to object
            self.par_list=np.append(self.par_list,name)
            dimlist_cp=self.dim_list.tolist()
            dimlist_cp.append(dimName.tolist())
            self.dim_list=np.array(dimlist_cp)
#            self.dim_list=np.append(self.dim_list,dimName)
            
            updateDim=True
        
        elif (not createDim.any()) & (not toCreate) :
            #append variable
            self.message(4,'Append data to variable '+name)
            value = np.ma.masked_array(np.append(self.__getattribute__(name).data,value.data),mask=np.append(self.__getattribute__(name).mask,value.mask))
#             cmd = 'self.'+name+'=np.ma.concatenate((self.'+name+',value),mask=np.append(self.dac))' #buggy command : concatenate do not always work..
#            cmd = 'self.{0}=np.ma.masked_array(np.append(self.{0}.data,value.data),mask=np.append(self.{0}.mask,value.mask))'.format(name)
        
#            exec('res = [self.{0}.mask.size,value.mask.size]'.format(name))
#                raise 'error!'
        
            updateDim=True
        
        elif createDim.any() & (not toCreate) :
            #Impossible case ?
            self.Error('Impossible case : create dimensions and variable {0} already existing'.format(name))
        
#        self.message(4,'exec : '+cmd)
        
        if name == 'nsct' :
            pass
        
        #Append dimensions to variable
        if not dimensions.has_key('_ndims') :
            dumDim={'_ndims':len(dimensions.keys())}
            dumDim.update(dimensions)
            dimensions=dumDim.copy()
        value.__setattr__('_dimensions',dimensions)
        
        try : self.__setattr__(name,value)
        except np.ma.core.MaskError : raise 'mask error' 
#           
#        try : exec(cmd)
#        except np.ma.core.MaskError :
#            raise 'mask error' 
#        exec(cmd)
         
        return updateDim
        
#        if toCreate :
#            
#            #This variable do not exists and needs to be created
#            if not createDim :
#                dumVar = np.ma.masked_array(np.zeros(curDim[1][curInd]),mask=np.repeat(True,curDim[1][curInd]))
#                value = np.ma.concatenate((dumVar,value))
#            cmd='self.'+name+'=value'
#            self.message(4,'exec : '+cmd)
#            exec(cmd)
#            
#
#        
#        #data should be append to existing data
#        else :
#            pass
#        
        
        
#             #Append data to existing variables
#            else :
#                #Use np.ma.concatenate if masked array
#                #and np.append for other cases
#                if isinstance(dum,np.ma.core.MaskedArray) : cmd = 'self.'+key+'=np.ma.concatenate((self.'+key+',dum))'
#                else : cmd='self.'+key+'=np.append(self.'+key+',dum)'
#            self.message(4,'exec : '+cmd)
            
    
    def get_currentDim(self):
        selfDim = self._dimensions.copy()
        if selfDim.has_key('_ndims') : nself = selfDim.pop('_ndims')
        else : 
            self.warning(1, 'self._dimensions does not have the _ndims key')
            nself = len(selfDim)
        curDim = [[key for key in selfDim.keys()],[selfDim[key] for key in selfDim.keys()]]
        return curDim, nself
    
    def update(self):
        self.message(0, "DEPRECATED function : -> hydro_data.update() called by "+inspect.stack()[2][3])
        
#        #Compute density if all variable exist
#        if hasattr(self,'psal') &  hasattr(self,'temp') & hasattr(self,'pres') : 
#            self.rho= gsw.rho(self.psal,self.temp,self.pres)
        
#        keys=self.__dict__.keys()
#        
#        for key in keys:
#            count=[np.size(0)]

    def message(self,MSG_LEVEL,str):
        """
         MESSAGE : print function wrapper. Print a message depending on the verbose level
         
         @param MSG_LEVEL {in}{required}{type=int} level of the message to be compared with self.verbose
         
         @example self.message(0,'This message will be shown for any verbose level') 
        """
        
        caller=get_caller()
        if MSG_LEVEL <= self.verbose :  print('[{0}.{1}()] {2}'.format(__name__,caller.co_name,str))
    
    def warning(self,MSG_LEVEL,str):
        """
         WARNING : Wrapper to the warning function. Returns a warning when verbose level is not 0. 
         
         @param MSG_LEVEL {in}{required}{type=int} level of the message to be compared with self.verbose
         
         @example self.waring(1,'Warning being issued) 
        """
        
        if self.verbose <= 1 : warn(str)
    
    def Error(self,ErrorMsg):    
        raise Exception(ErrorMsg)
    
    def copy(self,deep=True):
        return copy.deepcopy(self) if deep else copy.copy(self)
    
    def slice(self,param,range,surf=False,dimension=None):
        if np.size(range) == 2 :
            if not surf : fg = (self.__dict__[param] >= np.min(range)) & (self.__dict__[param] < np.max(range))
            else : fg = (self.__dict__[param+'_surf'] >= np.min(range)) & (self.__dict__[param+'_surf'] < np.max(range))
        elif np.size(range) == 1 :
            if not surf : fg = (self.__dict__[param] == range)
            else : fg = (self.__dict__[param+'_surf'] == range)
        else : self.Error('Range array must have 1 or 2 elements max')
        return fg

    def time_slice(self,timerange,surf=False):
        return self.slice('date',timerange,surf=surf)
#        if not surf : return (self.date >= timerange[0]) & (self.date < timerange[1])
#        else : return (self.date_surf >= timerange[0]) & (self.date_surf < timerange[1])

    def update_with_slice(self,flag):
        
        N=flag.sum()
        
        #Get object attributes to update
        varSize = np.array([np.size(self.__dict__[k]) for k in self.__dict__.keys()])
        par_list=np.array(self.__dict__.keys())[varSize == self._dimensions['time']].tolist()
        
        self._dimensions['time']=N
        
        for par in par_list :
            self.__setattr__(par,self.__dict__[par][flag])
            if (hasattr(self.__dict__[par], '_dimensions')) : self.__dict__[par]._dimensions['time']=self._dimensions['time']

    def get_file(self,pattern):
        flag=[fnmatch.fnmatch(l,pattern) for l in self.filelist]
        id=self.fid_list.compress(flag)
        flag = self.fileid == id
        return flag
    
#    def get_id(self,id):
#        return self.id == id
    
    def time_range(self,flag=None):
        if flag is None : return cnes_convert([self.date.min(),self.date.max()])
        else : return cnes_convert([self.date.compress(flag).min(),self.date.compress(flag).max()])
    
    def extension(self,flag=None,round=True):
        if flag is None : limit = [self.lat.min(),self.lon.min(),self.lat.max(),self.lon.max()]
        else : limit = [self.lat.compress(flag).min(),self.lon.compress(flag).min(),self.lat.compress(flag).max(),self.lon.compress(flag).max()]
        if round :
            limit[0]=np.floor(limit[0])
            limit[2]=np.floor(limit[2])
            limit[1]=np.ceil(limit[1])
            limit[3]=np.ceil(limit[3])
        return limit
    
    def get_id_list(self,flag=None):
        if flag is None : return np.unique(self.id)
        else : return np.unique(self.id.compress(flag))
    
    def get_stats(self,flag):
        par_list=self.par_list.compress([(not par.endswith('_surf')) & ( par != 'id') for par in self.par_list])
        valid= np.array([len(self.__dict__[par].compress(flag).compressed()) for par in par_list])
        per_valid = (100.0* valid) /float(np.sum(flag))
        
        fname = [self.filelist[i] for i in np.unique(self.fileid.compress(flag)).astype('i')]
        trange = self.time_range(flag)
        extent = self.extension(flag)
        N = np.sum(flag)
        avail_par = par_list.compress(per_valid > 0)
        avail_par_per = per_valid.compress(per_valid > 0)
        
        return par_list, valid, per_valid, fname, trange, extent, N, avail_par, avail_par_per
    
    def get_platform_stats(self,id):
        par_list, valid, per_valid, fname, trange, extent, N, avail_par, avail_par_per = self.get_stats(self.id == id)
        return (fname,trange,extent,N,avail_par,avail_par_per)

    def get_object_stats(self):
        return self.get_stats(np.ones(self.count,dtype='bool'))
        
    def platform_summary(self,id,col='.k'):
        stats = self.get_platform_stats(id)
        self.message(0, '\tPlatform {0} : '.format(id))
        self.message(0, '\t-> file : {0}'.format(stats[0]))
        self.message(0, '\t-> from : '+' - '.join(map(str,stats[1][0])))
        self.message(0, '\t-> extent : ['+', '.join(['{0:.1f}'.format(x) for x in stats[2]])+']')
        self.message(0, '\t-> size : {0} pts'.format(stats[3]))
        self.message(0, '\t-> variables :'+', '.join(['{0}({1:.0f} %)'.format(i[0],i[1]) for i in zip(*(stats[4],stats[5]))])+']')
        
    def map(self, flag=None,  fname=None, zoom=False, pmap=None, show=True, **kwargs):
                
        if zoom : limit = self.extension(flag)
        else : limit = self.limit
        
        if pmap is None : pmap=plot_map(0,0,0,limit=limit)
        p,=self.plot_track(pmap, flag,**kwargs)
        
        if show : 
            pmap.setup_map()
            pmap.show()

        return p,pmap

    def summary(self,all=False,fig=None,col='.k',legend=None,**kwargs):
        """
        SUMMARY : Returns a summary of all the data contained in the current object
        """
        
        par_list, valid, per_valid, fname, trange, extent, N, avail_par, avail_par_per = self.get_object_stats()
        
        #Print parameter list
        self.message(0, '\n\t  DATA SUMMARY')
        self.message(0, '\tObject type <class \'{0}\'> ({1})'.format(self.__class__.__name__,','.join(map(str,self.__class__.__bases__))))
        self.message(0, '\tLoaded from : '+self.dirname)#', '.join(map(str,[os.path.basename(file) for file in self.filelist])))
        self.message(0, '\tLimits = '+', '.join(map(str,self.limit)))       
        self.message(0, '\tDate range : '+' - '.join(map(str,trange[0])))
        self.message(0, '\tSpatial extension: ['+', '.join(['{0:.1f}'.format(x) for x in extent])+']')
        self.message(0, '\tRecords : [{0}]'.format(N))
        self.message(0, '\tAvailable parameters : '+', '.join(map(str,self.par_list))) 
        self.message(0, '\tParam summary :'+', '.join(['{0}({1:.0f} %)'.format(i[0],i[1]) for i in zip(*(avail_par,avail_par_per))])+']')

        
        if isinstance(fig,str) :
            ffig = fig+'dataset.png'
            show=False
        else :
            ffig = fig
            show=True
        
        
        if fig is not None :
            if self.__dict__.has_key('id_surf') : n=len(self.id_surf)
            else : n=N
#            p,pmap=self.map(np.ones(n,dtype='bool'),ffig,col=col,show=show,**kwargs)
            p,pmap=self.map(col=col,show=show,**kwargs)
            pmap.setup_map()
            if legend is not None : plt.legend([p],[legend])
            if not show :
                pmap.savefig(ffig)
                plt.clf()

        if all is not True :
            return
        self.message(0, '\t##### Details by float ID')
        
        #Get all platform ID's
        id_list=self.get_id_list().astype('a')
        for id in id_list : 
            self.platform_summary(id)
            if isinstance(fig,str) : ffig = fig+'{0}.png'.format(str(id))
            else : ffig = None
            if fig is not None : 
                if self.__dict__.has_key('id_surf') :
                    flag=self.id_surf == id
                    n=len(self.id_surf.compress(flag))
                else :
                    flag=self.id == id
                    n=len(self.id.compress(flag))
                p,pmap=self.map(flag=flag,col=col,show=show,**kwargs)
                pmap.setup_map()
                if legend is not None : plt.legend([p],['#{0}'.format(id)])
                if show is not True :
                    plt.savefig(ffig)
                    plt.clf()
            
#            for enum in enumerate(par_list.compress(per_valid > 0)) : print '\t{0} (valid:{1:4.1f}%)'.format(enum[1], per_valid.compress(per_valid > 0)[enum[0]])
#            per_valid.compress(per_valid > 0)
        
    def in_limits(self,limit=None):
        if limit is None : limit = self.limit
        flag=in_limits(self.lon, self.lat, limit)
        return flag

    def plot_track(self,pmap,flag=None,col='.k',endpoint='*r',endpoint_size=None,title=None,fontsize=8,textcolor='b',ms=5,linewidth=1,**kwargs):
        
        
        if self.__dict__.has_key('lon_surf') :
            if flag is None : flag=np.ones(self.lon_surf.size,dtype='bool')
            lon=self.lon_surf.compress(flag)
            lat=self.lat_surf.compress(flag)
            dat=self.date_surf.compress(flag)
            id=self.id_surf.compress(flag)
        else :
            if flag is None : flag=np.ones(self.lon.size,dtype='bool')
            lon=self.lon.compress(flag)
            lat=self.lat.compress(flag)
            dat=self.date.compress(flag)
            id=self.id.compress(flag)
        
        if endpoint_size is None : endpoint_size=2*ms
        
        id_list=np.unique(id)
        cnt=np.size(id_list)
        
#        #Go to next iteration when no data
#        if cnt == 0 :
#            continue
        
        for j in id_list :
            dumflag = (id == j)
            dumdat=dat.compress(dumflag)
            sort = dumdat.argsort()
            dumlon=lon.compress(dumflag)[sort]
            dumlat=lat.compress(dumflag)[sort]
            dumid=id.compress(dumflag)[sort]
            dumcnt=np.size(dumid)
            dumvalid=(~dumlon.mask & ~dumlat.mask & ~dumid.mask).sum()
            if dumvalid > 0:
                if isinstance(col,str) : p=pmap.plot(dumlon,dumlat,col,ms=ms,linewidth=linewidth,**kwargs)
                else : p=pmap.scatter(dumlon,dumlat,col.compress(dumflag)[sort],s=ms,linewidth=linewidth,**kwargs)
                pmap.plot(dumlon[dumcnt-1],dumlat[dumcnt-1],endpoint,ms=endpoint_size)
                pmap.text(dumlon[dumcnt-1],dumlat[dumcnt-1],'{0}'.format(j),fontsize=fontsize,color=textcolor)
                if title is not None : pmap.title(title)
#                 print j
#                 pmap.show()
            
        try : return p
        finally : return pmap.plot(-100,-370,'')
      
    def plot_track_old(self,*args,**kwargs):
        '''
        PLOT_TRACK: plot a surface map of sampling track
        
        @param map {in}{optional}{type=atools.pmap} Input plot_map or pyplot object. Must be passed as 1st argument
        @param lon {in}{required}{type=numeric} longitude
        @param lat {in}{required}{type=numeric} latitude
        @param var {in}{optional}{type=numeric} variable to be displayed (scatter plot)
        
        @keyword julian {in}{optional}{type=BOOLEAN} True if output date is in julian format
        
        @examples 1) Converts an array of calendar dates into julian dates<br />
        
        @author : Renaud DUSSURGET (RD), IFREMER LER/PAC, La Seyne/Mer
        @history
           - Created by RD (May 2012)
        '''
        
        #Get map object (plot_map or Basemap)
        flag_list=[inspect.ismodule(type(i)) | inspect.ismodule(type(i)) for i in args]
        
        if (np.array(flag_list).max()) :
            map_obj = args[0]
            obj_plot=map_obj
            args.pop(0)
        else :
            obj_plot=plt

        
        #Retrieve other params
        lon=args[0]
        lat=args[1]
        
        #override if passed directly through 3rd argument
        if len(args) == 3 :
            var = args[2]
            scatter=True
        else :
            if 'c' in kwargs :
                c = kwargs['c']
                del kwargs['c']
            else : c = '.k'
            scatter=False
        
        if scatter :obj_plot.scatter(lon,lat,var,**kwargs)
        else : obj_plot.plot(lon,lat,**kwargs)
    
    def plot_transect(self, x, z, var,
                      xrange=None, zrange=None,
                      vmin=None,vmax=None, #colormap limits
                      xstep=1,zstep=10,    #Interpolation step
                      s=10, edgecolor='none',
                      **kwargs):

        ax=plt.gca()
        if isinstance(var,np.ma.masked_array) :
            flag=np.reshape(var.mask ,var.size)
            values=var.compress(~flag)
            x=x.compress(~flag)
            z=z.compress(~flag)
        else : values = var
        
        return ax.scatter(x,z,c=values,s=s,edgecolor=edgecolor,vmin=vmin,vmax=vmax)
#        plt.show()
        
        pass
    
    def contour_transect(self, x, z, var, xrange=None,zrange=None, xstep=1,zstep=10,vmin=None,vmax=None,marker='.k', **kwargs):
        
        
        if xrange is None : xrange=(x.min(),x.max())
        if zrange is None : zrange=(z.min(),z.max())
        
        gx=np.arange(xrange[0],xrange[1],xstep)
        gz=np.arange(zrange[0],zrange[1],zstep)

        grid_x, grid_y = np.meshgrid(gx,gz) 
        
        #Remove masks
        flag=np.reshape(var.mask ,var.size)
        values=var.compress(~flag)
        x=x.compress(~flag)
        z=z.compress(~flag)
        
        npts=values.size        
        
#        plt.tricontourf(x,z,values) #VERY long!!
        
        points = zip(*(np.reshape(x,npts),np.reshape(z,npts)))
        Z = interpolate.griddata(points, values, (grid_x, grid_y), method='cubic')
#        plt.imshow(gdist,gdepth,Z,vmin=13.1,vmax=13.6)
        plt.pcolormesh(gx,gz,Z,vmin=vmin,vmax=vmax)
        if marker is not None : plt.plot(x,z,marker,ms=1)
#        plt.scatter(dist2d,depth,s=10,c=values,edgecolor='none',vmin=13.1,vmax=13.6)
#        plt.show()

        return Z, gx, gz

    def read_ArgoNC(self,filename,params=None,force=False,dephrange=None,timerange=None,**kwargs):
        """
        READ_ARGONC : Argo NetCDF reader
        
        @return: outStr {type:dict} Output data structure containing all recorded parameters as specificied by NetCDF file PARAMETER list.
        @author: Renaud Dussurget
        """
        
        #Open file
        self._filename = filename
        self._ncfile = ncfile(self._filename, "r")
     
        #Get list of recorded parameters:
        dum=self._ncfile.variables['STATION_PARAMETERS'][0,:]
        nparam=np.shape(dum)[0]
        par_list=np.array([''.join(self._ncfile.variables['STATION_PARAMETERS'][0,i,:].compressed()) for i in np.arange(nparam)])
        #remove empty items and update nparam
        par_list=par_list.compress([len(par) != 0 for par in par_list])
        nparam=par_list.size
        
        if params is not None :
            if force : par_list=[i.upper() for i in params]
            else :par_list=list(set(params).intersection(par_list))
        else : par_list=par_list.tolist()
        
        self.message(1,'Recorded parameters : '+str(nparam)+' -> '+str(par_list))
     
        lon = self.load_ncVar('LONGITUDE',**kwargs)
        lat = self.load_ncVar('LATITUDE',**kwargs)
        date = self.load_ncVar('JULD',**kwargs)
#        lon = self._ncfile.variables['LONGITUDE'][:]
#        lat = self._ncfile.variables['LATITUDE'][:]
        
        #Extract within limits
        ind, flag = in_limits(lon['data'],lat['data'],limit=self.limit)
        if timerange is not None :
            dflag = (date['data'] >= np.min(timerange)) & (date['data'] < np.max(timerange))
            flag = flag & dflag
            ind = np.where(flag)[0]
        dim_lon = lon['_dimensions']
        dim_date = date['_dimensions']
        lat['data'] = lat['data'].compress(flag)
        lon['data'] = lon['data'].compress(flag)
        date['data'] = date['data'].compress(flag)
        
        #Update dimensions
        lat['_dimensions']['N_PROF']=flag.sum()
        lon['_dimensions']['N_PROF']=flag.sum()
        date['_dimensions']['N_PROF']=flag.sum()
#        dist=cumulative_distance(lat['data'], lon['data'])
        
        sz=np.shape(lon)
        ndims=np.size(sz)

#        if dephrange is not None :
#            pres = self.load_ncVar('PRES',N_PROF=ind,**kwargs)
#            deph=gsw.z_from_p(pres['data'],lat[0])
#            depind=(np.abs(tutu-(-np.max(dephrange)))).argmin(1)
#        if timerange is not None :deph=gsw.z_from_p(pres['data'],lat[0])
            
        date = self.load_ncVar('JULD',N_PROF=ind,**kwargs)
        dimStr = date['_dimensions']
#        date=date['data']
        
#        #Create dimension structure
#        curDim = [str(dimname) for dimname in dimStr.keys()[1:]] #[str(dimname) for dimname in self._ncfile.variables['LONGITUDE'].dimensions]
#        curDimval = [dimStr[dim] for dim in curDim] #[len(self._ncfile.dimensions[dimname]) for dimname in curDim]
#   
#        outStr={'_dimensions':{'_ndims':ndims,'nbpoints':sz[0]},'lon':lon,'lat':lat,'date':date}
        outStr={'_dimensions':dimStr,'lon':lon,'lat':lat,'date':date}
        
        for param in par_list :
            dumVar = self.load_ncVar(param,N_PROF=ind,**kwargs) #Load variables
            dimStr=dumVar.get('_dimensions')
#            dimStr=dumVar.pop('_dimensions')
#            dumVar['data'].__setattr__('_dimensions',dimStr)
            
            #update current object dimensions with missing variable dimensions
            curDim = [str(dimname) for dimname in dimStr.keys()[1:]] #[str(dimname) for dimname in self._ncfile.variables['LONGITUDE'].dimensions]
            curDimval = [dimStr[dim] for dim in curDim] #[len(self._ncfile.dimensions[dimname]) for dimname in curDim]
#            curDim = [str(dimname) for dimname in self._ncfile.variables[param].dimensions]
#            curDimval = [len(self._ncfile.dimensions[dimname]) for dimname in curDim]
            flag = [(np.array(dimname) == outStr['_dimensions'].keys()).sum() == 0 for dimname in curDim] #find dimensions to update
            dimUpdate = np.array(curDim).compress(flag)
            for enum in enumerate(dimUpdate) : 
                self.message(2, 'Appending dimensions {0}:{1} to dataStructure'.format(enum[1],np.array(curDimval).compress(flag)[enum[0]]))
                outStr['_dimensions'].update({enum[1]:np.array(curDimval).compress(flag)[enum[0]]}) #Append new dimension
                outStr['_dimensions']['_ndims']+=1 #update dimension counts
            
            self.message(4, 'Loading {0} ({1})'.format(param.lower(),','.join(str(dimStr[key]) for key in dimStr.keys()[1:])))
            dumStr=dumVar.copy()
#            locals()[param.lower()] = {param.lower():dumVar['data']}
#            cmd = 'dumStr = {\''+param.lower()+'\':dumVar[\'data\']}'
#            self.message(4, 'exec : '+cmd)
#            exec(cmd)
            outStr.update({param.lower():dumStr})
        
#        id=np.array([''.join([num for num in j.compressed()]) for j in self._ncfile.variables['PLATFORM_NUMBER'][ind]]) #get platform id into array
        pid=self.load_ncVar('PLATFORM_NUMBER',N_PROF=ind,**kwargs)
        pid_var=np.ma.MaskedArray([''.join([num for num in j.compressed()]) for j in pid['data']],mask=lon['data'].mask)
        pid['data']=pid_var.copy()
        for key in pid['_dimensions'].keys() :
            if key.startswith('STRING') :
                pid['_dimensions'].pop(key)
                pid['_dimensions']['_ndims']-=1
        
        
        outStr.update({'id':pid})
     
        self._ncfile.close()
        return outStr

    def load_ncVar(self,varName,nc=None,**kwargs):
        return load_ncVar_v2(varName, nc=self._ncfile, **kwargs)

class buoy_data(hydro_data):
    
    def __init__(self,file_pattern,**kwargs):
        
        hydro_data.__init__(self,file_pattern,**kwargs)
        
       
#        self.count=np.size(self.lon)

    def read(self,filename,**kwargs):
        
        fname,extension = os.path.splitext(filename)
        
#        if extension == '.nc' : outStr=self.read_ArgoNC(filename,N_LEVELS=0,**kwargs)
        if extension == '.nc' : outStr=self.read_ArgoNC(filename,**kwargs)
        elif extension == '.dat' : outStr = self.read_txt(filename)
        elif extension == '.asc' :
            kwread=kwargs.get('lon_name',{})
            kwread=kwargs.get('lat_name',{})
            outStr = self.read_asc(filename,**kwread) #CLS dumps
        else : self.Error('Unknown formatting')
        
        
        #Rework dimensions
        ##################
        # Dimensions may be in 2D (more not compliant with current code)
        # Dataobject only accepts flattened (1D) variables
        dimStr=outStr.pop('_dimensions')
#        if dimStr['_ndims'] != 2 : self.Error('Variables with {0} dimensions can not be used as glider profiles (2 dimensions)'.format(dimStr['_ndims']))
        if not dimStr.has_key('N_LEVELS') or not dimStr.has_key('N_PROF') :
            self.Error('Data structure must contain dimensions N_LEVELS and N_PROF (only contains {0})'.format(dimStr.keys()[1:]))
        
        datalen = [np.size(outStr[key]) for key in outStr.keys()] #[np.shape(outStr[key]) for key in outStr.keys()] 
                
        nblevels = dimStr['N_LEVELS']
        nbprofiles = dimStr['N_PROF']
        nbpoints = nblevels*nbprofiles
        
        #Change dimension names towards new ones
        for key in outStr.keys() :
            if outStr[key].has_key('_dimensions'):
                outStr[key]['_dimensions'].pop('N_LEVELS',None)
                outStr[key]['_dimensions'].pop('N_PROF',None)
#                if outStr[key]['_dimensions'].has_key('N_PROF') : outStr[key]['_dimensions'].update({'nbprofiles':outStr[key]['_dimensions'].pop('N_PROF')})
                outStr[key]['_dimensions'].update({'nbpoints':nbpoints})  
        
        #Reform 2D variables from 1D
        ############################
        twoD_flag = []
        for key in outStr.keys() : twoD_flag.append(len(np.shape(outStr[key]['data'])) == 2) if isinstance(outStr[key],dict) else twoD_flag.append(len(np.shape(outStr[key])) == 2)
        twoD_flag = np.array(twoD_flag)
        oneD_vars = np.array(outStr.keys()).compress(~twoD_flag)
        twoD_vars = np.array(outStr.keys()).compress(twoD_flag)
        
        #Transform 2D variables into 1D
        # -> For each record (N_PROF), we only keep the valid depth bin (N_PROF)
        for key in twoD_vars:
            if isinstance(outStr[key],dict) : outStr[key]['data']=outStr[key]['data'][:,outStr[key]['data'].mask.argmin(1)].diagonal()
            else : outStr[key]=outStr[key][:,outStr[key].mask.argmin(1)].diagonal()
 
        #Flatten variables (force flattening in case...)
        for key in outStr.keys():
            if isinstance(outStr[key],dict) : outStr[key]['data']=outStr[key]['data'].flatten()
            else : outStr[key]=outStr[key].flatten()
        
        #Additonnal variables
        dst=outStr['lat'].copy()
        if isinstance(dst,dict) :
            dst['data']=calcul_distance(outStr['lat']['data'],outStr['lon']['data'])
            outStr.update({'dist':dst})
        else : outStr.update({'dist':calcul_distance(outStr['lat'],outStr['lon'])})
                        
        #Update dimensions
        newDim = {'_dimensions' : {'_ndims':1,'nbprofiles':nbprofiles}}
        outStr.update(newDim)
        
        self.update_fid_list(os.path.basename(filename),outStr['_dimensions']['nbprofiles'])
        
        
        return outStr
        
    def read_txt(self,filename):
        
        #Open file
        self._filename = filename
        data=np.genfromtxt(filename,skip_header=1)
        
        #Convert to numpy arrays
        lon=np.ma.masked_array(data[:,3])
        lat=np.ma.masked_array(data[:,2])
        import datetime
        date=np.ma.masked_array(data[:,1])-datetime.date(1951,1,2).toordinal()
        id=np.ma.masked_array(data[:,0])
        
        sz=np.shape(lon)
        ndims=np.size(sz)
        
        #Get SST, U & V data i available
        basename=os.path.basename(filename)
        dirname=os.path.dirname(filename)
        sstuv_name=dirname+'/tempuv_'+basename[basename.find('_')+1:len(basename)]

        data=np.genfromtxt(sstuv_name,skip_header=1)
        temp=np.ma.masked_array(data[:,0])
        temp=np.ma.masked_equal(temp, 9999.)
        u=np.ma.masked_array(data[:,1])
        u=np.ma.masked_equal(u, 9999.)
        v=np.ma.masked_array(data[:,2])
        v=np.ma.masked_equal(v, 9999.)
        
        
        
        return {'_dimensions':{'_ndims':ndims,'N_PROF':sz[0],'N_LEVELS':1},'lon':lon,'lat':lat,'date':date,'id':id,'temp':temp,'u':u,'v':v}

    def read_asc(self,filename,lon_name='LON_FILTRE',lat_name='LAT_FILTRE'):
        
        #Open file
        self._filename = filename
        
        asc=open(self._filename)
        
        #get header properties
        l=0
        par_list=[]
        par_dv=[]
        for line in asc.readlines():
            if not line.startswith('//') : break
            else :
                if line.startswith('//\tClip') :
                    par_list.append(line.split(' ')[-1].replace('\n',''))
                if line.startswith('//\tDefaut') :
                    par_dv.append(line.split('=')[-1].replace('\n',''))
                l+=1
        
        nPar=len(par_list)
        col_id=np.arange(nPar)+2
        
        #Get lon & lat
        par_list[par_list.index(lon_name)]='lon'
        par_list[par_list.index(lat_name)]='lat'
        
        data=np.genfromtxt(filename,skip_header=l)
        
        sz=data.shape[0]    
        mask=np.zeros(sz,dtype=bool)    
        
        #Construct output structure
        dimStr={'_dimensions':{'_ndims':2,'N_PROF':sz,'N_LEVELS':1}}
        
        outStr = OrderedDict()
        outStr.update(dimStr)
        
        outStr.update({'N_PROF':{'_dimensions':{'_ndims':1,'N_PROF':sz},'data':np.ma.array(np.arange(sz),mask=mask.copy())}})
        outStr.update({'N_LEVELS':{'_dimensions':{'_ndims':1,'N_PROF':sz},'data':np.ma.array(np.repeat(25.,sz),mask=mask.copy())}})
        
        #1st row is assigned to float ID, 2nd to date (julian seconds)
        outStr.update({'id':{'_dimensions':{'_ndims':1,'N_PROF':sz},'data':np.ma.masked_array(data[:,0],mask=mask.copy())}})
        outStr.update({'date':{'_dimensions':{'_ndims':1,'N_PROF':sz},'data':np.ma.masked_array(data[:,1]/86400.,mask=mask.copy())}})


        for i,p in enumerate(par_list):
            dumVar=data[:,col_id[i]]
            mask=dumVar==eval('np.{0}'.format(dumVar.dtype))(par_dv[i])
            outStr.update({p:{'_dimensions':{'_ndims':1,'N_PROF':sz},'data':np.ma.masked_array(dumVar,mask=mask.copy())}})
        
        #recompute mask based on lon,lat,date validity
        mask=outStr['lon']['data'].mask.copy() \
            | outStr['lat']['data'].mask.copy() \
            | outStr['date']['data'].mask.copy()
        
        
        #Reapply mask (we may check dimensions?)
        outStr['date']['data'].mask=mask.copy()
        for k in par_list :
            if outStr[k].has_key('data') : outStr[k]['data'].mask=mask.copy()
        
        return outStr
             

#    def plot_track(self,pmap,date):
#        
#        bflag=abs(self.date - date) <= 0.5
#        blon=self.lon.compress(bflag)
#        blat=self.lat.compress(bflag)
#        bdat=self.date.compress(bflag)
#        bid=self.id.compress(bflag)
#        
#        bid_list=np.unique(bid)
#        cnt=np.size(bid_list)
#        
##        #Go to next iteration when no data
##        if cnt == 0 :
##            continue
#        
#        for j in bid_list :
#            dumflag = (bid == j)
#            dumdat=bdat.compress(dumflag)
#            id = dumdat.argsort()
#            dumlon=blon.compress(dumflag)[id]
#            dumlat=blat.compress(dumflag)[id]
#            dumid=bid.compress(dumflag)[id]
#            dumcnt=np.size(dumid)
#            pmap.plot(dumlon,dumlat,'.k',ms=5)
#            pmap.plot(dumlon[dumcnt-1],dumlat[dumcnt-1],'*r',ms=10)
#            pmap.text(dumlon[dumcnt-1],dumlat[dumcnt-1],str(int(j)),fontsize=8,color='b')
##            pmap.title('IMEDIA day '+str(i))
##    #            plt.text(, y, s, fontsize=12)
##            pmap.savefig(FIG_DIR+"/IMEDIA_alti_buoys_"+str(i)+".png")

class argo_trajectory(hydro_data):
    
    def __init__(self,file_pattern,**kwargs):
        
        hydro_data.__init__(self,file_pattern,**kwargs)
        
       
#        self.count=np.size(self.lon)

    def read(self,filename,**kwargs):
        
        fname,extension = os.path.splitext(filename)
        
#        if extension == '.nc' : outStr=self.read_ArgoNC(filename,N_LEVELS=0,**kwargs)
        if extension == '.nc' : outStr=self.read_ArgoNC(filename,**kwargs)
        elif extension == '.dat' : outStr = self.read_txt(filename)
        else : self.Error('Unknown formatting')
        
        
        #Rework dimensions
        ##################
        # Dimensions may be in 2D (more not compliant with current code)
        # Dataobject only accepts flattened (1D) variables
        dimStr=outStr.pop('_dimensions')
#        if dimStr['_ndims'] != 2 : self.Error('Variables with {0} dimensions can not be used as glider profiles (2 dimensions)'.format(dimStr['_ndims']))
        if not dimStr.has_key('N_PROF') :
            self.Error('Data structure must contain dimensions N_LEVELS and N_PROF (only contains {0})'.format(dimStr.keys()[1:]))
        
        datalen = [np.size(outStr[key]) for key in outStr.keys()] #[np.shape(outStr[key]) for key in outStr.keys()] 
                
        nbprofiles = dimStr['N_PROF']
        nbpoints = nbprofiles
        
        #Reform 2D variables from 1D
        ############################
        twoD_flag = np.array([len(np.shape(outStr[key])) for key in outStr.keys()]) == 2 
        oneD_vars = np.array(outStr.keys()).compress(~twoD_flag)
        twoD_vars = np.array(outStr.keys()).compress(twoD_flag)
        
        #Transform 2D variables into 1D
        # -> For each record (N_PROF), we only keep the valid depth bin (N_PROF)
        for key in twoD_vars: outStr.update({key:outStr[key][:,outStr[key].mask.argmin(1)].diagonal()})
        
        #Flatten variables (force flattening in case...)
        for key in outStr.keys(): outStr.update({key:outStr[key].flatten()})
        
        #Additonnal variables
        outStr.update({'dist':calcul_distance(outStr['lat'],outStr['lon'])})
                
        #Update dimensions
        newDim = {'_dimensions' : {'_ndims':1,'nbpoints':nbprofiles}}
        outStr.update(newDim)
        
        self.update_fid_list(os.path.basename(filename),outStr['_dimensions']['nbpoints'])
        
        return outStr

class glider_data(hydro_data):
    
    def __init__(self,file_pattern,**kwargs):
        
        #Init hydro_data class (calling gilder_data.read() function)
        hydro_data.__init__(self,file_pattern,**kwargs)

    def read(self,filename,**kwargs):
        
        fname,extension = os.path.splitext(filename)
        outStr=self.read_ArgoNC(filename,**kwargs)
        
        #Rework dimensions
        ##################
        # Dimensions may be in 2D (more not compliant with current code)
        # Dataobject only accepts flattened (1D) variables
        dimStr=outStr.pop('_dimensions')
        if dimStr['_ndims'] != 2 : self.Error('Variables with {0} dimensions can not be used as glider profiles (2 dimensions)'.format(dimStr['_ndims']))
        if not dimStr.has_key('N_LEVELS') or not dimStr.has_key('N_PROF') :
            self.Error('Data structure must contain dimensions N_LEVELS and N_PROF (only contains {0})'.format(dimStr.keys()[1:]))
        
#        datalen = [np.size(outStr[key]) for key in outStr.keys()] #[np.shape(outStr[key]) for key in outStr.keys()] 
                
        nblevels = dimStr['N_LEVELS']
#        nblevels = nblevels.astype(nblevels.dtype)
        nbprofiles = dimStr['N_PROF']
        nbpoints = nblevels*nbprofiles
        
        #Change dimension names towards new ones
        for key in outStr.keys() :
            if outStr[key].has_key('_dimensions'):
                outStr[key]['_dimensions'].pop('N_LEVELS',None)
                outStr[key]['_dimensions'].pop('N_PROF',None)
#                if outStr[key]['_dimensions'].has_key('N_PROF') : outStr[key]['_dimensions'].update({'nbprofiles':outStr[key]['_dimensions'].pop('N_PROF')})
                outStr[key]['_dimensions'].update({'nbpoints':nbpoints})  
        
        #Reform 2D variables from 1D
        ############################
        twoD_flag = []
        for key in outStr.keys() : twoD_flag.append(len(np.shape(outStr[key]['data'])) == 2) if isinstance(outStr[key],dict) else twoD_flag.append(len(np.shape(outStr[key])) == 2)
        twoD_flag = np.array(twoD_flag)
        oneD_vars = np.array(outStr.keys()).compress(~twoD_flag)
        twoD_vars = np.array(outStr.keys()).compress(twoD_flag)
        
        #Transform 1D variables into 2D
        for key in oneD_vars:
            if isinstance(outStr[key],dict) : outStr[key]['data'] = np.reshape(np.repeat(outStr[key]['data'],nblevels),[nbprofiles,nblevels])
            else :  outStr[key]=np.reshape(np.repeat(outStr[key],nblevels),[nbprofiles,nblevels])
        
#        [{key:np.shape(outStr[key])} for key in outStr.keys()] #Check variables
        
        #Load surface data
        surfStr=self.read_ArgoNC(filename,N_LEVELS=0) #read surface data
        surfStr.pop('_dimensions')
        
        #Change dimension names towards new ones
        for key in surfStr.keys() :
            if surfStr[key].has_key('_dimensions'):
                surfStr[key]['_dimensions'].pop('N_LEVELS',None)
                if surfStr[key]['_dimensions'].has_key('N_PROF') : surfStr[key]['_dimensions'].update({'nbprofiles':surfStr[key]['_dimensions'].pop('N_PROF')})
                        
        #Flatten surface variables
        for key in surfStr.keys():
            if isinstance(surfStr[key],dict) : surfStr[key]['data'] = surfStr[key]['data'].flatten()
            else : surfStr[key] = surfStr[key].flatten()
            outStr.update({key+'_surf':surfStr[key]})
        
        #Flatten 2D variables (
        for key in outStr.keys():
            if isinstance(outStr[key],dict) : outStr[key]['data']=outStr[key]['data'].flatten()
            else : outStr[key]=outStr[key].flatten()
        
        #Additonnal variables
        dst_surf=outStr['lat_surf'].copy()
        dst=outStr['lat'].copy()
        if isinstance(dst_surf,dict) :
            dst_surf['data']=calcul_distance(outStr['lat_surf']['data'],outStr['lon_surf']['data'])
            dst['data']=np.repeat(dst_surf['data'],nblevels)
            outStr.update({'dist_surf':dst_surf})
            outStr.update({'dist':dst})
        else :
            outStr.update({'dist_surf':calcul_distance(outStr['lat_surf'],outStr['lon_surf'])})
            outStr.update({'dist':np.repeat(outStr['dist_surf'],nblevels)})

        
        ###TODO : Could be simplified? (if has_key(var) --> must have key var_surf
        if outStr.has_key('pres') : 
            deph = outStr['pres'].copy()
            deph_surf = outStr['pres_surf'].copy()
            try:
                deph['data']=gsw.z_from_p(outStr['pres']['data'],outStr['lat']['data'])
                deph_surf['data']=gsw.z_from_p(outStr['pres_surf']['data'],outStr['lat_surf']['data'])
            except :
                deph['data'][:]=deph['data'].fill_value
                deph_surf['data'][:]=deph['data'].fill_value
            outStr.update({'deph' : deph})
            outStr.update({'deph_surf' : deph_surf})
        if outStr.has_key('psal') and outStr.has_key('temp') and outStr.has_key('pres') :
            rho=outStr['psal'].copy()
            rho_surf=outStr['psal_surf'].copy()
            try:
                rho['data']=gsw.rho(outStr['psal']['data'],outStr['temp']['data'],outStr['pres']['data'])
                rho_surf['data']=gsw.rho(outStr['psal_surf']['data'],outStr['temp_surf']['data'],outStr['pres_surf']['data'])
            except :
                rho['data'][:]=deph['data'].fill_value
                rho_surf['data'][:]=deph['data'].fill_value
            outStr.update({'rho' : rho})
            outStr.update({'rho_surf' : rho_surf})
        
        #Update dimensions
        newDim = {'_dimensions' : {'_ndims':2,'nbpoints':nbpoints,'nbprofiles':nbprofiles}}
        outStr.update(newDim)
        
        #Update fid
        self.update_fid_list(os.path.basename(filename),outStr['_dimensions']['nbpoints'])
        
        return outStr
       
    
class TSG_data(hydro_data) :
    def __init__(self,file_pattern,**kwargs):
        
        hydro_data.__init__(self,file_pattern,**kwargs)
    
    def read(self,filename):
        
        fname,extension = os.path.splitext(filename)
        
        if extension == '.nc' :
            outStr=self.read_ArgoNC(filename,params=['TEM2','PSAL'])
            outStr.update({'temp':outStr.pop('tem2')}) #Use TEM2 as temperature field
#            outStr.update({'depth':outStr.pop('deph')}) #Use DEPH as depth field
            return outStr
        elif extension == '.dat' :
            return self.read_txt(filename)
        else : self.Error('Unknown formatting')
        
    def read_txt(self,filename):
         #Open file
        self._filename = filename
        data=np.genfromtxt(filename,skip_header=1)
        
        #Convert to numpy arrays
        lon=np.ma.masked_array(data[:,2])
        lat=np.ma.masked_array(data[:,1])
        
        import datetime
        date=np.ma.masked_array(data[:,0])+datetime.date(2012,1,1).toordinal()-datetime.date(1950,1,2).toordinal()
        
        temp=np.ma.masked_array(data[:,3])
        temp=np.ma.masked_greater(temp, 99.)
        psal=np.ma.masked_array(data[:,4])
        psal=np.ma.masked_less(psal, 35.)
        fluo=np.ma.masked_array(data[:,5])
        fluo=np.ma.masked_where((fluo == 0.0) | (fluo >= 99.),fluo)
        
        id=np.repeat('IMEDIA_TSG',len(psal))
        
        
        sz=np.shape(lon)
        ndims=np.size(sz)
                
        return {'_dimensions':{'_ndims':ndims,'nbpoints':sz[0]},'lon':lon,'lat':lat,'date':date,'id':id,'temp':temp,'psal':psal,'fluo':fluo}

class CTD_data(hydro_data):
    
    def __init__(self,file_pattern,**kwargs):
        
        #Init hydro_data class (calling gilder_data.read() function)
        hydro_data.__init__(self,file_pattern,**kwargs)

    def read(self,filename,**kwargs):
        
        fname,extension = os.path.splitext(filename)
        
        if extension == '.nc' :
            outStr=self.read_ArgoNC(filename,params=['TEM2','PSAL'])
            outStr.update({'temp':outStr.pop('tem2')}) #Use TEM2 as temperature field
#            outStr.update({'depth':outStr.pop('deph')}) #Use DEPH as depth field
            return outStr
        elif extension == '.dat' :
            return self.read_txt(filename)
        
        #Seabird CTD data
        elif extension == '.cnv' :
            return self.read_cnv(filename,**kwargs)
        elif extension == '.asc' :
            return self.read_asc(filename,**kwargs)
        else : self.Error('Unknown formatting')
    
    def read_asc(self,filename,**kwargs):
        
        self._filename = filename
        
        #Check file length
        nlines=0
        for line in open(filename): nlines+=1            
        
        if nlines > 1 :
            #Open file
            data=np.genfromtxt(filename,skip_header=1)
            
            #Convert to numpy arrays
            lon=kwargs['lon']
            lat=kwargs['lat']
            
            import datetime
            date=np.ma.masked_array(data[:,0])+datetime.date(2012,1,1).toordinal()-datetime.date(1950,1,2).toordinal()
            pres=np.ma.masked_array(data[:,1])
            temp=np.ma.masked_array(data[:,2])
            cond=np.ma.masked_array(data[:,3])
            obs1=np.ma.masked_array(data[:,4])
            obs2=np.ma.masked_array(data[:,5])
            descend_rate=np.ma.masked_array(data[:,6])
            scan=np.ma.masked_array(data[:,7])
            fluoro=np.ma.masked_array(data[:,8])
            depth=np.ma.masked_array(data[:,9])
            potemp=np.ma.masked_array(data[:,10])
            psal=np.ma.masked_array(data[:,11])
            dens=np.ma.masked_array(data[:,12])
            svCM=np.ma.masked_array(data[:,13])
            flag=np.ma.masked_array(data[:,14])
            
            reclen=len(pres)
            
            date.mask=np.zeros(reclen,dtype='bool')
            pres.mask=np.zeros(reclen,dtype='bool')
            temp.mask=np.zeros(reclen,dtype='bool')
            cond.mask=np.zeros(reclen,dtype='bool')
            obs1.mask=np.zeros(reclen,dtype='bool')
            obs2.mask=np.zeros(reclen,dtype='bool')
            descend_rate.mask=np.zeros(reclen,dtype='bool')
            scan.mask=np.zeros(reclen,dtype='bool')
            fluoro.mask=np.zeros(reclen,dtype='bool')
            depth.mask=np.zeros(reclen,dtype='bool')
            potemp.mask=np.zeros(reclen,dtype='bool')
            psal.mask=np.zeros(reclen,dtype='bool')
            dens.mask=np.zeros(reclen,dtype='bool')
            svCM.mask=np.zeros(reclen,dtype='bool')
            flag.mask=np.zeros(reclen,dtype='bool')
            
            
            
            id=np.repeat('{0}'.format(kwargs['stationid']),reclen)
            
            lon=np.ma.masked_array(np.repeat(lon,reclen),mask=np.zeros(reclen,dtype='bool'))
            lat=np.ma.masked_array(np.repeat(lat,reclen),mask=np.zeros(reclen,dtype='bool'))
            
    #        psal=csw.salt(cond/gsw.cte.C3515, temp, pres)
    #        depth=gsw.z_from_p(pres,lat)
    #        dens= gsw.rho(psal,temp,pres)
            
            
            sz=np.shape(lon)
            ndims=np.size(sz)
        else :
            ndims = 1.
            sz=(1,1)
            lon=np.array(np.NaN)
            lat=np.array(np.NaN)
            date=np.array(np.NaN)
            id=np.array(np.NaN)
            depth=np.array(np.NaN)
            pres=np.array(np.NaN)
            temp=np.array(np.NaN)
            psal=np.array(np.NaN)
            fluoro=np.array(np.NaN)
            dens=np.array(np.NaN)
            potemp=np.array(np.NaN)
            cond=np.array(np.NaN)
            obs1=np.array(np.NaN)
            obs2=np.array(np.NaN)
            svCM=np.array(np.NaN)
            descend_rate=np.array(np.NaN)
            flag=np.array(np.NaN)
            
        return {'_dimensions':{'_ndims':ndims,'nbpoints':sz[0]},'lon':lon,'lat':lat,'date':date,'id':id,'depth':depth,'pres':pres,'temp':temp, \
                'psal':psal,'fluoro':fluoro,'dens':dens,'potemp':potemp,'cond':cond,'obs1':obs1,'obs2':obs2,'svCM':svCM,'descend_rate':descend_rate,'flag':flag}    

def read_cnv(self,filename,**kwargs):
        #Open file
        data=np.genfromtxt(filename,skip_header=328)
        
        #Convert to numpy arrays
        lon=kwargs['lon']
        lat=kwargs['lat']
        
        import datetime
        date=np.ma.masked_array(data[:,0])+datetime.date(2012,1,1).toordinal()-datetime.date(1950,1,2).toordinal()
        
        pres=np.ma.masked_array(data[:,1])
        temp=np.ma.masked_array(data[:,2])
        cond=np.ma.masked_array(data[:,3])
        obs1=np.ma.masked_array(data[:,6])
        obs2=np.ma.masked_array(data[:,7])
        fluoro=np.ma.masked_array(data[:,])
        
        pres=np.ma.masked_greater(pres, 99999.)
        temp=np.ma.masked_greater(temp, 99.)
        cond=np.ma.masked_greater(cond, 99.)
        #fluo=np.ma.masked_where((fluo == 0.0) | (fluo >= 99.),fluo)
        
        reclen=len(pres)
        
        id=np.repeat('{0}'.format(kwargs['stationid']),reclen)
        
        lon=np.ma.masked_array(np.repeat(lon,reclen))
        lat=np.ma.masked_array(np.repeat(lat,reclen))
        
        try :
            psal=csw.salt(cond/gsw.cte.C3515, temp, pres)
            depth=gsw.z_from_p(pres,lat)
            dens= gsw.rho(psal,temp,pres)
        except :
            psal=np.ma.masked_array(np.repeat(lon.fill_value,reclen),mask=True)
            depth=np.ma.masked_array(np.repeat(lon.fill_value,reclen),mask=True)
            dens=np.ma.masked_array(np.repeat(lon.fill_value,reclen),mask=True)
        
        
        sz=np.shape(lon)
        ndims=np.size(sz)
                
        return {'_dimensions':{'_ndims':ndims,'nbpoints':sz[0]},'lon':lon,'lat':lat,'date':date,'id':id,'depth':depth,'pres':pres,'temp':temp,'psal':psal}  
   