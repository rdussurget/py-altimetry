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

def load_ncVar(varName, nc=None, **kwargs):
        
        if (nc is None) : raise Exception('No Netcdf file passed')
        
        #Load variable
        var = nc.variables[varName]
        
        var.set_auto_maskandscale(False)
        
        #Load dimensions
        varDim = [str(dim) for dim in var.dimensions]
        missDim=len(varDim) == 0
        if (missDim): warn('No dimension found')
        else : varDimval = [len(nc.dimensions[dimname]) for dimname in varDim]
        
        #Load Attributes
        attrStr=var.__dict__
        
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
                if fg.sum() == 0 :
                    raise IndexError('{0} {1} is not matching given dimensions {2}'.format(vn,(np.nanmin(nc.variables[vn][:]),np.nanmax(nc.variables[vn][:])),drange))
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
        
        #check index list
        sz=[len(i) for i in ind_list]
        
        dstr=','.join(dstr) #invert dimension list for numpy
#        dstr=','.join(dstr[::-1]) #invert dimension list for numpy
        if missDim : cmd = 'varOut = var[:]'
        else : cmd = 'varOut = var[{0}]'.format(dstr)
        exec(cmd)
        
        #find empty variables
#        if not (atools.where_list([0], shape)[0] == -1) : varOut = var[[0]][[]]
#        else : varOut = var[ind_list]
        
        #Mask it!
        if var.__dict__.has_key('_FillValue') :
            fill_value=var._FillValue
            mask = varOut == var._FillValue
        elif var.__dict__.has_key('missing_value') :
            fill_value=var._FillValue
            mask = varOut == var._FillValue
        else :
            fill_value = None
            mask = np.zeros(varOut.shape, dtype='bool')
        
        #Scale it
        #note : we do not use the *= or += operators to force casting to scaling attribute types
        if var.__dict__.has_key('scale') : varOut = varOut * var.scale
        elif var.__dict__.has_key('scale_factor') : varOut = varOut * var.scale_factor
        if var.__dict__.has_key('add_offset') : varOut = varOut + var.add_offset
        
        #Set masks properly
        if isinstance(varOut, np.ndarray) : varOut = np.ma.masked_array(varOut, mask=mask,dtype=varOut.dtype,fill_value=fill_value)
        elif isinstance(varOut, np.ma.masked_array) : var.mask = mask
        else : raise 'This data type {} has not been defined - code it!'.format(type(varOut))
        
        #Update masked data properly
        varOut.data[varOut.mask]=varOut.fill_value
        
        #Switch dimensions 
        if not missDim : varOut=np.transpose(varOut,tuple(range(len(dims.keys()[1:]))[::-1]))
        
        #Build up output structure
        dims.update({'_ndims':len(dims.keys()[1:])})
        outStr = {'_dimensions':dims, 'data':varOut}
        
        #Add variable attributes
        for A in var.__dict__.keys():
            outStr[A]=var.getncattr(A)
        
        return outStr
#            ind_list=[[]] 
    
