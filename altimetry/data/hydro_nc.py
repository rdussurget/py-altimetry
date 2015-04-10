import os
import glob
import inspect
import fnmatch

import numpy as np

import matplotlib.pyplot as plt
from netCDF4 import Dataset as ncfile
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


#Additional functions
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
                ind_list.append(xrange(varDimval[vid]))
                dims.update({enum[1]:varDimval[enum[0]]})
            
            #Data is indexed along current dimension
            else :
                dumind = kwargs[vn]
#                if not isinstance(dumind,list) : dumind=tuple(dumind)
                if isinstance(dumind,np.ndarray) : dumind=dumind.tolist() #Rq: tolist() can take a very long time to run on large arrays
                if type(dumind) is not list : dumind=[dumind] 
                ind_list.append(dumind)
#                ind_list=(ind_list,dumind)  if len(ind_list) != 0 else (dumind,)
                dims.update({vid:len(dumind)})
        
        #check index list
        sz=[len(i) for i in ind_list]
        
        #find empty dimensions
        if not (where_list([0],sz)[0] == -1 ) : varOut=var[[0]][[]] #np.array(where_list(sz,[0])) == -1
        else :
            varOut=var[ind_list]#.copy() #THIS IS LONG!!
            if var.shape == (1,1) : varOut=varOut.reshape(var.shape)
        
        #Mask it!
        if var.__dict__.has_key('_FillValue') : mask=varOut == var._FillValue
        elif var.__dict__.has_key('missing_value') : mask=varOut == var.missing_value
        else : mask=np.zeros(varOut.shape,dtype='bool')
        
        #Scale it
        #note : we do not use the *= or += operators to force casting to scaling attribute types
        if var.__dict__.has_key('scale') : varOut =varOut * var.scale
        elif var.__dict__.has_key('scale_factor') : varOut = varOut * var.scale_factor
        if var.__dict__.has_key('add_offset') : varOut = varOut + var.add_offset
        
        #Set masks properly
        if isinstance(varOut,np.ndarray) : varOut=np.ma.masked_array(varOut,mask=mask)
        elif isinstance(varOut,np.ma.masked_array) : var.mask=mask
        else : raise 'This data type {} has not been defined - code it!'.format(type(varOut))
        
        #Get attributes
        attrStr=var.__dict__
        attrStr.pop('_FillValue',None) #Remove this attributed as it is overidden
        
        #Append attributes to varOut
        varOut.__dict__.update(attrStr)
        
        #Build up output structure
        outStr={'_dimensions':dims,'data':varOut}
        dims.update({'_ndims':len(dims.keys()[1:])})
        
        return outStr
#            ind_list=[[]] 
    
