'''
CPT_READER : Read cpt palette and returns a segmented color dictionary for use in matplotlib

@author: David Huard (david.huard-Re5JQEeQqe8AvxtiuMwx3w@public.gmane.org)
         Renaud Dussurget, LER PAC/IFREMER (RD, renaud.dussurget@ifremer.fr)
@change: Created by DH in February 2006
    ported into PyValid on 21 janv. 2013 (RD)
'''
from scipy import zeros, linspace, shape, float32 as Float, concatenate
import numpy as np
import matplotlib.colors
import glob, os
from altimetry import defaults


def cpt2seg(file_name, sym=False, discrete=False):
    """Reads a .cpt palette and returns a segmented colormap.

    sym : If True, the returned colormap contains the palette and a mirrored copy.
    For example, a blue-red-green palette would return a blue-red-green-green-red-blue colormap.

    discrete : If true, the returned colormap has a fixed number of uniform colors.
    That is, colors are not interpolated to form a continuous range.

    Example :
    >>> _palette_data = cpt2seg('palette.cpt')
    >>> palette = matplotlib.colors.LinearSegmentedColormap('palette', _palette_data, 100)
    >>> imshow(X, cmap=palette)
    """
   
   
    dic = {}
#    io
#    f = scipy.io.open(file_name, 'r')
#    rgb = f.read_array(f)
    
    #Check flags:
#    with open(file_name) as f:
#        content = f.readlines()
#        header=np.where(np.array([c.startswith('#') for c in content]))[0][-1]
#        footer=np.where(np.array([c[0].isalpha() for c in content]))[0][0]
    
    rgb = np.genfromtxt(file_name,comments='#',invalid_raise=False)
#    rgb = np.genfromtxt(file_name)
    rgb = rgb/255.
    s = shape(rgb)
    colors = ['red', 'green', 'blue']
    for c in colors:
        i = colors.index(c)
        x = rgb[:, i+1]

        if discrete:
            if sym:
                dic[c] = zeros((2*s[0]+1, 3), dtype=Float)
                dic[c][:,0] = linspace(0,1,2*s[0]+1)
                vec = concatenate((x ,x[::-1]))
            else:
                dic[c] = zeros((s[0]+1, 3), dtype=Float)
                dic[c][:,0] = linspace(0,1,s[0]+1)
                vec = x
            dic[c][1:, 1] = vec
            dic[c][:-1,2] = vec
               
        else:
            if sym:
                dic[c] = zeros((2*s[0], 3), dtype=Float)
                dic[c][:,0] = linspace(0,1,2*s[0])
                vec = concatenate((x ,x[::-1]))
            else:
                dic[c] = zeros((s[0], 3), dtype=Float)
                dic[c][:,0] = linspace(0,1,s[0])
                vec = x
            dic[c][:, 1] = vec
            dic[c][:, 2] = vec
   
    return dic

def revert_cpt(cptdata):
    outStr=cptdata.copy()
    for k in outStr.keys() : outStr[k][:,1:]=outStr[k][::-1,1:]
    
    return outStr

def get_cmap(cmap_name,revert=False,N=256):
    
    cpt_dir=defaults.cptDir
    
    ls = glob.glob(cpt_dir + cmap_name + '.cpt')
    if len(ls) == 0 : raise Exception('file {0} do not exists'.format(os.path.basename(cpt_dir + cmap_name + '.cpt')))
    cptdata=cpt2seg(ls[0])
    if revert : cptdata=revert_cpt(cptdata)
    cmap = matplotlib.colors.LinearSegmentedColormap ('palette', cptdata, N)
    
    return cmap
