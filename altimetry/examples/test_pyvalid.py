#import altimetry
import numpy as np
from altimetry.tools import bytscl, map_tools
from altimetry.data import alti_data
from altimetry.tools.nctools import nc
import altimetry

print bytscl(np.arange(20.0))
obj=nc()

pmap=map_tools.plot_map(limit=[43,5,44,8],resolution='f')
#pmap.shadedrelief()
#pmap.bluemarble()
pmap.show()

print 'done'