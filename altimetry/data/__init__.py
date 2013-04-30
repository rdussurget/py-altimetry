# outer __init__.py
'''
altimetry.data module
@summary: The following classes and associated functions are involved with reading and processing altimetry, remote sensing and in-situ data.
@author: Renaud DUSSURGET, LER/PAC IFREMER.
@change: Create in December 2012 by RD.
'''
import altimetry.tools
from hydro import *
from alti_data import *

try:
 from rs_data import *
except ImportError, e:
 pass # module doesn't exist, deal with it.
