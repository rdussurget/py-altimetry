# outer __init__.py
'''
altimetry module
@summary: Contains several functions and routines to deal with altimetry, remote sensing and in-situ data.
@author: Renaud DUSSURGET, LER/PAC IFREMER.
@change: Create in December 2012 by RD.
'''

#Setup package defaults
from config import defaults as def_class
defaults = def_class()

from tools import *
from externals import *
from data import *

#from external-tools import *