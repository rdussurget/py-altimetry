#! /usr/bin/env python
# -*- coding: utf-8 -*- 
# ________________________________________________________________
#
# S. Theetten
#
# Created: 01/2011
# Last update: 
# ________________________________________________________________


from datetime import datetime, date, time
from vacumm.misc.atime import add

# find dates for diag

# today date
end_date = datetime.now()

# starting date (i.e. 7 days before)
start_date = add(end_date,-7,'days')

# update config.cfg

from ConfigParser import SafeConfigParser
from vacumm.misc.atime import strftime

config = SafeConfigParser()
config.read('config.cfg')

config.set('Time Period', 'andeb',strftime('%Y',start_date))
config.set('Time Period', 'mdeb',strftime('%m',start_date))
config.set('Time Period', 'jdeb',strftime('%d',start_date))

config.set('Time Period', 'anfin',strftime('%Y',end_date))
config.set('Time Period', 'mfin',strftime('%m',end_date))
config.set('Time Period', 'jfin',strftime('%d',end_date))

# On sauvegarde
fc = open('config.cfg', 'w')
config.write(fc)
fc.close()

