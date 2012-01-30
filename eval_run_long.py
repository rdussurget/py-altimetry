#! /usr/bin/env python
# -*- coding: utf-8 -*- 
# ________________________________________________________________
#
# G. Charria, S. Theetten, and J. Gatti 
#
# Created: 02/2011
# Last update: 
# ________________________________________________________________
#
import os,sys,shutil,ConfigParser, gc
# -- These two line allows running without display
# Warning: does not allow screen figure display
import matplotlib
matplotlib.use('Agg')
# --
import cdms2



import matplotlib.pyplot as P
import cdtime,  MV2
from pylab import figure,  show,  close, title, ylabel, ylim
from matplotlib.dates import datestr2num
import numpy as np

# Import de module specifiques a Vacumm
# Import de module specifiques a Vacumm
#from vacumm.data.satellite.old_nar import *
from vacumm.data.satellite.nar import Nar
from vacumm.data.satellite.seviri import Seviri
from vacumm.data.model.mars.get_ftp import  get_ftp_f1
from vacumm.data.model.mars.get_cp import  get_cp_f1
from vacumm.data.misc.handle_directories import make_directories
from vacumm.validator.valid.ValidXYT import ValidXYT
from vacumm.misc.plot import map,  curve, curve2,  savefigs
from vacumm.misc.atime import * # a detailler: are_same_units + ...
from vacumm.misc.axes import create_time,  set_order
from vacumm.misc.grid.regridding import regrid1d, regrid2d
from vacumm.misc.io import ncread_best_estimate
#from vacumm.misc.grid.regridding import regrid1d,  regrid2d,  interp2d
import glob
from vacumm.markup import html_tools
from global_stat import allstat, monthlystat, seasonalstat, regionalstat

# Permet la generation des Bounds automatiquement ... indispensable pour cdutil !
cdms2.setAutoBounds('on')


print 65*'-'
print ' Validation de la SST '
print 65*'-'

# -- Periode
#ctdeb=cdtime.comptime(1996,1,10,0,0,0)
#ctfin=cdtime.comptime(1996,2,10,0,0,0)
ctdeb=cdtime.comptime(2004,6,1,0,0,0)
ctfin=cdtime.comptime(2006,6,1,0,0,0)
#ctdeb=cdtime.comptime(2004,2,1,0,0,0)
#ctfin=cdtime.comptime(2005,2,1,0,0,0)

# Pour test car pas meme periode modele et sat ...
#ctdeb1=cdtime.comptime(2005,1,10,0,0,0)
#ctfin1=cdtime.comptime(2005,2,10,0,0,0)
ctdeb1=ctdeb
ctfin1=ctfin

# Level max !
#levmax = 30
levmax = 1
levmaxm1 = levmax-1


# -- Obs: Rapatriement et lecture des observations corespondantes a la periode consideree
# ----------------------------------------------------------------
print ' -- Lecture des observations satellites --    '
print 65*'-'
dir_obs = '/home11/caparmor/mars/VALID/DATA/MANGA/SST/SST_SEVIRI/TRAITEMENT/manga_extended/'
filename_obs = 'SEVIRI_MANGA_2004-2008_daily.nc'
obs = ncread_best_estimate('temp',os.path.join(dir_obs,filename_obs), (ctdeb1, ctfin1))



# -- Model : Rapatriement et lecture des sorties de modele corespondantes a la periode consideree
# ----------------------------------------------------------------
# ----------------------------------------------------------------
print ' -- Lecture des sorties de modeles --    '
print 65*'-'
#dir_model = '/work/theetten/MARS/RUN_MARS/MANGA_4000/MANGA_4000-03/rank_1/'
#filename_model = 'manga_4000_1_03.nc.all'

# Jerlov eau tres claire
#dir_model = '/work/theetten/MARS/RUN_MARS/MANGA_4000/MANGA_4000-JerlovTC/rank_1/'
#filename_model = 'sst_manga_4000_JerlovTC.nc.all'
#dir_fig = 'Fig_jerlovETC'
#dir_fig = 'Fig_jerlovETC2004'

# Jerlov eau tres trouble
dir_model = '/work/theetten/MARS/RUN_MARS/MANGA_4000/MANGA_4000-JerlovTB/rank_1/'
filename_model = 'sst_manga_4000_JerlovTB.nc.all'
dir_fig = 'Fig_jerlovETB'
#dir_fig = 'Fig_jerlovETB2004'


# Jerlov eau claire
#dir_model = '/work/theetten/MARS/RUN_MARS/MANGA_4000/MANGA_4000-Jerlov/rank_1/'
#filename_model = 'sst_manga_4000_Jerlov_Jerlov.nc.all'
#dir_fig = 'Fig_jerlovEC'
#dir_fig = 'Fig_jerlovEC2004'


# No Jerlov
#dir_model = '/work/theetten/MARS/RUN_MARS/MANGA_4000/MANGA_4000-No_Jerlov/rank_1/'
#filename_model = 'sst_manga_4000_No_Jerlov_No_Jerlov.nc.all'
#dir_fig = 'Fig_nojerlov'
#dir_fig = 'Fig_nojerlov2004'



model = ncread_best_estimate('TEMP',os.path.join(dir_model,filename_model), (ctdeb, ctfin),select=dict(level=slice(levmaxm1,levmax)))



# Creation des repertoires de resultats: Figures
if os.path.isdir(dir_fig)==False:
        os.mkdir(dir_fig)
SCRIPT_DIR = os.getcwd()
FIG_DIR = os.path.join(SCRIPT_DIR,dir_fig)


# Ouverture fichier config
config = ConfigParser.RawConfigParser()
config.read(os.path.join(SCRIPT_DIR,'config.cfg'))


# -----------------------------------------------------------------

if are_same_units(model.getTime().units,obs.getTime().units)==False:
    print '-- Conversion des unites de temps --'
    print 65*'-'
    ch_units(model, obs.getTime().units)

print '-- Interpolation temporelle --'
print 65*'-'
time_ref = model.getTime()
#modelregridontime = regrid1d(model, time_ref, 'linear', axis=0)
obsregridontime = regrid1d(obs, time_ref, 'remap', axis=0)

if config.get('Control', 'tseries') == 'True':
    print '-- Trace de controle (1) --'
    print 65*'-'
    x1 = datestr2num(strtime(obs.getTime()))
    x2 = datestr2num(strtime(obsregridontime.getTime()))
    itest = 40
    fig = figure()
    ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
    lo = obsregridontime.getLongitude().getValue()
    la = obsregridontime.getLatitude().getValue()
    title('%(#)3.2fE / %(##)4.2fN - Control temporal interpolation' %{'#':lo[itest], '##':la[itest]})
    ylabel('Temperature ($^{\circ}C$)')
    ax.plot_date(x1,obs.data[:,itest,itest],linestyle='-', marker='o', color='r')
    ax.plot_date(x2,obsregridontime.data[:,itest,itest],'bs')
    ax.legend(('Obs Before intepolation',  'Obs after interpolation'),  loc = 'upper left',  shadow = True,  fancybox=True)
    fig.autofmt_xdate()
    savefigs(FIG_DIR+'/control_interp_time')
    close()

print '-- Interpolation spatiale --'
print 65*'-'

cgrido = model.getGrid()
# !!!! Essayer de faire fonctionner 'remap' !!!!
obsregridspatial = regrid2d(obsregridontime,  cgrido,  method = 'nearest')

del obsregridontime

if config.get('Control', 'interp_map') == 'True':
    print '-- Trace de controle (2) --'
    print 65*'-'
    # Trace des resultats
    kwplot = dict(show=False, colorbar=True, colorbar_horizontal = True,  vmin=obsregridspatial.min(), vmax=obsregridspatial.max(),
        drawparallels_size=8, drawmeridians_size=8, drawmeridians_rotation=45.,  clabel_hide=True)
    if float(np.ma.count(obsregridspatial[-1,:,:])) / float(np.ma.count_masked(obsregridspatial[-1,:,:])) < 0.1:
        tt = 'Obs - Original\n (< 10% of non-masked values)'
    else:
        tt = 'Obs - Original\n (last time step)'
    m=map(obs[-1, :, :], title=tt, subplot=223, **kwplot)
    #m = map(obs.data[0, :, :], title='Obs - Original', subplot=223, **kwplot)
    #add_grid(cgridi, lw=.7, alpha=.3)
    if float(np.ma.count(obsregridspatial[-1,:,:])) / float(np.ma.count_masked(obsregridspatial[-1,:,:])) < 0.1:
        tt = 'Obs - after spatial interpolation\n (< 10% of non-masked values)'
    else:
        tt = 'Obs - after spatial interpolation\n (last time step)'
    map(obsregridspatial[-1, :, :], title=tt,  subplot=224,yhide=1, **kwplot)
    #add_grid(cgridi, lw=.7, alpha=.3)
    savefigs(FIG_DIR+'/control_interp_map')
    P.close()

del obs

print '-- Conservation du masque le plus restrictif  --'
print 65*'-'
model = model[:, 0, :, :](squeeze=1)
m3 = np.ma.mask_or(obsregridspatial.mask, model.mask)
model.mask = m3
obsregridspatial.mask = m3
#modelregridontime.mask = np.ma.getmask(obsregridspatial)

if config.get('Statistics', 'to_do') == 'True':
    allstat(model, obsregridspatial,FIG_DIR, SCRIPT_DIR)
    seasonalstat(model, obsregridspatial,FIG_DIR, SCRIPT_DIR)
    #monthlystat(model, obsregridspatial, FIG_DIR, SCRIPT_DIR)
    regionalstat(model, obsregridspatial, FIG_DIR, SCRIPT_DIR)
    print 65*'-'
    

if config.get('Report', 'rep_html') == 'True':
    print ' -- Generation du Rapport de Validation --    ' 
    print 65*'-'
    #title = "Test(1)"
    #header = ""
    #footer = ""
    
    #ttforintro = '<u>Observations</u>: %(#)s <br /> <u>Model</u>: %(##)s %(###)s <br /><hr>' %{'#':result.obs.long_name, '##':config.get('Model Description', 'name') , '###':result.model.long_name}
    ttforintro = '<u>File</u>: %(#)s <br /> <hr>' %{'#':filename_model}
    
    
    #if config.get('Statistics', 'extrema') == 'True' and config.get('Statistics', 'to_do') == 'True':
    #    mimamodel = '=> Modelled minimum and maximum: %(#)6.3f / %(##)6.3f %(u)s' %{'#':result.model.extrema[0], '##':result.model.extrema[1], 'u':result.model.units}
    #    mimaobs = '=> Observed minimum and maximum: %(###)6.3f / %(####)6.3f %(u)s'  %{'###':result.obs.extrema[0],'####':result.obs.extrema[1], 'u':result.obs.units}
    #else:
    #    mimamodel = ''
    #    mimaobs = ''
    mimamodel = ''
    mimaobs = ''
            
    images_control = glob.glob(os.path.join(FIG_DIR,'control_*.png'))
    for i,pat in enumerate(images_control):
        (pa,fil)=os.path.split(images_control[i])
        images_control[i] = fil   
    
    images_results = glob.glob(os.path.join(FIG_DIR,'*_mean*'))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_std*')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_statmean*')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_statstd*')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_mean_bias*')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_syst_bias*')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_min_bias*')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_max_bias*')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_evol_bias*')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_corr*')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*rms*')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_coverage*')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'delta_t*.png')))
    #images_results.extend(glob.glob(os.path.join(FIG_DIR,'result_*.png')))
    
    for i,pat in enumerate(images_results):
        (pa,fil)=os.path.split(images_results[i])
        images_results[i] = fil


    
    if images_control == []:
        images_control = None
        
    if images_results == []:
        images_results = None
    
    html_tools.simplereporthtml(images_results=images_results, images_control=images_control,  mimamodel=mimamodel, mimaobs=mimaobs,  intro=ttforintro,  
                                file_out=FIG_DIR+'/MANGA_SEVIRI.html')
    
    

    
    
print "Ok"
# -- Fin de la validation
# ________________________________________________________________
