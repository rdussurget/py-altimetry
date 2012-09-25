#! /usr/bin/env python
# -*- coding: utf-8 -*- 
# ________________________________________________________________
#
# S.Skrypnikov-Laviolle
#
# Created: 02/2012
# Last update: 05/2012
# Modified: 09/2012 (G. Charria)
# ________________________________________________________________
#
# Base sur la chaine d'evaluation des performances modele en SST
# developpe par P. Craneguy (Actimar) - Projet PRECOC (ANR)
# dans la version 1.0 sur 25/01/2007 
# ________________________________________________________________

import os,sys,shutil,ConfigParser
# -- These two line allows running without display
# Warning: does not allow screen figure display
import matplotlib
matplotlib.use('Agg')
# --
import matplotlib.pyplot as P
import cdtime,  MV2, cdms2
from pylab import figure,  show,  close, title, ylabel, ylim
from matplotlib.dates import datestr2num
import numpy as np

# Import de module specifiques a Vacumm
from vacumm.data.satellite.seviri import Seviri
from vacumm.data.model.mars.get_ftp import  get_ftp_f1
from vacumm.data.model.mars.get_cp import  get_cp_f1
from vacumm.data.misc.handle_directories import make_directories
from vacumm.validator.valid.ValidXYT import ValidXYT
from vacumm.misc.plot import map2 as map
from vacumm.misc.plot import curve, curve2,  savefigs
from vacumm.misc.atime import add, strtime, ch_units, are_same_units, comptime, strftime
from vacumm.misc.axes import  set_order
from vacumm.misc.grid.regridding import regrid1d, regrid2d
from vacumm.misc.io import ncread_best_estimate

import glob
from vacumm.markup import html_tools
from vacumm.misc.color import cmap_custom, StepsNorm, cmap_magic
from global_stat import allstat, regionalstat, detailedstat, seasonalstat, monthlystat



SCRIPT_DIR=os.getcwd()

config = ConfigParser.RawConfigParser()
config.read(os.path.join(SCRIPT_DIR,'config.cfg'))

os.chdir(SCRIPT_DIR+'/Figures')
rr = config.get('Observations', 'product')
rr = rr.replace(' ', '_')

FIG_DIR = os.path.join(SCRIPT_DIR,'Figures') 
FIG_DIR = os.path.join(FIG_DIR,rr)
	
andeb = config.getint('Time Period', 'andeb')
anfin = config.getint('Time Period', 'anfin')
mdeb = config.getint('Time Period', 'mdeb')
mfin = config.getint('Time Period', 'mfin')
jdeb = config.getint('Time Period', 'jdeb')
jfin = config.getint('Time Period', 'jfin')
hdeb = config.getint('Time Period', 'hdeb')
hfin = config.getint('Time Period', 'hfin')

ctdeb=cdtime.comptime(andeb,mdeb,jdeb,hdeb,0,0)

print 'La derniere heure consideree est forcee 23h du dernier jour !!!'
hfin = 23

ctfin=cdtime.comptime(anfin,mfin,jfin,hfin,0,0)
tag1= strftime('%Y%m%d',ctdeb)
tag2= strftime('%Y%m%d',ctfin)
tagforfiledate1 = '_'.join(tag1.split(' ')).lower()    
tagforfiledate2 = '_'.join(tag2.split(' ')).lower()
# ouverture, lecture des fichiers, extraction de la variable temperature de la couche de surface

# lecture de plusieurs fichiers entre ctdeb et ctfin et chargement de la variable TEMP de la couche de surface (dernier indice : 30)
ctfinplusuneheure=ctfin.add(1,cdtime.Days) # Sinon ne lit pas le dernier fichier correspondant a ctfin ...

if config.get('Model Description', 'download') == 'local_dir':
# Pour cette option, les sorties modeles sont directement lues dans "output_dir"
    dir_model = config.get('Model Description',  'output_dir')
else:
    dir_model = os.path.join(config.get('Env', 'workdir'), 'MODEL',config.get('Model Description', 'name').upper()) 
# F1 MANGA

if config.get('Model Description', 'name') == 'mars_manga':
    model = ncread_best_estimate('TEMP',os.path.join(dir_model,"PREVIMER_F1-MARS3D-MANGA4000_%Y%m%dT%H00Z.nc"), (ctdeb, ctfinplusuneheure),select=dict(level=slice(29,30)))
# modif J.Gatti : F1 MANGA devient F2 MENOR
if config.get('Model Description', 'name') == 'mars_menor':
    model = ncread_best_estimate('TEMP',dir_model+"PREVIMER_F2-MARS3D-MENOR_%Y%m%dT%H00Z.nc", (ctdeb, ctfinplusuneheure),select=dict(level=slice(29,30)))
# pour MFS
if config.get('Model Description', 'name') == 'mfs':
    model = ncread_best_estimate('votemper',dir_model+"INGV_MFSTEP_%Y%m%d_T.nc", (ctdeb, ctfinplusuneheure),select=dict(level=slice(0,1)))
    set_order(model, 'tzyx')   # pour que averager.py fonctionne
    lomin = float(config.get('Domain', 'lomin') ) 
    lomax = float(config.get('Domain', 'lomax') )
    lamin = float(config.get('Domain', 'lamin') )
    lamax = float(config.get('Domain', 'lamax') )   
    model = model(lon=(lomin, lomax), lat=(lamin, lamax))
# fin modif J.Gatti

# a faire : test sur le contenu du tableau ...
# Supprime la dimension à 1 (vertical level dans ce cas)




# ------------------------------------------------------------
# on a 2 objets "CDMS" :
# model : objet CDMS
# obs : n'est pas exactement un objet CDMS mais un objet contenant l'objet cdms data
# Conversion / renommage de obs.data en obs
#obs=ncread_best_estimate('temp',os.path.join(dir_obs,filename_obs), (ctdeb1, ctfin1))
os.chdir(SCRIPT_DIR)
obs=Seviri()
obs.read() # lecture des observations chargées dans init
obs = obs.data
# ------------------------------------------------------------



if are_same_units(model.getTime().units,obs.getTime().units)==False:
    print '-- Conversion des unites de temps --'
    print 65*'-'
    ch_units(model, obs.getTime().units)

print '-- Interpolation temporelle --'
print 65*'-'
time_ref = obs.getTime()
#modelregridontime = regrid1d(model, time_ref, 'linear', axis=0)
modelregridontime = regrid1d(model, time_ref, 'remap', axis=0)

if config.get('Control', 'tseries') == 'True':
    print '-- Trace de controle (1) --'
    print 65*'-'
    x1 = datestr2num(strtime(model.getTime()))
    x2 = datestr2num(strtime(modelregridontime.getTime()))
    itest = 40
    # Trace des resultats
    fig = figure()
    ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])   
    lo = model.getLongitude().getValue()
    la = model.getLatitude().getValue()
    #modif J.Gatti - conversion des variables lo, la 2D en 1D
    if config.get('Model Description', 'name') == 'mfs':
	lo=lo[0, :]
	la=la[:, 0]
    #fin modif J.Gatti    
    title('%(#)3.2fE / %(##)4.2fN - Control temporal interpolation' %{'#':lo[itest], '##':la[itest]})
    ylabel('Temperature ($^{\circ}C$)')
    ax.plot_date(x1,model.data[:,0,itest,itest],linestyle='-', marker='o', color='r')
    ax.plot_date(x2,modelregridontime.data[:,0,itest,itest],'bs')
    ax.legend(('Model Before intepolation',  'Model after interpolation'),  loc = 'upper left',  shadow = True,  fancybox=True)
    fig.autofmt_xdate()
    #show()
    # #P.show()
    savefigs(FIG_DIR+'/control_interp_time'+'_' +tagforfiledate1 +'_'+tagforfiledate2)
    close()

print '-- Interpolation spatiale --'
print 65*'-'


cgrido = modelregridontime.getGrid()
# !!!! Essayer de faire fonctionner 'remap' !!!!
obsregridspatial = regrid2d(obs,  cgrido,  method = 'nearest')
#cgrido = obs.getGrid()
#model = interp2d(model,  cgrido,  method = 'nat')

if config.get('Control', 'interp_map') == 'True':
    print '-- Trace de controle (2) --'
    print 65*'-'
    # Trace des resultats
    kwplot1 = dict(show=False, colorbar=False, vmin=obs.min(), vmax=obs.max(), 
	drawparallels_size=8, drawmeridians_size=8, drawmeridians_rotation=45.,  clabel_hide=True)
    kwplot = dict(show=False, colorbar=True, vmin=obs.min(), vmax=obs.max(), 
	    drawparallels_size=8, drawmeridians_size=8, drawmeridians_rotation=45.,  clabel_hide=True)
    m = map(model[-1, 0, :, :], title='Model - original\n (last time step)',  subplot=221, xhide=1,**kwplot1)
    #add_grid(cgrido, lw=.7, alpha=.3)
    map(modelregridontime[-1, 0, :, :], title='Model - after time interpolation\n (last time step)', subplot=222, xhide=1,yhide=1,m=m, **kwplot1)
    #add_grid(cgridi, lw=.7, alpha=.3)
    if float(np.ma.count(obs[-1,:,:])) / float(np.ma.count_masked(obs[-1,:,:])) < 0.1:
		tt = 'Obs - Original\n (< 10% of non-masked values)'
    else:
		tt = 'Obs - Original\n (last time step)'
    map(obs[-1, :, :], title=tt, subplot=223, m=m, **kwplot)
    #m = map(obs.data[0, :, :], title='Obs - Original', subplot=223, **kwplot)
    #add_grid(cgridi, lw=.7, alpha=.3)
    if float(np.ma.count(obsregridspatial[-1,:,:])) / float(np.ma.count_masked(obsregridspatial[-1,:,:])) < 0.1:
		tt = 'Obs - after spatial interpolation\n (< 10% of non-masked values)'
    else:
		tt = 'Obs - after spatial interpolation\n (last time step)'
    map(obsregridspatial[-1, :, :], title=tt,  subplot=224,yhide=1,m=m, **kwplot)
    #add_grid(cgridi, lw=.7, alpha=.3)
    savefigs(FIG_DIR+'/control_interp_map'+'_' +tagforfiledate1 +'_'+tagforfiledate2)
    P.close()

print 65*'-'
print '-- Conservation du masque le plus restrictif  --'
print 65*'-'
modelregridontime = modelregridontime[:, 0, :, :](squeeze=1)
m3 = np.ma.mask_or(obsregridspatial.mask, modelregridontime.mask)
modelregridontime.mask = m3
obsregridspatial.mask = m3
#modelregridontime.mask = np.ma.getmask(obsregridspatial)

if config.get('Statistics', 'to_do') == 'True':
    #print ' -- Validation XYT -- ' 
    print 65*'-'
    #Appel de global_stat /allstat ou /seasonalstat ou /monthlystat ou /regionalstat
    os.chdir(SCRIPT_DIR)
    if config.get('Statistics', 'allstat') == 'True':
		allstat(modelregridontime, obsregridspatial, FIG_DIR, SCRIPT_DIR) 
    if config.get('Statistics', 'regionalstat') == 'True':
		regionalstat(modelregridontime, obsregridspatial, FIG_DIR, SCRIPT_DIR)
    if config.get('Statistics', 'seasonalstat') == 'True':
		seasonalstat(modelregridontime, obsregridspatial, FIG_DIR, SCRIPT_DIR)  
    if config.get('Statistics', 'monthlystat') == 'True':
		monthlystat(modelregridontime, obsregridspatial, FIG_DIR, SCRIPT_DIR)



    
print "Ok"
# -- Fin de la validation
# ________________________________________________________________
