#! /usr/bin/env python
# -*- coding: utf-8 -*- 
# ________________________________________________________________
#
# G. Charria, S. Theetten, and J. Gatti 
#
# Created: 04/2010
# Last update: 01/2011
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
import cdtime,  MV2
from pylab import figure,  show,  close, title, ylabel, ylim
from matplotlib.dates import datestr2num
import numpy as np

# Import de module specifiques a Vacumm
#from vacumm.data.satellite.old_nar import *
from vacumm.data.satellite.nar import Nar
from vacumm.data.satellite.seviri import Seviri
from vacumm.data.model.mars.get_ftp import  get_ftp_f1
from vacumm.data.model.mars.get_cp import  get_cp_f1
from vacumm.data.misc.handle_directories import make_directories
from vacumm.validator.valid.ValidXYT import ValidXYT
from vacumm.misc.plot import map,  curve, curve2,  savefigs
from vacumm.misc.atime import * # a detailler: are_same_units + ...
from vacumm.misc.axes import create_time,  set_order
from vacumm.misc.grid.regridding import regrid1d, regrid2d
from vacumm.misc.io import ncread_best_estimate
#from vacumm.misc.grid.regridding import regrid1d,  regrid2d,  interp2d
import glob
from vacumm.markup import html_tools

# Remplacer les import * par des as ....

print 65*'-'
print ' Validation de la SST '
print 65*'-'

# -- Actions et Environnement de SST (Utilisation de la SST NAR)
# ----------------------------------------------------------------
print ' Definition des parametres d environnement '
print 65*'-'

# ---------------------------------------------------------
SCRIPT_DIR=os.getcwd()

config = ConfigParser.RawConfigParser()
config.read(os.path.join(SCRIPT_DIR,'config.cfg'))    
model_name = config.get('Model Description', 'name')

# Workdir
WK_DIR = config.get('Env', 'workdir')

# Creation des repertoires de resultats: Figures/nom_des_donnees/
if os.path.isdir(SCRIPT_DIR+'/Figures')==False:
        os.mkdir('Figures')
os.chdir(SCRIPT_DIR+'/Figures')
rr = config.get('Observations', 'product')
rr = rr.replace(' ', '_')
if os.path.isdir(SCRIPT_DIR+'/Figures/'+rr)==False:
        os.mkdir(rr)
FIG_DIR = os.path.join(SCRIPT_DIR,'Figures') 
FIG_DIR = os.path.join(FIG_DIR,rr) 
# Retour dans le dossier du Script
os.chdir(SCRIPT_DIR)
# Fin setup.py
# -----------------------------------------------------------------

# -- Obs: Rapatriement et lecture des observations corespondantes a la periode consideree
# ----------------------------------------------------------------
print ' -- Lecture des observations satellites --    '
print 65*'-'

# ----------------------------------------
# Tests pour l'utilisation des ConfigObj
# ----------------------------------------
#print 65*'++'
#test = Nar()
#print test.get_default_config()
#test.load_config('/export/home/gcharria/Work/VACUMM/vacumm/bin/test.cfg', nested=True)
#print test.get_config()
#print test
# ----------------------------------------
#print 65*'++'


if config.get('Observations', 'product') == 'NAR SST':
    obs=Nar()
if config.get('Observations', 'product') == 'OLD NAR SST':
    obs=Old_nar()
if config.get('Observations', 'product') == 'SEVIRI SST':
    obs=Seviri()

# -- Cette partie permet de charger une fichier de config qui viendra remplacer le .ini
# (dans le cas de l'utilisation des ConfigObj)
# --
#cfgfile = '/export/home/gcharria/Work/VACUMM/vacumm/bin/test.cfg'
#if os.path.isfile(cfgfile) == True :
#    obs.load_config(cfgfile, nested=True)
#    #print obs.get_config()
#else:
#    print 'Configuration file not found !'
#    sys.exit(0)
# -------------------------------------------------------------------------------------

# To update with .ini
if config.get('Observations', 'download') == 'True':
    obs.rappatrie()
# To update with .ini
if config.get('Observations', 'download') == 'assim_dir':
    obs.read_assim()
else:    
    obs.read()
    if config.get('Observations', 'clean_data') == 'True':
        obs.clean_data()

# -- Model : Rapatriement et lecture des sorties de modele corespondantes a la periode consideree
# ----------------------------------------------------------------
# ----------------------------------------------------------------
print ' -- Lecture des sorties de modeles --    '
print 65*'-'

# changement de repertoire et creation eventuelle de nouveaux repertoires
# les fichiers modeles sont rapatries dans vacumm/work/MODEL/MARS
dir_one ='MODEL'

#dir_two = 'MARS'
dir_two=config.get('Model Description', 'name').upper() 

# Cree les repertoires et se positionne dans le repertoire de rapatriment des fichiers
make_directories(SCRIPT_DIR, WK_DIR,  dir_one, dir_two)

# lecture de config.cfg
config = ConfigParser.RawConfigParser()
config.read(os.path.join(SCRIPT_DIR,'config.cfg'))	
andeb = config.getint('Time Period', 'andeb')
anfin = config.getint('Time Period', 'anfin')		
mdeb = config.getint('Time Period', 'mdeb')
mfin = config.getint('Time Period', 'mfin')
jdeb = config.getint('Time Period', 'jdeb')
jfin = config.getint('Time Period', 'jfin')
hdeb = config.getint('Time Period', 'hdeb')
hfin = config.getint('Time Period', 'hfin')
# fin lecture de config.cfg

ctdeb=cdtime.comptime(andeb,mdeb,jdeb,hdeb,0,0)

print 'La derniere heure consideree est forcee 23h du dernier jour !!!'
hfin = 23

ctfin=cdtime.comptime(anfin,mfin,jfin,hfin,0,0)

# rapatriement donnees F1 
# voir note readme_rapatriement_f1.py
dir_model = './' # Default / 
if config.get('Model Description', 'download') == 'ftp':
    # recuperation f1 via FTP
    get_ftp_f1(ctdeb, ctfin,andeb,SCRIPT_DIR)

if config.get('Model Description', 'download') == 'nfs':
    # recuperation f1 via NFS (cp command)
    dir_f1='/home/oo8/oo/modeles_previmer/f1/best_estimate/2011/'
    #dir_f1='/export/home/seb/vacumm/work/MODEL/MARS.old/'
    get_cp_f1(ctdeb, ctfin, dir_f1)

if config.get('Model Description', 'download') == 'opendap':
    print 'recup via opendap'
    # recuperation via OPENDAP (ncread_best_estimate command)
    # a faire ....
    # voir note dans vacumm/data/model/mars/
    # get_opendap.py, readme_ncread_best_estimate.py, 
    
if config.get('Model Description', 'download') == 'local_dir':
    # Pour cette option, les sorties modeles sont directement lues dans "output_dir"
    dir_model = config.get('Model Description',  'output_dir')
else:
    dir_model = os.getcwd()
print 'Read model output in ', dir_model 
print '---'
#--------------------------------------------------------------------
# ouverture, lecture des fichiers, extraction de la variable temperature de la couche de surface

# lecture de plusieurs fichiers entre ctdeb et ctfin et chargement de la variable TEMP de la couche de surface (dernier indice : 30)
ctfinplusuneheure=ctfin.add(1,cdtime.Days) # Sinon ne lit pas le dernier fichier correspondant a ctfin ...

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

# Repositionnement dans le repertoire 'bin' de lancement du script
os.chdir(SCRIPT_DIR)

# ------------------------------------------------------------
# on a 2 objets "CDMS" :
# model : objet CDMS
# obs : n'est pas exactement un objet CDMS mais un objet contenant l'objet cdms data
# Conversion / renommage de obs.data en obs
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
    #x1 = comptime(model.getTime())
    #x1 = create_time(x1)
    x2 = datestr2num(strtime(modelregridontime.getTime()))
    #x2 = comptime(modelregridontime.getTime())
    #x2 = create_time(x2)
    itest = 40
    # Trace des resultats
    fig = figure()
    # !! Ne fonctionne pas !!
    #    lp = curve2(model.data[:, 0, 40,40 ], color='g',  show=False,                       
    #        title=False, ylabel = 'Model [${^o}C$]',  
    #        ylabel_color='g', shadow=True, subplot=211)
    #    curve2(modelregridontime.data[:, 0, 40,40 ], color='r',  show=True,   
    #        marker = 'o', shadow=True, subplot=212, zorder = 100)
    #    savefigs(FIG_DIR+'/control_interp_time')
    #ax = fig.add_subplot(335)
    ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
    #ax.plot_date(x1,model.data[:,0,40,40], '-')    
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
    savefigs(FIG_DIR+'/control_interp_time')
    close()

print '-- Interpolation spatiale --'
print 65*'-'

#print modelregridontime.shape
#print obs.shape

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
    kwplot = dict(show=False, colorbar=True, colorbar_horizontal = True,  vmin=obs.min(), vmax=obs.max(), 
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
    savefigs(FIG_DIR+'/control_interp_map')
    P.close()

print '-- Conservation du masque le plus restrictif  --'
print 65*'-'
modelregridontime = modelregridontime[:, 0, :, :](squeeze=1)
m3 = np.ma.mask_or(obsregridspatial.mask, modelregridontime.mask)
modelregridontime.mask = m3
obsregridspatial.mask = m3
#modelregridontime.mask = np.ma.getmask(obsregridspatial)

if config.get('Statistics', 'to_do') == 'True':
    print ' -- Validation XYT --    ' 
    print 65*'-'
    # Le choix du Valid??? depend des champs lus ... ici SST en fonction de lon lat time => XYT
    result=ValidXYT(modelregridontime, obsregridspatial)
    
    if config.get('Statistics', 'temporal_mean') == 'True': 
        # -- it works !
        print ' - Calcul de la moyenne temporelle pour le modele et les observations - '
        result.temporal_mean()        
        kwplot = dict(vmin=result.obs.temp_mean.min(), vmax=result.obs.temp_mean.max(), nmax = 30, colorbar_extend='both')
        P.figure()
        map(result.model.temp_mean, title='Mean modelled Sea Surface Temperature',  show=False,  clabel_hide=True, **kwplot)
        savefigs(FIG_DIR+'/model_temporal_mean')
        P.close()
        P.figure()
        map(result.obs.temp_mean, title='Mean observed Sea Surface Temperature',  show=False,  clabel_hide=True, **kwplot)
        savefigs(FIG_DIR+'/obs_temporal_mean')
        P.close()
        
    if config.get('Statistics', 'spatial_mean') == 'True': 
        # -- it works !
        print ' - Calcul de la moyenne spatiale pour le modele et les observations - '
        result.spatial_mean()       
        
        # Trace des resultats
        fig = figure()
        ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        x = datestr2num(strtime(result.obs.spa_mean.getTime()))
        title('Spatial average over the domain')
        ylabel(obs.units)
        ax.plot_date(x,result.obs.spa_mean.data, marker='o', color='r',  markerfacecolor='r')
        ax.plot_date(x,result.model.spa_mean.data, marker='o', color='b',  markerfacecolor='b')  
        ax.legend(('Observations',  'Model'),  loc = 'upper right',  shadow = True,  fancybox=True)
        fig.autofmt_xdate()
        savefigs(FIG_DIR+'/result_spatial_statmean')
        close()
      
    
        
        
    if config.get('Statistics', 'temporal_std') == 'True':    
        # -- it works !
        print ' Calcul de la std en chaque point geographique'
        result.temporal_std()
        kwplot = dict(vmin=result.obs.temp_std.min(), vmax=result.obs.temp_std.max(), nmax = 30, colorbar_extend='max')
        P.figure()
        map(result.model.temp_std, title='Model - std(SST(t))',  show=False,  clabel_hide=True, **kwplot)
        savefigs(FIG_DIR+'/model_temporal_std')
        P.close()
        P.figure()
        map(result.obs.temp_std, title='Observation - std(SST(t))', show=False,  clabel_hide=True, **kwplot)
        savefigs(FIG_DIR+'/obs_temporal_std')
        P.close()
    
    if config.get('Statistics', 'spatial_std') == 'True': 
        # -- it works !
        print ' - Calcul de l ecart type pour le modele et les observations pour chaque pas de temps - '
        result.spatial_std()       
        
        # Trace des resultats
        fig = figure()
        ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        x = datestr2num(strtime(result.obs.spa_std.getTime()))
        title('Spatial STD over the domain')
        ylabel(obs.units)
        ax.plot_date(x,result.obs.spa_std.data, marker='o', color='r',  markerfacecolor='r')
        ax.plot_date(x,result.model.spa_std.data, marker='o', color='b',  markerfacecolor='b')  
        ax.legend(('Observations',  'Model'),  loc = 'upper right',  shadow = True,  fancybox=True)
        fig.autofmt_xdate()
        savefigs(FIG_DIR+'/result_spatial_statstd')
        close()
        
    if config.get('Statistics', 'evol_biais') == 'True':
        # -- it works !
        print ' - Calcul de l''evolution du biais moyen entre le modele et les observations - '
        if config.get('Statistics', 'spatial_mean') == 'False': # Pour ne pas faire le calcul 2 fois
            result.spatial_mean()       
        
        
        
        # Trace des resultats
        fig = figure()
        ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        x = datestr2num(strtime(result.obs.spa_mean.getTime()))
        title('Mean bias over the domain (<E> = <OBS> - <MODEL>)')
        ylabel(obs.units)
        ax.plot_date(x,result.obs.spa_mean.data-result.model.spa_mean.data, marker='o', markeredgecolor='b',  markerfacecolor='b')
        ax.plot_date(x,np.zeros(result.obs.spa_mean.data.shape), 'k:')
        
        fig.autofmt_xdate()
        savefigs(FIG_DIR+'/result_evol_bias')
        close()
        
    if config.get('Statistics', 'mean_biais') == 'True':
        # -- it works !
        print ' - Calcul du biais moyen entre le modele et les observations - '
        if config.get('Statistics', 'temporal_mean') == 'False': # Pour ne pas faire le calcul 2 fois
            result.temporal_mean()       
        
        interm = result.obs.temp_mean-result.model.temp_mean
        
        #map(result.model.temp_mean, res='c', show=False)        
        P.figure()
        kwplot = dict(vmin=-abs(interm).max(), vmax=abs(interm).max(), nmax=20)
        map(interm, title='Mean bias over the domain (<E> = <OBS> - <MODEL>)',  show=False,  clabel_hide=True,  **kwplot)
        savefigs(FIG_DIR+'/result_mean_bias')
        P.close()                
        
    if config.get('Statistics', 'extreme_bias') == 'True':        
        print ' - Calcul du biais minimum et maximum '        
        result.bias_extrema()
        
        P.figure()
        kwplot = dict(vmin=-abs(result.biasmin).max(), vmax=abs(result.biasmin).max())
        map(result.biasmin, title='Minimum bias',  show=False,  clabel_hide=True,  **kwplot)
        savefigs(FIG_DIR+'/result_min_bias')        
        P.close()       
        
        P.figure()
        kwplot = dict(vmin=-abs(result.biasmax).max(), vmax=abs(result.biasmax).max())
        map(result.biasmax, title='Maximum bias',  show=False,  clabel_hide=True,  **kwplot)
        savefigs(FIG_DIR+'/result_max_bias')        
        P.close()  
        
    if config.get('Statistics', 'systematic_bias') == 'True':     
        print ' - Calcul du biais systematique '
        result.systematic_bias()
        
        P.figure()
        kwplot = dict(vmin=-1, vmax=1,  levels=[-1, -0.1,0.1,  1], colorbar_ticks=[1,0,-1])
        map(result.biassyst, title='Systematic bias\n (+1: model underestimation, -1: model overestimation)',  show=False,  clabel_hide=True,  **kwplot)
        savefigs(FIG_DIR+'/result_syst_bias')        
        P.close()
        
        
    if config.get('Statistics', 'extrema') == 'True':     
        # -- it works !
        print ' - Extractions des minima et maxima - '
        result.extrema()
    
    if config.get('Statistics', 'obs_coverage') == 'True':    
        # -- it works !
        print ' - Calcul du taux de couverture des observations en chaque point geographique - '
        result.obs_coverage()        
        P.figure()
        map(result.obs_cov,  title='Observation coverage during the period', show=False,  clabel_hide=True)
        savefigs(FIG_DIR+'/result_obs_coverage')
        P.close()
    
    if config.get('Statistics', 'obs_spatialcoverage') == 'True':    
        # -- it works !
        print ' - Calcul du taux de couverture des observations en chaque pas de temps - '
        result.spatial_obs_coverage()        
        
        # Ne fonctionne pas !!!!
        #select = dict(color='g', shadow=True, ylabel=result.obs_spacov.units , show=False,  title='Observation coverage following the time period period')
        #curve2(result.obs_spacov, **select)
        
        # Trace des resultats
        fig = figure()
        ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        x = datestr2num(strtime(result.obs_spacov.getTime()))
        title('Observation coverage following the time period')
        ylabel('Coverage (%)')
        ax.plot_date(x,result.obs_spacov.data,linestyle='-', marker='o', color='r',  markerfacecolor='b')
        ylim(0, 100)
        fig.autofmt_xdate()
        savefigs(FIG_DIR+'/result_obs_spatial_coverage')
        close()
        
    
  
    
    if config.get('Statistics', 'temporal_rms') == 'True':    
        # -- it works !
        print ' - Calcul de la RMS temporelle - '
        result.temporal_rms()
        P.figure()
        map(result.temp_rms, title='RMS\n (uncentered and biased)', show=False,  clabel_hide=True)
        savefigs(FIG_DIR+'/result_temporal_rms')
        P.close()
        
    if config.get('Statistics', 'temporal_rmsc') == 'True':    
        # -- it works !
        print ' - Calcul de la RMS temporelle  centree- '
        result.temporal_rmsc()
        P.figure()
        map(result.temp_rmsc, title='RMS\n (centered and biased)', show=False,  clabel_hide=True)
        savefigs(FIG_DIR+'/result_temporal_rmsc')
        P.close()
        
    if config.get('Statistics', 'spatial_rms') == 'True': 
        # -- it works !
        print ' - Calcul de la RMS entre le modele et les observations pour chaque pas de temps - '
        result.spatial_rms()       
        
        # Trace des resultats
        fig = figure()
        ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        x = datestr2num(strtime(result.spa_rms.getTime()))
        title('Spatial RMS over the domain')
        ylabel(obs.units)
        ax.plot_date(x,result.spa_rms.data, marker='o', color='r',  markerfacecolor='r')
        fig.autofmt_xdate()
        savefigs(FIG_DIR+'/result_spatial_statrms')
        close()
        
    if config.get('Statistics', 'spatial_rmsc') == 'True': 
        # -- it works !
        print ' - Calcul de la RMS centree entre le modele et les observations pour chaque pas de temps - '
        result.spatial_rmsc()       
        
        # Trace des resultats
        fig = figure()
        ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        x = datestr2num(strtime(result.spa_rmsc.getTime()))
        title('Spatial centered RMS over the domain')
        ylabel(obs.units)
        ax.plot_date(x,result.spa_rmsc.data, marker='o', color='r',  markerfacecolor='r')
        fig.autofmt_xdate()
        savefigs(FIG_DIR+'/result_spatial_statrmsc')
        close()
    
    if config.get('Statistics', 'temporal_corr') == 'True':    
        # -- it works !
        print ' - Calcul de la Correlation temporelle - '
        result.temporal_corr()
        P.figure()
        map(result.temp_corr, title='Time series correlation', show=False,  clabel_hide=True)
        savefigs(FIG_DIR+'/result_temporal_corr')
        P.close()
    print 65*'-'
    

if config.get('Report', 'rep_html') == 'True':
    print ' -- Generation du Rapport de Validation --    ' 
    print 65*'-'
    #title = "Test(1)"
    #header = ""
    #footer = ""
    
    ttforintro = '<u>Observations</u>: %(#)s <br /> <u>Model</u>: %(##)s %(###)s <br /><hr>' %{'#':obs.long_name, '##':config.get('Model Description', 'name') , '###':model.long_name}
    
    
    if config.get('Statistics', 'extrema') == 'True' and config.get('Statistics', 'to_do') == 'True':
        mimamodel = '=> Modelled minimum and maximum: %(#)6.3f / %(##)6.3f %(u)s' %{'#':result.model.extrema[0], '##':result.model.extrema[1], 'u':model.units}
        mimaobs = '=> Observed minimum and maximum: %(###)6.3f / %(####)6.3f %(u)s'  %{'###':result.obs.extrema[0],'####':result.obs.extrema[1], 'u':obs.units}
    else:
        mimamodel = ''
        mimaobs = ''
            
    images_control = glob.glob(os.path.join(FIG_DIR,'control_*.png'))
    for i,pat in enumerate(images_control):
        (pa,fil)=os.path.split(images_control[i])
        images_control[i] = fil   
    
    images_results = glob.glob(os.path.join(FIG_DIR,'*_mean.png'))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_std.png')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_statmean.png')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_statstd.png')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_mean_bias.png')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_syst_bias.png')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_min_bias.png')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_max_bias.png')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_evol_bias.png')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_corr.png')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*rms*.png')))
    images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_coverage.png')))
    #images_results.extend(glob.glob(os.path.join(FIG_DIR,'result_*.png')))
    
    for i,pat in enumerate(images_results):
        (pa,fil)=os.path.split(images_results[i])
        images_results[i] = fil


    
    if images_control == []:
        images_control = None
        
    if images_results == []:
        images_results = None
    
    html_tools.simplereporthtml(images_results=images_results, images_control=images_control,  mimamodel=mimamodel, mimaobs=mimaobs,  intro=ttforintro,  
                                file_out=FIG_DIR+'/'+rr+'.html')
    
    

    
    
print "Ok"
# -- Fin de la validation
# ________________________________________________________________
