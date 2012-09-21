#! /usr/bin/env python
# -*- coding: utf-8 -*- 
# ________________________________________________________________
#
# S.Skrypnikov-Laviolle
#
# Created: 02/2012
# Last update: 05/2012
# Modified: 05/2012 (G. Charria)
# ________________________________________________________________
#
# Initialisation de la validation Modèle/Observation
# ________________________________________________________________

import os,sys,shutil,ConfigParser
# -- These two line allows running without display
# Warning: does not allow screen figure display
import matplotlib
matplotlib.use('Agg')
# --
import cdtime,  MV2
import numpy as np


# Import de module specifiques a Vacumm
#from vacumm.data.satellite.old_nar import *
from vacumm.data.satellite.seviri import Seviri
from vacumm.data.model.mars.get_ftp import  get_ftp_f1
from vacumm.data.model.mars.get_cp import  get_cp_f1
from vacumm.data.misc.handle_directories import make_directories
from vacumm.misc.atime import comptime
from ConfigParser import SafeConfigParser
from vacumm.misc.atime import strftime

import glob
from subprocess import call
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
print SCRIPT_DIR

config = ConfigParser.RawConfigParser()
config.read(os.path.join(SCRIPT_DIR,'config.cfg'))    
model_name = config.get('Model Description', 'name')

# Workdir
WK_DIR = config.get('Env', 'workdir')
print WK_DIR
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
#print 65*'++'add, strtime, ch_units, are_same_units, comptime
#test = Nar()
#print test.get_default_config()
#test.load_config('/export/home/gcharria/Work/VACUMM/vacumm/bin/test.cfg', nested=True)
#print test.get_config()
#print test
# ----------------------------------------
#print 65*'++'


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
    if config.get('Observations', 'clean_data') == 'True':
        obs.clean_data()

# -- Model : RapatrieSCRIPTment et lecture des sorties de modele corespondantes a la periode consideree
# ----------------------------------------------------------------
# ----------------------------------------------------------------
print ' -- Lecture des sorties de modeles --    ' 
tag = 'all'
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
flagYear=0
if (andeb==anfin and mdeb>mfin) or (andeb>anfin):
    print 'la validation concerne une plage de temps non valide (debut ulterieur a la fin)'
if andeb<anfin:
    print 'la validation est effectuée sur plusieurs années '
    flagYear=1
    
    #condition dans la lecture des model et obs pour changer de repertoires
ctdeb=cdtime.comptime(andeb,mdeb,jdeb,hdeb,0,0)

print 'La derniere heure consideree est forcee 23h du dernier jour !!!'
hfin = 23

ctfin=cdtime.comptime(anfin,mfin,jfin,hfin,0,0)

# rapatriement donnees F1 
# voir note readme_rapatriement_f1.py
dir_model = './' # Default 
dir_f1 = ' '
if config.get('Model Description', 'download') == 'ftp':
    # recuperation f1 via FTP
    get_ftp_f1(ctdeb, ctfin,andeb,SCRIPT_DIR)

if config.get('Model Description', 'download') == 'nfs':
    print 'nfs prise en compte'
    
    # recuperation f1 via NFS (cp command)
    if config.get('Model Description', 'name') == 'mars_manga':
	print 'modèle Mars_Manga4000 pris en compte'
	    
	if flagYear==1:
	    diff=anfin-andeb
	    ctfin=cdtime.comptime(andeb,12,31,23,0,0)
	    dir_f1=os.path.join(config.get('MARS F1','dir'),config.get('Time Period','andeb')+'/')
	    get_cp_f1(ctdeb, ctfin, dir_f1)
	    for i in np.arange(diff):
		
		ctdeb=cdtime.comptime(andeb+1,01,01,0,0,0)
		    
		if andeb+1==anfin:
		    ctfin=cdtime.comptime(anfin,mfin,jfin,hfin,0,0)
		    dir_f1=os.path.join(config.get('MARS F1','dir'),config.get('Time Period','anfin')+'/')
		    get_cp_f1(ctdeb, ctfin, dir_f1)
		    print dir_f1
		else :
		    ctfin=cdtime.comptime(andeb+1,12,31,23,0,0)
		    dir_f1=os.path.join(config.get('MARS F1','dir'),config.get('Time Period','andeb'+1)+'/')
		    get_cp_f1(ctdeb, ctfin, dir_f1)
		    print dir_f1
		andeb+=andeb
		
	else:
	    dir_f1= os.path.join(config.get('MARS F1','dir'),config.get('Time Period','andeb')+'/')
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
print 25*'-'
#--------------------------------------------------------------------
#Lancement de l'evaluation avec eval_run_sst_op_previmer.py

# Lancer Subprocess pour démarrer la routine batch

os.chdir(SCRIPT_DIR)
call(['vacumm.batch' , 'eval_run_sst_long.py'])


