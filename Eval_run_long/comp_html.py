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
import os,sys,shutil,ConfigParser, gc, glob
from vacumm.markup import html_tools

# Ouverture fichier config
config = ConfigParser.RawConfigParser()
config.read('config_html.cfg')

# FIG_DIR1
FIG_DIR1=config.get('Model Description','dir1')

# FIG_DIR2
FIG_DIR2=config.get('Model Description','dir2')

fileout=config.get('General','html_file')


# Deplacement des figures et renommage.
repout = config.get('General','output_main_dir')
if os.path.isdir(repout)==False:
    os.mkdir(repout)
    os.mkdir(repout+'/model1')
    os.mkdir(repout+'/model2')

os.chdir(repout)

for i,fln in enumerate(glob.glob(os.path.join(FIG_DIR1,'*.png'))):
    (pa,fil)=os.path.split(fln)
    shutil.copyfile(fln, os.path.join('model1/',fil))
for i,fln in enumerate(glob.glob(os.path.join(FIG_DIR2,'*.png'))):
    (pa,fil)=os.path.split(fln)
    shutil.copyfile(fln, os.path.join('model2/',fil))


FIG_DIR1 = os.path.join('model1')
FIG_DIR2 = os.path.join('model2')


# ------------------
print ' -- Generation du Rapport de Validation --    ' 
print 65*'-'
title = "Comparaison 2 runs"
#header = ""
#footer = ""

ttforintro = '<u>Model 1</u>: %(#)s <br /> <u>Model 2</u>: %(#2)s <br /><hr>' %{'#':config.get('Model Description', 'model1'), '#2':config.get('Model Description', 'model2')}
   
images_mean = glob.glob(os.path.join(FIG_DIR1,'model_temporal_mean_all.png'))
images_mean.extend(glob.glob(os.path.join(FIG_DIR2,'model_temporal_mean_all.png')))
images_mean.extend(glob.glob(os.path.join(FIG_DIR2,'obs_temporal_mean_all.png')))

images_mean.extend(glob.glob(os.path.join(FIG_DIR1,'model_temporal_mean_djf.png')))
images_mean.extend(glob.glob(os.path.join(FIG_DIR2,'model_temporal_mean_djf.png')))
images_mean.extend(glob.glob(os.path.join(FIG_DIR2,'obs_temporal_mean_djf.png')))

images_mean.extend(glob.glob(os.path.join(FIG_DIR1,'model_temporal_mean_mam.png')))
images_mean.extend(glob.glob(os.path.join(FIG_DIR2,'model_temporal_mean_mam.png')))
images_mean.extend(glob.glob(os.path.join(FIG_DIR2,'obs_temporal_mean_mam.png')))

images_mean.extend(glob.glob(os.path.join(FIG_DIR1,'model_temporal_mean_jja.png')))
images_mean.extend(glob.glob(os.path.join(FIG_DIR2,'model_temporal_mean_jja.png')))
images_mean.extend(glob.glob(os.path.join(FIG_DIR2,'obs_temporal_mean_jja.png')))

images_mean.extend(glob.glob(os.path.join(FIG_DIR1,'model_temporal_mean_son.png')))
images_mean.extend(glob.glob(os.path.join(FIG_DIR2,'model_temporal_mean_son.png')))
images_mean.extend(glob.glob(os.path.join(FIG_DIR2,'obs_temporal_mean_son.png')))

images_cov = glob.glob(os.path.join(FIG_DIR1,'result_obs_coverage_all.png'))

# ------------------------------

images_std =  glob.glob(os.path.join(FIG_DIR1,'model_temporal_std_all.png'))
images_std.extend(glob.glob(os.path.join(FIG_DIR2,'model_temporal_std_all.png')))
images_std.extend(glob.glob(os.path.join(FIG_DIR2,'obs_temporal_std_all.png')))

images_std.extend(glob.glob(os.path.join(FIG_DIR1,'model_temporal_std_djf.png')))
images_std.extend(glob.glob(os.path.join(FIG_DIR2,'model_temporal_std_djf.png')))
images_std.extend(glob.glob(os.path.join(FIG_DIR2,'obs_temporal_std_djf.png')))

images_std.extend(glob.glob(os.path.join(FIG_DIR1,'model_temporal_std_mam.png')))
images_std.extend(glob.glob(os.path.join(FIG_DIR2,'model_temporal_std_mam.png')))
images_std.extend(glob.glob(os.path.join(FIG_DIR2,'obs_temporal_std_mam.png')))

images_std.extend(glob.glob(os.path.join(FIG_DIR1,'model_temporal_std_jja.png')))
images_std.extend(glob.glob(os.path.join(FIG_DIR2,'model_temporal_std_jja.png')))
images_std.extend(glob.glob(os.path.join(FIG_DIR2,'obs_temporal_std_jja.png')))

images_std.extend(glob.glob(os.path.join(FIG_DIR1,'model_temporal_std_son.png')))
images_std.extend(glob.glob(os.path.join(FIG_DIR2,'model_temporal_std_son.png')))
images_std.extend(glob.glob(os.path.join(FIG_DIR2,'obs_temporal_std_son.png')))



# -----------------------------------------

images_biasmap =  glob.glob(os.path.join(FIG_DIR1,'result_mean_bias_all.png'))
images_biasmap.extend(glob.glob(os.path.join(FIG_DIR2,'result_mean_bias_all.png')))

images_biasmap.extend(glob.glob(os.path.join(FIG_DIR1,'result_mean_bias_djf.png')))
images_biasmap.extend(glob.glob(os.path.join(FIG_DIR2,'result_mean_bias_djf.png')))

images_biasmap.extend(glob.glob(os.path.join(FIG_DIR1,'result_mean_bias_mam.png')))
images_biasmap.extend(glob.glob(os.path.join(FIG_DIR2,'result_mean_bias_mam.png')))

images_biasmap.extend(glob.glob(os.path.join(FIG_DIR1,'result_mean_bias_jja.png')))
images_biasmap.extend(glob.glob(os.path.join(FIG_DIR2,'result_mean_bias_jja.png')))

images_biasmap.extend(glob.glob(os.path.join(FIG_DIR1,'result_mean_bias_son.png')))
images_biasmap.extend(glob.glob(os.path.join(FIG_DIR2,'result_mean_bias_son.png')))

# -----------------------------------------

images_bias1d =  glob.glob(os.path.join(FIG_DIR1,'result_evol_bias_all.png'))
images_bias1d.extend(glob.glob(os.path.join(FIG_DIR2,'result_evol_bias_all.png')))

images_bias1d.extend(glob.glob(os.path.join(FIG_DIR1,'result_evol_bias_manche_ouest.png')))
images_bias1d.extend(glob.glob(os.path.join(FIG_DIR2,'result_evol_bias_manche_ouest.png')))
images_bias1d.extend(glob.glob(os.path.join(FIG_DIR1,'result_evol_bias_manche_est.png')))
images_bias1d.extend(glob.glob(os.path.join(FIG_DIR2,'result_evol_bias_manche_est.png')))
images_bias1d.extend(glob.glob(os.path.join(FIG_DIR1,'result_evol_bias_bretagne_sud.png')))
images_bias1d.extend(glob.glob(os.path.join(FIG_DIR2,'result_evol_bias_bretagne_sud.png')))
images_bias1d.extend(glob.glob(os.path.join(FIG_DIR1,'result_evol_bias_vendee.png')))
images_bias1d.extend(glob.glob(os.path.join(FIG_DIR2,'result_evol_bias_vendee.png')))
images_bias1d.extend(glob.glob(os.path.join(FIG_DIR1,'result_evol_bias_basque.png')))
images_bias1d.extend(glob.glob(os.path.join(FIG_DIR2,'result_evol_bias_basque.png')))

html_tools.simplecomparisonhtml(title = title, images_cov = images_cov, images_mean = images_mean,  images_std = images_std,  images_biasmap = images_biasmap, images_bias1d = images_bias1d, intro=ttforintro, 
    file_out=fileout)
    
    

    
    
print "Ok"
# -- Fin de la validation
# ________________________________________________________________
