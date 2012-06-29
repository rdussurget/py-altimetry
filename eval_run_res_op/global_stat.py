# -*- coding: utf8 -*-

"""
Global statistics for pyvalid project 

Author: S. Skrypnikov-Laviolle (03/2012)
Last update 06/2012 (G. Charria)


   Statistics on the whole dataset/run ... 
"""

import os,sys,shutil,ConfigParser, gc
#import cdms2



import matplotlib.pyplot as P
import cdtime,  MV2, cdms2
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
from vacumm.misc.plot import map2 as map
from vacumm.misc.plot import curve2,  savefigs
from vacumm.misc.atime import add, strtime, ch_units, are_same_units, comptime, strftime
from vacumm.misc.axes import create_time,  set_order, create_dep, create_lat, create_lon
from vacumm.misc.grid.regridding import regrid1d, regrid2d
from vacumm.misc.io import ncread_best_estimate
#from vacumm.misc.grid.regridding import regrid1d,  regrid2d,  interp2d
import glob, time
from vacumm.markup import html_tools
from vacumm.misc.color import cmap_custom, StepsNorm, cmap_magic
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import inspect
global SCRIPT_DIR
SCRIPT_DIR=os.getcwd()
print SCRIPT_DIR


    # Pour produire du NetCDF3
cdms2.setNetcdfShuffleFlag(0); cdms2.setNetcdfDeflateFlag(0); cdms2.setNetcdfDeflateLevelFlag(0)
def write_nc_cf(filename,dual=False, varin1=None,varin2=None):
    config = ConfigParser.RawConfigParser()
    config.read(os.path.join(SCRIPT_DIR,'config.cfg'))
    rep_ecriture= config.get('Output','rep_ecriture')
    filename= os.path.join(rep_ecriture,filename)
    f = cdms2.open(filename, 'w')
    f.write(varin1) # ecriture d'une variable
    if dual==True:
	f.write(varin2) #ecriture
    creation_date = time.strftime('%Y-%m-%dT%H:%M:%SZ')
    f.creation_date = creation_date
    f.title = config.get('Output','title')
    #print f
    f.close() # fermeture

    #-- Probleme a voir ... comment ne pas ecrire les bounds de maniere concise.    


__all__ = ['allstat','monthlystat','detailedstat','regionalstat']
__all__.sort()

def allstat(model, obs, FIG_DIR, SCRIPT_DIR):
    """
    
    """

    from global_stat import detailedstat
    print ' -- Validation XYT (Global) --    ' 
    print 65*'-'
    result=ValidXYT(model, obs)
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


    ctdeb=cdtime.comptime(andeb,mdeb,jdeb,hdeb,0,0)
    hfin = 23
    ctfin=cdtime.comptime(anfin,mfin,jfin,hfin,0,0)
    
    tag = 'all'
    tag1= strftime('%Y%m%d',ctdeb)
    tag2= strftime('%Y%m%d',ctfin)
    print 65*'-'
    print  'La validation concerne les modeles et observations compris entre '+tag1+' jusqu au ' +tag2
    print 65*'-'

    detailedstat(result,tag,tag1,tag2,SCRIPT_DIR,FIG_DIR)
    del result
    gc.collect()

def monthlystat(model, obs, FIG_DIR, SCRIPT_DIR):
    """
    
    """
    import vacumm.misc.atime as atime
    from global_stat import detailedstat

    model = atime.monthly(model)
    obs = atime.monthly(obs) 
    gc.collect()

    print ' -- Validation XYT (Monthly) --    ' 
    print 65*'-'

    tags=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']

    for i, tag in enumerate(tags):
        model2 = model()
        obs2 = obs()
        result=ValidXYT(model2,obs2)

        detailedstat(result,tag,tag1,tag2,SCRIPT_DIR,FIG_DIR)
        del result
        del model2
        del obs2
        gc.collect()

    del model 
    del obs
    gc.collect()

def regionalstat(model, obs, FIG_DIR, SCRIPT_DIR):
    """
    - SST statistics over 5 regions in the Bay of Biscay / Channel -
    Manche Est: 3.5W/5E 48N/52.8N
    Manche Ouest: 3.5W/7W 48N/52.8N
    Bretagne Sud: 1W/5.2W 47N/48N
    Vendee: 1W/4.7W 45N/47N
    Basque: 1W/3W 43N/45N
    """

    lo_min=[-3.5, -7, -5.2, -4.7, -3]
    la_min=[48, 48, 47, 45, 43]
    lo_max=[5, -3.5, -1, -1, -1]
    la_max=[52.8, 52.8, 48, 47, 45]
    tags=['Manche Est','Manche Ouest','Bretagne Sud','Vendee','Basque']
    
    result=ValidXYT(model, obs)
    
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


    ctdeb=cdtime.comptime(andeb,mdeb,jdeb,hdeb,0,0)
    hfin = 23
    ctfin=cdtime.comptime(anfin,mfin,jfin,hfin,0,0)
    
    tag = 'regionalstat'
    tag1= strftime('%Y%m%d',ctdeb)
    tag2= strftime('%Y%m%d',ctfin)
    for i, tag in enumerate(tags):
	
        model2 = model(lon=[lo_min[i],lo_max[i]],lat=[la_min[i],la_max[i]])
        obs2 = obs(lon=[lo_min[i],lo_max[i]],lat=[la_min[i],la_max[i]])
        result=ValidXYT(model2,obs2)
        
        detailedstat(result,tag,tag1,tag2,SCRIPT_DIR,FIG_DIR)
        del result
        del model2
        del obs2
        gc.collect()
    
    

def seasonalstat(model, obs, FIG_DIR, SCRIPT_DIR):
    """
    """
    import cdutil
    from global_stat import detailedstat

    print ' -- Validation XYT (Seasonal) -- '
    
    model2 = cdutil.DJF(model)
    model2.units = model.units
    obs2 = cdutil.DJF(obs)
    obs2.units = result.obs.units
    result=ValidXYT(model2,obs2)
    tag = 'DJF'
    detailedstat(result,tag,tag1,tag2,SCRIPT_DIR,FIG_DIR)
    del result
    del model2
    del obs2
    gc.collect()

    model2 = cdutil.MAM(model)
    model2.units = model.units
    obs2 = cdutil.MAM(obs)
    obs2.units = result.obs.units
    result=ValidXYT(model2,obs2)
    tag = 'MAM'
    detailedstat(result,tag,SCRIPT_DIR,FIG_DIR)
    del result
    del model2
    del obs2
    gc.collect()

    model2 = cdutil.JJA(model)
    model2.units = model.units
    obs2 = cdutil.JJA(obs)
    obs2.units = result.obs.units
    result=ValidXYT(model2,obs2)
    tag = 'JJA'
    detailedstat(result,tag,SCRIPT_DIR,FIG_DIR)
    del result
    del model2
    del obs2
    gc.collect()

    model2 = cdutil.DJF(model)
    model2.units = model.units
    obs2 = cdutil.DJF(obs)
    obs2.units = result.obs.units
    result=ValidXYT(model2,obs2)
    tag = 'SON'
    detailedstat(result,tag,SCRIPT_DIR,FIG_DIR)
    del result
    del model2
    del obs2
    gc.collect()




    del model
    del obs
    gc.collect()



def detailedstat(result,tag,tag1,tag2,SCRIPT_DIR,FIG_DIR):
    """
    """

    tagforfilename = '_'.join(tag.split(' ')).lower()
    tagforfiledate1 = '_'.join(tag1.split(' ')).lower()    
    tagforfiledate2 = '_'.join(tag2.split(' ')).lower()
    #lon = create_lon( id='ni', lonid='lon') 

    #LONGITUDE = create_lon(loi, id='ni', attributes=dict(long_name='Longitude of each location',standard_name='longitude',units='degrees_east',valid_min='-180.',valid_max='180.',axis='X'))
    #LATITUDE = create_lat(lai, id='nj', attributes=dict(long_name='Latitude of each location',standard_name='latitude',units='degrees_north',valid_min='-90.',valid_max='90.',axis='Y'))
    
    
    # Ouverture fichier config
    config = ConfigParser.RawConfigParser()
    if tag == 'all':
        config.read(os.path.join(SCRIPT_DIR,'config.cfg'))
    elif tag == 'DJF' or tag == 'MAM' or tag == 'SON' or tag == 'JJA' or tag == 'monthly':
        config.read(os.path.join(SCRIPT_DIR,'config_s.cfg'))
    else:
        config.read(os.path.join(SCRIPT_DIR,'config_r.cfg'))    
    rep_ecriture= config.get('Output','rep_ecriture')   

    if config.get('Statistics', 'temporal_mean') == 'True': 
        # -- it works !
        
        print ' - Calcul de la moyenne temporelle pour le modele et les observations - '
	print 25*'-'

        result.temporal_mean()  
        
        post = np.array([2,3,4,5,6,7,8,9,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,20.5,21,22,23,24])        
        #pos = (post-post[0])/(post[-1]-post[0])
        
        
        # - normalisation
        norm = StepsNorm(post)        
        
        #pos=np.linspace(0,1,len(post))
        
        #r1=np.array([212,191,165,153,97,83,65,48,34,12,8,24,42,56,75,90,105,122,132,134,175,216,250,252,252,252,252,245,231,213,200,185,167,152,136])       
        r1=np.array([212,191,165,153,97,83,65,48,34,12,8,24,42,56,75,90,105,122,132,134,175,216,250,252,252,252,252,245,231,213,200,185,167,152])       
        r1 = r1/255.
        
        #v1=np.array([154,117,95,66,29,46,61,78,94,114,134,150,165,182,199,216,231,244,250,250,250,250,247,231,207,174,141,118,101,87,70,55,40,22,6])  
        v1=np.array([154,117,95,66,29,46,61,78,94,114,134,150,165,182,199,216,231,244,250,250,250,250,247,231,207,174,141,118,101,87,70,55,40,22])  
        v1 = v1/255.
        
        #b1=np.array([229,213,208,209,159,173,119,205,222,242,252,252,252,252,252,252,252,252,228,155,132,132,126,105,82,49,15,4,4,4,4,4,4,4,4])      
        b1=np.array([229,213,208,209,159,173,119,205,222,242,252,252,252,252,252,252,252,252,228,155,132,132,126,105,82,49,15,4,4,4,4,4,4,4])      
        b1 = b1/255.
        
        colors = []
        for icolor in np.arange(len(r1)):
            colors.append(((r1[icolor],v1[icolor],b1[icolor]),norm.positions[icolor]))
        
        cmap_previmer = cmap_custom(colors,name='cmap_previmer',ncol=len(colors)-1)
        
	
	
	kwplot = dict(vmin=result.obs.temp_mean.min(), vmax=result.obs.temp_mean.max(), nmax = 30, colorbar_extend='both')
	#kwplot = dict(cmap=cmap_previmer, vmin=2, vmax=24, levels=post, norm=norm, colorbar_extend='both')
        P.figure()
        map(result.model.temp_mean, title='Mean modelled Sea Surface Temperature',  show=False,  clabel_hide=True, **kwplot)
        savefigs(FIG_DIR+'/model_temporal_mean_'+ tagforfilename +'_'+tagforfiledate1+'_' +tagforfiledate2)
        P.close()
        P.figure()
        map(result.obs.temp_mean, title='Mean observed Sea Surface Temperature',  show=False,  clabel_hide=True, **kwplot)
        savefigs(FIG_DIR+'/obs_temporal_mean_' + tagforfilename+'_' +tagforfiledate1 +'_'+tagforfiledate2)
        P.close()
        
	#Ecriture du fichier netcdf
	#print result.model.temp_mean
	result.model.temp_mean = cdms2.createVariable(result.model.temp_mean, typecode='f',id='TEMP', attributes=dict(long_name='Sea Surface Temp',standard_name='sst',units='degres celsius',valid_min='-20.',valid_max='30.'))
	result.obs.temp_mean=cdms2.createVariable(result.obs.temp_mean, typecode='f',id='TEMP', attributes=dict(long_name='Sea Surface Temp',standard_name='sst',units='degres celsius',valid_min='-20.',valid_max='30.'))
	write_nc_cf( 'model_obs_'+tagforfilename+'_temp_mean_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',True,result.model.temp_mean,result.obs.temp_mean)
	print 65*'-'
	print 'Temperature mean stat written to'+rep_ecriture+'/'+tagforfiledate1 +'_'+tagforfiledate2+'.nc'

	
    if config.get('Statistics', 'spatial_mean') == 'True': 
        # -- it works !
        print 65*'-'
        print ' - Calcul de la moyenne spatiale pour le modele et les observations - '
        result.spatial_mean()
        
        print 65*'-'
        print ' - Calcul de l ecart type pour le modele et les observations pour chaque pas de temps - '
        result.spatial_std()
        # Masque moyennes egales a zero.
        #result.model.spa_mean = MV2.masked_equal(result.model.spa_mean,0)
        result.model.spa_mean[result.model.spa_mean.data<0.1]=np.NaN
        result.obs.spa_mean[result.obs.spa_mean.data<0.1]=np.NaN
        
        # Trace des resultats
        fig = figure()
        # Creation d'un axe
        ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        # Vecteur temps pour l'axe des abscisses
        x = datestr2num(strtime(result.obs.spa_mean.getTime()))
        # Titre
        title('Spatial average over the domain')
        # 
        ylabel('Temperature (degC)')
        ax.plot_date(x,result.obs.spa_mean.data, linestyle='-', color='r', marker='None',xdate=True)
        ax.plot_date(x,result.model.spa_mean.data, linestyle='-', color='b', marker='None',xdate=True)
        
        # ATTENTION NaN ... Faire par bloc de donnees non NaN ...        
        # Identification des blocs
        index = np.where(~np.isnan(result.obs.spa_mean.data))
        di = np.diff(index[0][:])
        ifin=index[0][di>1]
        ideb=index[0][di>1]+di[di>1]
        # On entame le processing de l'image seulement si les données ne sont pas vides.

	if len(ideb)!=0:
	      if (index[0][0]!=ideb[0]) or (index[0][-1]!=ideb[-1]):
		  ideb = np.append(index[0][0],ideb)
		  ifin = np.append(ifin,index[0][-1])
		  ifin[ifin<len(result.obs.spa_mean.data)-1]=ifin[ifin<len(result.obs.spa_mean.data)-1]+1
	    ## +1 pour les appels dans les tableaux
	else :
	    ideb= np.append(ideb,0)
	    ifin= np.append(ifin,len(di))
	    
	index = np.where(~np.isnan(result.model.spa_mean.data))
	di = np.diff(index[0][:])
	mifin=index[0][di>1]
	mideb=index[0][di>1]+di[di>1]
	if len(mideb)!=0: 
	    if (index[0][0]!=mideb[0]) or (index[0][-1]!=mideb[-1]) :
		mideb = np.append(index[0][0],mideb)
		mifin = np.append(mifin,index[0][-1])
		mifin[mifin<len(result.model.spa_mean.data)-1]=mifin[mifin<len(result.model.spa_mean.data)-1]+1    
	    ## +1 pour les appels dans les tableaux
	else:
	    mideb= np.append(mideb,0)
	    mifin= np.append(mifin,len(di))                        
	    ## Boucle sur les differents blocs
	for ibloc in np.arange(len(ideb)):

	    #if len(ifin)==len(ideb):  
	    if ifin[ibloc]-ideb[ibloc] > 1:
				    
		    ## Generation du patch
		zobs=np.zeros(2*len(result.obs.spa_mean.data[ideb[ibloc]:ifin[ibloc]])+1)
		interm1=result.obs.spa_mean.data[ideb[ibloc]:ifin[ibloc]]
		interm2=result.obs.spa_std.data[ideb[ibloc]:ifin[ibloc]]
		    
		zobs[:-1]=np.concatenate([result.obs.spa_mean.data[ideb[ibloc]:ifin[ibloc]]-result.obs.spa_std.data[ideb[ibloc]:ifin[ibloc]],
		interm1[::-1]+interm2[::-1]])
		    ## Dernier point = premier point
		zobs[-1]=result.obs.spa_mean.data[ideb[ibloc]]-result.obs.spa_std.data[ideb[ibloc]]
		    
		    
		    ## Axe des x
		xx=np.zeros(2*len(result.obs.spa_mean.data[ideb[ibloc]:ifin[ibloc]])+1)
		intermt=x[ideb[ibloc]:ifin[ibloc]]
		xx[:-1]=np.concatenate([x[ideb[ibloc]:ifin[ibloc]],intermt[::-1]])
		xx[-1]=x[ideb[ibloc]]
		    
		    ## Creation du tuple
		Path = mpath.Path
		pathdata = []
		for i,ii in enumerate(zobs):
		    if i==0:
			pathdata.append((Path.MOVETO,(tuple([xx[i],zobs[i]]))))
		    elif i==(len(zobs)-1):
			pathdata.append((Path.CLOSEPOLY,(tuple([xx[i],zobs[i]]))))             
		    else:
			pathdata.append((Path.LINETO,(tuple([xx[i],zobs[i]]))))
			
		    
		ocodes, overts = zip(*pathdata)
		opath = mpath.Path(overts, ocodes)
		opatch = mpatches.PathPatch(opath, facecolor='red', edgecolor='red', alpha=0.3)
		ax.add_patch(opatch)
		    
		  
	    
	for ibloc in np.arange(len(mideb)):
	    #if len(ifin)==len(ideb):
		    
	    if mifin[ibloc]-mideb[ibloc] > 1:
		
		    
		    ## Generation du patch
		zmodel=np.zeros(2*len(result.model.spa_mean.data[mideb[ibloc]:mifin[ibloc]])+1)
		interm1=result.model.spa_mean.data[mideb[ibloc]:mifin[ibloc]]
		interm2=result.model.spa_std.data[mideb[ibloc]:mifin[ibloc]]
		zmodel[:-1]=np.concatenate([result.model.spa_mean.data[mideb[ibloc]:mifin[ibloc]]-result.model.spa_std.data[mideb[ibloc]:mifin[ibloc]],
					    interm1[::-1]+interm2[::-1]])
		    ## Dernier point = premier point
		zmodel[-1]=result.model.spa_mean.data[mideb[ibloc]]-result.model.spa_std.data[mideb[ibloc]]
		    
		    ## Axe des x
		xx=np.zeros(2*len(result.model.spa_mean.data[mideb[ibloc]:mifin[ibloc]])+1)
		intermt=x[mideb[ibloc]:mifin[ibloc]]
		xx[:-1]=np.concatenate([x[mideb[ibloc]:mifin[ibloc]],intermt[::-1]])
		xx[-1]=x[mideb[ibloc]]
		    
		    ## Creation du tuple
		Path = mpath.Path
		pathmodel = []
		for i,ii in enumerate(zmodel):
		    if i==0:
			pathmodel.append((Path.MOVETO,(tuple([xx[i],zmodel[i]]))))
		    elif i==(len(zmodel)-1):
			pathmodel.append((Path.CLOSEPOLY,(tuple([xx[i],zmodel[i]]))))                
		    else:
			pathmodel.append((Path.LINETO,(tuple([xx[i],zmodel[i]]))))
			
		    
		mocodes, moverts = zip(*pathmodel)
		mopath = mpath.Path(moverts, mocodes)
		mopatch = mpatches.PathPatch(mopath, facecolor='blue', edgecolor='blue', alpha=0.3)
		ax.add_patch(mopatch)
	    
 
	ax.grid()
	ylim( (2,25) )
	ax.legend(('Observations',  'Model'),  loc = 'upper right',  shadow = True,  fancybox=True)
	fig.autofmt_xdate()
	savefigs(FIG_DIR+'/result_spatial_statmean_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)
	close()
	
	  #Ecriture du fichier netcdf
	result.model.spa_mean = cdms2.createVariable(result.model.spa_mean, typecode='f',id='SST_model', attributes=dict(long_name='Sea Surface Temp',standard_name='sst',units='degres celsius',valid_min='-20.',valid_max='30.'))
	result.obs.spa_mean=cdms2.createVariable(result.obs.spa_mean, typecode='f',id='SST_obs', attributes=dict(long_name='Sea Surface Temp',standard_name='sst',units='degres celsius',valid_min='-20.',valid_max='30.'))
	write_nc_cf('model_obs_'+tagforfilename+'_spa_mean_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',True,result.model.spa_mean,result.obs.spa_mean)
	print 65*'-'
	print 'Spatial mean stat written to'+rep_ecriture+'/'+tagforfilename+tagforfiledate1 +'_'+tagforfiledate2+'.nc'
	 # Nettoyage
	gc.collect()
        
    if config.get('Statistics', 'temporal_std') == 'True':    
        # -- it works !
        print 65*'-'
        print ' Calcul de la std en chaque point geographique'
        result.temporal_std()
        print 65*'-'
     
       
	kwplot = dict(vmin=result.obs.temp_std.min(), vmax=result.obs.temp_std.max(), nmax = 30, colorbar_extend='max')
	P.figure()
	map(result.model.temp_std, title='Model - std(SST(t))',  show=False,  clabel_hide=True, **kwplot)
	savefigs(FIG_DIR+'/model_temporal_std_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)
	P.close()
	P.figure()
	map(result.obs.temp_std, title='Observation - std(SST(t))', show=False,  clabel_hide=True, **kwplot)
	savefigs(FIG_DIR+'/obs_temporal_std_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)
	P.close()
	
	  #Ecriture du fichier netcdf
	result.model.temp_std = cdms2.createVariable(result.model.temp_std, typecode='f',id='SST_model', attributes=dict(long_name='Temporal Deviation of Modelled Sea Surface Temperature',standard_name='spa_std_model_sst',units='degres celisus',valid_min='-30.',valid_max='30.'))
	result.obs.temp_std=cdms2.createVariable(result.obs.temp_std, typecode='f',id='SST_obs', attributes=dict(long_name='Temporal Deviation of Observed Sea Surface Temperature',standard_name='spa_std_obs_sst',units='degres celsius',valid_min='-30.',valid_max='30.'))
	write_nc_cf( 'model_obs_'+tagforfilename+'_temp_std_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',True,result.model.temp_std,result.obs.temp_std)
	print 65*'-'
	print 'Temperature std stat written to'+rep_ecriture+'/'+tagforfilename+tagforfiledate1 +'_'+tagforfiledate2+'.nc'
	print 65*'-'
	  # Nettoyage
        result.obs.temp_std=[]
        result.model.temp_std=[]
        gc.collect()


	
    
    if config.get('Statistics', 'spatial_std') == 'True': 
        # -- it works !
        print 25*'-'
        print ' - Calcul de l ecart type pour le modele et les observations pour chaque pas de temps - '
        print 25*'-'
        result.spatial_std()       
        
        fig = figure()
	ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
	x = datestr2num(strtime(result.obs.spa_std.getTime()))
	title('Spatial STD over the domain')
	ylabel(result.obs.units)
	ax.plot_date(x,result.obs.spa_std.data, marker='o', color='r',  markerfacecolor='r',xdate=True)
	ax.plot_date(x,result.model.spa_std.data, marker='o', color='b',  markerfacecolor='b',xdate=True)  
	ax.legend(('Observations',  'Model'),  loc = 'upper right',  shadow = True,  fancybox=True)
	fig.autofmt_xdate()
	savefigs(FIG_DIR+'/result_spatial_statstd_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)
	close()
        
        
	    #Ecriture du fichier netcdf
        result.model.spa_std = cdms2.createVariable(result.model.spa_std, typecode='f',id='SST_model', attributes=dict(long_name='Standard Deviation of Modelled Sea Surface Temperature',standard_name='std_model_sst',units='D_C',valid_min='-30.',valid_max='30.'))
	result.obs.spa_std=cdms2.createVariable(result.obs.spa_std, typecode='f',id='SST_obs', attributes=dict(long_name='Standard Deviation of Observed Sea Surface Temperature',standard_name='std_obs_sst',units='D_C',valid_min='-30.',valid_max='30.'))
	write_nc_cf('model_obs_'+tagforfilename+'_spa_std_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',True,result.model.spa_std,result.obs.spa_std)
        print 'Temperature spatial std stat written to'+rep_ecriture+'/'+tagforfilename+tagforfiledate1 +'_'+tagforfiledate2+'.nc'
                
	    # Nettoyage
        result.obs.spa_std=[]
        result.model.spa_std=[]
        gc.collect()
        
    if config.get('Statistics', 'evol_biais') == 'True':
        # -- it works !
        print 25*'-'
        print ' - Calcul de l''evolution du biais moyen entre le modele et les observations - '  
        print 25*'-'
        if config.get('Statistics', 'spatial_mean') == 'False': # Pour ne pas faire le calcul 2 fois
	    result.obs.spatial_mean()
	    result.model.spatial_mean()
        
        
	    # Trace des resultats
	fig = figure()
	ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
	x = datestr2num(strtime(result.obs.spa_mean.getTime()))
	title('Mean bias over the domain (<E> = <OBS> - <MODEL>)')
	ylabel(result.obs.units)
	ax.plot_date(x,result.obs.spa_mean.data-result.model.spa_mean.data, marker='o', markeredgecolor='b',  markerfacecolor='b',xdate=True)
	ax.plot_date(x,np.zeros(result.obs.spa_mean.data.shape), 'k:',xdate=True)
        
	fig.autofmt_xdate()
	savefigs(FIG_DIR+'/result_evol_bias_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)
	close()
	
	evol_biais=cdms2.createVariable(result.obs.spa_mean.data-result.model.spa_mean.data ,typecode='f',id='SST_bias', attributes=dict(long_name='Temporal Evolution of the Bias between Observed and Modelled Sea Surface Temperature',standard_name='bias_sst',units='D_C',valid_min='-30.',valid_max='30.'))
	write_nc_cf('model_obs_'+tagforfilename+'_evol_biais_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',varin1=evol_biais )
        print 'Evolution of Bias stat written to'+rep_ecriture+'/'+tagforfilename+tagforfiledate1 +'_'+tagforfiledate2+'.nc'

        # Nettoyage
        result.obs.spa_mean=[]
        result.model.spa_mean=[]
        gc.collect()
        
    if config.get('Statistics', 'mean_biais') == 'True':
        # -- it works !
        print 25*'-'
        print ' - Calcul du biais moyen entre le modele et les observations - '
        print 25*'-'
        if config.get('Statistics', 'temporal_mean') == 'False': # Pour ne pas faire le calcul 2 fois
	    result.model.temporal_mean()       
	    result.obs.temporal_mean()
        interm = result.obs.temp_mean-result.model.temp_mean

        bias_list=np.array([-10,-5,-4,-3,-2,-1,-0.5,0.5,1,2,3,4,5,10])        
        # - normalisation
        norm = StepsNorm(bias_list)        
        
        #map(result.model.temp_mean, res='c', show=False)        
        P.figure()
        #kwplot = dict(vmin=-abs(interm).max(), vmax=abs(interm).max(), nmax=20)
        kwplot = dict(cmap=cmap_magic(bias_list,anomaly=True), vmin=bias_list.min(), vmax=bias_list.max(), levels=bias_list, norm=norm, colorbar_extend='both')
        map(interm, title='Mean bias over the domain (<E> = <OBS> - <MODEL>)',  show=False,  clabel_hide=True,  **kwplot)
        savefigs(FIG_DIR+'/result_mean_bias_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)
        P.close()  
        
        #Ecriture du fichier netcdf   
	mean_biais=cdms2.createVariable(interm ,typecode='f',id='SST_bias', attributes=dict(long_name='Mean Bias between Observed and Modelled (<Obs>-<Mod>) Sea Surface Temperature',standard_name='bias_sst',units='D_C',valid_min='-30.',valid_max='30.'))
	write_nc_cf('model_obs_'+tagforfilename+'_mean_biais_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',varin1=mean_biais)
        print 'Mean Bias stat written to'+rep_ecriture+'/'+tagforfilename+tagforfiledate1 +'_'+tagforfiledate2+'.nc'
	
        # Nettoyage
        result.obs.temp_mean=[]
        result.model.temp_mean=[]   
        gc.collect()
        
    if config.get('Statistics', 'extreme_bias') == 'True':        
        print 25*'-'
        print ' - Calcul du biais minimum et maximum '        
        print 25*'-'
        result.bias_extrema()
        
	P.figure()
	kwplot = dict(vmin=-abs(result.biasmin).max(), vmax=abs(result.biasmin).max())
	map(result.biasmin, title='Minimum bias',  show=False,  clabel_hide=True,  **kwplot)
	savefigs(FIG_DIR+'/result_min_bias_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)        
	P.close()       
        
	P.figure()
	kwplot = dict(vmin=-abs(result.biasmax).max(), vmax=abs(result.biasmax).max())
	map(result.biasmax, title='Maximum bias',  show=False,  clabel_hide=True,  **kwplot)
	savefigs(FIG_DIR+'/result_max_bias_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)        
	P.close()  
	
	#Ecriture du fichier netcdf	ne marche pas
	#print result.bias_extrema
        result.biasmax=cdms2.createVariable(result.biasmax,typecode='f',id='SST_biasmax', attributes=dict(long_name='Maximum Bias between Observed and Modelled (<Obs>-<Mod>) Sea Surface Temperature',standard_name='max_bias',units='D_C',valid_min='-30.',valid_max='30.'))
        result.biasmin=cdms2.createVariable(result.biasmin,typecode='f',id='SST_biasmin', attributes=dict(long_name='Minimum Bias between Observed and Modelled (<Obs>-<Mod>) Sea Surface Temperature',standard_name='min_bias',units='D_C',valid_min='-30.',valid_max='30.'))
	write_nc_cf('model_obs_'+tagforfilename+'_extreme_biais_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',True,result.biasmin,result.biasmax)
        print 'Extreme Bias stat written to'+rep_ecriture+'/'+tagforfilename+tagforfiledate1 +'_'+tagforfiledate2+'.nc'
	
        # Nettoyage
        result.biasmin = []
        result.biasmax = []
        gc.collect()
 
    if config.get('Statistics', 'systematic_bias') == 'True':     
        print 25*'-'
        print ' - Calcul du biais systematique '
        print 25*'-'
        result.systematic_bias()

        P.figure()
	kwplot = dict(vmin=-1, vmax=1,  levels=[-1, -0.1,0.1,  1], colorbar_ticks=[1,0,-1])
	map(result.biassyst, title='Systematic bias\n (+1: model underestimation, -1: model overestimation)',  show=False,  clabel_hide=True,  **kwplot)
	savefigs(FIG_DIR+'/result_syst_bias_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)        
	P.close()
        
        #Ecriture du fichier netcdf
        
	sys_bias=cdms2.createVariable(result.biassyst, typecode='f',id='SYS_BIAS', attributes=dict(long_name='Systemic bias between Observed and Modelled Sea Surface Temperature',standard_name='bias_systematic',units='D_C',valid_min='-100.',valid_max='100.'))
	write_nc_cf('model_obs_'+tagforfilename+'_systematic_bias_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',varin1=sys_bias )
	print 'Systematic bias stat written to'+rep_ecriture+'/'+tagforfilename+tagforfiledate1 +'_'+tagforfiledate2+'.nc'
        # Nettoyage
        result.biassyst = []
        sys_bias=[]
        gc.collect()

        
    if config.get('Statistics', 'extrema') == 'True':     
        # -- it works !
        print 25*'-'
        print ' - Extractions des minima et maxima - '
        print 25*'-'
        result.extrema()
        
	
    if config.get('Statistics', 'obs_coverage') == 'True':    
        # -- it works !
        print 25*'-'
        print ' - Calcul du taux de couverture des observations en chaque point geographique - '
        print 25*'-'
        result.obs_coverage()        

	P.figure()
	map(result.obs_cov,  title='Observation coverage during the period', show=False,  clabel_hide=True)
	savefigs(FIG_DIR+'/result_obs_coverage_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)
	P.close()
    
	#Ecriture du fichier netcdf
	result.obs_cov=cdms2.createVariable(result.obs_cov, typecode='f',id='Obs_Cov', attributes=dict(long_name='Observation coverage',standard_name='observation_coverage',units='%',valid_min='0.',valid_max='100.'))
	write_nc_cf('obs_'+tagforfilename+'_obs_coverage_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',varin1= result.obs_cov )
	print 'Observation coverage stat written to'+rep_ecriture+'/'+tagforfilename+tagforfiledate1 +'_'+tagforfiledate2+'.nc'
	
        # Nettoyage
        result.obs_cov = []
        gc.collect()

    
    if config.get('Statistics', 'obs_spatialcoverage') == 'True':    
        # -- it works !
        print 25*'-'
        print ' - Calcul du taux de couverture des observations en chaque pas de temps - '
        print 25*'-'
        result.spatial_obs_coverage()        
        
       
        fig = figure()
	ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
	x = datestr2num(strtime(result.obs_spacov.getTime()))
	title('Observation coverage following the time period')
	ylabel('Coverage (nb)')
	ax.plot_date(x,result.obs_spacov.data,linestyle='-', marker='o', color='r',  markerfacecolor='b',xdate=True)
	ylim()
	fig.autofmt_xdate()
	savefigs(FIG_DIR+'/result_obs_spatial_coverage_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)
	close()
        
        #Ecriture du fichier netcdf
        result.obs_spacov=cdms2.createVariable(result.obs_spacov, typecode='f',id='Obs_Cov', attributes=dict(long_name='Spatial observation coverage',standard_name='spatial_observation_coverage',units='%',valid_min='0.',valid_max='100.'))
	write_nc_cf('obs_'+tagforfilename+'_spatial_obs_coverage_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',varin1=result.obs_spacov )
	print 'Spatial observation coverage stat written to'+rep_ecriture+'/'+tagforfilename+tagforfiledate1 +'_'+tagforfiledate2+'.nc'
        # Nettoyage
        result.obs_spacov = []
        gc.collect()
     
    
    
    if config.get('Statistics', 'temporal_rms') == 'True':    
        # -- it works !
        print 25*'-'
        print ' - Calcul de la RMS temporelle - '
        print 25*'-'
        result.temporal_rms()
         
        P.figure()
	map(result.temp_rms, title='RMS\n (uncentered and biased)', show=False,  clabel_hide=True)
	savefigs(FIG_DIR+'/result_temporal_rms_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)
	P.close()
        
        #Ecriture du fichier netcdf
        result.temp_rms=cdms2.createVariable(result.temp_rms, typecode='f',id='Temporal_Rms', attributes=dict(long_name='Temporal RMS',standard_name='temporal_rms',units='D_C',valid_min='-100.',valid_max='100.'))
	write_nc_cf('model_obs_'+tagforfilename+'_temporal_rms_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',varin1=result.temp_rms )
	print 'Temporal Rms stat written to'+rep_ecriture+'/'+tagforfilename+tagforfiledate1 +'_'+tagforfiledate2+'.nc'
        # Nettoyage
        result.temp_rms = []
        gc.collect()
    
        
    if config.get('Statistics', 'temporal_rmsc') == 'True':    
        # -- it works !
        print 25*'-'
        print ' - Calcul de la RMS temporelle  centree- '
        print 25*'-'
        result.temporal_rmsc()
       
        P.figure()
	map(result.temp_rmsc, title='RMS\n (centered and biased)', show=False,  clabel_hide=True)
	savefigs(FIG_DIR+'/result_temporal_rmsc_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)
	P.close()
       
        #Ecriture du fichier netcdf
        result.temp_rmsc=cdms2.createVariable(result.temp_rmsc, typecode='f',id='Temporal_Rmsc', attributes=dict(long_name='Temporal RMSC',standard_name='temporal_rmsc',units='D_C',valid_min='-100.',valid_max='100.'))
	write_nc_cf('model_obs_'+tagforfilename+'_temporal_rmsc_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',varin1=result.temp_rmsc )
	print 'Temporal Rmsc stat written to'+rep_ecriture+'/'+tagforfilename+tagforfiledate1 +'_'+tagforfiledate2+'.nc'
        # Nettoyage
        result.temp_rmsc = []
        gc.collect()
    
 
    if config.get('Statistics', 'spatial_rms') == 'True': 
        # -- it works !
        print 25*'-'
        print ' - Calcul de la RMS entre le modele et les observations pour chaque pas de temps - '
        print 25*'-'
        result.spatial_rms()       
        
	fig = figure()
	ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
	x = datestr2num(strtime(result.spa_rms.getTime()))
	title('Spatial RMS over the domain')
	ylabel(result.obs.units)
	ax.plot_date(x,result.spa_rms.data, marker='o', color='r',  markerfacecolor='r',xdate=True)
	fig.autofmt_xdate()
	savefigs(FIG_DIR+'/result_spatial_stat_rms_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)
	close()
        
        #Ecriture du fichier netcdf
        result.spa_rms=cdms2.createVariable(result.spa_rms, typecode='f',id='Spatial_Rms', attributes=dict(long_name='Spatial RMS',standard_name='spatial_rms',units='D_C',valid_min='-100.',valid_max='100.'))
	write_nc_cf('model_obs_'+tagforfilename+'_spatial_rms_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',varin1=result.spa_rms )
	print 'Spatial Rms stat written to'+rep_ecriture+'/'+tagforfilename+tagforfiledate1 +'_'+tagforfiledate2+'.nc'        
        # Nettoyage
        result.spa_rms = []
        gc.collect()
    
        
    if config.get('Statistics', 'spatial_rmsc') == 'True': 
        # -- it works !
        print 25*'-'
        print ' - Calcul de la RMS centree entre le modele et les observations pour chaque pas de temps - '
	print 25*'-' 
        result.spatial_rmsc()       
        
          
        fig = figure()
	ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
	x = datestr2num(strtime(result.spa_rmsc.getTime()))
	title('Spatial centered RMS over the domain')
	ylabel(result.obs.units)
	ax.plot_date(x,result.spa_rmsc.data, marker='o', color='r',  markerfacecolor='r',xdate=True)
	fig.autofmt_xdate()
	savefigs(FIG_DIR+'/result_spatial_statrmsc_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)
	close()
        
        #Ecriture du fichier netcdf
        result.spa_rmsc=cdms2.createVariable(result.spa_rmsc, typecode='f',id='Spatial_Rmsc', attributes=dict(long_name='Spatial RMSC',standard_name='spatial_rmsc',units='D_C',valid_min='0.',valid_max='30.'))
	write_nc_cf('model_obs_'+tagforfilename+'_spatial_rmsc_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',varin1=result.spa_rmsc)
	print 'Spatial Rmsc stat written to'+rep_ecriture+'/'+tagforfilename+tagforfiledate1 +'_'+tagforfiledate2+'.nc'            
        # Nettoyage
        result.spa_rmsc = []
        gc.collect()
    
    
    if config.get('Statistics', 'temporal_corr') == 'True':    
        # -- it works !
        print 10*'-'
        print ' - Calcul de la Correlation temporelle - '
        result.temporal_corr()

	P.figure()
	map(result.temp_corr, title='Time series correlation', show=False,  clabel_hide=True)
	savefigs(FIG_DIR+'/result_temporal_corr_'+ tagforfilename+'_'+tagforfiledate1+'_' +tagforfiledate2)
	P.close()
	
	#Ecriture du fichier netcdf
        result.temp_corr=cdms2.createVariable(result.temp_corr, typecode='f',id='Temp_corr', attributes=dict(long_name='Temporal_correlation',standard_name='temporal_corr',units='%',valid_min='0.',valid_max='100.'))
	write_nc_cf('model_Obs_'+tagforfilename+'_temp_corr_' +tagforfiledate1 +'_'+tagforfiledate2+'.nc',varin1=result.temp_corr)
	print 'Temporal Correlation stat written to'+rep_ecriture+'/'+tagforfilename+tagforfiledate1 +'_'+tagforfiledate2+'.nc'       	
        # Nettoyage
        result.temp_corr = []
        gc.collect()
    
    if config.get('Report', 'rep_html') == 'True':
	print 86*'*'
	print ' -- Generation du Rapport de Validation --    ' 
	print 86*'*'
	#title = "Test(1)"
	#header = ""
	#footer = ""
	config2 = ConfigParser.RawConfigParser()
	config2.read(os.path.join(SCRIPT_DIR,'config.cfg'))
	rr = config2.get('Observations', 'product')
        rr = rr.replace(' ', '_')    
	ttforintro = '<u>Observations</u>: %(#)s <br /> <u>Model</u>: %(##)s %(###)s <br /><hr>' %{'#':result.obs.long_name, '##':config2.get('Model Description', 'name') , '###':result.model.long_name}
	  
	
	if config2.get('Statistics', 'extrema') == 'True' and config2.get('Statistics', 'to_do') == 'True':
	    mimamodel = '=> Modelled minimum and maximum: %(#)6.3f / %(##)6.3f %(u)s' %{'#':result.model.extrema[0], '##':result.model.extrema[1], 'u':result.model.units}
	    mimaobs = '=> Observed minimum and maximum: %(###)6.3f / %(####)6.3f %(u)s'  %{'###':result.obs.extrema[0],'####':result.obs.extrema[1], 'u':result.obs.units}
	else:
	    mimamodel = ''
	    mimaobs = ''
		
	images_control = glob.glob(os.path.join(FIG_DIR,'control_*'+tagforfiledate1 +'_'+tagforfiledate2+'.png'))
	for i,pat in enumerate(images_control):
	    (pa,fil)=os.path.split(images_control[i])
	    images_control[i] = fil   
	
	images_results = glob.glob(os.path.join(FIG_DIR,'*_mean'+tagforfilename+'_'+tagforfiledate1 +'_'+tagforfiledate2+'.png'))
	images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_std'+tagforfilename+'_'+tagforfiledate1 +'_'+tagforfiledate2+'.png')))
	images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_statmean'+tagforfilename+'_'+tagforfiledate1 +'_'+tagforfiledate2+'.png')))
	images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_statstd'+tagforfilename+'_'+tagforfiledate1 +'_'+tagforfiledate2+'.png')))
	images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_mean_bias'+tagforfilename+'_'+tagforfiledate1 +'_'+tagforfiledate2+'.png')))
	images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_syst_bias'+tagforfilename+'_'+tagforfiledate1 +'_'+tagforfiledate2+'.png')))
	images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_min_bias'+tagforfilename+'_'+tagforfiledate1 +'_'+tagforfiledate2+'.png')))
	images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_max_bias'+tagforfilename+'_'+tagforfiledate1 +'_'+tagforfiledate2+'.png')))
	images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_evol_bias'+tagforfilename+'_'+tagforfiledate1 +'_'+tagforfiledate2+'.png')))
	images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_corr'+tagforfilename+'_'+tagforfiledate1 +'_'+tagforfiledate2+'.png')))
	images_results.extend(glob.glob(os.path.join(FIG_DIR,'*rms'+tagforfilename+'_'+tagforfiledate1 +'_'+tagforfiledate2+'.png')))
	images_results.extend(glob.glob(os.path.join(FIG_DIR,'*_coverage'+tagforfilename+'_'+tagforfiledate1 +'_'+tagforfiledate2+'.png')))
	#images_results.extend(glob.glob(os.path.join(FIG_DIR,'result_*.png')))

	for i,pat in enumerate(images_results):
	    (pa,fil)=os.path.split(images_results[i])
	    images_results[i] = fil


	
	if images_control == []:
	    images_control = None
	    
	if images_results == []:
	    images_results = None
	
	html_tools.simplereporthtml(images_results=images_results, images_control=images_control,  mimamodel=mimamodel, mimaobs=mimaobs,  intro=ttforintro,  
				    file_out=FIG_DIR+'/'+rr+'_'+tagforfilename+'.html')
	
	
    