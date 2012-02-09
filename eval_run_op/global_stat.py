# -*- coding: utf8 -*-

"""
Global statistics

   Statistics on the whole dataset/run ... 
"""

import os,sys,shutil,ConfigParser, gc
#import cdms2



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
from vacumm.misc.atime import add, strtime, ch_units, are_same_units, comptime
from vacumm.misc.axes import create_time,  set_order
from vacumm.misc.grid.regridding import regrid1d, regrid2d
from vacumm.misc.io import ncread_best_estimate
#from vacumm.misc.grid.regridding import regrid1d,  regrid2d,  interp2d
import glob
from vacumm.markup import html_tools
from vacumm.misc.color import cmap_custom, StepsNorm, cmap_magic

__all__ = ['allstat','monthlystat','detailedstat','regionalstat']
__all__.sort()

def allstat(model, obs, FIG_DIR, SCRIPT_DIR):
    """
    
    """

    from global_stat import detailedstat
    print ' -- Validation XYT (Global) --    ' 
    print 65*'-'
    result=ValidXYT(model, obs)
    #del model 
    #del obs
    #gc.collect()

    tag = 'all'

    detailedstat(result,tag, SCRIPT_DIR, FIG_DIR)
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

        detailedstat(result,tag,SCRIPT_DIR,FIG_DIR)
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

    for i, tag in enumerate(tags):
        model2 = model(lon=[lo_min[i],lo_max[i]],lat=[la_min[i],la_max[i]])
        obs2 = obs(lon=[lo_min[i],lo_max[i]],lat=[la_min[i],la_max[i]])
        result=ValidXYT(model2,obs2)
        
        detailedstat(result,tag,SCRIPT_DIR,FIG_DIR)
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
    print 65*'-'
    
    model2 = cdutil.DJF(model)
    model2.units = model.units
    obs2 = cdutil.DJF(obs)
    obs2.units = result.obs.units
    result=ValidXYT(model2,obs2)
    tag = 'DJF'
    detailedstat(result,tag,SCRIPT_DIR,FIG_DIR)
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



def detailedstat(result,tag,SCRIPT_DIR,FIG_DIR):
    """
    """

    tagforfilename = '_'.join(tag.split(' ')).lower()

    # Ouverture fichier config
    config = ConfigParser.RawConfigParser()
    if tag == 'all':
        config.read(os.path.join(SCRIPT_DIR,'config.cfg'))
    elif tag == 'DJF' or tag == 'MAM' or tag == 'SON' or tag == 'JJA' or tag == 'monthly':
        config.read(os.path.join(SCRIPT_DIR,'config_s.cfg'))
    else:
        config.read(os.path.join(SCRIPT_DIR,'config_r.cfg'))

    if config.get('Statistics', 'delta_t') == 'True':
        # Quantification du delta de temperature de surface entre le talus et le plateau: T(plateau) - T(talus)
        # ==> Evolution temporelle
        
        # Premier couple: i=235, j=168 (5.5W/47N) --- i=267, j=168 (3.7W/47N)
        # (-1) car tableaux commencent a zero.
        i1=266
        j1=167
        i2=234
        j2=167
        serie1o = result.obs[:,i1,j1]-result.obs[:,i2,j2] 
        serie1m = result.model[:,i1,j1]-result.model[:,i2,j2]

        lo = result.obs.getLongitude().getValue()
        la = result.obs.getLatitude().getValue()


        # Trace des resultats
        kwploto=dict(marker = 'o', color='r', ls=' ')
        kwplotm=dict(marker = 'o', color='b', ls=' ')
        P.figure()
        curve2(serie1o,**kwploto)
        curve2(serie1m,**kwplotm)

        #fig = figure()
        #ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        #x = datestr2num(strtime(result.obs.getTime()))
        title('SST ( %(#)3.2fE / %(##)4.2fN ) - SST ( %(###)3.2fE / %(####)4.2fN )' %{'#':lo[i1], '##':la[j1], '###':lo[i2], '####':la[j2]})
        #ylabel(result.obs.units)
        #ax.plot_date(x,serie1o.data, color='r', ls='-', marker=' ')
        #ax.plot_date(x,serie1m.data,  color='b', ls='-', marker=' ')
        #ax.legend(('Observations',  'Model'),  loc = 'upper right',  shadow = True,  fancybox=True)
        #fig.autofmt_xdate()
        savefigs(FIG_DIR+'/delta_t_serie1_'+tagforfilename)
        #close()
        P.close()

        # 2eme couple: i=282, j=111 (2.9W/45N) --- i=303, j=111 (1.8W/45N)
        # (-1) car tableaux commencent a zero.
        i1=302
        j1=110
        i2=281
        j2=110
        serie2o = result.obs[:,i1,j1]-result.obs[:,i2,j2]
        serie2m = result.model[:,i1,j1]-result.model[:,i2,j2]

        lo = result.obs.getLongitude().getValue()
        la = result.obs.getLatitude().getValue()


        # Trace des resultats
        kwploto=dict(marker = 'o', color='r', ls=' ')
        kwplotm=dict(marker = 'o', color='b', ls=' ')
        P.figure()
        curve2(serie2o,**kwploto)
        curve2(serie2m,**kwplotm)

        #fig = figure()
        #ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        #x = datestr2num(strtime(result.obs.getTime()))
        title('SST ( %(#)3.2fE / %(##)4.2fN ) - SST ( %(###)3.2fE / %(####)4.2fN )' %{'#':lo[i1], '##':la[j1], '###':lo[i2], '####':la[j2]})
        #ylabel(result.obs.units)
        #ax.plot_date(x,serie2o.data, color='r', ls='-', marker=' ')
        #ax.plot_date(x,serie2m.data,  color='b', ls='-', marker=' ')
        #ax.legend(('Observations',  'Model'),  loc = 'upper right',  shadow = True,  fancybox=True)
        #fig.autofmt_xdate()
        savefigs(FIG_DIR+'/delta_t_serie2_'+tagforfilename)
        #close()
        P.close()

        
        
        
    if config.get('Statistics', 'temporal_mean') == 'True': 
        # -- it works !
        
        print ' - Calcul de la moyenne temporelle pour le modele et les observations - '
	print 25*'-'
	#Vieille version du code
        #result.temporal_mean()        
        #kwplot = dict(vmin=result.obs.temp_mean.min(), vmax=result.obs.temp_mean.max(), nmax = 30, colorbar_extend='both')
        #P.figure()
        #map(result.model.temp_mean, is !title='Mean modelled Sea Surface Temperature - '+tag,  show=False,  clabel_hide=True, **kwplot)
        #savefigs(FIG_DIR+'/model_temporal_mean_'+tagforfilename)
        #P.close()
        #P.figure()
        #map(result.obs.temp_mean, title='Mean observed Sea Surface Temperature - '+tag,  show=False,  clabel_hide=True, **kwplot)
        #savefigs(FIG_DIR+'/obs_temporal_mean_'+tagforfilename)
        #P.close()

        ## Nettoyage
        #result.obs.temp_mean=[]
        #result.model.temp_mean=[]
        #gc.collect()

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
        
	
	
	kwplot = dict(cmap=cmap_previmer, vmin=2, vmax=24, levels=post, norm=norm, colorbar_extend='both')
        P.figure()
        map(result.model.temp_mean, title='Mean modelled Sea Surface Temperature',  show=False,  clabel_hide=True, **kwplot)
        savefigs(FIG_DIR+'/model_temporal_mean')
        P.close()
        P.figure()
        map(result.obs.temp_mean, title='Mean observed Sea Surface Temperature',  show=False,  clabel_hide=True, **kwplot)
        savefigs(FIG_DIR+'/obs_temporal_mean')
        P.close()
        
        #kwplot = dict(vmin=result.obs.temp_mean.min(), vmax=result.obs.temp_mean.max(), nmax = 30, colorbar_extend='both')
        #P.figure()
        #map(result.model.temp_mean, title='Mean modelled Sea Surface Temperature - '+tag,  show=False,  clabel_hide=True, **kwplot)
        #savefigs(FIG_DIR+'/model_temporal_mean_'+tagforfilename)
        #P.close()
        #P.figure()
        #map(result.obs.temp_mean, title='Mean observed Sea Surface Temperature - '+tag,  show=False,  clabel_hide=True, **kwplot)
        #savefigs(FIG_DIR+'/obs_temporal_mean_'+tagforfilename)
        

    if config.get('Statistics', 'spatial_mean') == 'True': 
        # -- it works !
        print 25*'-'
        print ' - Calcul de la moyenne spatiale pour le modele et les observations - '

        result.spatial_mean()
	    #Vieille version du code
        #result.obs.spa_mean = MV2.masked_equal(result.obs.spa_mean,0)
        #result.model.spa_mean = MV2.masked_equal(result.model.spa_mean,0)


        ## Trace des resultats
        #kwploto=dict(marker = ' ', color='r', ls='-', label='Observations')
        #kwplotm=dict(marker = ' ', color='b', ls='--', label='Model')
        #P.figure()
        #curve2(result.obs.spa_mean,**kwploto)
        #curve2(result.model.spa_mean,**kwplotm)

        #title('Spatial average over the domain - '+tag)
        #savefigs(FIG_DIR+'/result_spatial_statmean_'+tagforfilename)
        #P.close()
        #close()      
        ## Nettoyage
        #result.obs.spa_mean=[]
        #result.model.spa_mean=[]
        #gc.collect()
    
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
        if ideb == [] :
	    ideb=0
	    ifin=len(index)
	elif len(ideb)!=0:
	      if index[0][0]!=ideb[0]:
		  ideb = np.append(index[0][0],ideb)
	      if index[0][-1]!=ideb[-1]:
		  ifin = np.append(ifin,index[0][-1])
	    
	    ## +1 pour les appels dans les tableaux
	ifin[ifin<len(result.obs.spa_mean.data)-1]=ifin[ifin<len(result.obs.spa_mean.data)-1]+1
	    
	index = np.where(~np.isnan(result.model.spa_mean.data))
	di = np.diff(index[0][:])
	mifin=index[0][di>1]
	mideb=index[0][di>1]+di[di>1]
	if mideb ==[]:
	    mideb=0
	    mifin=len(index)
	elif len(ideb)!=0: 
	    if index[0][0]!=mideb[0]:
		mideb = np.append(index[0][0],mideb)
	    if index[0][-1]!=mideb[-1]:
		mifin = np.append(mifin,index[0][-1])
	    
	    ## +1 pour les appels dans les tableaux
	mifin[mifin<len(result.model.spa_mean.data)-1]=mifin[mifin<len(result.model.spa_mean.data)-1]+1                        
	    
	  
	    
	    ## Boucle sur les differents blocs
	for ibloc in np.arange(len(ideb)):
		
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
	savefigs(FIG_DIR+'/result_spatial_statmean')
	close()
	      # Nettoyage
	gc.collect()
        
    if config.get('Statistics', 'temporal_std') == 'True':    
        # -- it works !
        print 25*'-'
        print ' Calcul de la std en chaque point geographique'
        result.model.temporal_std()
        print 25*'-'
        result.obs.temporal_std()
        
        #kwplot = dict(vmin=0, vmax=4, nmax = 30, colorbar_extend='max')
        #P.figure()
        #map(result.model.temp_std, title='Model - std(SST(t)) - '+tag,  show=False,  clabel_hide=True, **kwplot)
        #savefigs(FIG_DIR+'/model_temporal_std_'+tagforfilename)
        #P.close()
        #P.figure()
        #map(result.obs.temp_std, title='Observation - std(SST(t)) - '+tag, show=False,  clabel_hide=True, **kwplot)
        #savefigs(FIG_DIR+'/obs_temporal_std_'+tagforfilename)
        #P.close()

	kwplot = dict(vmin=result.obs.temp_std.min(), vmax=result.obs.temp_std.max(), nmax = 30, colorbar_extend='max')
	P.figure()
	map(result.model.temp_std, title='Model - std(SST(t))',  show=False,  clabel_hide=True, **kwplot)
	savefigs(FIG_DIR+'/model_temporal_std')
	P.close()
	P.figure()
	map(result.obs.temp_std, title='Observation - std(SST(t))', show=False,  clabel_hide=True, **kwplot)
	savefigs(FIG_DIR+'/obs_temporal_std')
	P.close()


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
        # Trace des resultats
        #fig = figure()
        #ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        #x = datestr2num(strtime(result.obs.spa_std.getTime()))
        #title('Spatial STD over the domain - '+tag)
        #ylabel(result.obs.units)
        #ax.plot_date(x,result.obs.spa_std.data, color='r', marker=' ', ls='-')
        #ax.plot_date(x,result.model.spa_std.data, color='b', marker=' ', ls='--')  
        #ax.legend(('Observations',  'Model'),  loc = 'upper right',  shadow = True,  fancybox=True)
        #fig.autofmt_xdate()
        #savefigs(FIG_DIR+'/result_spatial_statstd_'+tagforfilename)
        #close()
        
        fig = figure()
	ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
	x = datestr2num(strtime(result.obs.spa_std.getTime()))
	title('Spatial STD over the domain')
	ylabel(result.obs.units)
	ax.plot_date(x,result.obs.spa_std.data, marker='o', color='r',  markerfacecolor='r',xdate=True)
	ax.plot_date(x,result.model.spa_std.data, marker='o', color='b',  markerfacecolor='b',xdate=True)  
	ax.legend(('Observations',  'Model'),  loc = 'upper right',  shadow = True,  fancybox=True)
	fig.autofmt_xdate()
	savefigs(FIG_DIR+'/result_spatial_statstd')
	close()
        
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
	savefigs(FIG_DIR+'/result_evol_bias')
	close()
        
        ## Trace des resultats
        #fig = figure()
        #ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        #x = datestr2num(strtime(result.obs.spa_mean.getTime()))
        #title('Mean bias over the domain (<E> = <OBS> - <MODEL>) - '+tag)
        #ylabel(result.obs.units)
        #ax.plot_date(x,result.obs.spa_mean.data-result.model.spa_mean.data, marker=' ', ls='-')
        #ax.plot_date(x,np.zeros(result.obs.spa_mean.data.shape), 'k:')
        
        #fig.autofmt_xdate()
        #savefigs(FIG_DIR+'/result_evol_bias_'+tagforfilename)
        #close()

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
        
        ##map(result.model.temp_mean, res='c', show=False)        
        #P.figure()
        #kwplot = dict(vmin=-2, vmax=2, nmax=20)
        #map(interm, title='Mean bias over the domain (<E> = <OBS> - <MODEL>) - '+tag,  show=False,  clabel_hide=True,  **kwplot)
        #savefigs(FIG_DIR+'/result_mean_bias_'+tagforfilename)
        #P.close()             
        
        bias_list=np.array([-10,-5,-4,-3,-2,-1,-0.5,0.5,1,2,3,4,5,10])        
        # - normalisation
        norm = StepsNorm(bias_list)        
        
        #map(result.model.temp_mean, res='c', show=False)        
        P.figure()
        #kwplot = dict(vmin=-abs(interm).max(), vmax=abs(interm).max(), nmax=20)
        kwplot = dict(cmap=cmap_magic(bias_list,anomaly=True), vmin=bias_list.min(), vmax=bias_list.max(), levels=bias_list, norm=norm, colorbar_extend='both')
        map(interm, title='Mean bias over the domain (<E> = <OBS> - <MODEL>)',  show=False,  clabel_hide=True,  **kwplot)
        savefigs(FIG_DIR+'/result_mean_bias')
        P.close()  
        
        # Nettoyage
        result.obs.temp_mean=[]
        result.model.temp_mean=[]   
        gc.collect()
        
    if config.get('Statistics', 'extreme_bias') == 'True':        
        print 25*'-'
        print ' - Calcul du biais minimum et maximum '        
        print 25*'-'
        result.model.bias_extrema()
        result.obs.bias_extrema()
        #P.figure()
        #kwplot = dict(vmin=-abs(result.biasmin).max(), vmax=abs(result.biasmin).max())
        #map(result.biasmin, title='Minimum bias - '+tag,  show=False,  clabel_hide=True,  **kwplot)
        #savefigs(FIG_DIR+'/result_min_bias_'+tagforfilename)        
        #P.close()       
        
        #P.figure()
        #kwplot = dict(vmin=-abs(result.biasmax).max(), vmax=abs(result.biasmax).max())
        #map(result.biasmax, title='Maximum bias - '+tag,  show=False,  clabel_hide=True,  **kwplot)
        #savefigs(FIG_DIR+'/result_max_bias_'+tagforfilename)        
        #P.close()  
	
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
	
       
        # Nettoyage
        result.biasmin = []
        result.biasmax = []
        gc.collect()
 
    if config.get('Statistics', 'systematic_bias') == 'True':     
        print 25*'-'
        print ' - Calcul du biais systematique '
        print 25*'-'
        result.model.systematic_bias()
        result.obs.systematic_bias()
        #P.figure()
        #kwplot = dict(vmin=-1, vmax=1,  levels=[-1, -0.1,0.1,  1], colorbar_ticks=[1,0,-1])
        #map(result.biassyst, title='Systematic bias - '+tag+'\n (+1: model underestimation, -1: model overestimation)',  show=False,  clabel_hide=True,  **kwplot)
        #savefigs(FIG_DIR+'/result_syst_bias_'+tagforfilename)        
        #P.close()
        
        P.figure()
	kwplot = dict(vmin=-1, vmax=1,  levels=[-1, -0.1,0.1,  1], colorbar_ticks=[1,0,-1])
	map(result.biassyst, title='Systematic bias\n (+1: model underestimation, -1: model overestimation)',  show=False,  clabel_hide=True,  **kwplot)
	savefigs(FIG_DIR+'/result_syst_bias')        
	P.close()
        
        
        # Nettoyage
        result.biassyst = []
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
        #P.figure()
        #map(result.obs_cov,  title='Observation coverage during the period - '+tag, show=False,  clabel_hide=True, nmax=10)
        #savefigs(FIG_DIR+'/result_obs_coverage_'+tagforfilename)
        #P.close()

	P.figure()
	map(result.obs_cov,  title='Observation coverage during the period', show=False,  clabel_hide=True)
	savefigs(FIG_DIR+'/result_obs_coverage')
	P.close()
    

        # Nettoyage
        result.obs_cov = []
        gc.collect()

    
    if config.get('Statistics', 'obs_spatialcoverage') == 'True':    
        # -- it works !
        print 25*'-'
        print ' - Calcul du taux de couverture des observations en chaque pas de temps - '
        print 25*'-'
        result.spatial_obs_coverage()        
        
        # Ne fonctionne pas !!!!
        #select = dict(color='g', shadow=True, ylabel=result.obs_spacov.units , show=False,  title='Observation coverage following the time period period')
        #curve2(result.obs_spacov, **select)
        
        # Trace des resultats
        #fig = figure()
        #ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        #x = datestr2num(strtime(result.obs_spacov.getTime()))
        #title('Observation coverage following the time period - '+tag)
        #ylabel('Number of observations')
        #ax.plot_date(x,result.obs_spacov.data,linestyle='-', color='r', marker=' ')
        ##ylim(0, 100)
        #fig.autofmt_xdate()
        #savefigs(FIG_DIR+'/result_obs_spatial_coverage_'+tagforfilename)
        #close()
        
        fig = figure()
	ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
	x = datestr2num(strtime(result.obs_spacov.getTime()))
	title('Observation coverage following the time period')
	ylabel('Coverage (nb)')
	ax.plot_date(x,result.obs_spacov.data,linestyle='-', marker='o', color='r',  markerfacecolor='b',xdate=True)
	ylim()
	fig.autofmt_xdate()
	savefigs(FIG_DIR+'/result_obs_spatial_coverage')
	close()
        
        # Nettoyage
        result.obs_spacov = []
        gc.collect()
     
    
    
    if config.get('Statistics', 'temporal_rms') == 'True':    
        # -- it works !
        print 25*'-'
        print ' - Calcul de la RMS temporelle - '
        print 25*'-'
        result.temporal_rms()
        #P.figure()
        #map(result.temp_rms, title='RMS - '+tag+'\n (uncentered and biased)', show=False,  clabel_hide=True, vmin=0., vmax=2.)
        #savefigs(FIG_DIR+'/result_temporal_rms_'+tagforfilename)
        #P.close()
         
        P.figure()
	map(result.temp_rms, title='RMS\n (uncentered and biased)', show=False,  clabel_hide=True)
	savefigs(FIG_DIR+'/result_temporal_rms')
	P.close()
         
        # Nettoyage
        result.temp_rms = []
        gc.collect()
    
        
    if config.get('Statistics', 'temporal_rmsc') == 'True':    
        # -- it works !
        print 25*'-'
        print ' - Calcul de la RMS temporelle  centree- '
        print 25*'-'
        result.temporal_rmsc()
        #P.figure()
        #map(result.temp_rmsc, title='RMS - '+tag+'\n (centered and biased)', show=False,  clabel_hide=True, vmin=0., vmax=2.)
        #savefigs(FIG_DIR+'/result_temporal_rmsc_'+tagforfilename)
        #P.close()
       
        P.figure()
	map(result.temp_rmsc, title='RMS\n (centered and biased)', show=False,  clabel_hide=True)
	savefigs(FIG_DIR+'/result_temporal_rmsc')
	P.close()
       
        # Nettoyage
        result.temp_rmsc = []
        gc.collect()
    
 
    if config.get('Statistics', 'spatial_rms') == 'True': 
        # -- it works !
        print 25*'-'
        print ' - Calcul de la RMS entre le modele et les observations pour chaque pas de temps - '
        print 25*'-'
        result.spatial_rms()       
        
        ## Trace des resultats
        #fig = figure()
        #ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        #x = datestr2num(strtime(result.spa_rms.getTime()))
        #title('Spatial RMS over the domain - '+tag)
        #ylabel(result.obs.units)
        #ax.plot_date(x,result.spa_rms.data, marker='o', color='r',  markerfacecolor='r', ls='-')
        #fig.autofmt_xdate()
        #savefigs(FIG_DIR+'/result_spatial_statrms_'+tagforfilename)
        #close()
      
	fig = figure()
	ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
	x = datestr2num(strtime(result.spa_rms.getTime()))
	title('Spatial RMS over the domain')
	ylabel(result.obs.units)
	ax.plot_date(x,result.spa_rms.data, marker='o', color='r',  markerfacecolor='r',xdate=True)
	fig.autofmt_xdate()
	savefigs(FIG_DIR+'/result_spatial_statrms')
	close()
      
        # Nettoyage
        result.spa_rms = []
        gc.collect()
    
        
    if config.get('Statistics', 'spatial_rmsc') == 'True': 
        # -- it works !
        print 25*'-'
        print ' - Calcul de la RMS centree entre le modele et les observations pour chaque pas de temps - '
	print 25*'-' 
        result.spatial_rmsc()       
        
        ## Trace des resultats
        #fig = figure()
        #ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
        #x = datestr2num(strtime(result.spa_rmsc.getTime()))
        #title('Spatial centered RMS over the domain - '+tag)
        #ylabel(result.obs.units)
        #ax.plot_date(x,result.spa_rmsc.data, marker='o', color='r',  markerfacecolor='r', ls='-')
        #fig.autofmt_xdate()
        #savefigs(FIG_DIR+'/result_spatial_statrmsc_'+tagforfilename)
        #close()
          
        fig = figure()
	ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
	x = datestr2num(strtime(result.spa_rmsc.getTime()))
	title('Spatial centered RMS over the domain')
	ylabel(result.obs.units)
	ax.plot_date(x,result.spa_rmsc.data, marker='o', color='r',  markerfacecolor='r',xdate=True)
	fig.autofmt_xdate()
	savefigs(FIG_DIR+'/result_spatial_statrmsc')
	close()
        
        # Nettoyage
        result.spa_rmsc = []
        gc.collect()
    
    
    if config.get('Statistics', 'temporal_corr') == 'True':    
        # -- it works !
        print 10*'-'
        print ' - Calcul de la Correlation temporelle - '
        result.temporal_corr()
        #P.figure()
        #map(result.temp_corr, title='Time series correlation - '+tag, show=False,  clabel_hide=True)
        #savefigs(FIG_DIR+'/result_temporal_corr_'+tagforfilename)
        #P.close()

	P.figure()
	map(result.temp_corr, title='Time series correlation', show=False,  clabel_hide=True)
	savefigs(FIG_DIR+'/result_temporal_corr')
	P.close()
	

        # Nettoyage
        result.temp_corr = []
        gc.collect()
    

