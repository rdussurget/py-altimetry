# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 17:46:24 2012

@author: sskrypni

Last update: 03/2012
"""
import os

# --------------------------------------------------------------------
# Définition des classes utiles pour l'analyse
class Profile() :
    
    glider_begin = "EXGL"
    
    glider_type = "GLIDER"
    recopeca_type = "RECOPESCA"
    
    #cle des differentes tables
    time_key = "time"
    lon_key = "lon"
    lat_key = "lat"
    depth_key = "depth"
    temp_key = "temp"
    sal_key = "sal"
    std_temp_key="std_temp"
    std_sal_key="std_sal"
    """ Profile vertical """
    def __init__(self, time, platform_name, lon, lat, cfg):   
        
        self.name = "Vertical Profile"
        self.shortname = "VertProf"   
        
        self.time = time
        self.lon = lon
        self.lat = lat
        
        #tableau numpy
        self.depth = []
        self.temp = []
        self.sal = []
        self.std_temp = []
        self.std_sal = []
        #...
        # AJOUTER NOUVELLES VARIABLES ICI
        #...
                
        self.platform_name = platform_name
        if self.platform_name.startswith(self.glider_begin):
            self.platform_type = self.glider_type
        else:
            self.platform_type = self.recopeca_type


    # initialisation : creation des list du profil
    def add_time_step(self, depth=None, temp=None, sal=None,std_temp=None, std_sal=None):
        self.depth.append(depth)
        self.temp.append(temp)
        self.sal.append(sal)
        self.std_temp.append(std_temp)
	self.std_sal.append(std_sal)

        #...
        # AJOUTER NOUVELLES VARIABLES ICI
        #...

    
    #fin initialistion :creation table numpy
    def convert_tables(self):
        from numpy import array
        self.depth = array(self.depth, dtype='float')
        self.temp = array(self.temp, dtype='float')
        self.sal = array(self.sal, dtype='float')        
        self.std_temp =array(self.std_temp, dtype='float')
        self.std_sal = array(self.std_sal, dtype='float')
        #...
        # AJOUTER NOUVELLES VARIABLES ICI
        #...

    
    # creation de tables pour tracer les profils
    # retourne une map de ndarray        
    def create_tabs(self):
        from numpy import ones_like
        from matplotlib.dates import date2num
        
        #duplication sur le nombre de profondeur (utile pour tracer le profil) 
        time_array = ones_like(self.depth) * date2num(self.time)
        lon_array = ones_like(self.depth) * float(self.lon)
        lat_array = ones_like(self.depth) * float(self.lat)        
        
        map_return = {self.time_key : time_array,
                      self.lon_key :lon_array,
                      self.lat_key :lat_array
                      }
        
        return map_return

    
    # retourne les infos date,lon et lat du profil
    def get_position(self):
        return [self.time, self.lon, self.lat]




class Prof_insitu():
    qc_ok = '0111111'
    ''' Donnees RECOPESCA '''
    def __init__(self, cfg):
        import cdtime
        SCRIPT_DIR = os.getcwd()
        self.SCRIPT_DIR = SCRIPT_DIR
        andeb = cfg['Time Period']['andeb']        
        anfin = cfg['Time Period']['anfin']		
        mdeb = cfg['Time Period']['mdeb']
        mfin = cfg['Time Period']['mfin']
        jdeb = cfg['Time Period']['jdeb']
        jfin = cfg['Time Period']['jfin']
        hdeb = cfg['Time Period']['hdeb']
        hfin = cfg['Time Period']['hfin']

        self.WORKDIR = cfg['Env']['workdir']
        self.ctdeb = cdtime.comptime(andeb, mdeb, jdeb, hdeb, 0, 0)
        self.ctfin = cdtime.comptime(anfin, mfin, jfin, hfin, 0, 0)
       
        #print self.ctfin
 
        DIR_REC = os.path.join(self.WORKDIR, 'RECOPESCA')    # repertoire de stockage des donnees     
        if os.path.isdir(DIR_REC) == False:
            os.chdir(self.WORKDIR)
            os.mkdir('RECOPESCA')
        self.WORKDIR = os.path.join(self.WORKDIR, 'RECOPESCA')

        self.name = "RECOPESCA data"
        self.shortname = "Recopesca"

        #dictionnaire contient les profils (cle : date)
        self.map_profiles = {} 
        
    def rappatrie(self, cfg):
        """ Rappatrie par ftp les donnees RECOPESCA """
        import subprocess, cdtime
        
        os.chdir(self.WORKDIR)  # on se place dans le repertoire de travail
        #----------------------------------------------------
        print '---------- RECUPERATION FICHIERS RECOPESCA ----------'
        print 65 * '-'
        #----------------------------------------------------
    
    
        #-------------------------------------------------------------
        #------- recuperation des donnees : generalites (site ftp, etc)
        URL_CDOCO = cfg['RECOPESCA']['url_cdoco']
        DATA_DIR = cfg['RECOPESCA']['data_dir']
        URL_REC_DATA = os.path.join(URL_CDOCO, DATA_DIR)
        EXT = ".csv"         
        usr = cfg['RECOPESCA']['user']
        pwd = cfg['RECOPESCA']['pwd']
    
        os.chdir(self.WORKDIR)

        #-- recuperation des donnees par FTP anonymous
        #----------------------------------------------------
        ctest = self.ctdeb

        # prevoir un test pour les cas ou le fichier de donnees n'existe pas !!! 
        while ctest.month <= self.ctfin.month and ctest.year <= self.ctfin.year:  
            Prof_insitu.recup_ftp_recopesca(self, URL_CDOCO, DATA_DIR, ctest.year, ctest.month, EXT, usr, pwd, '9')                  
            # On incremente le temps "test" de 1 mois (va au premier du mois suivant)
            print ctest
            ctest = ctest.add(1, cdtime.Months) 
       
        os.chdir(self.WORKDIR)
        #-- Fin de recuperation des donnees
        #----------------------------------------------------
    def recup_ftp_recopesca(self, URL_DATA, DATA_DIR, YYYY, MM, ext, usr, pwd, ntry):
        """ Recuperation par FTP d un fichier RECOPESCA mensuel""" 
        from ftplib import FTP
        
        if MM < 10:
          MM = "0%(MM)s" % vars()
        
        # connection au CDOCO
        ftp = FTP(host=URL_DATA, user=usr, passwd=pwd)
        
        # positionnement sur le repertoire "best_estimate"
        ftp.cwd(DATA_DIR)
        
        # utile pour ftp.retrbinary
        def handleDownload(block):
            file.write(block)
        
        filename = '*%(#)s%(##)s*%(###)s' % {'#':YYYY, '##':MM, '###':ext}    
        #print filename
        
        list_file = ftp.nlst(filename)
        #print list_file
        
        for file_to_read in list_file:
            # recuperation fichier si non present dans work/MODEL/MARS
            if os.path.isfile(file_to_read) == False:
                # rajouter test sur date du fichier ... telecharge que si plus recent ...
                
                file = open(file_to_read, 'wb')
                
                ftp.retrbinary('RETR ' + file_to_read , handleDownload)
            # executer un test pour voir si le fichier existe dans la base de donnee
            # ....
            # ....
        
        ftp.quit()
        
        
        #-- Fin de ftp
        #----------------------------------------------------
    
    def read(self,cfg):
        """ Lecture des fichiers csv de RECOPESCA """
        import glob
        #from vacumm.misc.atime import strptime
        from vacumm.misc.atime import numtime, strtime
        from matplotlib.dates import num2date
        import cdtime
        
        # - convertisseur pour le temps
        def convtime(s):
            return numtime(s)
        
        #----------------------------------------------------
        ctest = self.ctdeb
        ext = ".csv"
        platform_name = list()

        #print self.ctfin.month

        while ctest.month <= self.ctfin.month and ctest.year <= self.ctfin.year:  

            #print ctest.month 
            
            if ctest.month < 10:
                MM = "0%(#)s" % {'#':ctest.month}
            else:
                MM = "%(#)s" % {'#':ctest.month}
                
            # On incremente le temps "test" de 1 mois
            filename = '*%(#)s%(##)s*%(###)s' % {'#':ctest.year, '##':MM, '###':ext}  
            #print glob.glob(filename)
            
            #creaton map temporaire
            map_profiles_tmp = {}
            
            #parcours des fichiers
            for file in glob.glob(filename):
                
                #lecture fichier
                file_data = open(file).read()
                #suppression ligne 1
                first_line, file_data = file_data.split('\n', 1)
                
                #extraction des profil du fichier (ligne a ligne)
                for line in file_data.splitlines():
                    
                    line = line.strip().split(',')
                    
                    platform_name = line[0]
                    time = line[1]
                    lat = line[2]
                    lon = line[3]
                    depth = line[4]
                    temp = line[5]
                    sal = line[6]
                    qc = line[7]
                    
                    #test qc (si != 0111111 non traite)
                    if qc != self.qc_ok:
                        continue
                    
                    time = num2date(convtime(time))
                    key = strtime(time)
                    
                    #si profil existe
                    if key in map_profiles_tmp:
                        profile = map_profiles_tmp.get(key)                
                        profile.add_time_step(depth, temp, sal)
                    
                    #sinon creation profil    
                    else:
                        profile = Profile(time, platform_name, lon, lat,cfg)     
                        profile.add_time_step(depth, temp, sal)  
                        map_profiles_tmp[key] = profile
                
                #creation tables numpy
                for key in map_profiles_tmp:
                    profile=map_profiles_tmp[key]
                    profile.convert_tables()
                
                #copie des profil du fichier
                self.map_profiles.update(map_profiles_tmp)
                
                #reinitialisation
                map_profiles_tmp.clear()
            
            ctest = ctest.add(1, cdtime.Months) # Passe au premier du mois suivant
    
    def get_profiles_position(self):
        from matplotlib.dates import date2num
        from numpy import array
        time_list = []
        lon_list = []
        lat_list = []
        
        #creation list
        for key in self.map_profiles:
            profile=self.map_profiles[key]
            time, lon, lat = profile.get_position()
            time_list.append(date2num(time))
            lat_list.append(lat)
            lon_list.append(lon)
        
        #conversion list table numpy
        time_list = array(time_list, dtype='float')    
        lon_list = array(lon_list, dtype='float')    
        lat_list = array(lat_list, dtype='float')
        
        return [time_list, lon_list, lat_list]
