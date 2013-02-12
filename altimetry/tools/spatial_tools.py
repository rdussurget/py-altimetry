# -*- coding: utf-8 -*-
import numpy as np
from altimetry.tools.others import get_zero_element
if __debug__ : import matplotlib.pyplot as plt

def in_limits(lon, lat, limit):
    
    limit = recale_limits(limit)
    
    lats = limit[0]
    lons = limit[1]
    late = limit[2]
    lone = limit[3]

    ind=[]
    
    lon = recale(lon,degrees=True)

    # find out lon between S and E (start and end)
    if lons < lone:
        lon_flag = (lon >= lons) & (lon <= lone)
        
    # find out lon outside S and E (between E and S)
    else:
        lon_flag = (lon >= lons) | (lon <= lone)

    # nothing remains
    if not np.sometrue(lon_flag):
        print "* WARNING : No lon [", lons, ",", lone, "]"
        return ind, lon_flag

    # find lat constraints
    lat_flag = (lat >= lats) & (lat <= late)
 
    # nothing remains
    if not np.sometrue(lat_flag):
        print "* WARNING : No Lat [", lats, ",", lats, "]  with lon  [", lons, ",", lone, "]"
        return ind, lat_flag
    
    #Construct flag and index arrays
    if (len(lon_flag) == len(lat_flag)) :
        flag = lon_flag & lat_flag
        ind = np.arange(len(lon)).compress(flag)
    else :
        flag = (lon_flag,lat_flag)
        ind = (np.arange(len(lon)).compress(flag[0]),np.arange(len(lat)).compress(flag[1]))
    

    # compress indexes
#    ind = ind.compress(flag)
    
    return ind, flag

def recale_limits(limitin, zero_2pi=True, minpi_pi=None): 
    
    #defaults
    if (minpi_pi is not None) : zero_2pi != minpi_pi
    
    limitout=limitin
    if isinstance(limitin,np.ndarray) : limitout[[1,3]]=recale(limitin[[1,3]], zero_2pi=zero_2pi, minpi_pi=minpi_pi,degrees=True) #Always in degrees!
    else :
        limitout[1]=recale(limitin[1], zero_2pi=zero_2pi, minpi_pi=minpi_pi,degrees=True)
        limitout[3]=recale(limitin[3]-1e2*np.MachAr().resolution, zero_2pi=zero_2pi, minpi_pi=minpi_pi,degrees=True)

    return limitout

def recale(angle, degrees=False, zero_2pi=True, minpi_pi=None) :

    #defaults
    if (minpi_pi is not None) : zero_2pi = not minpi_pi
    radians = not degrees
  
    modulus= 360.0 if degrees else 2.0 * np.pi
  
    limits=np.array([0.0, 2.0*np.pi])
    if minpi_pi : limits-=np.pi
    if degrees : limits = np.rad2deg(limits)
    
    angle = np.mod(angle + limits[0],limits[0]+limits[1]) - limits[0]
    
#    over=angle >= limits[1]
#    under=angle < limits[0]
#    angle+=under*modulus
#    angle-=over*modulus
  
    return angle


def rgb_int2tuple(rgbint):
    return (rgbint // 256 // 256 % 256, rgbint // 256 % 256, rgbint % 256)


def calcul_distance(*args):
    """
    ;+
    ; CALCUL_DISTANCE : fonction permettant de calculer la distance entre deux points sur la sphere
    ;
    ; @Author : Mathilde FAILLOT (LEGOS/CTOH), 2005
    ; @History :
    ;    - 2005 : First release, source:<br />
    ;             Weisstein, Eric W. "Great Circle." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/GreatCircle.html<br /> 
    ;    - 2008 : Renaud DUSSURGET, optimization + vectors + Mercurial<br />
    ;    - Feb. 2009 : RD, debug on repeated positions removal<br />
    ;    - March 2010 : RD, force coordinates to FLOATS if Floating-Point operand errors<br />
    ;    - Nov. 2010 : 2 params option for calculating distance vector
    ;    - Apr. 2012 : Converted to Python
    ;-
    """
#def calcul_distance(norepeat=False, *args):
    
    lon_a_in = args[1]
    lat_a_in = args[0]
    
    if np.size(lon_a_in) != np.size(lat_a_in) : raise 'Error : arguments are not the same size'
    
    if len(args) == 2 :
        lon_b_in = lon_a_in
        lat_b_in = lat_a_in
        lon_a_in = get_zero_element(lon_a_in)
        lat_a_in = get_zero_element(lat_a_in)
    elif len(args) == 4 :
        lat_b_in = args[2]
        lon_b_in = args[3]
        if np.size(lon_b_in) != np.size(lat_b_in) : raise 'Error : arguments are not the same size'
    else :
        print "ERROR"
        return -1

    #   Define constants
    
    #Degree to radians conversion
    #deg_rad = np.pi/360.0
    
    #Earth radius (km)
    rt = 6378.137
    
    #;Remove repeated positions
    #IF (KEYWORD_SET(NOREPEAT)) THEN id = WHERE(~((lon_b_in EQ lon_a_in[0]) * (lat_b_in EQ lat_a_in[0])),okCnt) ELSE okcnt = 1 
    #IF (okCnt EQ 0) THEN RETURN lon_a_in, DOUBLE(0.)

    #Degree to radians conversion
    lon_a = np.deg2rad(lon_a_in)
    lat_a = np.deg2rad(lat_a_in)
    lon_b = np.deg2rad(lon_b_in)
    lat_b = np.deg2rad(lat_b_in)


    #;   calcul de la distance entre les deux points consideres
#    interm = (np.cos(lon_b) * np.cos(lat_b)) * (np.cos(lon_a) * np.cos(lat_a)) + \
#      (np.sin(lon_a) * np.cos(lat_a)) * (np.sin(lon_b) * np.cos(lat_b)) + \
#      np.sin(lat_a) * np.sin(lat_b)
    
#    interm = (np.cos(lon_b) * np.cos(lat_b)) * (np.cos(lon_a) * np.cos(lat_a)) + \
#      (np.sin(lon_a) * np.cos(lat_a)) * (np.sin(lon_b) * np.cos(lat_b)) + \
#      np.sin(lat_a) * np.sin(lat_b)

    interm = np.cos(lat_a) * np.cos(lat_b) * np.cos(lon_a - lon_b) + np.sin(lat_a) * np.sin(lat_b) #Wolfram Math World definition
    
    dist = rt * np.arccos(interm)
    
    #Remove computation errors
    fgerr=(lon_b == lon_a) & (lat_b == lat_a)
    if fgerr.sum() > 0 : dist[fgerr] = 0.
    
    return dist

def cumulative_distance(lat, lon,dist_int=None):
    '''
    ;+
    ; CUMULATIVE_DISTANCE : permet de calculer la distance le long d'une ligne.
    ;
    ; @Author : Renaud DUSSURGET, LEGOS/CTOH
    ; @History :
    ;    - Feb. 2009 : First release (adapted from calcul_distance)
    ;-
    '''

    #   Define constants
    
    #Degree to radians conversion
    #deg_rad = np.pi/360.0
    
    #Earth radius (km)
    rt = 6378.137
    
    #;Remove repeated positions
    #IF (KEYWORD_SET(NOREPEAT)) THEN id = WHERE(~((lon_b_in EQ lon_a_in[0]) * (lat_b_in EQ lat_a_in[0])),okCnt) ELSE okcnt = 1 
    #IF (okCnt EQ 0) THEN RETURN lon_a_in, DOUBLE(0.)

    nelts=lon.size

    #Degree to radians conversion
    lon_a = np.deg2rad(lon[0:nelts - 1])
    lon_b = np.deg2rad(lon[1:nelts])
    lat_a = np.deg2rad(lat[0:nelts - 1])
    lat_b = np.deg2rad(lat[1:nelts])
    
    interm = np.cos(lat_a) * np.cos(lat_b) * np.cos(lon_a - lon_b) + np.sin(lat_a) * np.sin(lat_b) #Wolfram Math World definition
    
    dist_int=np.append(0,rt*np.arccos(interm))
    
    return dist_int.cumsum()
