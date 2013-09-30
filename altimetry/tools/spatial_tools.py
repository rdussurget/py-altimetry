# -*- coding: utf-8 -*-
import numpy as np
from altimetry.tools.others import get_zero_element
if __debug__ : import matplotlib.pyplot as plt

def in_limits(lon, lat, limit):
    
    limitin = recale_limits(limit)
    
    lats = limitin[0]
    lons = limitin[1]
    late = limitin[2]
    lone = limitin[3]

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

def recale_limits(limitin, zero_2pi=True, minpi_pi=None,type=None): 
    
    if type is not None : limitin=[type(l) for l in limitin]
    
    #defaults
    if (minpi_pi is not None) : zero_2pi != minpi_pi
    if (zero_2pi is not None) & (zero_2pi is False) : minpi_pi = True
    
    limitout=limitin[:]
    if isinstance(limitin,np.ndarray) : limitout[[1,3]]=recale(limitin[[1,3]], zero_2pi=zero_2pi, minpi_pi=minpi_pi,degrees=True) #Always in degrees!
    else :
        limitout[1]=recale(limitin[1], zero_2pi=zero_2pi, minpi_pi=minpi_pi,degrees=True)
        limitout[3]=recale(limitin[3], zero_2pi=zero_2pi, minpi_pi=minpi_pi,degrees=True)

    if type is not None : limitin=[type(l) for l in limitout]
    return limitout

def recale(angle, degrees=False, zero_2pi=True, minpi_pi=None) :

    #defaults
    if (minpi_pi is not None) : zero_2pi = not minpi_pi
    
    modulus = 360.0 if degrees else 2.0*np.pi
    
    limits=np.array([0.0, 2.0*np.pi])
    if minpi_pi : limits-=np.pi
    if degrees : limits = np.rad2deg(limits)
#    if minpi_pi : angle+=limits[0]
    
#    angle = np.mod(angle + limits[0],limits[0]+limits[1]) - limits[0]
    
#     res = angle
#     if angle > limits[1] or angle == limits[0]:
#         angle = limits[0] + abs(angle + limits[1]) % (abs(limits[0]) + abs(limits[1]))
#     if angle < limits[0] or angle == limits[1]:
#         angle = limits[1] - abs(limits[1] - lower) % (abs(limits[0]) + abs(limits[1]))
    
    over=angle > limits[1]
    under=angle < limits[0]
#    over=np.mod(angle,limits[1])
#    under=angle < limits[0]
    angle+=under*modulus#*np.floor(np.abs(angle)/modulus)
    angle-=over*modulus#*np.floor(np.abs(angle)/modulus)
  
    
#    if minpi_pi : angle+=limits[0]
    
#    over=angle >= limits[1]
#    under=angle < limits[0]
#    angle = over *(limits[0] + abs(angle[over] + limits[1]) % (abs(limits[0]) + abs(limits[1])))
#    angle= under * (limits[1] - abs(angle[under] - limits[0]) % (abs(limits[0]) + abs(limits[1])))

  
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
    if fgerr.sum() > 0 :
        if np.size(dist) == 1 : dist = 0. 
        else : dist[fgerr] = 0.
    
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



# From http://www.movable-type.co.uk/scripts/LatLongVincenty.html:
#   The most accurate and widely used globally-applicable model for the earth
#   ellipsoid is WGS-84, used in this script. Other ellipsoids offering a
#   better fit to the local geoid include Airy (1830) in the UK, International
#   1924 in much of Europe, Clarke (1880) in Africa, and GRS-67 in South
#   America. America (NAD83) and Australia (GDA) use GRS-80, functionally
#   equivalent to the WGS-84 ellipsoid.
ELLIPSOIDS = {
    # model           major (km)   minor (km)     flattening
    'WGS-84':        (6378.137,    6356.7523142,  1 / 298.257223563),
    'GRS-80':        (6378.137,    6356.7523141,  1 / 298.257222101),
    'Airy (1830)':   (6377.563396, 6356.256909,   1 / 299.3249646),
    'Intl 1924':     (6378.388,    6356.911946,   1 / 297.0),
    'Clarke (1880)': (6378.249145, 6356.51486955, 1 / 293.465),
    'GRS-67':        (6378.1600,   6356.774719,   1 / 298.25),}


def VincentyDistance(lon1_in,lat1_in,lon2_in,lat2_in,ELLIPSOID='GRS-80'):
    """
    Calculate the geodesic distance between two points using the formula
    devised by Thaddeus Vincenty, with an accurate ellipsoidal model of the
    earth.
   
    The class attribute `ELLIPSOID` indicates which ellipsoidal model of the
    earth to use. If it is a string, it is looked up in the `ELLIPSOIDS`
    dictionary to obtain the major and minor semiaxes and the flattening.
    Otherwise, it should be a tuple with those values. The most globally
    accurate model is WGS-84. See the comments above the `ELLIPSOIDS`
    dictionary for more information.
   
    """
    lng1 = np.deg2rad(lon1_in)
    lat1 = np.deg2rad(lat1_in)
    lng2 = np.deg2rad(lon2_in)
    lat2 = np.deg2rad(lat2_in)

    major, minor, f = ELLIPSOIDS[ELLIPSOID]

    delta_lng = lng2 - lng1

    reduced_lat1 = np.arctan((1 - f) * np.tan(lat1))
    reduced_lat2 = np.arctan((1 - f) * np.tan(lat2))

    sin_reduced1, cos_reduced1 = np.sin(reduced_lat1), np.cos(reduced_lat1)
    sin_reduced2, cos_reduced2 = np.sin(reduced_lat2), np.cos(reduced_lat2)

    lambda_lng = delta_lng
    lambda_prime = 2 * np.pi

    iter_limit = 20

    while abs(lambda_lng - lambda_prime) > 10e-12 and iter_limit > 0:
        sin_lambda_lng, cos_lambda_lng = np.sin(lambda_lng), np.cos(lambda_lng)

        sin_sigma = np.sqrt( (cos_reduced2 * sin_lambda_lng) ** 2 +
                             (cos_reduced1 * sin_reduced2 -
                              sin_reduced1 * cos_reduced2 * cos_lambda_lng) ** 2
                            )

        if sin_sigma == 0:
            return 0 # Coincident points

        cos_sigma = (
            sin_reduced1 * sin_reduced2 +
            cos_reduced1 * cos_reduced2 * cos_lambda_lng
        )

        sigma = np.arctan2(sin_sigma, cos_sigma)

        sin_alpha = (
            cos_reduced1 * cos_reduced2 * sin_lambda_lng / sin_sigma
        )
        cos_sq_alpha = 1 - sin_alpha ** 2

        if cos_sq_alpha != 0:
            cos2_sigma_m = cos_sigma - 2 * (
                sin_reduced1 * sin_reduced2 / cos_sq_alpha
            )
        else:
            cos2_sigma_m = 0.0 # Equatorial line

        C = f / 16. * cos_sq_alpha * (4 + f * (4 - 3 * cos_sq_alpha))

        lambda_prime = lambda_lng
        lambda_lng = (
            delta_lng + (1 - C) * f * sin_alpha * (
                sigma + C * sin_sigma * (
                    cos2_sigma_m + C * cos_sigma * (
                        -1 + 2 * cos2_sigma_m ** 2
                    )
                )
            )
        )
        iter_limit -= 1

    if iter_limit == 0:
        raise ValueError("Vincenty formula failed to converge!")

    u_sq = cos_sq_alpha * (major ** 2 - minor ** 2) / minor ** 2

    A = 1 + u_sq / 16384. * (
        4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq))
    )

    B = u_sq / 1024. * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))

    delta_sigma = (
        B * sin_sigma * (
            cos2_sigma_m + B / 4. * (
                cos_sigma * (
                    -1 + 2 * cos2_sigma_m ** 2
                ) - B / 6. * cos2_sigma_m * (
                    -3 + 4 * sin_sigma ** 2
                ) * (
                    -3 + 4 * cos2_sigma_m ** 2
                )
            )
        )
    )

    s = minor * A * (sigma - delta_sigma)
    return s
