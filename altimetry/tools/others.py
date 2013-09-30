# -*- coding: utf-8 -*-
import numpy as np
from Tkinter import Tk
import inspect
import getpass, socket #User and Host names
from warnings import warn

def get_zero_element(array):
    try : array = array.flat.next()
    except StopIteration : pass
    return array

def bytscl(array, maximum=None , minimum=None , nan=0, top=255 ):
    """
    see http://star.pst.qub.ac.uk/idl/BYTSCL.html
    note that IDL uses slightly different formulae for bytscaling floats and ints. 
    here we apply only the FLOAT formula...
    """
    if maximum is None: maximum = np.nanmax(array)
    if minimum is None: minimum = np.nanmin(array)
    return np.maximum(np.minimum(((top+1.0)*(array-minimum)/(maximum-minimum)).astype(np.int16), top),0)

def scale_vector(x,bottom,top):
    return (x - bottom) / (top-bottom)

#Screen size widgets
def get_screen_size(units='p'):
    units='1'+units
    width = Tk().winfo_fpixels(str(Tk().winfo_screenwidth())+'p')/Tk().winfo_fpixels(units)
    height = Tk().winfo_fpixels(str(Tk().winfo_screenheight())+'p')/Tk().winfo_fpixels(units)
    return np.array([width,height])

def get_screen_dpi(units='1i'):
    return Tk().winfo_fpixels(units)

def where_list(list1,list2):
    index=[]
    for i in list1 :
        try:
            index.append(list2.index(i))
        except ValueError:
            index.append(-1)

    return index

##from scipy.interpolate import lagrange, BarycentricInterpolator, barycentric_interpolate
#def deprecated_deriv(*args):
#    if len(args) == 1 :
#        y = args[0]
#        x = np.arange(len(y))
#    if len(args) == 2 :
#        x = args[0]
#        y = args[1]
#    
#    dx = x[1:] - x[:-1]
#    dy = y[1:] - y[:-1]
#    
##    d2 = lagrange(np.arange(len(x)),np.arange(len(dx))+0.5,dy/dx)
#
#    dx = np.append(np.array(2*dx[1] - dx[0]),dx)
#    dx = np.append(dx,np.array(2*dx[-1] - dx[-2]))
#    
#    dy = np.append(np.array(2*dy[1] - dy[0]),dy)
#    dy = np.append(dy,np.array(2*dy[-1] - dy[-2]))
#    
##    dxint=lagrange(np.arange(nx),np.arange(nx-1)+0.5,dx)
##    dyint = lagrange(np.arange(nx),np.arange(nx-1)+0.5,dy)
#    
##    d = BarycentricInterpolator(dx)
##    dxint = lagrange(np.arange(len(dx)-1) + 0.5,np.arange(len(dx)),dx)
##    dyint = lagrange(np.arange(len(dx)-1) + 0.5,np.arange(len(dx)),dy)
#    d = lagrange(np.arange(len(x)) + 0.5,np.arange(len(dx)),dy/dx)
#    
##    p = barycentric_interpolate(np.arange(len(dx))+0.5,dy/dx,np.arange(len(dx)))
##    p(np.arange(len(dx)))
##    p(np.arange(len(dx)))
#    
#    return d


def deriv(*args):
    """
    ; Copyright (c) 1984-2009, ITT Visual Information Solutions. All
    ;       rights reserved. Unauthorized reproduction is prohibited.
    ;
    
    ;+
    ; NAME:
    ;    DERIV
    ;
    ; PURPOSE:
    ;    Perform numerical differentiation using 3-point, Lagrangian 
    ;    interpolation.
    ;
    ; CATEGORY:
    ;    Numerical analysis.
    ;
    ; CALLING SEQUENCE:
    ;    Dy = Deriv(Y)         ;Dy(i)/di, point spacing = 1.
    ;    Dy = Deriv(X, Y)    ;Dy/Dx, unequal point spacing.
    ;
    ; INPUTS:
    ;    Y:  Variable to be differentiated.
    ;    X:  Variable to differentiate with respect to.  If omitted, unit 
    ;        spacing for Y (i.e., X(i) = i) is assumed.
    ;
    ; OPTIONAL INPUT PARAMETERS:
    ;    As above.
    ;
    ; OUTPUTS:
    ;    Returns the derivative.
    ;
    ; COMMON BLOCKS:
    ;    None.
    ;
    ; SIDE EFFECTS:
    ;    None.
    ;
    ; RESTRICTIONS:
    ;    None.
    ;
    ; PROCEDURE:
    ;    See Hildebrand, Introduction to Numerical Analysis, Mc Graw
    ;    Hill, 1956.  Page 82.
    ;
    ; MODIFICATION HISTORY:
    ;    Written, DMS, Aug, 1984.
    ;    Corrected formula for points with unequal spacing.  DMS, Nov, 1999.
    ;-
    ;
    ; on_error,2              ;Return to caller if an error occurs
    """
    x = args[0]
    n = x.size
    if n < 3 : raise Exception('Parameters must have at least 3 points')


    if (len(args) == 2) :
        y=args[1]
        if n != y.size : raise 'Vectors must have same size'
    
        #;df/dx = y0*(2x-x1-x2)/(x01*x02)+y1*(2x-x0-x2)/(x10*x12)+y2*(2x-x0-x1)/(x20*x21)
        #; Where: x01 = x0-x1, x02 = x0-x2, x12 = x1-x2, etc.
    
        if isinstance(x,np.ma.masked_array) :  x = x.data                # Convert masked arrays to classic arrays
        if not isinstance(x,np.float) : x.astype(np.float) #;If not floating type, ensure floating... 
        
        x12 = x - np.roll(x,-1)                                      #;x1 - x2
        x01 = np.roll(x,1) - x                                       #;x0 - x1
        x02 = np.roll(x,1) - np.roll(x,-1)                           #;x0 - x2
    
        d = np.roll(y,1) * (x12 / (x01*x02)) \
            + y * (1./x12 - 1./x01) \
            - np.roll(y,-1) * (x01 / (x02 * x12))                    #Middle points
        
        
        #Formulae for the first and last points:
        d[0] = y[0] * (x01[1]+x02[1])/(x01[1]*x02[1]) \
            - y[1] * x02[1]/(x01[1]*x12[1]) \
            + y[2] * x01[1]/(x02[1]*x12[1])                          #;First point
        n2 = n-2
        d[n-1] = -y[n-3] * x12[n2]/(x01[n2]*x02[n2]) \
            + y[n-2] * x02[n2]/(x01[n2]*x12[n2]) \
            - y[n-1] * (x02[n2]+x12[n2]) / (x02[n2]*x12[n2])             #;Last point

    #Equally spaced point case
    else :
        d = (np.roll(x,-1) - np.roll(x,1))/2.
        d[0] = (-3.0*x[0] + 4.0*x[1] - x[2])/2.
        d[n-1] = (3.*x[n-1] - 4.*x[n-2] + x[n-3])/2.

    return d


def mask2NaN(array):
    n=array.size
    if array.mask.size != n :
        array.mask = np.zeros(n,dtype=bool)
#        raise np.ma.MaskError("[mask2NaN]Error : mask length is not consistent with data")
#        array.mask = np.ones(n,dtype=bool)
    array.data[np.arange(n).compress(array.mask)] = np.NaN
    return array


def histogram_indices(hist,R):
    ind = []
    for k in np.arange(len(hist)) : ind.append(R[R[k] : R[k+1]]) if hist[k] > 0 else ind.append([])
    return ind
#    for k in notempty : ind.append(R[R[k] : R[k+1]])



def nanargmin(array,axis=None):
    shape=array.shape
    nx=shape[0]
    ny=shape[1]
    if axis is None : return np.nanargmin()

def nearest(t, x):
    adiff=np.abs(t-np.float(x))
    i=np.argmin(adiff)
    return i

def cart2polar(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return r, theta

def polar2cart(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

def rad2geo(alpha):
    theta=(np.rad2deg(alpha)*-1)+360+90
    theta=np.mod(theta,360)
    return theta

def cart2geo(u,v):
    spd,dir=cart2polar(u, v)
    dir=rad2geo(dir)
    return spd,dir

def rms(array):
    """
    ;+
    ; RMS : Returns root-mean-squared deviation of an array
    ; 
    ; @param array {in}{required}{type:NUMERIC} 2-dimensionnal array. Time dimension<br />
    ;   should be the last dimension.
    ;   
    ; @returns RMS array (1 value per time serie)
    ; 
    ; @author Renaud DUSSURGET, LEGOS/CTOH
    ;-
    """
    #;time should be the last dimension
  
    n=array.size
    
#    IF sz[ndims-1] EQ 1 THEN RETURN, 0.
    
    nans = np.where(~np.isnan(array))
    cnt=len(nans)
    
    if (cnt == 0) : return np.NaN
    
    nval=np.nansum(~np.isnan(array))
    
    mn = np.nansum(array) / nval

    return np.sqrt(np.nansum((array - mn)**2.0)/nval)
#    RETURN, SQRT((TOTAL((array - REBIN(mn,sz))^2D,ndims,/DOUBLE,/NAN) / nval))
#;  RETURN, SQRT((TOTAL((array - REBIN(mn,sz))^2D,ndims,/DOUBLE,/NAN) / sz[ndims -1]))

def get_caller(level=2):
        frame=inspect.currentframe()
        for l in np.arange(level) : frame=frame.f_back
        code=frame.f_code
        return code

def username():
    return getpass.getuser()

def hostname():
    return socket.gethostbyaddr(socket.gethostname())[0]

def message(MSG_LEVEL,str,verbose=1):
    """
     MESSAGE : print function wrapper. Print a message depending on the verbose level
     
     @param MSG_LEVEL {in}{required}{type=int} level of the message to be compared with self.verbose
     
     @example self.message(0,'This message will be shown for any verbose level') 
    """
    
    caller=get_caller()
    if MSG_LEVEL <= verbose :  print('[{0}.{1}()] {2}'.format(__name__,caller.co_name,str))

def warning(MSG_LEVEL,str,verbose=1):
    """
     WARNING : Wrapper to the warning function. Returns a warning when verbose level is not 0. 
     
     @param MSG_LEVEL {in}{required}{type=int} level of the message to be compared with self.verbose
     
     @example self.waring(1,'Warning being issued) 
    """
    
    if verbose <= 1 : warn(str)
