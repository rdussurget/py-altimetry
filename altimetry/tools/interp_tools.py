# -*- coding: utf-8 -*-
import numpy as np
import scipy
import threading
import Queue
import bisect
import operator

def lagrange(x, x_values, y_values):
    def _basis(j):
        p = [(x - x_values[m])/(x_values[j] - x_values[m]) for m in xrange(k) if m != j]
        return reduce(operator.mul, p)
    assert len(x_values) != 0 and (len(x_values) == len(y_values)), 'x and y cannot be empty and must have the same length'
    k = len(x_values)
    return sum(_basis(j)*y_values[j] for j in xrange(k))

def bilinear_interpolation(x, y, points):
    '''Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

        >>> bilinear_interpolation(12, 5.5,
        ...                        [(10, 4, 100),
        ...                         (20, 4, 200),
        ...                         (10, 6, 150),
        ...                         (20, 6, 300)])
        165.0
        
        @note : source -> http://stackoverflow.com/questions/8661537/how-to-perform-bilinear-interpolation-in-python
    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    points = sorted(points)               # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

#    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
#        raise ValueError('points do not form a rectangle')
#    if not x1 <= x <= x2 or not y1 <= y <= y2:
#        raise ValueError('(x, y) not within the rectangle')

    try :
        if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
            raise ValueError('points do not form a rectangle')
        if not x1 <= x <= x2 or not y1 <= y <= y2:
            raise ValueError('(x, y) not within the rectangle')
        return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
           ) / ((x2 - x1) * (y2 - y1) + 0.0)
    except ValueError : return np.nan

def interp1d(x,Z,xout,spline=False,kind='linear',**kwargs):
    """    
    INTERP1D : Interpolate values from a 1D vector at given positions
    
    @param x: 1st dimension vector of size NX
    
    @author: Renaud DUSSURGET, LER/PAC, Ifremer La Seyne
    """

    linear = not spline
    
    
    nx=len(x)
    
    if linear :
        
        try :
            f = scipy.interpolate.interp1d(x, Z, kind=kind,bounds_error=False)
            Zout = f(xout)
        except RuntimeError : Zout = np.repeat(np.NaN,nx)

    else :
        tck = scipy.interpolate.splrep(x,Z,s=0)
        try : Zout = scipy.interpolate.splev(xout,tck,der=0)
        except RuntimeError : Zout = np.repeat(np.NaN,nx)

    return Zout

def interp2d1d(x,y,Z,xout,yout,**kwargs):
    """    
    INTERP2D1d : Interpolate values from a 2D matrix along a 1D trajectory
    
    @param x: 1st dimension vector of size NX
    @param y: 2nd dimension vector of size NY
    @param Z: Array to interpolate (NXxNY)
    @param xout: 1st dimension vector of size Nout
    @param yout: 2nd dimension vector of size Nout
    
    @return: Interpolated array (Nout)
    
    @author: Renaud DUSSURGET, LER/PAC, Ifremer La Seyne
    """

    gx = np.reshape(np.repeat(x,y.size),(x.size,y.size))
    gy = np.reshape(np.repeat(y,x.size),(y.size,x.size)).transpose((1,0))

    gxout = xout
    gyout = yout
    gz = Z.flatten()
     
    points = zip(*(gx.flatten(),gy.flatten())) 
    xi = zip(*(gxout,gyout))
    
    try : Zout = scipy.interpolate.griddata(points, gz, xi, **kwargs)
    except RuntimeError : Zout = np.NaN
#    Zout = sc.interpolate.griddata(points, gz, xi, **kwargs)
    return Zout      

def interp2d2d(x,y,Z,xout,yout,split_factor=2,**kwargs):
    """    
    INTERP2D2D : Interpolate a 2D matrix into another 2D matrix
    
    @param x: 1st dimension vector of size NX
    @param y: 2nd dimension vector of size NY
    @param Z: Array to interpolate (NXxNY)
    @param xout: 1st dimension vector of size NXout
    @param yout: 2nd dimension vector of size NYout
    
    @keyword split_factor:Nummber of times to split arrays. Nb of threads is equal to split_factor**2. 
    
    @return: Interpolated array (NXoutxNYout)
    
    @author: Renaud DUSSURGET, LER/PAC, Ifremer La Seyne
    """
    
    #Queued thread
    ##############
    class ThreadClass(threading.Thread):
        
        #We override the __init__ method
        def __init__(self, input_q,indices_q,result_q):
            threading.Thread.__init__(self)
            self.Z = input_q
#            self.indices = indices_q
            self.result = result_q
        
        def task(self,NaN=True):
            
            #grabs host from queue
            inargs = self.input.get()
            x = inargs[0]
            y = inargs[1]
            Z = inargs[2]
            xout = inargs[3]
            yout = inargs[4]
            
            #Get dimensions to split matrix
            nxin=x.size
            nyin=y.size
            nxout=xout.size
            nyout=yout.size
            
#            if verbose > 0 : print "%s started at time: %s" % (self.getName(), datetime.datetime.now())
            
            gx = np.reshape(np.repeat(x,nyin),(nxin,nyin))
            gy = np.reshape(np.repeat(y,nxin),(nyin,nxin)).transpose((1,0))
        
            gxout = np.reshape(np.repeat(xout,nyout),(nxout,nyout)).flatten()
            gyout = np.reshape(np.repeat(yout,nxout),(nyout,nxout)).transpose((1,0)).flatten()
        
            gz = Z.flatten()
             
            points = zip(*(gx.flatten(),gy.flatten())) 
            xi = zip(*(gxout,gyout))
            
            try : Zout = scipy.interpolate.griddata(points, gz, xi, **kwargs)
            except RuntimeError : Zout = np.ones((xout.size,yout.size))*np.NaN
        #    Zout = sc.interpolate.griddata(points, gz, xi, **kwargs)
            return Zout.reshape((xout.size,yout.size))   
                   
#            self.indices.put(h)
            self.result.put(Zout)
    
            del Zout
            
#            if verbose > 0 : print "%s ended at time: %s" % (self.getName(), datetime.datetime.now())
        
        def run(self):
            
                #Starts the queue
                self.task()
                #signals to queue job is done
                self.input.task_done()
       

    #Setup input and output queues
    input_q = Queue.Queue()
    indices_q = Queue.Queue()
    result_q = Queue.Queue()
    
    #Map the data along X axis
    N_threads=split_factor**2
#    over=np.ceil(cut/dx) #Overlay between each time series processed in parallel
    
    #Get dimensions to split matrix
    nxin=x.size
    nyin=y.size
    nxout=xout.size
    nyout=yout.size
    
    gx = np.reshape(np.repeat(x,nyin),(nxin,nyin))
    gy = np.reshape(np.repeat(y,nxin),(nyin,nxin)).transpose((1,0))

    gxout = np.reshape(np.repeat(xout,nyout),(nxout,nyout))
    gyout = np.reshape(np.repeat(yout,nxout),(nyout,nxout)).transpose((1,0))
    
    #Map output coordinates
    xsplit=[]
    ysplit=[]
    for i in np.arange(split_factor) : 
        for j in np.arange(split_factor) :
            xsplit.append(i*(nxout/float(split_factor)))
            xsplit.append((i+1)*(nxout/float(split_factor)))
            ysplit.append(j*(nyout/float(split_factor)))
            ysplit.append((j+1)*(nyout/float(split_factor)))
    
    #Round
    xsplit=np.round(xsplit).astype(int)
    ysplit=np.round(ysplit).astype(int)
    
    N_threads = len(xsplit)/2
    
    
    
#    th_xout=gxout[0:nxout/split_factor,0:2]
#    th_yout=

#    gz = Z.flatten()
     
#    points = zip(*(gx.flatten(),gy.flatten())) 
    
    #spawn a pool of threads, and pass them queue instance 
    for i in np.arange(N_threads):
        t = ThreadClass(input_q,indices_q,result_q)
        t.setDaemon(True)
        t.start()
    
    #Feed threads with data
    for i in range(N_threads):
        xoind=np.array(xsplit)[[i*2,(i*2)+1]]
        yoind=np.array(ysplit)[[i*2,(i*2)+1]]
        xiind=x[xout[xoind].astype(int)].astype(int)
        yiind=y[yout[yoind].astype(int)].astype(int)
        th_x=gx[xiind[0]:xiind[1],yiind[0]:yiind[1]]
        th_y=gy[xiind[0]:xiind[1],yiind[0]:yiind[1]]
        th_xo=gx[xoind[0]:xoind[1],yoind[0]:yoind[1]]
        th_yo=gy[xoind[0]:xoind[1],yoind[0]:yoind[1]]
        th_z=Z[xiind[0]:xiind[1],yiind[0]:yiind[1]]
        print i
#        input_q.put((i,th_x,th_y,Z[indices[i]],x[indices[i]],cut,ind_out[i]))

#    #wait on the queue until everything has been processed     
#    input_q.join()
#    
#    #Sort threads
#    tnb=[]
#    for i in range(N_threads) : tnb.append(indices_q.get(i))
#    tsort=np.argsort(tnb)
#    
#    #Get back the results for each thread in a list of results
#    for i in np.arange(N_thread) :
#        r=result_q.get(i)
#        if i == 0 : dum = [r]
#        else : dum.append(r)
#
#    #Reorder data from each thread into output matrix
#    for i in tsort :
#        if i == tsort[0] : outmat = dum[i]
#        else : outmat=np.ma.concatenate((outmat,dum[i]),0) if isinstance(h,np.ma.masked_array) else np.concatenate((outmat,dum[i]),0)
#    
#    if len(outmat) != len(h) :
#        raise '[ERROR]Output array is not coherent with input array - check array reconstruction'
#    return(outmat)

def interp3d(x,y,t,Z,xout,yout,tout,**kwargs):
    
    nx=x.size
    ny=y.size
    nt=t.size
    
    #Turn around matrix to (nx,ny,nt)
    order=(Z.shape.index(nx),Z.shape.index(ny),Z.shape.index(nt))
    Z=Z.transpose(order)
    
    N=xout.size
    
    for enum in enumerate(zip(*(xout,yout,tout))):
        time=enum[1][2]
        lon=enum[1][0]
        lat=enum[1][1]
        idx = np.arange(bisect.bisect_left(x,lon)-1,bisect.bisect_left(x,lon)+1)
        idy = np.arange(bisect.bisect_left(y,lat)-1,bisect.bisect_left(y,lat)+1)
        idt = np.arange(bisect.bisect_left(t,time)-1,bisect.bisect_left(t,time)+1)
        
        idx=idx.compress((idx >= 0) & (idx < nx))
        idy=idy.compress((idy >= 0) & (idy < ny))
        idt=idt.compress((idt >= 0) & (idt < nt))
        nxid=idx.size
        nyid=idy.size
        ntid=idt.size
        
#        print enum[0]
        
        if ntid == 0 :
            res = np.NaN
        elif ((ntid == 1) & (nxid == 2) & (nyid == 2)):
            zi = np.squeeze(Z[idx,:,:][:,idy,:][:,:,idt])
#            res = interp2d(x[idx], y[idy], zi, [lon], [lat])
            res = bilinear_interpolation(lon, lat, zip(*(x[idx[[0,0,1,1]]],y[idy[[0,1,0,1]]], np.reshape(zi,4))))
        elif ((ntid == 2) & (nxid == 2) & (nyid == 2)) :
            zi = Z[idx,:,:][:,idy,:][:,:,idt]
#            a=interp2d(x[idx], y[idy], zi[:,:,0], [lon], [lat])
#            b=interp2d(x[idx], y[idy], zi[:,:,1], [lon], [lat])
            a = bilinear_interpolation(lon, lat, zip(*(x[idx[[0,0,1,1]]],y[idy[[0,1,0,1]]], np.reshape(zi[:,:,0],4))))
            b = bilinear_interpolation(lon, lat, zip(*(x[idx[[0,0,1,1]]],y[idy[[0,1,0,1]]], np.reshape(zi[:,:,1],4))))
            res = ((b-a)/(t[idt[1]] - t[idt[0]]))*(time - t[idt[0]]) + a #linear interpolation
#            res = interp3d_core(x, y, t[idt], zi, [enum[1][0]], [enum[1][1]], [time])
        else : res = np.NaN
        if enum[0] == 0 : Zout=res
        else : Zout=np.append(Zout,res)  
        
    return Zout

def interp3d_core(x,y,t,Z,xout,yout,tout,**kwargs):
    """    
    INTERP3D : Interpolate values from a 3D matrix along a 1D trajectory
    
    @param x: 1st dimension vector of size NX
    @param y: 2nd dimension vector of size NY
    @param t: 3rd dimension vector of size NT
    
    @author: Renaud DUSSURGET, LER/PAC, Ifremer La Seyne
    """
    
    #this below can take a very LOOOOOONG time
    gx = np.reshape(np.repeat(x,y.size*t.size),(x.size,y.size,t.size))
    gy = np.reshape(np.repeat(y,x.size*t.size),(y.size,x.size,t.size)).transpose((1,0,2))
    gt = np.reshape(np.repeat(t,x.size*y.size),(t.size,x.size,y.size)).transpose((1,2,0))

    gz = Z.flatten()
     
    points = zip(*(gx.flatten(),gy.flatten(),gt.flatten())) 
    xi = zip(*(xout,yout,tout))
    
    Zout = scipy.interpolate.griddata(points, gz, xi, **kwargs)
    return Zout    

