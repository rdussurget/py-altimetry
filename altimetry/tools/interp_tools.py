# -*- coding: utf-8 -*-
import numpy as np
import scipy.interpolate
import threading
import Queue
import bisect
import operator
import datetime
if __debug__ : import matplotlib.pyplot as plt

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

def extrap1d(Z,mask):
    """    
    EXTRAP1D : Extrapolate values from a 1D vector at its beginning and end using reversal of this array
    
    @note : gaps in vector Z should be filled first as values from the vector are replicated out of edges
    @note : if isinstance(Z,np.ma.masked_array) : mask = Z.mask
    
    @param Z: 1D vector to extrapolate
    @param mask: mask flag of the vector 
    
    @author: Renaud DUSSURGET, LER/PAC, Ifremer La Seyne
    """
    Zout=Z.copy()
    N=len(Zout)
    ind=np.arange(N)
    xout=ind[mask]
    hist=(~mask).astype(int)
    dhist=hist[1:]-hist[:-1]
    st=ind[dhist==1]+1
    if len(st) > 0 :
        st=st[0]
        Zout[:st]=Z[st] - (np.roll(Z,-st-1)[:st][::-1] - Z[st])
    en=ind[dhist==-1]
    if len(en) > 0 :
        en=en[-1]
        Zout[en+1:]=Z[en] - (np.roll(Z,-en)[::-1][:N-en-1] - Z[en])
    return Zout
    

def interp1d(x,Z,xout,spline=False,kind='linear',fill_value=np.NaN,**kwargs):
    """    
    INTERP1D : Interpolate values from a 1D vector at given positions
    
    @param x: 1st dimension vector of size NX
    
    @author: Renaud DUSSURGET, LER/PAC, Ifremer La Seyne
    """

    linear = not spline
    
    
    nx=len(x)
    
    if linear :
        
        try :
            f = scipy.interpolate.interp1d(x, Z, kind=kind,bounds_error=False,fill_value=fill_value,**kwargs)
            Zout = f(xout)
        except RuntimeError : Zout = np.repeat(np.NaN,nx)

    else :
        tck = scipy.interpolate.splrep(x,Z,s=0)
        try : Zout = scipy.interpolate.splev(xout,tck,der=0,**kwargs)
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

    gx = np.reshape(np.repeat(x,y.size),(x.size,y.size)).transpose((1,0))
    gy = np.reshape(np.repeat(y,x.size),(y.size,x.size))

    gxout = xout
    gyout = yout
    gz = Z.flatten()
    
    Zout=np.empty(Z.shape,dtype=Z.dtype)
    
#    if isinstance(Z,np.ma.masked_array):
#        mask=Z.mask
#        Zout=np.ma.array(Zout,mask=True)
#    else : mask=np.ones(gz.shape,dtype=bool)
     
#    points_masked = zip(*(gx[mask].flatten(),gy[mask].flatten()))
    points = zip(*(gx.flatten(),gy.flatten()))
#    gz=gz[mask.flatten()]
    xi = zip(*(gxout,gyout))
    
    try : Zout = scipy.interpolate.griddata(points, gz, xi, **kwargs)
    except RuntimeError : Zout = np.NaN
#    Zout = sc.interpolate.griddata(points, gz, xi, **kwargs)
    return Zout      


def interp2d1d_parallel(x,y,Z,xout,yout,split_factor=4,**kwargs):
    #Queued thread
    ##############
    class ThreadClass(threading.Thread):
        
        #We override the __init__ method
        def __init__(self, input_q,indices_q,result_q):#,points,gz,xout,yout):
            threading.Thread.__init__(self)
            
            self.input = input_q
            self.indices = indices_q
            self.result = result_q
            
#            self.points=points
#            self.gz=gz
#            
#            self.xout=xout
#            self.yout=yout
            
            self.verbose=True
            
            
        def task(self,NaN=True):
            
            
            
            #grabs host from queue
            input = self.input.get()
            id = input[0]
            #index = input[1]
            points=input[1]
            gz=input[2]
            xout=input[3]
            yout=input[4]
            
            #Get dimensions to split matrix
            
            
            xi = zip(*(xout,yout))
#            xi = zip(*(xout[index],yout[index]))
            
            if self.verbose: print "%s started at time: %s" % (self.getName(), datetime.datetime.now())
            try : Zout = scipy.interpolate.griddata(points, gz, xi, **kwargs) 
#            try : Zout = scipy.interpolate.griddata(self.points, self.gz, xi, **kwargs)
            except RuntimeError : Zout = np.NaN
            #    Zout = sc.interpolate.griddata(points, gz, xi, **kwargs)
            if self.verbose: print "%s ended at time: %s" % (self.getName(), datetime.datetime.now())
            
            self.indices.put(id)
            self.result.put(Zout)
            del Zout
            
#            return Zout.reshape((xout.size,yout.size))
                        
            
        
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
    N_threads=split_factor
    
#    over=np.ceil(cut/dx) #Overlay between each time series processed in parallel
    
    #Get dimensions to split matrix
    nxin=x.size
    nyin=y.size
    N=xout.size
    
    nxin=x.size
    nyin=y.size
    gx = np.reshape(np.repeat(x,nyin),(nxin,nyin)).transpose((1,0)).flatten()
    gy = np.repeat(y,nxin)
    gz = Z.flatten()
    
    points=zip(*(gx,gy))

    #Map output coordinates
    ind=[]
    
    for i in np.arange(split_factor) :
        ind.append(i*(N/float(split_factor)))
        ind.append((i+1)*(N/float(split_factor))-1)
    
    ind=np.array(ind)
    
    #Round
    ind[1::2]=np.ceil(ind[1::2]).astype(int)
    ind[0::2]=np.floor(ind[0::2]).astype(int) 
    
    N_threads = len(ind)/2    
    
    
#    th_xout=gxout[0:nxout/split_factor,0:2]
#    th_yout=

#    gz = Z.flatten()
     
#    points = zip(*(gx.flatten(),gy.flatten())) 
    
    #spawn a pool of threads, and pass them queue instance 
    for i in np.arange(N_threads):
        #t = ThreadClass(input_q,indices_q,result_q,points,gz,xout,yout)
        t = ThreadClass(input_q,indices_q,result_q)
        t.setDaemon(True)
        t.start()
    
    iind=[]
    outvar=np.ma.array(np.empty(N),mask=True,dtype=Z.dtype)
    outvar.data[:]=outvar.fill_value
    
    #Feed threads with data
    for i in range(N_threads):
        oind=np.array(ind)[[i*2,(i*2)+1]].astype(int)
        iind.append(np.arange(oind.min(),oind.max(),dtype=int))
        print i,oind,iind[i].max()
    
    for i in range(N_threads):
        print "%s launched time: %s" % (i, datetime.datetime.now())
        input_q.put((i,np.copy(points),gz.copy(),xout[iind[i]].copy(),yout[iind[i]].copy()))

    #wait on the queue until everything has been processed     
    input_q.join()
    
    #Sort threads
    tnb=[]
    for i in range(N_threads) : tnb.append(indices_q.get(i))
    tsort=np.argsort(tnb)

    #Get back the results for each thread in a list of results
    for i in np.arange(N_threads) :
        r=result_q.get(i)
        outvar[iind[tsort[i]]]=r.astype(outvar.dtype)

    if len(outvar) != len(xout) :
        raise Exception('[ERROR]Output matrix is not coherent with input data - check matrix reconstruction')

    return(outvar)
    


def interp2d2d(x,y,Z,xout,yout,split_factor=1,**kwargs):
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
            self.input = input_q
            self.indices = indices_q
            self.result = result_q
        
        def task(self,NaN=True):
            
            #grabs host from queue
            inargs = self.input.get()
            ind = inargs[0] #thread index
            x = inargs[1]
            y = inargs[2]
            Z = inargs[3]
            xout = inargs[4]
            yout = inargs[5]
            
            Zout=_interp2d2d(x,y,Z,xout,yout,**kwargs)
            
        #    Zout = sc.interpolate.griddata(points, gz, xi, **kwargs)
            
            self.indices.put(ind)
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
    ind=[]
    
    #Map output coordinates
    xsplit=[]
    ysplit=[]
    for i in np.arange(split_factor) : 
        for j in np.arange(split_factor) :
            xsplit.append(i*(nxout/float(split_factor)))
            xsplit.append((i+1)*(nxout/float(split_factor))-1)
            ysplit.append(j*(nyout/float(split_factor)))
            ysplit.append((j+1)*(nyout/float(split_factor))-1)
    
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
        xoind=xsplit[[i*2,(i*2)+1]]
        yoind=ysplit[[i*2,(i*2)+1]]
        xiind=x[xout[xoind].astype(int)].astype(int)
        yiind=y[yout[yoind].astype(int)].astype(int)
        th_x=x[xiind[0]:xiind[1]+1]
        th_y=y[yiind[0]:yiind[1]+1]
        th_xo=xout[xoind[0]:xoind[1]+1]
        th_yo=yout[yoind[0]:yoind[1]+1]
        th_z=Z[xiind[0]:xiind[1]+1,yiind[0]:yiind[1]+1]
#        input_q.put((i,np.copy(points),gz.copy(),xout[iind[i]].copy(),yout[iind[i]].copy()))
        input_q.put((i,th_x,th_y,th_z,th_xo,th_yo))

   
#    outvar=Z.copy()
#    outvar.data[:]=outvar.fill_value
    
#    for i in range(N_threads):
#        print "%s launched time: %s" % (i, datetime.datetime.now())
#        input_q.put((i,np.copy(points),gz.copy(),xout[iind[i]].copy(),yout[iind[i]].copy()))

    

    #wait on the queue until everything has been processed     
    input_q.join()
    
    #Sort threads
    tnb=[]
    for i in range(N_threads) : tnb.append(indices_q.get(i))
    tsort=np.argsort(tnb)
    
    #Get back the results for each thread in a list of results
    for i in np.arange(N_threads) :
        r=result_q.get(i)
        if i == 0 : dum = [r]
        else : dum.append(r)

    #Reorder data from each thread into output matrix
    for i in tsort :
        if i == tsort[0] : outmat = dum[i]
        else : outmat=np.ma.concatenate((outmat,dum[i]),0) if isinstance(outmat,np.ma.masked_array) else np.concatenate((outmat,dum[i]),0)
    
    if len(outmat) != len(outmat) :
        raise '[ERROR]Output array is not coherent with input array - check array reconstruction'
    return(outmat)

def _interp2d2d(x,y,Z,xout,yout,**kwargs):
            
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
            if isinstance(Z,np.ma.masked_array): gmask=gz.mask
             
            points = zip(*(gx.flatten(),gy.flatten())) 
            xi = zip(*(gxout,gyout))
            
            try : Zout = scipy.interpolate.griddata(points, gz, xi, **kwargs)
            except RuntimeError : Zout = np.ones((xout.size,yout.size))*np.NaN
            
            Zout=Zout.reshape((xout.size,yout.size))
            
            if isinstance(Z,np.ma.masked_array) :
                try : maskout = scipy.interpolate.griddata(points, gmask, xi, **kwargs).astype(bool)
                except RuntimeError : Zout = np.ones((xout.size,yout.size))*np.NaN
                Zout = np.ma.array(Zout,mask=maskout.reshape(xout.size,yout.size))
            
            return Zout

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

