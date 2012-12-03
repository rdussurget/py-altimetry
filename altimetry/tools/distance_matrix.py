# -*- coding: utf-8 -*-
import threading
import datetime
import numpy as np
import Queue

#import alti_tools as atools


##Global variables (rq: global keyword should be called by a function to set/modify global vars
#global x; x=10000L
#global y; y=2500L
#
#global A; A=np.matrix(np.repeat(np.NaN,x*y).reshape((x,y)))
#global map; map = np.arange(x)


#Init global variables
#dist_matrix=np.empty((0,0),dtype=np.float64)
rt = 6378.137 #Earth radius (km)
verbose = 0

#Computation thread
####################
#Queued thread
##############
class ThreadClass(threading.Thread):
    
    #We override the __init__ method
    def __init__(self, input_q,indices_q,result_q, alat, alon, blat, blon):
        threading.Thread.__init__(self)
        self.input = input_q
        self.indices = indices_q
        self.result = result_q
        self.alat=alat
        self.alon=alon
        self.blat=blat
        self.blon=blon
        del alat,alon,blat,blon #clean memory
        self.na = len(self.alon)
        self.nb = len(self.blon)
    
    def task(self):
        
        #grabs host from queue
        input = self.input.get()
        index = input[1]
        id = input[0]
        
        if verbose > 0 : print "%s started at time: %s" % (self.getName(), datetime.datetime.now())

        #Map the first lon/lat pairs
        alat = self.alat[index]
        alon = self.alon[index]
        
        na = len(alat)
        nb = len(self.blat)
        
        dum = np.tile(np.deg2rad(alat), (nb, 1))
        dumdst = np.cos(dum)
        interm = np.sin(dum)
        
        del dum
        dum =  np.tile(np.deg2rad(self.blat).transpose(), (na, 1)).transpose()
        dumdst *= np.cos(dum)
        interm *= np.sin(dum)
        
        del dum
        dum = np.tile(np.deg2rad(alon), (nb, 1))
        dum -= np.tile(np.deg2rad(self.blon).transpose(), (na, 1)).transpose()
        dumdst *= np.cos(dum)
        
        del dum
        
        dumdst += interm
        dumdst=rt * np.arccos(dumdst)
        
        del interm
        
        self.indices.put(id)
        self.result.put(dumdst)
#        lock = threading.Lock()
#        lock.acquire()
##        dist_matrix[:,index] = dumdst
#        lock.release()
        del dumdst
#        return(dumdst)
        
        if verbose > 0 : print "%s ended at time: %s" % (self.getName(), datetime.datetime.now())
    
    def run(self):
        
            #Starts the queue
            self.task()
            #signals to queue job is done
            self.input.task_done()
       

def distance_matrix_parallel(alat,alon,blat,blon,N_thread=4):
    
    
    #Setup input and output queues
    input_q = Queue.Queue()
    indices_q = Queue.Queue()
    result_q = Queue.Queue()
    
    #Get dimensions to split matrix
    x=len(alat)
    y=len(blat)
    
    
    #Map the data along X axis
    step = np.ceil(np.float(x) / N_thread)
    if N_thread *step < x : raise '[ERROR] splitting of time series is not valid'
    
    #If the serie is not long enough, go back to sequential version
    if step <= N_thread :
        if verbose > 1 : print "[WARNING] series too short. back to sequential code"
        return distance_matrix(alat,alon,blat,blon)
    
    indices = []
    for i in np.arange(N_thread) : indices.append(np.arange(i*step,i*step+step))
    
    indices = []
    for i in np.arange(N_thread) :
        dumind = np.arange(i*step,i*step+step,dtype=np.long)
        dumind=dumind.compress((dumind >= 0) & (dumind <= x - 1))
        indices.append(dumind)
        
    
    #redefine output matrix
#    global dist_matrix
#    dist_matrix = np.zeros((x,y))    
    
    #spawn a pool of threads, and pass them queue instance 
    for i in range(N_thread):
        t = ThreadClass(input_q,indices_q,result_q,alat,alon,blat,blon)
        t.setDaemon(True)
        t.start()
    
    #Feed threads with data
    for i in range(N_thread):
        input_q.put((i,indices[i]))

    #wait on the queue until everything has been processed     
    input_q.join()
    
    #Sort threads
    tnb=[]
    for i in range(N_thread) : tnb.append(indices_q.get(i))
    tsort=np.argsort(tnb)
    
    #Get back the results for each thread in a list of results
    for i in np.arange(N_thread) :
        r=result_q.get(i)
        if i == 0 : dum = [r]
        else : dum.append(r)

    #Reorder data from each thread into output matrix
    for i in tsort :
        if i == tsort[0] : dist_matrix = dum[i]
        else : dist_matrix=np.append(dist_matrix,dum[i],1)
    
    if dist_matrix.shape != (y,x) :
        raise Exception('[ERROR]Output matrix is not coherent with input data - check matrix reconstruction')

    return(dist_matrix)

def distance_matrix_inline(lat_a, lon_a, lat_b, lon_b):
    
    na = np.size(lat_a)
    nb = np.size(lat_b)
    
    #Earth radius (km)
    rt = 6378.137
    
#    lon_a_tile = np.tile(np.deg2rad(lon_a), (nb, 1))
#    lat_a_tile = np.tile(np.deg2rad(lat_a), (nb, 1))
#    lon_b_tile = np.tile(np.deg2rad(lon_b).transpose(), (na, 1)).transpose()
#    lat_b_tile = np.tile(np.deg2rad(lat_b).transpose(), (na, 1)).transpose()

#    plt.figure()
#    plt.imshow(lon_a_tile,interpolation='NEAREST',aspect='auto')
#    plt.show()
#    
#    plt.figure()
#    plt.imshow(lon_b_tile,interpolation='NEAREST',aspect='auto')
#    plt.show()
    
#    interm = np.cos(lat_a_tile) * np.cos(lat_b_tile) * np.cos(lon_a_tile - lon_b_tile) + np.sin(lat_a_tile) * np.sin(lat_b_tile)
    
    #lower memory usage for huge matrices
    dum = np.tile(np.deg2rad(lat_a), (nb, 1))
    dist_matrix = np.cos(dum)
    interm = np.sin(dum)
    
    del dum
    dum =  np.tile(np.deg2rad(lat_b).transpose(), (na, 1)).transpose()
    dist_matrix *= np.cos(dum)
    interm *= np.sin(dum)
    
    del dum
    dum = np.tile(np.deg2rad(lon_a), (nb, 1))
    dum -= np.tile(np.deg2rad(lon_b).transpose(), (na, 1)).transpose()
    dist_matrix *= np.cos(dum)
    
    del dum
    
    dist_matrix += interm
    
    del interm
        
    dist_matrix = rt * np.arccos(dist_matrix)
    
    return dist_matrix

def distance_matrix(alat,alon,blat,blon,parallel=True,**kwargs):
    if parallel : return distance_matrix_inline(alat, alon, blat, blon,**kwargs)
    else : return distance_matrix_inline(alat, alon, blat, blon)