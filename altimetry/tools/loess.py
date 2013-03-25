# -*- coding: utf-8 -*-
import threading
import datetime
import numpy as np
import Queue
from altimetry.tools.others import deriv, mask2NaN

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
    def __init__(self, input_q,indices_q,result_q):
        threading.Thread.__init__(self)
        self.input = input_q
        self.indices = indices_q
        self.result = result_q
    
    def task(self,NaN=True):
        
        #grabs host from queue
        inargs = self.input.get()
        id = inargs[0]
        h = inargs[1]
        x = inargs[2]
        cut = inargs[3]
        edges = inargs[4]
        
        
        if verbose > 0 : print "%s started at time: %s" % (self.getName(), datetime.datetime.now())
        
        try :
            res=loess(h, x, cut)
        except (ValueError):
            res=loess(h, x, cut) #Retry if problem
               
        self.indices.put(id)
        self.result.put(res[edges[0]:edges[1]])

        del res
        
        if verbose > 0 : print "%s ended at time: %s" % (self.getName(), datetime.datetime.now())
    
    def run(self):
        
            #Starts the queue
            self.task()
            #signals to queue job is done
            self.input.task_done()
       

def loess_parallel(h,x,cut,N_thread=4):
    
    
    #Setup input and output queues
    input_q = Queue.Queue()
    indices_q = Queue.Queue()
    result_q = Queue.Queue()
    
    #Get dimensions to split matrix
    n=len(h)
    
    #Map the data along X axis
    step = np.ceil(np.float(n) / N_thread)
    if N_thread *step < n : raise '[ERROR] splitting of time series is not valid'
    dx=np.median(deriv(x))
    over=np.ceil(cut/dx) #Overlay between each time series processed in parallel
    
    #If the serie is not long enough, go back to sequential version
    if step <= N_thread :
        if verbose > 1 : print "[WARNING] series too short. back to sequential code"
        return loess(h, x, cut)
    
    if ~np.isfinite(over) :
        if verbose > 1 : print "[WARNING] not valid overlay. returning raw data"
        return(h) 
    
    #Check data consistency
    
    
    indices = []
    ind_out = []
    for i in np.arange(N_thread) :
        try :
            dumind = np.arange(i*step - over,i*step+step+over,dtype=np.long)
        except (ValueError):
            print 'pb'
        dumind=dumind.compress((dumind >= 0) & (dumind <= n - 1))
        indices.append(dumind)
        dumind=np.arange(len(dumind))[(dumind >= i*step) & (dumind <= (i*step + step - 1))]
        ind_out.append([dumind.min(),dumind.max()+1]) #Pb with NP indexation... Add + 1 to get the last element

    #Check output indexes with initial length of data
    if np.sum(np.transpose(ind_out)[1] - np.transpose(ind_out)[0]) != n :
        raise '[ERROR] array reconstruction is not valid'

    #spawn a pool of threads, and pass them queue instance 
    for i in range(N_thread):
        t = ThreadClass(input_q,indices_q,result_q)
        t.setDaemon(True)
        t.start()
    
    #Feed threads with data
    for i in range(N_thread):
        input_q.put((i,h[indices[i]],x[indices[i]],cut,ind_out[i]))

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
        if i == tsort[0] : outmat = dum[i]
        else : outmat=np.ma.concatenate((outmat,dum[i]),0) if isinstance(h,np.ma.masked_array) else np.concatenate((outmat,dum[i]),0)
    
    if len(outmat) != len(h) :
        raise '[ERROR]Output array is not coherent with input array - check array reconstruction'
    return(outmat)

def loess(h,x,cut):
    """
    #===========================================================================
    # ;+
    # ; LOESS : Filter a signal H along X with a cutoff frequency CUT<br /><br />
    # ; 
    # ; Reference : Schlax, M. G., et D. B. Chelton (1992), Frequency Domain Diagnostics<br />
    # ;   for Linear Smoothers, Journal of the American Statistical Association, <br />
    # ;   87(420), 1070-1081. 
    # ; 
    # ; @param h {in}{required}{type:NUMERIC} Array of height anomalies
    # ; @param x {in}{required}{type:NUMERIC} Array of corresponding distances
    # ; @param cut {in}{required}{type:NUMERIC} Cutoff distance in units of X.
    # ; @param nan {in}{optional}{type:LOGICAL}{default:FALSE} Use NaNs<br />
    # ;          
    # ; @returns Hfilt (H filtered at cutoff CUT)
    # ; 
    # ; 
    # ; @uses exist
    # ; 
    # ; @author Renaud DUSSURGET, LEGOS/CTOH
    # ; @history Created Jan. 2010 by R.DUSSURGET from loess.c develloped by Mog2d team<br />
    # ;          May 2010 : Take NaNs into account for computation<br />
    # ;                     Checked the filters formulation (Schlax and Chelton, 1992) <br /><br />
    # ;
    # ;  LEGACY :<br />
    # ;  void dloess1d_irregular(int n, double l_c, double *h, double *x, double mask,double *lf,int *rstatus)<br /><br /> 
    # ;  
    # ;  ./varReduction.cpp:          dloess1d_irregular(n,cut,h,t,(double) satData.mask,lf,&status);<br />
    # ;  ./functions.extraction.cpp:  dloess1d_irregular(n,cut,h2,x,mask,lf,&status);<br /><br />
    # ;  
    # ;  rq: n-> nb points dans s�rie<br />
    # ;  l_c -> cut<br />
    # ;  *h -> sla<br />
    # ;  *x -> dimensional array<br />
    # ;  mask -> valeur du masque<br />
    # ;  *lf -> partie basse freq<br />
    # ;  *status -> retour<br />
    # ;-
    #===========================================================================
    """
          
    n=np.size(h)

    #If masked_array --> Get data and replace masked values by NaNs
    if isinstance(h,np.ma.masked_array): h=mask2NaN(h)
    else : h=np.ma.masked_array(h,mask=np.zeros(n,dtype=bool))

    #Start filtering if there are more than 1 point
    if n == 1 : lf = h[0]
    else : pass
    
    #Get the half-width of filter window
    l_c=cut/2.0          
    
    #Initiate weights and outputs
    w=np.zeros(n)
    lf = np.repeat(np.NaN,n)
    
    #Get valid points in serie
    flag = ~np.isnan(h)
    fcnt = flag.sum()
    fnt = np.arange(n).compress(flag)
    
    x = np.ma.masked_array(x,mask=(np.isnan(h) + h.mask) )
#        WHERE(FINITE(h),fcnt)

    #Vectorized form
    q = np.reshape(np.repeat(x,n),(n,n))
    q = (np.abs(q - q.transpose())/l_c).reshape(n*n)
    
    s = (1.0 - q*q*q)
    
#    outOfFilter_flag = q > 1.0
    outOfFilter_flag = np.greater(q,1)
    outCnt = np.nansum(outOfFilter_flag)
    outOfFilter = np.arange(n*n).compress(outOfFilter_flag)
    if (outCnt > 0) :
        s[outOfFilter]=0.0
        s.mask[outOfFilter]=True
    
    w=(s*s*s).reshape((n,n))  
    sumvar=np.nansum(w,1)

    #Compute current height (heights times weigths, divided by sum of weights)
    try :
        lf=np.nansum(w*np.repeat(h,n).reshape((n,n)),0)/sumvar
    except (ValueError):
        lf=np.nansum(w*np.repeat(h,n).reshape((n,n)),0)/sumvar #Retry

    return lf


def loess_inline(h,x,cut,nan=True):
    """
    #This is a raw, inline version of the loess filtering function. Runs much more slowlier.
    """
    
    n=np.size(h)

    #Start filtering if there are more than 1 point
    if n == 1 : lf = h[0]
    else : pass
    
    #Get the half-width of filter window
    l_c=cut/2.0          
    
    #Initiate weights and outputs
    w=np.zeros(n)
    lf = np.repeat(np.NaN,n)
    
    #Get valid points in serie
    flag = ~np.isnan(h)
    fcnt = flag.sum()
    fnt = np.arange(n).compress(flag)
#        WHERE(FINITE(h),fcnt)

    #Loop from 1st to last valid point
    for i in fnt :

        icur=i
      
        #Get Q (distance from central point divided by half-window width)
        q=(np.abs((x-x[icur])/l_c))
            
        #Compute S from Q�
        s = (1.0 - q*q*q)
        
     
        #Set all values out of filter window to 0
        outOfFilter_flag = q > 1.0
        outCnt = outOfFilter_flag.sum() 
        outOfFilter = np.arange(n).compress(outOfFilter_flag)
        if (outCnt > 0) : s[outOfFilter]=0.0
      
        #Compute weights (S� - Tricube)
        w=s*s*s  
        sumvar=np.nansum(w)

        #Compute current height (heights times weigths, divided by sum of weights)
        lf[icur]=np.nansum(w*h)/sumvar

    return lf
