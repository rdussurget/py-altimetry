import numpy as np
import alti_tools as atools

def date2seas( date, seas=['DJF','MAM','JJA','SON'], shift=1):
    """
    #===============================================================================
    # ;+
    # ;
    # ; CDO_date2seas : provides the equivalent CDO season to a julian date array
    # ;
    # ; @param date {in}{required}{type=NUMERIC} julian date (CNES days)
    # ; @keyword seas {in}{optional}{type=STRING} seasons vector !! CAUTION !! This <br />
    # ;   must be accompanied by the SHIFT variable
    # ; @keyword shift {in}{optional}{type=NUMERIC} Backward in time offset of the 
    # ;   first element of SEAS vector in number of months from January <br />
    # ;   (e.g. Decembre -> 1, Novembre -> 2, etc... )
    # ; 
    # ;
    # ;-
    #===============================================================================
    """
        

    corr_months=np.roll(np.arange(12).reshape((4,3))+1,shift)
    
    nt = len(date)
    outvec=(atools.cnes_convert(date))[1]
    month=np.array([d.month for d in outvec])
    outseas_ind=-np.ones(nt)
#    outseas=np.repeat('   ',nt)
    
    for i,m in enumerate(month) : outseas_ind[i] = int(np.where(m == corr_months)[0])
    outseas =[seas[int(i)] for i in outseas_ind]
    
    return outseas, outseas_ind