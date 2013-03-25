# -*- coding: utf-8 -*-
'''
Created on 30 janv. 2013

@author: rdussurg
'''
import numpy as np
from altimetry.tools import interp1d, calcul_distance, deriv, gravity, coriolis
if __debug__ : import matplotlib.pyplot as plt

def uvgrid(*args,**kwargs) : #;lon, lat, time, sla, STRICT=True):
    
    lon=args[0]
    lat=args[1]
    
    if len(args) == 3:
        sla=args[2]
        time=np.arange(1)
    else :
        time=args[2]
        sla=args[3]
    
    strict=kwargs.get('strict',False)
    
    nx=len(lon)
    ny=len(lat)
    nt=len(time)
    
    sla=sla.reshape((nt,ny,nx))
    
    nxout= nx - 1 if strict else nx
    nyout= ny - 1 if strict else ny
    
    dx=np.median(deriv(lon))
    dy=np.median(deriv(lat))
    
    
    #Compute spatial distance gradients (convert to meters)
    xgrad=np.repeat([calcul_distance(l,0.0,l,dx)*1e3 for l in np.float64(lat)],nx).reshape((ny,nx))
    ygrad=np.repeat(calcul_distance(0,0,dy,0)*1e3,nx*ny).reshape((ny,nx)) #meridional distance gradient is constant everywhere
  
    lonout = interp1d(np.arange(nx),lon,np.arange(nxout)+0.5) if strict else lon
    latout = interp1d(np.arange(ny),lat,np.arange(nyout)+0.5) if strict else lat
    
    glon,glat=np.meshgrid(lonout, latout)
    
    g = gravity(glat)
    f = coriolis(glat)
  
    #Init gradients
    dhx=np.ma.array(np.zeros((nt,nyout,nxout)),mask=True)
    dhx.data[:]=dhx.fill_value
    dhy=dhx.copy();dhdx=dhx.copy(); dhdy=dhx.copy(); u=dhx.copy(); v=dhx.copy();
    
    #Loop over time
    for i in np.arange(nt) :
        
        #2 points differenciation
        if strict :
            dhx[i,:,:] = np.diff(sla[i,:,:], axis=1) 
            dhy[i,:,:] = np.diff(sla[i,:,:], axis=0)
        else :
            #3 points differenciation
            dhy[i,:,:], dhx[i,:,:] = np.gradient(sla[i,:,:])
    
        dhdx[i,:,:] = dhx[i,:,:] / xgrad
        dhdy[i,:,:] = dhy[i,:,:] / ygrad
    
        u[i,:,:]=-(g*dhdy[i,:,:]) / (f) #m.s-1
        v[i,:,:]=(g*dhdx[i,:,:]) / (f) #m.s-1
    
    u=np.squeeze(u)
    v=np.squeeze(v)
    
    
    return u,v
