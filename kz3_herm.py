import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from config import *
import new_dd as dd
from ev_diags import *
from spectra import *
from nlt_diags import *
#from landau_tests import *
import os
from sys import argv

def diags_cp_cm (x,y):
    dd.read_parameters()
    diagdir = '/scratch/01658/drhatch/dna_out'
    par['diagdir']="\'"+diagdir+"\'"

    time = dd.get_time_from_gout()
    kx,ky,kz,hermiteNumbers = dd.get_grids()

    g_in = dd.read_time_step_g(len(time)-1)
    g_in=np.reshape(g_in,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')

    x = x
    y = y
    dist = g_in[x,y,:,:] #distribution function: 2d slice of 4d g_in

    #make a new 2d array with the same dimensions as dist_
    shape = np.shape(dist)
    dist_t = np.empty(shape)

    for i in range(len(kz)):
        for j in range(len(hermiteNumbers)):
    		dist_t[i,j] = (1j*np.sign(kz[i]))**hermiteNumbers[j] * dist[i,j]

    dist_tp = np.empty(shape)
    dist_tm = np.empty(shape)

    for i in range(len(kz)):
	for j in range(len(hermiteNumbers)+1):
 		if j < (len(hermiteNumbers)-1):
    			dist_tp[i,j]= (dist_t[i,j]+dist_t[i,j+1])/2
    			dist_tm[i,j]= (dist_t[i,j]-dist_t[i,j+1])/2

    c_p=np.empty(np.shape(dist_tp))
    c_m=np.empty(np.shape(dist_tm))

    for i in range(len(kz)):
    	for j in range(len(kz)):
            	c_p[i,j]= dist_tp[i,j] * np.conjugate(dist_tp[i,j])
            	c_m[i,j]= dist_tm[i,j] * np.conjugate(dist_tm[i,j])

    c_ps = np.sum(c_p,axis=0)
    c_ms = np.sum(c_m,axis=0)
    xf = kx[x]
    yf = ky[y]
    return c_ps,c_ms, xf, yf

def diags_kz(x,y,z):
     dd.read_parameters()
     diagdir = '/work/01658/drhatch/prod_nu05_omt7.5_zf0b'
     par['diagdir']="\'"+diagdir+"\'"

     time = dd.get_time_from_gout()
     kx,ky,kz,hermiteNumbers = dd.get_grids()

     g_in = dd.read_time_step_g(len(time)-1)
     g_in=np.reshape(g_in,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')

     x = x
     y = y
     dist = g_in[x,y,:,:] #distribution function: 2d slice of 4d g_in

             #make a new 2d array with the same dimensions as dist_
     shape = np.shape(dist)
     dist_t = np.empty(shape)

     for i in range(len(kz)):
         for j in range(len(hermiteNumbers)):
             dist_t[i,j] = (1j*np.sign(kz[i]))**hermiteNumbers[j] * dist[i,j]

     dist_tp = np.empty(shape)
     dist_tm = np.empty(shape)

     for i in range(len(kz)):
         for j in range(len(hermiteNumbers)+1):
             if j < (len(hermiteNumbers)-1):
                 dist_tp[i,j]= (dist_t[i,j]+dist_t[i,j+1])/2
                 dist_tm[i,j]= (dist_t[i,j]-dist_t[i,j+1])/2

     c_p=np.empty(np.shape(dist_tp))
     c_m=np.empty(np.shape(dist_tm))

     for i in range(len(kz)):
         for j in range(len(kz)):
             c_p[i,j]= dist_tp[i,j] * np.conjugate(dist_tp[i,j])
             c_m[i,j]= dist_tm[i,j] * np.conjugate(dist_tm[i,j])

     z = z
     cps = c_p[z,:]
     cms = c_m[z,:]
     xf = kx[x]
     yf = ky[y]
     zf = kz[z]
     return cps,cms,hermiteNumbers, xf,yf,zf
