import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from config import *
import new_dd as dd
from ev_diags import *
import time as t2
from spectra import *
from nlt_diags import *
#from landau_tests import *
import os
from sys import argv
import kz3_herm as d2

dd.read_parameters()
diagdir = '/scratch/01658/drhatch/dna_out'
par['diagdir']="\'"+diagdir+"\'"

time = dd.get_time_from_gout()
kx,ky,kzgrid,herm_grid = dd.get_grids()

start_time=time[100]

end_time=time[len(time)-1]
if start_time >= end_time:
    stop

include_kz0=True

prefactor=np.empty(par['nv0'])
prefactor[:]=1.0
plabel='Entropy'

istart=np.argmin(abs(time-start_time))
iend=np.argmin(abs(time-end_time))/16
ntime=iend-istart+1
#iend => 1/16

entn_sum=np.zeros((par['nv0'],11),dtype='float')
entnp_sum=np.zeros((par['nv0'],11),dtype='float')
entnm_sum=np.zeros((par['nv0'],11),dtype='float')
for i in range(istart,iend+1):
    print 'time=',time[i],' of ',time[iend]
    gt0=dd.read_time_step_g(i)
    gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
    gt1= gt0
    ##################################

    dist = gt1[:,:,:,:] #distribution function: 2d slice of 4d g_in

    #make a new 2d array with the same dimensions as dist_
    shape = np.shape(dist)
    g_t = np.empty(shape)
    for v in range(len(kx)):
        for w in range(len(ky)):
            for b in range(len(kzgrid)):
                for j in range(len(herm_grid)):
                    g_t[v,w,b,j] = (1j*np.sign(kzgrid[b]))**herm_grid[j] * dist[v,w,b,j]

## hazeltines equation
## (- 1j*np.sign(kzgrid[b]))**herm_grid * np.exp((herm_grid[j]**(.5)*D/(np.absolute(nu)))* (Npl)**(-np.sign(kzgrid[b]*(herm_grid[j]+2*nu**(-2)-.5))/(herm_grid**(.25)*D**(.5))
    g_tp = np.empty(shape)
    g_tm = np.empty(shape)
    for b in range(len(kzgrid)):
        for j in range(len(herm_grid)):
            #D = (herm_grid[j]*nu**2/4 + 1)**(.5) , nu =.001/kz_grid[b] , npl = D + herm_grid[j]**(.5)*nu/2
            #CapG = (- 1j*np.sign(kzgrid[b]))**herm_grid[j] * np.exp((herm_grid[j]**(.5)*D/(np.absolute(nu)))* (Npl)**(-np.sign(kzgrid[b]*(herm_grid[j]+2*nu**(-2)-.5))/(herm_grid[j]**(.25)*D**(.5))
    for v in range(len(kx)):
        for w in range(len(ky)):
            for b in range(len(kzgrid)):
                for j in range(len(herm_grid)+1):
                    if j < (len(herm_grid)-1):
                        g_tp[v,w,b,j]= (g_t[v,w,b,j]+g_t[v,w,b,j+1])/2
                        g_tm[v,w,b,j]= (g_t[v,w,b,j]-g_t[v,w,b,j+1])/2

    ##################################
    #Entropy
    entn_sum[:,10]=entn_sum[:,10]+get_entropy_hermite(gt0,kzind=-1,include_kz0=include_kz0)
    for k in range(10):
        kzindex=k*par['nkz0']/20
        #print 'kzindex',kzindex
        entn_sum[:,k]=entn_sum[:,k]+dd.get_entropy_hermite(gt0,kzind=kzindex)
        entnp_sum[:,k]= entnp_sum[:,k]+dd.get_entropy_hermite(g_tp,kzind=kzindex)
        entnm_sum[:,k]= entnm_sum[:,k]+dd.get_entropy_hermite(g_tm,kzind=kzindex)

entn_sum=entn_sum/float(ntime)
entnp_sum=entnp_sum/float(ntime)
entnm_sum=entnm_sum/float(ntime)

plt.loglog(herm_grid,prefactor*entn_sum[:,10],basex=10,basey=10)

temp=prefactor*entn_sum[20,10]
temp=temp/(herm_grid**(-1))[20]
plt.loglog(herm_grid,2.0*temp*herm_grid**(-1),'--',basex=10,basey=10,label=str(-1))
temp=prefactor*entn_sum[20,10]
temp=temp/(herm_grid**(-1.5))[20]
plt.loglog(herm_grid,2.0*temp*herm_grid**(-1.5),'--',basex=10,basey=10,label=str(-1.5))

plt.xlabel('Hermite n')
plt.ylabel(plabel)
plt.title(plabel+'(kz sum)')
plt.legend(loc='lower left')
plt.show()

#plots = np.shape(entn_sum)
print "entn:", entn_sum
for k in range(10):
    #print prefactor
    #print prefactor*entn_sum[:,k]
    kz0=kzgrid[k*par['nkz0']/20]
    plt.loglog(herm_grid,prefactor*entn_sum[:,k],basex=10,basey=10,label=\
              plabel+' (k_z='+str(kzgrid[k*par['nkz0']/20])+')')
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.legend(loc='lower left')
plt.show()
print 'plus sum:', entnp_sum
print 'minus sum', entnm_sum
for k in range(10):
    #print prefactor
    #print prefactor*entn_sum[:,k]
    kz0=kzgrid[k*par['nkz0']/20]
    #plots.append(prefactor*entn_sum[:,k])
    plt.loglog(herm_grid,prefactor*entnp_sum[:,k],basex=10,basey=10,label=\
              plabel+' (k_z='+str(kzgrid[k*par['nkz0']/20])+')')
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.title('Entropy+')
    plt.legend(loc='lower left')
temp=prefactor*entn_sum[20,10]
temp=temp/(herm_grid**(-1))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-1),'--',basex=10,basey=10,label=str(-1))
temp=prefactor*entn_sum[20,10]
temp=temp/(herm_grid**(-1.5))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-3.5),'--',basex=10,basey=10,label=str(-3.5))
plt.loglog(herm_grid, (10)*herm_grid**(-2),'--',basex=10,basey=10,label='n^(-2)')
plt.loglog(herm_grid, (10**-2.5)*herm_grid**(-.5),'--',basex=10,basey=10,label='n^(-.5)')
plt.legend(loc='lower left')
plt.show()

##################
#try log(value * herm_grid for intercept)
#####################
for k in range(10):
    #print prefactor
    #print prefactor*entn_sum[:,k]
    kz0=kzgrid[k*par['nkz0']/20]
    #plots.append(prefactor*entn_sum[:,k])
    plt.loglog(herm_grid,prefactor*entnm_sum[:,k],basex=10,basey=10,label=\
              plabel+' (k_z='+str(kzgrid[k*par['nkz0']/20])+')')
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.title('Entropy-')
    plt.legend(loc='lower left')
temp=prefactor*entn_sum[20,10]
temp=temp/(herm_grid**(-1))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-1),'--',basex=10,basey=10,label=str(-1))
temp=prefactor*entn_sum[20,10]
temp=temp/(herm_grid**(-1.5))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-3.5),'--',basex=10,basey=10,label=str(-3.5))
plt.loglog(herm_grid, 10*herm_grid**(-3/2),'--',basex=10,basey=10,label='n^(-.5)')
plt.loglog(herm_grid, 10**(-2.7)*herm_grid**(-3/2),'--',basex=10,basey=10,label='n^(-.5)')
plt.legend(loc='lower left')
plt.show()

 #split entropy into plus and minus, entropy is g^2
