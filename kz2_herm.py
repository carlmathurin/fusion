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

dd.read_parameters()
diagdir = '/scratch/01658/drhatch/dna_out'
par['diagdir']="\'"+diagdir+"\'"

time = dd.get_time_from_gout()
kx,ky,kzgrid,herm_grid = dd.get_grids()

start_time=time[0]
end_time=time[len(time)-1]
if start_time >= end_time:
    stop

include_kz0=True

prefactor=np.empty(par['nv0'])
prefactor[:]=1.0
plabel='Entropy'

istart=np.argmin(abs(time-start_time))
iend=np.argmin(abs(time-end_time))
ntime=iend-istart+1


entn_sum=np.zeros((par['nv0'],11),dtype='float')
for i in range(istart,iend+1):
    print 'time=',time[i],' of ',time[iend]
    gt0=dd.read_time_step_g(i)
    gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
    #Entropy
    entn_sum[:,10]=entn_sum[:,10]+get_entropy_hermite(gt0,kzind=-1,include_kz0=include_kz0)
    for k in range(10):
        kzindex=k*par['nkz0']/20
        #print 'kzindex',kzindex
        entn_sum[:,k]=entn_sum[:,k]+dd.get_entropy_hermite(gt0,kzind=kzindex)

entn_sum=entn_sum/float(ntime)

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

plot_list = np.empty(10)
for k in range(10):
    #print prefactor
    #print prefactor*entn_sum[:,k]
    kz0=kzgrid[k*par['nkz0']/20]
    #plt.loglog(herm_grid,prefactor*entn_sum[:,k],basex=10,basey=10,label=\
              #plabel+' (k_z='+str(kzgrid[k*par['nkz0']/20])+')')
    plot = prefactor*entn_sum[:,k]

    #########
    shape = np.shape(plot)
    plot_t = np.empty(shape)

    for i in range(len(kzgrid)):
        for j in range(len(hermiteNumbers)):
            plot_t[i,j] = (1j*np.sign(kzgrid[i]))**herm_grid[j] * plot[i,j]

    plot_tp = np.empty(shape)
    plot_tm = np.empty(shape)

    for i in range(len(kzgrid)):
        for j in range(len(herm_grid)+1):
            if j < (len(herm_grid)-1):
                plot_tp[i,j]= (plot_t[i,j]+plot_t[i,j+1])/2
                plot_tm[i,j]= (plot_t[i,j]-plot_t[i,j+1])/2

    c_p=np.empty(np.shape(plot_tp))
    c_m=np.empty(np.shape(plot_tm))

    for i in range(len(kzgrid)):
        for j in range(len(kzgrid)):
                c_p[i,j]= plot_tp[i,j] * np.conjugate(plot_tp[i,j])
                c_m[i,j]= plot_tm[i,j] * np.conjugate(plot_tm[i,j])

    c_ps = np.sum(c_p,axis=0)
    c_ms = np.sum(c_m,axis=0)
    ########

    plt.loglog(herm_grid,c_ps,basex=10,basey=10,label=\
              plabel+' (k_z+='+str(kzgrid[k*par['nkz0']/20])+')')

    plt.loglog(herm_grid,c_ms,basex=10,basey=10,label=\
              plabel+' (k_z-='+str(kzgrid[k*par['nkz0']/20])+')')

    #plt.xlabel('Hermite n')
    #plt.ylabel(plabel)
    #plt.legend(loc='lower left')
plt.xlabel('Hermite n')
plt.ylabel(plabel)
plt.legend(loc='lower left')
plt.show()
