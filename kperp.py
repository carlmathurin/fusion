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

k_bin=np.empty(15)
counter = np.zeros(15)

xmax = par['nkx0']
kzind = 2*par['nkz0']/20


end_time=time[len(time)-1]
if start_time >= end_time:
    stop

istart=np.argmin(abs(time-start_time))
iend=np.argmin(abs(time-end_time))/16
ntime=iend-istart+1
#iend => 1/16
#testing time =>1/27

entn_sum=np.zeros((par['nv0'],11,16),dtype='float')
entnp_sum=np.zeros((par['nv0'],11,11),dtype='float')
entnm_sum=np.zeros((par['nv0'],11,11),dtype='float')


for i in range(istart,iend+1):
    print 'time=',time[i],' of ',time[iend]
    gt0=dd.read_time_step_g(i)
    gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')


for i in range(15):
    k_bin[i]= .15*i




for i in range(xmax):
    for j in range(xmax):
        t = 12
        k_perp = sqrt(kx[i]**2 + ky[j]**2)
        if  k_perp < k_bin[1]:
            t = 0
            counter[t] = counter[t] + 1
        elif k_bin[1] <= k_perp < k_bin[2]:
            t= 1
            counter[t] = counter[t] + 1
        elif k_bin[2] <= k_perp < k_bin[3]:
            t= 2
            counter[t] = counter[t] + 1
        elif k_bin[3] <= k_perp < k_bin[4]:
            t= 3
            counter[t] = counter[t] + 1
        elif k_bin[4] <= k_perp < k_bin[5]:
            t= 4
            counter[t] = counter[t] + 1
        elif k_bin[5] <= k_perp < k_bin[6]:
            t= 5
            counter[t] = counter[t] + 1
        elif k_bin[6] <= k_perp < k_bin[7]:
            t= 6
            counter[t] = counter[t] + 1
        elif k_bin[7] <= k_perp < k_bin[8]:
            t= 7
            counter[t] = counter[t] + 1
        elif k_bin[8] <= k_perp < k_bin[9]:
            t= 8
            counter[t] = counter[t] + 1
        elif k_bin[9] <= k_perp  < k_bin[10]:
            t= 9
            counter[t] = counter[t] + 1
        elif k_bin[10] <= k_perp  < k_bin[11]:
            t= 10
            counter[t] = counter[t] + 1
        elif k_bin[11] <= k_perp  < k_bin[12]:
            t= 11
            counter[t] = counter[t] + 1
        elif k_bin[12] <= k_perp  < k_bin[13]:
            t= 12
            counter[t] = counter[t] + 1
        elif k_bin[13] <= k_perp  < k_bin[14]:
            t= 13
            counter[t] = counter[t] + 1
        elif k_bin[14] <= k_perp  :
            t= 14
            counter[t] = counter[t] + 1


        entn_sum[:,kzind,t] = entn_sum[:,kzind,t] + dd.get_entropy_hermite2(gt0,i,j,kzind)

print 'counter: ',counter
print np.shape(entn_sum)
print 'entn sum:' , entn_sum

prefactor=np.empty(par['nv0'])
prefactor[:]=1.0
plabel='Entropy'

entn_sum=entn_sum/float(ntime)

plt.loglog(herm_grid,prefactor*entn_sum[:,kzind,13],basex=10,basey=10)

plt.xlabel('Hermite n')
plt.ylabel(plabel)
plt.title(plabel+'(k_perp sum)')
plt.legend(loc='lower left')
plt.show()


for j in range(15):
#print prefactor
#print prefactor*entn_sum[:,kzind,]
#kz0=kzgrid[k*par['nkz0']/20]
    plt.loglog(herm_grid,prefactor*entn_sum[:,kzind,j],basex=10,basey=10,label=\
        plabel+' (k_perp=' + str(k_perp[j])+')')
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.legend(loc='lower left')
plt.show()
