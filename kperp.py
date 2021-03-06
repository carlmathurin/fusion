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
kzind = 2#(*par['nkz0']/20)
print kzgrid


end_time=time[len(time)-1]
if start_time >= end_time:
    stop

istart=np.argmin(abs(time-start_time))
iend=np.argmin(abs(time-end_time))/16
ntime=iend-istart+1
#iend => 1/16
#testing time =>1/27

entn_sum=np.zeros((par['nv0'],11,16),dtype='float')
entnp_sum=np.zeros((par['nv0'],11,16),dtype='float')
entnm_sum=np.zeros((par['nv0'],11,16),dtype='float')


for i in range(istart,iend+1):
    print 'time=',time[i],' of ',time[iend]
    gt0=dd.read_time_step_g(i)
    gt0=np.reshape(gt0,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
    dist=gt0[:,:,:,:]
    shape = np.shape(dist)
    g_t = np.empty(shape)
    for v in range(len(kx)):
        for w in range(len(ky)):
            for b in range(len(kzgrid)):
                for j in range(len(herm_grid)):
                    g_t[v,w,b,j] = (1j*np.sign(kzgrid[b]))**herm_grid[j] * dist[v,w,b,j]
    g_tp = np.empty(shape)
    g_tm = np.empty(shape)
    for v in range(len(kx)):
        for w in range(len(ky)):
            for b in range(len(kzgrid)):
                for j in range(len(herm_grid)+1):
                    if j < (len(herm_grid)-1):
                        g_tp[v,w,b,j]= (g_t[v,w,b,j]+g_t[v,w,b,j+1])/2
                        g_tm[v,w,b,j]= (g_t[v,w,b,j]-g_t[v,w,b,j+1])/2

for i in range(15):
    k_bin[i]= .15*i




for i in range(xmax):
    for j in range(xmax):
        t = 12
        k_perp= sqrt(kx[i]**2 + ky[j]**2)
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
        entnp_sum[:,kzind,t] = entn_sump[:,kzind,t] + dd.get_entropy_hermite2(g_tp,i,j,kzind)
        entnm_sum[:,kzind,t] = entn_summ[:,kzind,t] + dd.get_entropy_hermite2(g_tm,i,j,kzind)

print 'counter: ',counter

prefactor=np.empty(par['nv0'])
prefactor[:]=1.0
plabel='Entropy'

entn_sum=entn_sum/float(ntime)
entnp_sum=entnp_sum/float(ntime)
entnm_sum=entnm_sum/float(ntime)

plt.loglog(herm_grid,prefactor*entn_sum[:,kzind,13],basex=10,basey=10)

plt.xlabel('Hermite n')
plt.ylabel(plabel)
plt.title(plabel+ '(k_perp sum (kz =' + str(kzgrid[kzind]) + '))')
plt.legend(loc='lower left')
plt.show()


for j in range(15):
#print prefactor
#print prefactor*entn_sum[:,kzind,]
#kz0=kzgrid[k*par['nkz0']/20]
    plt.loglog(herm_grid,prefactor*entn_sum[:,kzind,j],basex=10,basey=10,label=\
        plabel+' (k_perp=' + str(k_bin[j])+')')
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.legend(loc='lower left')
plt.loglog(herm_grid, 10*herm_grid**(-1.5),'--',basex=10,basey=10,label='n^(-1.5)')
plt.loglog(herm_grid, 10**(-4)*herm_grid**(-1.5),'--',basex=10,basey=10,label='n^(-1.5)')
plt.legend(loc='lower left')
plt.title(plabel+'(k_perp sum [kz = ' + str(kzgrid[kzind])+'])')
plt.show()

#plus and minus terms
#~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
plt.loglog(herm_grid, (10)*herm_grid**(-1.5),'--',basex=10,basey=10,label='n^(-2)')
plt.loglog(herm_grid, (10**-2.5)*herm_grid**(-.5),'--',basex=10,basey=10,label='n^(-.5)')
plt.legend(loc='lower left')
plt.show()

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
plt.loglog(herm_grid, 10*herm_grid**(-3/2),'--',basex=10,basey=10,label='n^(-1.5)')
plt.loglog(herm_grid, 10**(-2.7)*herm_grid**(-3/2),'--',basex=10,basey=10,label='n^(-1.5)')
plt.legend(loc='lower left')
plt.show()
