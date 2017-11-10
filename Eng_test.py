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
iend=np.argmin(abs(time-end_time))/27
ntime=iend-istart+1
#iend => 1/16
#testing time =>1/27

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
    #for b in range(len(kzgrid)):
        #for j in range(len(herm_grid)):
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
              plabel+'+ (k_z='+str(kzgrid[k*par['nkz0']/20])+')')
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.title('Energy comparison')
    #print prefactor
    #print prefactor*entn_sum[:,k]
    kz0=kzgrid[k*par['nkz0']/20]
    #plots.append(prefactor*entn_sum[:,k])
    plt.loglog(herm_grid,prefactor*entnm_sum[:,k],basex=10,basey=10,label=\
              plabel+'- (k_z='+str(kzgrid[k*par['nkz0']/20])+')')
    plt.legend(loc='lower left')

    npr =[1+i for i in range(1000)]
    kzpr = [.2 + i*.2 for i in range(1000)]
    nuupr =.002

    sizepr = np.shape(npr)
    CapG = np.empty([1000,1000])
    Cap2G = np.empty([1000,1000])

    for j in range(1000):
        for i in range(1000):
            nupr = .01 #nuupr/kzpr[j]
            Dpr = ((npr[i]*nupr**2)/4 + 1)**(.5)
            npl = Dpr + npr[i]**(.5)*nupr/2
            CapG[j,i] = (np.exp((npr[i]**(.5))*\
                Dpr/(np.absolute(nupr))) \
                * (npl)**(-np.sign(kzpr[j])*(npr[i]+2*nupr**(-2)-.5))\
                /(npr[i]**(.25)*Dpr**(.5)))\
                #*(-1j*np.sign(kz[j]))**n[i]
            Cap2G[j,i] = CapG[j,i]* np.conjugate(CapG[j,i])
        plt.loglog(npr,Cap2G[j,:],label= 'G' + 'kz ('+ str(kzpr[j])+')')
        plt.title('G ')
        plt.xlabel('Hermite n')
        plt.ylabel('CapG')
    plt.show()
