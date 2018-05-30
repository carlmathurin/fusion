import numpy as np
import math
import matplotlib.pyplot as plt
import scipy as sc
from scipy.stats import linregress
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
k_bin=np.empty(16)
counter = np.zeros(16)
end_time=time[len(time)-1]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xmax = par['nkx0']
kzind = 2#(*par['nkz0']/20)
print 'kzgrid', np.shape(kzgrid), kzgrid
print 'kx' , np.shape(kx), kx
print 'ky',np.shape(ky), ky

if start_time >= end_time:
    stop

istart=np.argmin(abs(time-start_time))
iend=np.argmin(abs(time-end_time))/27
ntime=iend-istart+1
#iend => 1/16
#testing time =>1/27

entn_sum=np.zeros((par['nv0'],96,16),dtype='float')
entnp_sum=np.zeros((par['nv0'],96,16),dtype='float')
entnm_sum=np.zeros((par['nv0'],96,16),dtype='float')

print "herm_grid",np.shape(herm_grid)

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

for i in range(16):
    k_bin[i]= .15*i

for i in range(xmax):
    for j in range(xmax):
        t = 19
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
# add new loop for k_z here and sum over every k_z for each x-y combination.
        for b in range(len(kzgrid)):
             entn_sum[:,b,t] = entn_sum[:,b,t] + dd.get_entropy_hermite2(gt0,i,j,b)
             entnp_sum[:,b,t] = entnp_sum[:,b,t] + dd.get_entropy_hermite2(g_tp,i,j,b)
             entnm_sum[:,b,t] = entnm_sum[:,b,t] + dd.get_entropy_hermite2(g_tm,i,j,b)
"""
        entn_sum[:,kzind,t] = entn_sum[:,kzind,t] + dd.get_entropy_hermite2(gt0,i,j,kzind)
        entnp_sum[:,kzind,t] = entnp_sum[:,kzind,t] + dd.get_entropy_hermite2(g_tp,i,j,kzind)
        entnm_sum[:,kzind,t] = entnm_sum[:,kzind,t] + dd.get_entropy_hermite2(g_tm,i,j,kzind)
"""
print 'counter: ',counter,'kperp',np.shape(k_bin), k_bin

prefactor=np.empty(par['nv0'])
prefactor[:]=1.0
plabel='Entropy'
entn_sum=entn_sum/float(ntime)
entnp_sum=entnp_sum/float(ntime)
entnm_sum=entnm_sum/float(ntime)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TESTING BLOCK~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
#polyfit
hermy = []
enm = []
m = np.zeros(20)
b = np.zeros(20)
counter = 0

for j in range(20):
    for i in range(60):
        if herm_grid[i] > 0:
            kzs = j*par['nkz0']/20
            if entn_sum[i,kzs,6] > 0:
                lomein = np.log(entn_sum[i,kzs,6])
                friedrice = np.log(herm_grid[i])
                hermy.append(friedrice)
                enm.append(lomein)

    m[j],b[j] = polyfit(hermy,enm,1)

print np.size(m) ,'slopes:', m
m1 = 0
for j in range(20):
    m1 = m1 + m[j]

m1 = m1/ (20)
print 'm1 =', m1

fig1 = plt.figure()
ax1 = plt.subplot(111)
box = ax1.get_position()
ax1.set_position([box.x0, box.y0 , box.width* .8, box.height])

for j in range(20):
#kz0=kzgrid[k*par['nkz0']/20]
    k=j*par['nkz0']/20
    plt.loglog(herm_grid,prefactor*entn_sum[:,k,6],basex=10,basey=10,label=\
        '(k_z='+str(kzgrid[j*par['nkz0']/20]))
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.legend(loc='lower left')
plt.loglog(herm_grid, 10*herm_grid**(-1.5),'--',basex=10,basey=10,label='n^(-1.5)')
plt.loglog(herm_grid, 10**(-4)*herm_grid**(-1.5),'--',basex=10,basey=10,label='n^(-1.5)')
plt.loglog(herm_grid, (10**-1)*herm_grid**(m1),'--',basex=10,basey=10,label=('n^(%.4f)'% m1))
plt.legend(loc='center left', bbox_to_anchor=(1 ,.5) )
plt.title(plabel+'(kz sum [k_perp = ' + str(k_bin[6])+'])')
plt.show()
"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TESTING BLOCK~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TESTING BLOCK~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
#polyfit

hermy = []
enm = []
m = np.zeros(15)
b = np.zeros(15)
counter = 0

for j in range(15):
    for i in range(60):
        if herm_grid[i] > 0:
            if entnp_sum[i,67,j] > 0:
                lomein = np.log(entnp_sum[i,67,j])
                friedrice = np.log(herm_grid[i])
                hermy.append(friedrice)
                enm.append(lomein)
            if entnp_sum[i,67,j] == 0:
                counter = counter + 1

    m[j],b[j] = polyfit(hermy,enm,1)

m1 = 0
b1 = 0
for j in range(15):
    m1 = m1 + m[j]
    b1 = b1 + b[j]

m1 = m1/ (15 - counter)
b1 = b1/(15 - counter)
m2 = m[4]
b2 = b[4]

print 'm1 =', m1, 'm2 =', m2

for j in range(15):
    #kz0=kzgrid[k*par['nkz0']/20]
    plt.loglog(herm_grid,prefactor*entnp_sum[:,67,j],basex=10,basey=10)#,label=\
              #plabel+' (k_z='+str(kzgrid[k*par['nkz0']/20])+')')
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.title('Entropy+')
    plt.legend(loc='lower left')
#temp=prefactor*entn_sum[20,10]
#temp=temp/(herm_grid**(-1))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-1),'--',basex=10,basey=10,label=str(-1))
#temp=prefactor*entn_sum[20,10]
#temp=temp/(herm_grid**(-1.5))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-3.5),'--',basex=10,basey=10,label=str(-3.5))
plt.loglog(herm_grid, (10**-3)*herm_grid**(m1),'--',basex=10,basey=10,label=('n^(%.4f)'% m1))
plt.loglog(herm_grid, (10**-3)*herm_grid**(m2),'--',basex=10,basey=10,label=('n^(%.4f)'% m2))
plt.legend(loc='lower left')
plt.show()
## slope for E+ ~1.75
"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TESTING BLOCK~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
#poly fit~~~~~~~~~~~~~~~~~~~~~~~
print np.shape(entnm_sum)
print "entnm_sum:" ,entnm_sum[:,kzind,10]
print "herm", herm_grid
hermy = []
enm = []

for i in range(60):
     if herm_grid[i] > 0:
          if entnm_sum[i,kzind,10] > 0:
               lomein = np.log(entnm_sum[i,kzind,10])
               friedrice = np.log(herm_grid[i])
               hermy.append(friedrice)
               enm.append(lomein)
print "hermy" , hermy
print "enm" , enm

m,b = polyfit(hermy,enm,1)

print 'm = ',m

for j in range(15):
    #print prefactor
    #print prefactor*entn_sum[:,k]
    #kz0=kzgrid[k*par['nkz0']/20]
    #plots.append(prefactor*entn_sum[:,k])
    plt.loglog(herm_grid,prefactor*entnm_sum[:,kzind,j],basex=10,basey=10,)#label=\
              #plabel+' (k_z='+str(kzgrid[k*par['nkz0']/20])+')')
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.title('Entropy-')
    plt.legend(loc='lower left')
#temp=prefactor*entn_sum[20,10]
#temp=temp/(herm_grid**(-1))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-1),'--',basex=10,basey=10,label=str(-1))
#temp=prefactor*entn_sum[20,10]
#temp=temp/(herm_grid**(-1.5))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-3.5),'--',basex=10,basey=10,label=str(-3.5))
plt.loglog(herm_grid, (10**-1.25)*herm_grid**(-1.75),'--',basex=10,basey=10,label='n^(-1.75)')
plt.loglog(herm_grid, 10**(-2.25)*herm_grid**(m),'--',basex=10,basey=10,label= ('n^(%.4f)'% m))
plt.legend(loc='lower left')
plt.show()

print np.shape(entnm_sum)
print "entnm_sum:" ,entnm_sum[:,kzind,10]
print "herm", herm_grid
hermy = []
enm = []
"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TESTING BLOCK~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
entp_sm = 0
entm_sm = 0
dummyx = [0,1]

nmax=len(herm_grid)

for i in range(nmax):
     entp_sm = entp_sm + entnp_sum[i,kzind,13]
     entm_sm = entm_sm + entnm_sum[i,kzind,13]

ent_tot = entp_sm + entm_sm

fig, ax = plt.subplots()

pper = (entp_sm/ent_tot) * 100
mper = (entm_sm/ent_tot) * 100

labels = 'Ent +', 'Ent -'
fracs = [pper, mper]


# Make square figures and axes
the_grid = GridSpec(1, 1)
plt.subplot(the_grid[0, 0], aspect=1)

textstr = 'Ent + = %.5f, Ent - = %.5f' %(entp_sm,entm_sm)

plt.pie(fracs, labels=labels, autopct='%1.1f%%', shadow=True)
plt.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top')
plt.title('Ent Total='+str(ent_tot)+' k_perp ='+str(k_bin[13]))
plt.show()

"""
entp_sm = 0
entm_sm = 0
dummyx = [0,1]

nmax=len(herm_grid)

for i in range(nmax):
     entp_sm = entp_sm + entnp_sum[i,kzind,13]
     entm_sm = entm_sm + entnm_sum[i,kzind,13]

ent_tot = entp_sm + entm_sm

fig, ax = plt.subplots()

pper = (entp_sm/ent_tot) * 100
mper = (entm_sm/ent_tot) * 100

ax.fill_between(dummyx,0,entp_sm, facecolor='blue',label=('ent_p=%.4f' % entp_sm, ', %.2f %% of total ' % pper))
ax.fill_between(dummyx,entp_sm,ent_tot, facecolor='green',label=('ent_m= %.4f'  % entm_sm, ', %.2f %% of total '%mper))
ax.plot(dummyx,[0,0],label=('ent_tot = %.4f' % ent_tot, ';kz = %.4f'% kzind, ';k_perp = %.4f' % k_bin[13]))
ax.set_ylabel('Ent (units)')
box = ax.get_position()
#~~~~~~~~~~~~~~legend on the left
#ax.set_position([box.x0, box.y0, box.width * .8, box.height])
#plt.legend(loc='center left',bbox_to_anchor=(1,.5))
ax.set_position([box.x0, box.y0 + box.height*.1, box.width, box.height*.9])
ax.legend(loc='upper center', bbox_to_anchor=(.5,-.05))
plt.title('Total + and - Entropy')
plt.show()
"""
