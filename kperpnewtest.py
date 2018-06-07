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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TESTING BLOCK (fixed k_perp)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
#First k_perp
plt.figure(1)
#~~~~~~~~~~~~~~~~~~~~~~ 1) k_perp +
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
            if entn_sum[i,kzs,11] > 0:
                lomein = np.log(entn_sum[i,kzs,11])
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

ax1 = plt.subplot(121)
box = ax1.get_position()
ax1.set_position([box.x0, box.y0 , box.width* .8, box.height])

for j in range(20):
#kz0=kzgrid[k*par['nkz0']/20]
    k=j*par['nkz0']/20
    plt.loglog(herm_grid,prefactor*entnp_sum[:,k,11],basex=10,basey=10,label=\
        '(k_z='+str(kzgrid[j*par['nkz0']/20]))
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.legend(loc='lower left')
plt.loglog(herm_grid, 10*herm_grid**(-1.5),'--',basex=10,basey=10,label='n^(-1.5)')
plt.loglog(herm_grid, 10**(-4)*herm_grid**(-1.5),'--',basex=10,basey=10,label='n^(-1.5)')
plt.loglog(herm_grid, (10**-1)*herm_grid**(m1),'--',basex=10,basey=10,label=('n^(%.4f)'% m1))
plt.legend(loc='center left', bbox_to_anchor=(1 ,.5) )
plt.title(plabel+'(kz(+) sum [k_perp = ' + str(k_bin[11])+'])')
#plt.show()
#~~~~~~~~~~~~~~~~~~~~~~ 1) k_perp
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
            if entn_sum[i,kzs,11] > 0:
                lomein = np.log(entn_sum[i,kzs,11])
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

ax1 = plt.subplot(122)
box = ax1.get_position()
ax1.set_position([box.x0, box.y0 , box.width* .8, box.height])

for j in range(20):
#kz0=kzgrid[k*par['nkz0']/20]
    k=j*par['nkz0']/20
    plt.loglog(herm_grid,prefactor*entnm_sum[:,k,11],basex=10,basey=10,label=\
        '(k_z='+str(kzgrid[j*par['nkz0']/20]))
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.legend(loc='lower left')
plt.loglog(herm_grid, 10*herm_grid**(-1.5),'--',basex=10,basey=10,label='n^(-1.5)')
plt.loglog(herm_grid, 10**(-4)*herm_grid**(-1.5),'--',basex=10,basey=10,label='n^(-1.5)')
plt.loglog(herm_grid, (10**-1)*herm_grid**(m1),'--',basex=10,basey=10,label=('n^(%.4f)'% m1))
plt.legend(loc='center left', bbox_to_anchor=(1 ,.5) )
plt.title(plabel+'(kz(-) sum [k_perp = ' + str(k_bin[11])+'])')
plt.show()
"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TESTING BLOCK (fixed K_perp [+])~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### plot 1
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

plt.figure(1)
ax1 = plt.subplot(121)

for j in range(15):
    #kz0=kzgrid[k*par['nkz0']/20]
    plt.loglog(herm_grid,prefactor*entnp_sum[:,67,j],basex=10,basey=10,label= '(k_perp ='+str(k_bin[j])+')')
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.title(plabel+'+ (k_perp(+) sum [kz = ' + str(kzgrid[67])+'])')
    plt.legend(loc='lower left')
#temp=prefactor*entn_sum[20,10]
#temp=temp/(herm_grid**(-1))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-1),'--',basex=10,basey=10,label=str(-1))
#temp=prefactor*entn_sum[20,10]
#temp=temp/(herm_grid**(-1.5))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-3.5),'--',basex=10,basey=10,label=str(-3.5))
plt.loglog(herm_grid, (10**-3.5)*herm_grid**(m1),'--',basex=10,basey=10,label=('n^(%.4f)'% m1))
plt.loglog(herm_grid, (10**-3.5)*herm_grid**(m2),'--',basex=10,basey=10,label=('n^(%.4f)'% m2))
plt.legend(loc='lower left')

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 , box.width* .8, box.height])
plt.legend(loc='center left', bbox_to_anchor=(1 ,.5) )
#plt.show()

########## plot 2
#polyfit

hermy = []
enm = []
m = np.zeros(15)
b = np.zeros(15)
counter = 0

for j in range(15):
    for i in range(60):
        if herm_grid[i] > 0:
            if entnp_sum[i,40,j] > 0:
                lomein = np.log(entnp_sum[i,40,j])
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

plt.figure(2)
ax1=plt.subplot(121)

for j in range(15):
    #kz0=kzgrid[k*par['nkz0']/20]
    plt.loglog(herm_grid,prefactor*entnp_sum[:,40,j],basex=10,basey=10,label='(k_perp='+str(k_bin[j])+')')
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.title(plabel+'+ (k_perp(+) sum [kz = ' + str(kzgrid[40])+'])')
    plt.legend(loc='lower left')
#temp=prefactor*entn_sum[20,10]
#temp=temp/(herm_grid**(-1))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-1),'--',basex=10,basey=10,label=str(-1))
#temp=prefactor*entn_sum[20,10]
#temp=temp/(herm_grid**(-1.5))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-3.5),'--',basex=10,basey=10,label=str(-3.5))
plt.loglog(herm_grid, (10**-4.5)*herm_grid**(m1),'--',basex=10,basey=10,label=('n^(%.4f)'% m1))
plt.loglog(herm_grid, (10**-4.5)*herm_grid**(m2),'--',basex=10,basey=10,label=('n^(%.4f)'% m2))

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 , box.width* .8, box.height])
plt.legend(loc='center left', bbox_to_anchor=(1 ,.5) )
## slope for E+ ~1.75

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TESTING BLOCK (fixed K_z [-])~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###### plot 1
#poly fit~~~~~~~~~~~~~~~~~~~~~~~
#print np.shape(entnm_sum)
#print "entnm_sum:" ,entnm_sum[:,67,10]
#print "herm", herm_grid
hermy = []
enm = []

for i in range(60):
     if herm_grid[i] > 0:
          if entnm_sum[i,67,10] > 0:
               lomein = np.log(entnm_sum[i,67,10])
               friedrice = np.log(herm_grid[i])
               hermy.append(friedrice)
               enm.append(lomein)
print "hermy" , hermy
print "enm" , enm

m,b = polyfit(hermy,enm,1)

print 'm = ',m
plt.figure(1)
ax1= plt.subplot(122)

for j in range(15):
    #print prefactor
    #print prefactor*entn_sum[:,k]
    #kz0=kzgrid[k*par['nkz0']/20]
    #plots.append(prefactor*entn_sum[:,k])
    plt.loglog(herm_grid,prefactor*entnm_sum[:,67,j],basex=10,basey=10,label= '(k_perp='+str(k_bin[j])+')')
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.title(plabel+'- (k_perp(-) sum [kz = ' + str(kzgrid[67])+'])')
    plt.legend(loc='lower left')
#temp=prefactor*entn_sum[20,10]
#temp=temp/(herm_grid**(-1))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-1),'--',basex=10,basey=10,label=str(-1))
#temp=prefactor*entn_sum[20,10]
#temp=temp/(herm_grid**(-1.5))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-3.5),'--',basex=10,basey=10,label=str(-3.5))
plt.loglog(herm_grid, (10**-3.5)*herm_grid**(-1.75),'--',basex=10,basey=10,label='n^(-1.75)')
plt.loglog(herm_grid, 10**(-3.5)*herm_grid**(m),'--',basex=10,basey=10,label= ('n^(%.4f)'% m))
plt.legend(loc='lower left')

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 , box.width* .8, box.height])
plt.legend(loc='center left', bbox_to_anchor=(1 ,.5) )
##plt.show()

###### plot 2
#poly fit~~~~~~~~~~~~~~~~~~~~~~~
#print np.shape(entnm_sum)
#print "entnm_sum:" ,entnm_sum[:,kzind,10]
hermy = []
enm = []

for i in range(60):
     if herm_grid[i] > 0:
          if entnm_sum[i,kzind,10] > 0:
               lomein = np.log(entnm_sum[i,40,10])
               friedrice = np.log(herm_grid[i])
               hermy.append(friedrice)
               enm.append(lomein)
print "hermy" , hermy
print "enm" , enm

m,b = polyfit(hermy,enm,1)

print 'm = ',m
plt.figure(2)
ax1 = plt.subplot(122)

for j in range(15):
    #print prefactor
    #print prefactor*entn_sum[:,k]
    #kz0=kzgrid[k*par['nkz0']/20]
    #plots.append(prefactor*entn_sum[:,k])
    plt.loglog(herm_grid,prefactor*entnm_sum[:,40,j],basex=10,basey=10,label= '(k_perp='+str(k_bin[j])+')')
    plt.xlabel('Hermite n')
    plt.ylabel(plabel)
    plt.title(plabel+'- (k_perp(-) sum [kz = ' + str(kzgrid[40])+'])')
#temp=prefactor*entn_sum[20,10]
#temp=temp/(herm_grid**(-1))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-1),'--',basex=10,basey=10,label=str(-1))
#temp=prefactor*entn_sum[20,10]
#temp=temp/(herm_grid**(-1.5))[20]
#plt.loglog(herm_grid,2.0*temp*herm_grid**(-3.5),'--',basex=10,basey=10,label=str(-3.5))
plt.loglog(herm_grid, (10**-4.5)*herm_grid**(-1.75),'--',basex=10,basey=10,label='n^(-1.75)')
plt.loglog(herm_grid, 10**(-4.5)*herm_grid**(m),'--',basex=10,basey=10,label= ('n^(%.4f)'% m))
plt.legend(loc='lower left')

box = ax1.get_position()
ax1.set_position([box.x0, box.y0 , box.width* .8, box.height])
plt.legend(loc='center left', bbox_to_anchor=(1 ,.5) )

plt.show()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~TESTING BLOCK (pie graphs)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
the_grid = GridSpec(3, 1)

entp_sm = 0
entm_sm = 0
dummyx = [0,1]
nmax=len(herm_grid)

for i in range(nmax):
     entp_sm = entp_sm + entnp_sum[i,11,13]
     entm_sm = entm_sm + entnm_sum[i,11,13]

ent_tot = entp_sm + entm_sm
pper = (entp_sm/ent_tot) * 100
mper = (entm_sm/ent_tot) * 100
labels = 'Ent +', 'Ent -'
fracs = [pper, mper]
# Make square figures and axes
plt.subplot(the_grid[0, 0], aspect=1)

textstr = 'Ent + = %.5f\n Ent - = %.5f\n Ent tot =%.5f' %(entp_sm,entm_sm,ent_tot)

plt.pie(fracs, labels=labels, autopct='%1.1f%%', shadow=True)
plt.text(0.9, 0.95, textstr, fontsize=14,verticalalignment='top')
plt.title('Ent Total, k_perp ='+str(k_bin[13])+', kz ='+str(kzgrid[11]))
######plt.show()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~pie 2~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
entp_sm = 0
entm_sm = 0
dummyx = [0,1]
nmax=len(herm_grid)

for i in range(nmax):
     entp_sm = entp_sm + entnp_sum[i,32,13]
     entm_sm = entm_sm + entnm_sum[i,32,13]

ent_tot = entp_sm + entm_sm
pper = (entp_sm/ent_tot) * 100
mper = (entm_sm/ent_tot) * 100
labels = 'Ent +', 'Ent -'
fracs = [pper, mper]
# Make square figures and axes
plt.subplot(the_grid[1, 0], aspect=1)

textstr = 'Ent + = %.5f\n Ent - = %.5f\n Ent tot =%.5f' %(entp_sm,entm_sm,ent_tot)

plt.pie(fracs, labels=labels, autopct='%1.1f%%', shadow=True)
plt.text(0.9, 0.95, textstr, fontsize=14,verticalalignment='top')
plt.title('Ent Total, k_perp ='+str(k_bin[13])+', kz ='+str(kzgrid[32]))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~pie 3~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
entp_sm = 0
entm_sm = 0
dummyx = [0,1]
nmax=len(herm_grid)

for i in range(nmax):
     entp_sm = entp_sm + entnp_sum[i,80,13]
     entm_sm = entm_sm + entnm_sum[i,80,13]

ent_tot = entp_sm + entm_sm
pper = (entp_sm/ent_tot) * 100
mper = (entm_sm/ent_tot) * 100
labels = 'Ent +', 'Ent -'
fracs = [pper, mper]
# Make square figures and axes
plt.subplot(the_grid[2, 0], aspect=1)

textstr = 'Ent + = %.5f\n Ent - = %.5f\n Ent tot =%.5f' %(entp_sm,entm_sm,ent_tot)

plt.pie(fracs, labels=labels, autopct='%1.1f%%', shadow=True)
plt.text(0.9, 0.95, textstr, fontsize=14,verticalalignment='top')
plt.title('Ent Total, k_perp ='+str(k_bin[13])+', kz ='+str(kzgrid[80]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~pie 4~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
entp_sm = 0
entm_sm = 0
dummyx = [0,1]
nmax=len(herm_grid)

for i in range(nmax):
     entp_sm = entp_sm + entnp_sum[i,41,11]
     entm_sm = entm_sm + entnm_sum[i,41,11]

ent_tot = entp_sm + entm_sm
pper = (entp_sm/ent_tot) * 100
mper = (entm_sm/ent_tot) * 100
labels = 'Ent +', 'Ent -'
fracs = [pper, mper]
# Make square figures and axes
plt.subplot(the_grid[0, 0], aspect=1)

textstr = 'Ent + = %.5f\n Ent - = %.5f\n Ent tot =%.5f' %(entp_sm,entm_sm,ent_tot)

plt.pie(fracs, labels=labels, autopct='%1.1f%%', shadow=True)
plt.text(0.9, 0.95, textstr, fontsize=14,verticalalignment='top')
plt.title('Ent Total, k_perp ='+str(k_bin[11])+', kz ='+str(kzgrid[41]))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~pie 5~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
entp_sm = 0
entm_sm = 0
dummyx = [0,1]
nmax=len(herm_grid)

for i in range(nmax):
     entp_sm = entp_sm + entnp_sum[i,41,7]
     entm_sm = entm_sm + entnm_sum[i,41,7]

ent_tot = entp_sm + entm_sm
pper = (entp_sm/ent_tot) * 100
mper = (entm_sm/ent_tot) * 100
labels = 'Ent +', 'Ent -'
fracs = [pper, mper]
# Make square figures and axes
plt.subplot(the_grid[1, 0], aspect=1)

textstr = 'Ent + = %.5f\n Ent - = %.5f\n Ent tot =%.5f' %(entp_sm,entm_sm,ent_tot)

plt.pie(fracs, labels=labels, autopct='%1.1f%%', shadow=True)
plt.text(0.9, 0.95, textstr, fontsize=14,verticalalignment='top')
plt.title('Ent Total, k_perp ='+str(k_bin[7])+', kz ='+str(kzgrid[41]))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~pie 6~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
entp_sm = 0
entm_sm = 0
dummyx = [0,1]
nmax=len(herm_grid)

for i in range(nmax):
     entp_sm = entp_sm + entnp_sum[i,41,3]
     entm_sm = entm_sm + entnm_sum[i,41,3]

ent_tot = entp_sm + entm_sm
pper = (entp_sm/ent_tot) * 100
mper = (entm_sm/ent_tot) * 100
labels = 'Ent +', 'Ent -'
fracs = [pper, mper]
# Make square figures and axes
plt.subplot(the_grid[2, 0], aspect=1)

textstr = 'Ent + = %.5f\n Ent - = %.5f\n Ent tot =%.5f' %(entp_sm,entm_sm,ent_tot)

plt.pie(fracs, labels=labels, autopct='%1.1f%%', shadow=True)
plt.text(0.9, 0.95, textstr, fontsize=14,verticalalignment='top')
plt.title('Ent Total, k_perp ='+str(k_bin[3])+', kz ='+str(kzgrid[41]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plt.show()
"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Testing Block (Countour) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
kzbin = np.zeros(20)
for j in range(20):
    kzs = j*par['nkz0']/20
    kzbin[j] = kzgrid[kzs]

shape2 = (20,15)
entp_sm = np.zeros(shape2)
entm_sm = np.zeros(shape2)
ent_tot = np.zeros(shape2)
nmax=len(herm_grid)

for j in range(20):
    for i in range(15):
        for k in range(nmax):
            kzs = j*par['nkz0']/20
            kzbin[j] = kzgrid[kzs]

            entp_sm[j,i] = entp_sm[j,i] + entnp_sum[k,kzs,i]
            entm_sm[j,i] = entm_sm[j,i] + entnm_sum[k,kzs,i]

            ent_tot[j,i] = entp_sm[j,i] + entm_sm[j,i]

plt.figure(2)

X1,Y1 = np.meshgrid(k_bin,kzbin)



plt.contour(X1[:,1:],Y1[:,1:],ent_tot)
plt.xlabel('K_z')
plt.ylabel('K_perp')
plt.show()
"""
