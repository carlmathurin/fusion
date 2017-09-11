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
import kz3_herm as d2

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



ap = ''
bp = ''
cp = ''
ep = ''
am = ''
bm = ''
cm = ''
em = ''
a, a1, ax, ay  = d2.diags_cp_cm(6,16)
b, b1, bx, by = d2.diags_cp_cm(10,50)
c, c1, cx, cy  = d2.diags_cp_cm(30,63)
f, f1, fx, fy = d2.diags_cp_cm(25,38)
e, e1, n, ex, ey, ez  = d2.diags_kz(16,5,10)
g, g1, n, gx, gy, gz = d2.diags_kz(5,45,30)


ss = ""
ap = 'C+ (',str(ax),',',str(ay), ')'
ap = ss.join(ap)

bp = 'C+ (',str(bx),',',str(by), ')'
bp = ss.join(bp)

cp = 'C+ (',str(cx),',',str(cy), ')'
cp = ss.join(cp)

fp = 'C+ (',str(fx),',',str(fy), ')'
fp = ss.join(fp)

ep = 'C+ (',str(ex),',',str(ey),',',str(ez), ')'
ep = ss.join(ep)

gp = 'C+ (',str(gx),',',str(gy),',',str(gz), ')'
gp = ss.join(gp)

am = 'C- (',str(ax),',',str(ay), ')'
am = ss.join(am)

bm = 'C- (',str(bx),',',str(by), ')'
bm = ss.join(bm)

cm = 'C- (',str(cx),',',str(cy), ')'
cm = ss.join(cm)

fm = 'C+ (',str(fx),',',str(fy), ')'
fm = ss.join(fm)

em = 'C- (',str(ex),',',str(ey),',',str(ez), ')'
em = ss.join(em)

gm = 'C- (',str(gx),',',str(gy),',',str(gz), ')'
gm = ss.join(gm)


cp1, = plt.loglog(n,a,)
cp2, = plt.loglog(n,b,)
cp3, = plt.loglog(n,c,)
cp4, = plt.loglog(n,e,)
cp5, = plt.loglog(n,a1,)
cp6, = plt.loglog(n,b1,)
cp7, = plt.loglog(n,c1,)
cp8, = plt.loglog(n,e1,)
cp11, = plt.loglog(n,g)
cp12, = plt.loglog(n,g1)
plt.autoscale(enable=True,axis='both',tight=False)
plt.legend([cp1,cp2,cp3,cp4,cp5,cp6,cp7,cp8,cp11,cp12], [ap,bp,cp,ep,am,bm,cm,em,gp,gm], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.xlabel('Hermite Number')
plt.show()
