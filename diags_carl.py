import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from config import *
from dna_diags import *
from ev_diags import *
from spectra import *
from nlt_diags import *
#from landau_tests import *
import os


print "Hello world"
read_parameters()
diagdir = '/work/01658/drhatch/prod_nu05_omt7.5_zf0b'
par['diagdir']="\'"+diagdir+"\'"

time = get_time_from_gout()
time
kx,ky,kz,n = get_grids()
print "kx",kx
print "ky",ky
print "kz",kz
print "n",n
g_in = read_time_step_g(len(time)-1)
g_in=np.reshape(g_in,(par['nkx0'],par['nky0'],par['nkz0'],par['nv0']),order='F')
print np.shape(g_in)
gkx5ky5 = g_in[5,5,:,:]
print np.shape(gkx5ky5)

ra=(np.shape(gkx5ky5)[0],np.shape(gkx5ky5)[1])
print 'ra: ', ra
gkx5ky5t = np.empty(ra)
for i in range(len(kz)):
    for j in range(len(n)):
	if kz[i] == 0:
           	sgnkz = 0
        else:
            	sgnkz = kz[i]/abs(kz[i])
        gkx5ky5t[i,j] = (1j*sgnkz)**n[j]*gkx5ky5[i,j]
gkx5ky5tp = np.empty(ra)
gkx5ky5tm = np.empty(ra)
for i in range(len(kz)):
    for j in range(len(n)+1):
	if j < (len(n)-1):
		gkx5ky5tp[i,j]= (gkx5ky5t[i,j]+gkx5ky5t[i,j+1])/2
		gkx5ky5tm[i,j]= (gkx5ky5t[i,j]-gkx5ky5t[i,j+1])/2
print np.shape(gkx5ky5tm)
print np.shape(gkx5ky5tp)
c_p=np.empty(np.shape(gkx5ky5tp))
c_m=np.empty(np.shape(gkx5ky5tm))
for i in range(len(kz)):
	for j in range(len(kz)):
        	c_p[i,j]= gkx5ky5tp[i,j] * np.conjugate(gkx5ky5tp[i,j])
        	c_m[i,j]= gkx5ky5tm[i,j] * np.conjugate(gkx5ky5tm[i,j])
print np.shape(c_p)
print np.shape(c_m)
print 'time#=', np.shape(time)

c_ps = np.sum(c_p,axis=0)
plt.loglog(n,c_ps)
c_ms = np.sum(c_m,axis=0)
plt.loglog(n,c_ms)
plt.show()

#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*##*#*#*#*#

nuu =.002





#D = (n[j]*nu**2/4 + 1)**(.5)
#npl = D + n[j]**(.5)*nu/2

size = np.shape(n)
CapG = np.empty([64])
Cap2G = np.empty([64])

#for j in range(len(kz)):
    #print 'Kz =', kz[j]
j=10
for i in range(len(kz)):
        #print 'N =', n[i]
    nu = nuu/kz[j] #.01
    D = ((n[i]*nu**2)/4 + 1)**(.5)
    print 'D = ', D
    npl = D + n[i]**(.5)*nu/2
    result = (np.exp((n[i]**(.5))*\
        D/(np.absolute(nu))) \
        * (npl)**((-np.sign(kz[j]))*(n[i]+2*nu**(-2)-.5))
         #(npl)**((1)*(n[i]+2*nu**(-2)-.5))\
        #/(n[i]**(.25)*D**(.5)))\
        *(1j*np.sign(kz[j]))**n[i]
    print 'result: ', result
    Cap2G[i] = result * np.conjugate(result) #CapG[j,i]* np.conjugate(CapG[j,i])
    #print (-1j*np.sign(kz[j]))**n[i]
    #print (-1*np.sign(kz[j]))**n[i]
plt.loglog(n,np.real(Cap2G[:]))
plt.title('G ')
plt.xlabel('Hermite n')
plt.ylabel('CapG')

plt.show()

#*#*#*#*#*#*#*#*##*#*#*#*#*#*#*#*#*##*#*#*#*#*#

#change to take specific inputs of kx ky
#--Construct Eq 4.27 from Kanekar paper
#Construct C+ and C- by calculating g+*cc(g+) and sum over kz
#np.sum(C+,axis=0)
#plt.loglog(n,C+)
#plt.show()
