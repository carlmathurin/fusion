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
import matrix_maker as mat

dd.read_parameters()

diagdir = '/scratch/01658/drhatch/dna_out'
par['diagdir']="\'"+diagdir+"\'"
time = dd.get_time_from_gout()
kx,ky,kz,herm_grid = dd.get_grids()

par['omn'] = 0
par['omt'] = 0
par['hyp_x'] = 0

Gamma_0 = mat.get_gamma0()
print 'gamma: ', np.shape(Gamma_0)
print 'nu: ', par['nu']
print 'kx = ', kx[5],'ky = ', ky[5], 'kz = ', kz[5] , 'Gamma_0(index/value) [5,5]/',Gamma_0[5,5]
us_matrix = mat.matrix(kx[5],ky[5],kz[5],Gamma_0[5,5],par['nu'])
omega, freq , growth, evec = mat.get_spectrum(kx[5],ky[5],kz[5],Gamma_0[5,5],par['nu'])
print 'freq: ',np.shape(freq), freq, ', growth: ', np.shape(growth),', evec: ',np.shape(evec), evec[:,0]
print 'eval: ', omega

g_0 = evec[:,0]
print 'g1 = ', g_0[1]

#0 = (nu_bar*n - 1j*w/kz)*g(n) + 1j*(sqrt(n+1)*g(n+1) + sqrt(n)*g(n-1))
nmax = 48
nu_bar = par['nu']/kz[5]
g_calc = np.empty(48,dtype=complex)

for n in range(48):
    if n == 0 or n == 1 :
        g_calc[0] = g_0[0]
        g_calc[1] = g_0[1]
    else:
        g_calc[n] = ( (omega[0]/kz[5] + 1j*nu_bar*(n-1))*g_0[n-1] - (n-1)**(.5)*g_0[n-2])/ ((n)**(.5))

print 'g2 matrix vs calc: vs error '
error = np.empty(48,dtype=complex)
error2 = np.empty(48,dtype=complex)
for i in range(48):
    error[i] = abs(g_0[i] - g_calc[i])/ (abs(g_0[n]*g_calc[i])**.5)
    error2[i] = abs(g_0[i] - g_calc[i])/ (abs(g_0[n]))
    print 'n =', i,'  ', g_0[i], 'vs', g_calc[i] , '  error =', error[i]

print 'herm ',np.shape(herm_grid[0:48]), 'g_calc ', np.shape(g_calc)
print 'herm :', herm_grid[0:48]

plt.plot(herm_grid[0:48],g_calc,'b',label = 'calc')
# try plotting absolute value
plt.plot(herm_grid[0:48],g_0,'r',label= 'matrix')
plt.xlabel('herm #')
plt.ylabel('g')
plt.title('first g eigen vector [linear]')
plt.legend()
#label= 'k_p ='+str(k_bin[j])
plt.show()

plt.loglog(herm_grid[0:48],abs(g_calc),'b',label = 'calc')
# try plotting absolute value
plt.loglog(herm_grid[0:48],abs(g_0),'r',label= 'matrix')
plt.xlabel('herm #')
plt.ylabel('g')
plt.title('first g eigen vector [logbase(10)]')
plt.legend()
#label= 'k_p ='+str(k_bin[j])
plt.show()



# maybe look for most unstable eigenvalue


print 'error ', error
plt.plot(herm_grid[0:48],error,label= 'err1 [|g_r - g_m| / |g_r*g_m|^.5 ]')
plt.plot(herm_grid[0:48],error2,label= 'err2 [|g_r - g_m| / |g_m|]')
plt.legend(loc='top left')
plt.xlabel('herm #')
plt.ylabel('error')
plt.show()
# plt.xlabel('Hermite n')
# plt.ylabel('g')
# plt.title('spectra')

"""
f = open('g_data.txt','w+')
print g_calc[2].real, g_calc[2].imag
f.write('n      |g_0             |G_calc\n')
for i in range(48):
    f.write('%d     |%e + i(%e)      |%e + i(%e) \n' %( i,  g_0[i].real, g_0[i].imag, g_calc[i].real, g_calc[i].imag ))
f.close()
# w is freq/eigenvalues
"""
print 'Matrix: ', np.shape(us_matrix)
print 'here ~~', par['hyp_x']
