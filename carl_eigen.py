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
par['nu'] = .0015
print par['nu']

par['omn'] = 0
par['omt'] = 0
par['hyp_x'] = 0

Gamma_0 = mat.get_gamma0()
print 'gamma: ', np.shape(Gamma_0), Gamma_0
print 'nu: ', par['nu']
print 'kx = ', kx[5],'ky = ', ky[5], 'kz = ', kz[5] , 'Gamma_0(index/value) [5,5]/',Gamma_0[5,5]
us_matrix = mat.matrix(kx[5],ky[5],kz[5],Gamma_0[5,5],par['nu'])
omega, freq , growth, evec = mat.get_spectrum(kx[5],ky[5],kz[5],Gamma_0[5,5],par['nu'])
print 'freq: ',np.shape(freq), freq, ', growth: ', growth,', evec: ',np.shape(evec), evec[:,0]
print 'eval: ', omega

max_g = -10
max_g_i = 0
g_repeat = 0

for i in range(len(growth)):
    if growth[i] > max_g:
        max_g = growth[i]
        max_g_i = i
    elif growth[i] == max_g:
        g_repeat = i


print 'max Growth =', max_g, '  with index:', max_g_i, 'and repeats at', g_repeat



g_0 = evec[:,max_g_i]
print 'g = ', g_0

#0 = (nu_bar*n - 1j*w/kz)*g(n) + 1j*(sqrt(n+1)*g(n+1) + sqrt(n)*g(n-1))
nmax = get_nmax()
print 'nmax = ', nmax
nu_bar = par['nu']/kz[5]
g_calc = np.empty(nmax,dtype=complex)

for n in range(nmax):
    if n == 0 or n == 1 :
        g_calc[0] = g_0[0]
        g_calc[1] = g_0[1]
    else:
        g_calc[n] = ( (omega[0]/kz[5] + 1j*nu_bar*(n-1))*g_0[n-1] - (n-1)**(.5)*g_0[n-2])/ ((n)**(.5))

#print 'g2 matrix vs calc: vs error '
error = np.empty(nmax,dtype=complex)
error2 = np.empty(nmax,dtype=complex)
for i in range(nmax):
    error[i] = abs(g_0[i] - g_calc[i])/ (abs(g_0[n]*g_calc[i])**.5)
    error2[i] = abs(g_0[i] - g_calc[i])/ (abs(g_0[n]))
    #print 'n =', i,'  ', g_0[i], 'vs', g_calc[i] , '  error =', error[i]

print 'herm ',np.shape(herm_grid[0:nmax]), 'g_calc ', np.shape(g_calc)
print 'herm :', herm_grid[0:nmax]

plt.plot(growth,freq,'b*')
plt.grid() # color='blue')
plt.xlabel('real [growth]')
plt.ylabel('imagainary [freq]')
plt.title('complex eigenvalues')
plt.show()





#mat.plot_spectrum(kx[5],ky[5],kz[5],Gamma_0[5,5],par['nu'])

plt.plot(herm_grid[0:nmax],g_calc,'b',label = 'calc')
# try plotting absolute value
plt.plot(herm_grid[0:nmax],g_0,'r',label= 'matrix')
plt.xlabel('herm #')
plt.ylabel('g')
plt.title('first g eigen vector [linear]')
plt.legend()
#label= 'k_p ='+str(k_bin[j])
plt.show()

plt.loglog(herm_grid[0:nmax],abs(g_calc),'b',label = 'calc')
# try plotting absolute value
plt.loglog(herm_grid[0:nmax],abs(g_0),'r',label= 'matrix')
plt.xlabel('herm #')
plt.ylabel('g')
plt.title('first g eigen vector [logbase(10)]')
plt.legend()
#label= 'k_p ='+str(k_bin[j])
plt.show()



# maybe look for most unstable eigenvalue


print 'error ', error
plt.plot(herm_grid[0:nmax],error,label= 'err1 [|g_r - g_m| / |g_r*g_m|^.5 ]')
plt.plot(herm_grid[0:nmax],error2,label= 'err2 [|g_r - g_m| / |g_m|]')
plt.legend(loc='top left')
plt.xlabel('herm #')
plt.ylabel('error')
plt.show()
# plt.xlabel('Hermite n')
# plt.ylabel('g')
# plt.title('spectra')


# w is freq/eigenvalues
