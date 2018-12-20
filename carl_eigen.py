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
nmax = 40
mat.set_nmax(nmax)

par['nu'] = 0.05
print 'nu = '  , par['nu']
print 'kz = ', kz


par['omn'] = 0
par['omt'] = 0
par['hyp_x'] = 0
par['hyp_y'] = 0
par['hyp_v'] = 0
par['nuno_closure'] = False

kz_t = .3

Gamma_0 = mat.get_gamma0()
#print 'gamma: ', np.shape(Gamma_0), Gamma_0
print 'nu: ', par['nu']
print 'kx = ', kx[5],'ky = ', ky[5], 'kz = ', kz[5] , 'Gamma_0(index/value) [5,5]/',Gamma_0[5,5]
us_matrix = mat.matrix(0,0,kz_t,Gamma_0[0,0],par['nu'])

omega, freq , growth, evec = mat.get_spectrum(0.0,0.0,kz_t,Gamma_0[0,0],par['nu'])

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
Chi = np.exp(0.0)/(par['Ti0Te'] + 1.0 - Gamma_0[0,0])
nu_bar = par['nu']/kz_t
g_calc = np.empty(nmax,dtype=complex)

#omega = 1j*omega
for n in range(nmax):
    if n == 0 or n == 1 or n == 2 :
        g_calc[0] = g_0[0]
        g_calc[1] = omega[max_g_i]/kz[5]*g_calc[0]
        g_calc[2] = -(Chi + 1.0)/np.sqrt(2.0)*g_calc[0] + (omega[max_g_i]/np.sqrt(2.0)/kz_t + 1j*nu_bar/np.sqrt(2.0))*g_calc[1]
    else:
        g_calc[n] = ( (omega[max_g_i]/kz_t + 1j*nu_bar*(n-1))*g_calc[n-1] - (n-1)**(.5)*g_calc[n-2])/ ((n)**(.5))

#print 'g2 matrix vs calc: vs error '
error = np.empty(nmax,dtype=complex)
error2 = np.empty(nmax,dtype=complex)
for i in range(nmax):
    error[i] = abs(g_0[i] - g_calc[i])/ (abs(g_0[i]*g_calc[i])**.5)
    error2[i] = abs(g_0[i] - g_calc[i])/ (abs(g_0[i]))
    #print 'n =', i,'  ', g_0[i], 'vs', g_calc[i] , '  error =', error[i]

print 'herm ',np.shape(herm_grid[0:nmax]), 'g_calc ', np.shape(g_calc)
print 'herm :', herm_grid[0:nmax]

print 'eigenvector check: ', np.linalg.norm(np.matmul(us_matrix, g_0) - omega[max_g_i]*g_0)

print 'eigen vector check 2', np.divide(np.matmul(us_matrix, g_0) ,g_0)
print 'eigen value: ', omega[max_g_i]

plt.plot(growth,freq,'b*')
plt.plot(np.real(omega), np.imag(omega), 'rx')
plt.grid() # color='blue')
plt.xlabel('real [growth]')
plt.ylabel('imagainary [freq]')
plt.title('complex eigenvalues')
plt.show()





#mat.plot_spectrum(kx[5],ky[5],kz[5],Gamma_0[5,5],par['nu'])
fig, (f1,f2) = plt.subplots(1,2)
f1.plot(herm_grid[0:nmax],abs(g_calc)**2,'bo-',label = 'calc')
# try plotting absolute value
f1.plot(herm_grid[0:nmax],abs(g_0)**2,'rx-',label= 'matrix')
plt.xlabel('herm #')
plt.ylabel('g')
plt.title('first g eigen vector [linear]')
plt.legend()
#label= 'k_p ='+str(k_bin[j])
#plt.show()

f2.loglog(herm_grid[0:nmax],abs(g_calc)**2,'bo-',label = 'calc')
# try plotting absolute value
f2.loglog(herm_grid[0:nmax],abs(g_0)**2,'rx-',label= 'matrix')
plt.xlabel('herm #')
plt.ylabel('g')
plt.title('eigen vector [logbase(10)]')
plt.legend()
#label= 'k_p ='+str(k_bin[j])
plt.show()

plt.semilogy(herm_grid[0:nmax],abs(g_calc)**2,'bo-',label = 'calc')
plt.semilogy(herm_grid[0:nmax],abs(g_0)**2,'rx-',label = 'matrix')
plt.show()

print 'g_calc from 0-3: ', g_calc[0:3]
print 'g_matrix from 0-3: ', g_0[0:3]


# maybe look for most unstable eigenvalue


print 'error ', error
plt.plot(herm_grid[0:nmax],error,label= 'err1 [|g_r - g_m| / |g_r*g_m|^.5 ]')
#plt.plot(herm_grid[0:nmax],error2,label= 'err2 [|g_r - g_m| / |g_m|]')
plt.legend(loc='top left')
plt.xlabel('herm #')
plt.ylabel('error')
plt.show()
# plt.xlabel('Hermite n')
# plt.ylabel('g')
# plt.title('spectra')


# w is freq/eigenvalues
