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
nmax = 36
mat.set_nmax(nmax)
Gam0 = mat.get_gamma0()

par['nu'] = 0.002222
print 'nu = '  , par['nu']


par['omn'] = 0
par['omt'] = 0
par['hyp_x'] = 0

omega, freq , growth, evec = mat.get_spectrum(kx[5],ky[5],kz[5],Gam0[5,5],par['nu'])

max_g = -10
max_g_i = 0

for i in range(len(growth)):
    if growth[i] > max_g:
        max_g = growth[i]
        max_g_i = i


g_0 = evec[:,max_g_i]

nu_bar = par['nu']/kz[5]
g_calc = np.empty(nmax,dtype=complex)

for n in range(nmax):
    if n == 0 or n == 1 :
        g_calc[0] = g_0[0]
        g_calc[1] = g_0[1]
    else:
        g_calc[n] = ( (omega[max_g_i]/kz[5] + 1j*nu_bar*(n-1))*g_0[n-1] - (n-1)**(.5)*g_0[n-2])/ ((n)**(.5))

error = np.empty(nmax,dtype=complex)

for i in range(nmax):
    error[i] = abs(g_0[i] - g_calc[i])/ (abs(g_0[i]*g_calc[i])**.5)
    error2[i] = abs(g_0[i] - g_calc[i])/ (abs(g_0[i]))

def find_error(g_eigen,g_calc):
    for i in range(nmax):
        error[i] = abs(g_eigen[i] - g_calc[i])/ (abs(g_eigen[i]*g_calc[i])**.5)

"""
plt.plot(growth,freq,'b*')
plt.grid() # color='blue')
plt.xlabel('real [growth]')
plt.ylabel('imagainary [freq]')
plt.title('complex eigenvalues')
plt.show()
"""
