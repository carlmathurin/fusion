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

print 'max Growth =', max_g, '  with index:', max_g_i, 
"""
plt.plot(growth,freq,'b*')
plt.grid() # color='blue')
plt.xlabel('real [growth]')
plt.ylabel('imagainary [freq]')
plt.title('complex eigenvalues')
plt.show()
"""
