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

Gamma_0 = mat.get_gamma0()
print 'gamma: ', np.shape(Gamma_0)
print 'nu: ', par['nu']
print 'kx = ', kx[5],'ky = ', ky[5], 'kz = ', kz[5] , 'Gamma_0(index/value) [5,5]/',Gamma_0[5,5]
us_matrix = mat.matrix(kx[5],ky[5],kz[5],Gamma_0[5,5],par['nu'])
freq , growth, evec = mat.get_spectrum(kx[5],ky[5],kz[5],Gamma_0[5,5],par['nu'])
print 'freq: ',np.shape(freq), freq ,', growth: ', growth,', evec: ',np.shape(evec), evec

g0 = evec[:,0]
print 'g2 = ', g[2]

#0 = (nu_bar*n - 1j*w/kz)*g(n) + 1j*(sqrt(n+1)*g(n+1) + sqrt(n)*g(n-1))
n = 2
nu_bar = par['nu']/kz[5]
g3 = ( (1j*freq[0]/kz[5] - nu_bar*n)*g0[n] - 1j*n**(.5)*g0[n-1])/ (1j*( n+1)**(.5))
print 'g3 matrix vs calc: ', g[3], ' vs.', g3
# w is freq/eigenvalues
print 'Matrix: ', np.shape(us_matrix)
