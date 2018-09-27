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

Gamma_0 = mat.get_gamma0()
print 'gamma: ', np.shape(Gamma_0), Gamma_0
print 'nu: ', par['nu']
print 'kx = ', kx[20],'ky = ', ky[20], 'kz = ', kz[20]
us_matrix = mat.matrix(kx[20],ky[20],kz[20],Gamma_0[20],par['nu'])

print us_matrix
