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
print Gamma_0
us_matrix = mat.matrix(par[nkx0],par[nky0],par[nkz0],Gamma_0,par[nu])

print us_matrix
