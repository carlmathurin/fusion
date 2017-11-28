import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
from config import *
import new_dd as dd
from ev_diags import *
import time as t2
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

start_time=time[100]
k_perp = [.15 + i*.15 for i in range(10)]
kxind=np.empty(10)
kyind=np.empty(10)
for i in range(10):
    kxind[i]=i*par['nkx0']/20
    kyind[i]=i*par['nky0']/20
print 'kxindex' , kxind
print 'kyindex' , kyind
