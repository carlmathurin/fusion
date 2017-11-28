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

k_bin=np.empty(10)

xmax = par['nkx0']
kzind = 2*par['nkz0']/20

entn_sum=np.zeros((par['nv0'],11,11),dtype='float')
entnp_sum=np.zeros((par['nv0'],11,11),dtype='float')
entnm_sum=np.zeros((par['nv0'],11,11),dtype='float')

for i in range(10):
    k_bin[i]= .15*i
for i in range(xmax):
    for j in range(xmax):
        t = 12
        k_perp = sqrt(kx[i]**2 + ky[j]**2)
        if  k_perp < k_bin[1]:
            t = 0
        elif k_bin[1] <= k_perp < k_bin[2]:
            t= 1
        elif k_bin[2] <= k_perp < k_bin[3]:
            t= 2
        elif k_bin[3] <= k_perp < k_bin[4]:
            t= 3
        elif k_bin[4] <= k_perp < k_bin[5]:
            t= 4
        elif k_bin[5] <= k_perp < k_bin[6]:
            t= 5
        elif k_bin[6] <= k_perp < k_bin[7]:
            t= 6
        elif k_bin[7] <= k_perp < k_bin[8]:
            t= 7
        elif k_bin[8] <= k_perp < k_bin[9]:
            t= 8
        elif k_bin[9] <= k_perp :
            t= 9

        entn_sum[:,kzind,t] = entn_sum[:,kzind,t] + dd.get_entropy_hermite2(g_t0,i,j)

print np.shape(entn_sum)
print 'entn sum:' , entn_sum
