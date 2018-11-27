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


def max_g(growth):
    max_g = -10
    max_g_i = 0
    for i in range(len(growth)):
        if growth[i] > max_g:
            max_g = growth[i]
            max_g_i = i
    return max_g, max_g_i

g_0 = evec[:,max_g_i]

def find_error(g_eigen,g_calc,nmax):
    """"finds the error between g_calc and the eigenvector g"""
    error = np.empty(nmax,dtype=complex)

    for i in range(nmax):
        error[i] = abs(g_eigen[i] - g_calc[i])/ (abs(g_eigen[i]*g_calc[i])**.5)

    return error

def plot_error(error,herm_grid,nmax):
    plt.plot(herm_grid[0:nmax],error,label= 'err1 [|g_r - g_m| / |g_r*g_m|^.5 ]')
    plt.legend(loc='top left')
    plt.xlabel('herm #')
    plt.ylabel('error')
    plt.show()
    return

def plot_eig_vec(g_calc,g_eigen,herm_grid,nmax):
    plt.plot(herm_grid[0:nmax],abs(g_calc)**2,'b',label = 'calc')
    plt.plot(herm_grid[0:nmax],abs(g_eigen)**2,'r',label= 'matrix')
    plt.xlabel('herm #')
    plt.ylabel('g')
    plt.title('eigen vector [linear]')
    plt.legend()
    plt.show()

    plt.loglog(herm_grid[0:nmax],abs(g_calc)**2,'b',label = 'calc')
    # try plotting absolute value
    plt.loglog(herm_grid[0:nmax],abs(g_0)**2,'r',label= 'matrix')
    plt.xlabel('herm #')
    plt.ylabel('g')
    plt.title('eigen vector [logbase(10)]')
    plt.legend()
    #label= 'k_p ='+str(k_bin[j])
    plt.show()
    return

def calculate_g(kx,ky,kz,nmax,nu,gam0):
    """calculates g using eigenvector g"""

    omega, freq , growth, evec = mat.get_spectrum(kx,ky,kz,Gam0,nu)
    g_max, gi_max = max_g(growth)
    g_0 = evec[:,gi_max]

    nu_bar = nu/kz
    g_calc = np.empty(nmax,dtype=complex)

    for n in range(nmax):
        if n == 0 or n == 1 :
            g_calc[0] = g_0[0]
            g_calc[1] = g_0[1]
        else:
            g_calc[n] = ( (omega[gi_max]/kz + 1j*nu_bar*(n-1))*g_0[n-1] - (n-1)**(.5)*g_0[n-2])/ ((n)**(.5))
    err = find_error(g_0,g_calc,nmax)
    plot_eig_vec(g_calc,g_0,herm_grid,nmax)
    plot_error(err,herm_grid,nmax)

    return g_calc, g_0


"""
plt.plot(growth,freq,'b*')
plt.grid() # color='blue')
plt.xlabel('real [growth]')
plt.ylabel('imagainary [freq]')
plt.title('complex eigenvalues')
plt.show()
"""
