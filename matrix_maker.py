import re
import math
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
import scipy as sc
import scipy.special
import matplotlib
import scipy.optimize
import scipy.linalg as lin
from scipy import interpolate



nu = 0.00222;

par = {'omt':0, 'omn':0, 'Ti0Te':1.0, 'kxmin':0.05, 'kxmax0':1.55, 'kymin':0.05, 'kymax0':1.55, 'kzmin':0.1, 'kzmax0':3.1, 'nkx0':32, 'nky0':64, 'nkz0':64,
       'nv0':48, 'nh0':1, 'nspec':1, 'hyp_x':0, 'hyp_y':0, 'hyp_z':0.0, 'hypx_order':0, 'hypy_order':0, 'hypz_order':0, 'hyp_v':0, 'hypv_order':0, 'hyp_conv':0,
       'num_k_hyp_conv':0, 'hyp_conv_ky':False, 'np_herm':24, 'np_kz':1, 'np_hank':1, 'np_spec':1, 'hyp_nu':0, 'nuno_closure':True, 'em_conserve':True}

#changed omn from 1.0 to 0 and omt from 10.0 to 0

def get_grids():
    """Returns kx,ky,kz,Hermite grids in the same form as used in the code \n
    kxgrid = 0, kxmin, . . . kxmax \n
    kygrid = 0, kymin, . . . kymax, kymax+kymin, -kymax, . . . -kymin"""
    kxgrid=np.arange((par['nkx0']))
    kxgrid=kxgrid*par['kxmin']
    kygrid=np.empty(par['nky0'])
    kzgrid=np.empty(par['nkz0'])
    herm_grid=np.arange(par['nv0'])
    herm_grid=1.0*herm_grid
    for i in range(par['nky0']/2):
        kygrid[par['nky0']-1-i]=-float(i+1)*par['kymin']
        kygrid[i]=float(i)*par['kymin']
    kygrid[par['nky0']/2]=par['nky0']/2*par['kymin']
    for i in range(par['nkz0']/2):
        kzgrid[par['nkz0']-1-i]=-float(i+1)*par['kzmin']
        kzgrid[i]=float(i)*par['kzmin']
    kzgrid[par['nkz0']/2]=par['nkz0']/2*par['kzmin']
    return kxgrid,kygrid,kzgrid,herm_grid



def get_gamma0():
    """Returns Gamma0 as a function of kx,ky"""
    kx,ky,kz,herm=get_grids()
    Gam0=np.empty((par['nkx0'],par['nky0']))
    for i in range(par['nkx0']):
        for j in range(par['nky0']):
            Gam0[i,j]=sc.special.ive(0,kx[i]**2+ky[j]**2)
    return Gam0



def matrix(kx,ky,kz,Gam0,nu):
    """My function that constructs the linear matrix for eigenvalue computation. This version takes as input the indices of the k components instead it their values."""

    k_perp2 = kx**2 + ky**2
    factor = math.exp(-k_perp2)/(1 + par['Ti0Te'] - Gam0)


    mat = np.zeros((par['nv0'],par['nv0']), dtype='complex64')
    phase_mix = np.zeros((par['nv0'],par['nv0']), dtype='complex64')
    diag = np.zeros((par['nv0'],par['nv0']), dtype='complex64')

    mat[0,0] =  1j*ky*factor*( par['omt']*k_perp2/2.0 - par['omn'] ) #drive terms in matrix
    mat[2,0] = -1j*par['omt']*ky*factor/np.sqrt(2.0) #drive terms
    mat[1,0] = -1j*kz*factor #landau damping

    for i in range(par['nv0']-1):       # be careful here since we miss the last element in the lower right corner
        if ( ( np.abs(ky)-par['kymin']*par['num_k_hyp_conv'] < 1.0e-6 and np.abs(ky) > 1.0e-6 and par['hyp_conv_ky'] ) or ( np.abs(kz)-par['kzmin']*par['num_k_hyp_conv'] < 1.0e-6 ) ):
            diag[i,i] = nu*float(i) + par['hyp_v']*(float(i)/float(par['nv0']-1))**par['hypv_order']\
                      + par['hyp_nu']*par['hyp_conv']\
                      + par['hyp_x']*(kx/par['kxmax0'])**par['hypx_order']\
                      + par['hyp_y']*(ky/par['kymax0'])**par['hypy_order']\
                      + par['hyp_z']*(kz/par['kzmax0'])**par['hypz_order']
        else:
            diag[i,i] = nu*float(i) + par['hyp_v']*(float(i)/float(par['nv0']-1))**par['hypv_order']\
                      + par['hyp_x']*(kx/par['kxmax0'])**par['hypx_order']\
                      + par['hyp_y']*(ky/par['kymax0'])**par['hypy_order']\
                      + par['hyp_z']*(kz/par['kzmax0'])**par['hypz_order']
        phase_mix[i,i+1] = np.sqrt(float(i+1)) #phase mixing term4
        phase_mix[i+1,i] = np.sqrt(float(i+1)) #phase mixing term

    if ( ( np.abs(ky)-par['kymin']*par['num_k_hyp_conv'] < 1.0e-6 and np.abs(ky) > 1.0e-6 and par['hyp_conv_ky'] ) or ( np.abs(kz)-par['kzmin']*par['num_k_hyp_conv'] < 1.0e-6 ) ):
        diag[par['nv0']-1,par['nv0']-1] = nu*float(par['nv0']-1) + par['hyp_v']\
                                        + par['hyp_nu']*par['hyp_conv']\
                                        + par['hyp_x']*(kx/par['kxmax0'])**par['hypx_order']\
                                        + par['hyp_y']*(ky/par['kymax0'])**par['hypy_order']\
                                        + par['hyp_z']*(kz/par['kzmax0'])**par['hypz_order']
    else:
        diag[par['nv0']-1,par['nv0']-1] = nu*float(par['nv0']-1) + par['hyp_v']\
                                        + par['hyp_x']*(kx/par['kxmax0'])**par['hypx_order']\
                                        + par['hyp_y']*(ky/par['kymax0'])**par['hypy_order']\
                                        + par['hyp_z']*(kz/par['kzmax0'])**par['hypz_order']

    if ( par['em_conserve'] == True ):
        diag[1,1] -= nu*1.0
        diag[2,2] -= nu*2.0

#    if (par['nuno_closure']):
#        diag[par['nv0']-1,par['nv0']-1] = -kz**2*par['nv0']/( (nu*par['nv0'] + par['hyp_v']*(float(par['nv0'])/float(par['nv0']-1)))**par['hypv_order'] )

    mat = mat - 1j*kz*phase_mix - diag
    print 'HERE I AM~~~~~~~~~~~~~~~~~~~~~~~~~~~~~', par['hyp_x']

    return mat



def get_spectrum(kx,ky,kz,Gam0,nu):
    """My function to compute the eigenvalues."""

    mat = matrix(kx,ky,kz,Gam0,nu)

    omega,evec = lin.eig(mat,right=True)
    #lin.eig gets eigen vectors, but for general arrays, might need different value for symmetric or hermitain matrices
    freq = np.imag(omega)
    growth = np.real(omega)

    return omega,freq, growth, evec



def plot_spectrum(kx,ky,kz,Gam0,nu):
    """Plots the spectrum obtained by the function get_spectrum_VB()"""

    omega, evec,freq, gam = get_spectrum(kx,ky,kz,Gam0,nu)

    plt.plot(freq, gam, 'b*')
    plt.grid()
    plt.show()



def ev_scan():
    """This scans over all eigenvalues for all wave numbers (kx,ky,kz) and finds the least damped/most unstable eigenvalue for each wave number.
        Hence, that way for each set of parameters, i.e., gradient lengths, we have only one eigenmode with the corresponding wave number."""

    kxset, kyset, kzset, hnset = get_grids()

    Gam0 = get_gamma0()

    k_vec = np.zeros(3)    # will hold the position in k-space of the most unstable eigenvalue for each temperature gradient
    gam = -100.0
    count = 0

    nx = len(kxset)
    ny = len(kyset)
    nz = len(kzset)

    time_1 = time.time()

    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                w, g = get_spectrum(kxset[ix],kyset[iy],kzset[iz],Gam0[ix,iy])
                i_max = np.argmax(g)
                if ( g[i_max] > gam ):
                    gam = g[i_max]
                    freq = w[i_max]
                    k_vec = [kxset[ix], kyset[iy], kzset[iz]]
                if ( g[i_max] > 0.0 ):
                    count = count + 1

    time_2 = time.time()

    print 'eigenvalue scan over kx, ky and kz took ', time_2 - time_1, 'sec'

    print 'for the termperature gradient of ', par['omt']
    print 'the most unstable eigenvalue is W =', freq, '+ i*', gam
    print 'it is located at k = (kx,ky,kz) = ', k_vec
    print
    print "there are ", count, "unstable modes out of total of ", par['nkx0']*par['nky0']*par['nkz0']
    print
