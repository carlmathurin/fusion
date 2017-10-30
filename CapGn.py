import numpy as np
import matplotlib.pyplot as plt
import scipy as sc

import os
from sys import argv



n =[1+i for i in range(10)]
kz = [.2 + i*.2 for i in range(10)]
nuu =.002





#D = (n[j]*nu**2/4 + 1)**(.5)
#npl = D + n[j]**(.5)*nu/2

size = np.shape(n)
CapG = np.empty([10,10])

for j in range(10):
    print 'Kz =', kz[j]
    for i in range(10):
        print 'N =', n[i]
        nu = .01 #nuu/kz[i]
        D = ((n[i]*nu**2)/4 + 1)**(.5)
        npl = D + n[i]**(.5)*nu/2
        CapG[j,i] = (-1*np.sign(kz[j]))**n[i]\
        * np.exp((n[i]**(.5))*\
        D/(np.absolute(nu))) \
        * (npl)**(-np.sign(kz[j])*(n[i]+2*nu**(-2)-.5))\
        /(n[i]**(.25)*D**(.5))

    plt.loglog(n,CapG[j,:],basex=10,basey=10,label=\
              'CapG' +' (k_z='+str(kz)+')')
    plt.xlabel('Hermite n')
    plt.ylabel('CapG')
    plt.legend(loc='lower left')

#CapG = (- 1j*np.sign(kz))**n[j] * np.exp((n[j]**(.5)*D/(np.absolute(nu)))* \
#(Npl)**(-np.sign(kz*(n[j]+2*nu**(-2)-.5))/(n[j]**(.25)*D**(.5))
