import numpy as np
import matplotlib.pyplot as plt
import scipy as sc

import os
from sys import argv



n =[10**i for i in range(10)]
kz = .4
nu =.001/kz




#D = (n[j]*nu**2/4 + 1)**(.5)
#npl = D + n[j]**(.5)*nu/2

size = np.shape(n)
CapG = np.zeros(size)

for j in range(10):
    D = (n[j]*nu**2/4 + 1)**(.5)
    npl = D + n[j]**(.5)*nu/2
    CapG[j] = (-1*np.sign(kz))**n[j]* np.exp((n[j]**(.5))*D/(np.absolute(nu))) \
    * (npl)**(-np.sign(kz)*(n[j]+2*nu**(-2)-.5))/(n[j]**(.25)*D**(.5))
#CapG = (- 1j*np.sign(kz))**n[j] * np.exp((n[j]**(.5)*D/(np.absolute(nu)))* \
#(Npl)**(-np.sign(kz*(n[j]+2*nu**(-2)-.5))/(n[j]**(.25)*D**(.5))
