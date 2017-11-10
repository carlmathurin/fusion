import numpy as np
import matplotlib.pyplot as plt
import scipy as sc

import os
from sys import argv



n =[1+i for i in range(100)]
kz = [.2 + i*.2 for i in range(100)]
nuu =.002





#D = (n[j]*nu**2/4 + 1)**(.5)
#npl = D + n[j]**(.5)*nu/2

size = np.shape(n)
CapG = np.empty([100,100])
Cap2G = np.empty([100,100])

for j in range(100):
    print 'Kz =', kz[j]
    for i in range(100):
        print 'N =', n[i]
        nu = nuu/kz[j]
        D = ((n[i]*nu**2)/4 + 1)**(.5)
        npl = D + n[i]**(.5)*nu/2
        CapG[j,i] = (np.exp((n[i]**(.5))*\
            D/(np.absolute(nu))) \
            * (npl)**(-np.sign(kz[j])*(n[i]+2*nu**(-2)-.5))\
            /(n[i]**(.25)*D**(.5)))\
            #*(-1j*np.sign(kz[j]))**n[i]
        Cap2G[j,i] = CapG[j,i]* np.conjugate(CapG[j,i])


    plt.loglog(n,Cap2G[j,:],label= 'G' + 'kz ('+ str(kz[j])+')')
    plt.title('G ')
    plt.xlabel('Hermite n')
    plt.ylabel('CapG')
print CapG
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.show()



#(-1*np.cos(n[i]*np.pi/2)*np.sign(kz[j]))**n[i]\

#CapG = (- 1j*np.sign(kz))**n[j] * np.exp((n[j]**(.5)*D/(np.absolute(nu)))* \
#(Npl)**(-np.sign(kz*(n[j]+2*nu**(-2)-.5))/(n[j]**(.25)*D**(.5))
