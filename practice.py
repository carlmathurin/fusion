import new_diag2 as d2
import matplotlib.pyplot as plt

a, a1 = d2.diags_cp_cm(6,6)
b, b1 = d2.diags_cp_cm(10,10)
c, c1 = d2.diags_cp_cm(30,70)
e, e1, n = d2.diags_kz(10,10,10)

cp1, = plt.loglog(n,a,label = 'C+ 5,5' )
cp2, = plt.loglog(n,b,label = 'C+ 9,9' )
cp3, = plt.loglog(n,c,label = 'C+ 31,71' )
cp4, = plt.loglog(n,e,label = 'C+ 11,11,11')
cp1, = plt.loglog(n,a1,label = 'C+ 5,5' )
cp2, = plt.loglog(n,b1,label = 'C+ 9,9' )
cp3, = plt.loglog(n,c1,label = 'C+ 31,71' )
cp4, = plt.loglog(n,e1,label = 'C+ 11,11,11')
plt.legend([cp1,cp2,cp3])
plt.show()
