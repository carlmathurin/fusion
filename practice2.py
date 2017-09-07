import kz3_herm as d2
import matplotlib.pyplot as plt

ap = ''
bp = ''
cp = ''
ep = ''
am = ''
bm = ''
cm = ''
em = ''
a, a1, ax, ay  = d2.diags_cp_cm(6,16)
b, b1, bx, by = d2.diags_cp_cm(10,50)
c, c1, cx, cy  = d2.diags_cp_cm(30,63)
f, f1, fx, fy = d2.diags_cp_cm(25,38)
e, e1, n, ex, ey, ez  = d2.diags_kz(16,5,10)
g, g1, n, gx, gy, gz = d2.diags_kz(5,45,30)


ss = ""
ap = 'C+ (',str(ax),',',str(ay), ')'
ap = ss.join(ap)

bp = 'C+ (',str(bx),',',str(by), ')'
bp = ss.join(bp)

cp = 'C+ (',str(cx),',',str(cy), ')'
cp = ss.join(cp)

fp = 'C+ (',str(fx),',',str(fy), ')'
fp = ss.join(fp)

ep = 'C+ (',str(ex),',',str(ey),',',str(ez), ')'
ep = ss.join(ep)

gp = 'C+ (',str(gx),',',str(gy),',',str(gz), ')'
gp = ss.join(gp)

am = 'C- (',str(ax),',',str(ay), ')'
am = ss.join(am)

bm = 'C- (',str(bx),',',str(by), ')'
bm = ss.join(bm)

cm = 'C- (',str(cx),',',str(cy), ')'
cm = ss.join(cm)

fm = 'C+ (',str(fx),',',str(fy), ')'
fm = ss.join(fm)

em = 'C- (',str(ex),',',str(ey),',',str(ez), ')'
em = ss.join(em)

gm = 'C- (',str(gx),',',str(gy),',',str(gz), ')'
gm = ss.join(gm)


cp1, = plt.loglog(n,a,)
cp2, = plt.loglog(n,b,)
cp3, = plt.loglog(n,c,)
cp4, = plt.loglog(n,e,) 
cp5, = plt.loglog(n,a1,)
cp6, = plt.loglog(n,b1,)
cp7, = plt.loglog(n,c1,)
cp8, = plt.loglog(n,e1,)
cp11, = plt.loglog(n,g)
cp12, = plt.loglog(n,g1)
plt.legend([cp1,cp2,cp3,cp4,cp5,cp6,cp7,cp8,cp11,cp12], [ap,bp,cp,ep,am,bm,cm,em,gp,gm], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.xlabel('Hermite Number')
plt.show()
