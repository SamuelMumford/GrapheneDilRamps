# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 09:58:11 2021

@author: KGB
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

def prefilter(f0, b):
    keep = (f0 > 1000)
    f0 = f0[keep]
    b = b[keep]
    keep2 = np.insert((np.diff(f0) < .05), 0, True)
    f0 = f0[keep2]
    b = b[keep2]
    return f0, b

pt0V = '/home/sam/Documents/DriftTests/BRamp0V_6_21_slow.txt'
pt2V = '/home/sam/Documents/DriftTests/BRamp2V_6_21_slow.txt'
ptm2V = '/home/sam/Documents/DriftTests/BRamp_min2V_6_7_21.txt'


F0V = np.loadtxt(pt0V, delimiter="\t")
B0V = np.array(F0V[:, 0])
f0V = np.array(F0V[:, 1])
f0V, B0V = prefilter(f0V, B0V)

F2V = np.loadtxt(pt2V, delimiter="\t")
B2V = np.array(F2V[:, 0])
f2V = np.array(F2V[:, 1])
f2V, B2V = prefilter(f2V, B2V)

Fm2V = np.loadtxt(ptm2V, delimiter="\t")
Bm2V = np.array(Fm2V[:, 0])
fm2V = np.array(Fm2V[:, 1])
fm2V, Bm2V = prefilter(fm2V, Bm2V)

tol = .0005
m0V0B = (np.abs(B0V) < tol)
m2V0B = (np.abs(B2V) < tol)
mm2V0B = (np.abs(Bm2V) < tol)
f0V_zeroB = f0V[m0V0B]
f2V_zeroB = f2V[m2V0B]
fm2V_zeroB = fm2V[mm2V0B]

#Fm2V = np.loadtxt(ptm2V, delimiter="\t")
#Bm2V = np.array(Fm2V[:, 0])
#fm2V = np.array(Fm2V[:, 1])
#fm2V, B2V = prefilter(fm2V, Bm2V)

tol_final =1.1
sm = .1

plt.scatter(B0V, f0V, s=10, linewidth = .5,  zorder=1, marker = '.', c = 'black', label = 'O V')
plt.scatter(B2V, f2V, s=10, linewidth = .5,  zorder=1, marker = '.', c = 'red', label = '2 V')
plt.scatter(Bm2V, fm2V, s=10, linewidth = .5,  zorder=1, marker = '.', c = 'green', label = '-2 V')
plt.grid()
plt.xlabel('B (T)')
plt.ylabel('f0 (Hz)')
plt.legend(title = 'Gate Voltage')
plt.savefig('/home/sam/Documents/DriftTests/SweepBs_f0_odd.pdf', bbox_inches="tight")
plt.show()

def cubicC(x, a, b):
    return a + b*x**2

def makeZeroB(Bs, f0s, tol, smWin):
    Bwin = (np.abs(Bs - tol) < smWin)
    Bfit = Bs[Bwin]
    f0fit = f0s[Bwin]
    popt, _ = curve_fit(cubicC, Bfit, f0fit, p0=[0,0])
    return popt
    
p0V = makeZeroB(B0V, f0V-np.mean(f0V_zeroB), tol_final, sm)
p2V = makeZeroB(B2V, f2V-np.mean(f2V_zeroB), tol_final, sm)
pm2V = makeZeroB(Bm2V, fm2V-np.mean(fm2V_zeroB), tol_final, sm)

plotP = np.linspace(0, tol_final)
plotf0 = cubicC(plotP, 0, p0V[1])
plotf2 = cubicC(plotP, 0, p2V[1])
plotfm2 = cubicC(plotP, 0, pm2V[1])

maski = (np.abs(B0V) > tol_final)
maskii = (np.abs(B2V) > tol_final)
maskiii = (np.abs(Bm2V) > tol_final)

al0 = .1

plt.scatter(np.abs(B0V)[maski], (f0V-np.mean(f0V_zeroB)-p0V[0])[maski], s=10, linewidth = .5,  zorder=1, marker = '.', c = 'black', label = 'O V')
#plt.scatter(np.abs(B0V), (f0V-np.mean(f0V_zeroB)-p0V[0]), s=10, linewidth = .5,  zorder=1, marker = '.', c = 'black', alpha=al0)
plt.plot(plotP, plotf0, color = 'black')
plt.scatter(np.abs(B2V)[maskii], (f2V-np.mean(f2V_zeroB)-p2V[0])[maskii], s=10, linewidth = .5,  zorder=1, marker = '.', c = 'red', label = '2 V')
#plt.scatter(np.abs(B2V), (f2V-np.mean(f2V_zeroB)-p2V[0]), s=10, linewidth = .5,  zorder=1, marker = '.', c = 'red', alpha=al0)
plt.plot(plotP, plotf2, color = 'red')
plt.scatter(np.abs(Bm2V)[maskiii], (fm2V-np.mean(fm2V_zeroB)-pm2V[0])[maskiii], s=10, linewidth = .5,  zorder=1, marker = '.', c = 'green', label = '-2 V')
#plt.scatter(np.abs(Bm2V), (fm2V-np.mean(fm2V_zeroB)-pm2V[0]), s=10, linewidth = .5,  zorder=1, marker = '.', c = 'green', alpha=al0)
plt.plot(plotP, plotfm2, color = 'green')
plt.grid()
plt.axvline(tol_final, alpha=.5)
#plt.axvline(tol_final+sm, color = 'blue')
#plt.axvline(tol_final-sm, color = 'blue')
plt.xlabel('|B| (T)')
plt.ylabel('f0 (Hz)')
plt.legend(title = 'Gate Voltage')
plt.savefig('/home/sam/Documents/DriftTests/SweepBs_f0.pdf', bbox_inches="tight")
plt.show()

offy = 0
offy2 = 0
Bdata = np.array(B0V[maski])
fdata = np.array((f0V-np.mean(f0V_zeroB)-p0V[0])[maski] - offy)

def FF(x, a, b, d, e):
    c=1.3
    q1 = (2*a + 1)/(a)
    l1 = 1/(2*a)
    polyP = b*x*(q1/np.tanh(q1*c*a*x) - l1/(np.tanh(c*a*x*l1)))
    polyP += d*np.sign(x)*x
    polyP += e
    #polyP += g*np.abs(x)
    return polyP

ff = FF(Bdata,  2, 2, -2, 1.5)
p0 = [2, 2, -2, 1.5]
popt, pcov = curve_fit(FF, Bdata, fdata, p0=p0, maxfev=5000)
ff = FF(Bdata, popt[0],popt[1], popt[2], popt[3])
print(popt)
print(np.sqrt(np.diag(pcov)))

plt.plot(Bdata, ff, color = 'black', label = 'fit')
plt.plot(Bdata, fdata, color = 'red', label = 'data')
plt.ylabel(r'$f_{0}$ (Hz)')
plt.xlabel('B (T)')
plt.legend(loc = 'best')
plt.show()

bd0 = Bdata
fd0 = (fdata - ff)/Bdata

Bdata = np.array(Bm2V[maskiii])
fdata = np.array((fm2V-np.mean(fm2V_zeroB)-pm2V[0])[maskiii])

ff = FF(Bdata,  2, 2, -2, 1.5)
p0 = [2, 2, -2, 1.5]
popt, pcov = curve_fit(FF, Bdata, fdata, p0=p0, maxfev=5000)
ff = FF(Bdata, popt[0],popt[1], popt[2], popt[3])
print(popt)
print(np.sqrt(np.diag(pcov)))

bdm2 = Bdata
fdm2 = (fdata - ff)/Bdata

Bdata = np.array(B2V[maskii])
fdata = np.array((f2V-np.mean(f2V_zeroB)-p2V[0])[maskii])

ff = FF(Bdata,  2, 2, -2, 1.5)
p0 = [2, 2, -2, 1.5]
popt, pcov = curve_fit(FF, Bdata, fdata, p0=p0, maxfev=5000)
ff = FF(Bdata, popt[0],popt[1], popt[2], popt[3])
print(popt)
print(np.sqrt(np.diag(pcov)))

bd2 = Bdata
fd2 = (fdata - ff)/Bdata

cf = 5E7
plt.scatter(1/bdm2, fdm2*cf, s=1, linewidth = .5,  zorder=1, color = 'green', label = '-2 V')
L = sorted(zip(1/bdm2,fdm2*cf))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'green', alpha = .25)
plt.scatter(1/bd0, fd0*cf, s=1, linewidth = .5,  zorder=1, color = 'black', label = '0 V')
L = sorted(zip(1/bd0,fd0*cf))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'black', alpha = .25)
plt.scatter(1/bd2, fd2*cf, s=1, linewidth = .5,  zorder=1, color = 'red', label = '2 V')
L = sorted(zip(1/bd2,fd2*cf))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'red', alpha = .25)
plt.ylabel(r'm ($\mu_{B})$')
plt.xlabel(r'1/B (T$^{-1}$)')
plt.legend(loc = 'best')
plt.show()


plt.scatter(1/bdm2, fdm2*bdm2, s=1, linewidth = .5,  zorder=1, color = 'green', label = '-2 V')
L = sorted(zip(1/bdm2,fdm2*bdm2))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'green', alpha = .25)
plt.scatter(1/bd0, fd0*bd0, s=1, linewidth = .5,  zorder=1, color = 'black', label = '0 V')
L = sorted(zip(1/bd0,fd0*bd0))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'black', alpha = .25)
plt.scatter(1/bd2, fd2*bd2, s=1, linewidth = .5,  zorder=1, color = 'red', label = '2 V')
L = sorted(zip(1/bd2,fd2*bd2))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'red', alpha = .25)
plt.ylabel(r'$\delta f_{0}$')
plt.xlabel(r'1/B (T$^{-1}$)')
plt.legend(loc = 'best')
plt.show()
#plt.scatter(np.abs(B0V), f0V-np.mean(f0V_zeroB)-p0V[0], s=10, linewidth = .5,  zorder=1, marker = '.', c = 'black', label = 'O V')
#plt.scatter(np.abs(B2V), f2V-np.mean(f2V_zeroB)-p2V[0], s=10, linewidth = .5,  zorder=1, marker = '.', c = 'red', label = '2 V')
#plt.scatter(np.abs(Bm2V), fm2V-np.mean(fm2V_zeroB)-pm2V[0], s=10, linewidth = .5,  zorder=1, marker = '.', c = 'green', label = '-2 V')
#plt.grid()
#plt.axvline(tol_final, color = 'black')
#plt.axvline(tol_final+sm, color = 'blue')
#plt.axvline(tol_final-sm, color = 'blue')
#plt.xlim(0, 1)
#plt.xlabel('|B| (T)')
#plt.ylabel('f0 (Hz)')
#plt.legend(title = 'Gate Voltage')
#plt.savefig('/home/sam/Documents/DriftTests/SweepBs_f0.pdf', bbox_inches="tight")
#plt.show()

mask0V_high = (np.abs(B0V) >= tol_final)
mask2V_high = (np.abs(B2V) >= tol_final)
maskm2V_high = (np.abs(Bm2V) >= tol_final)

f0V_mask = f0V[mask0V_high]
f2V_mask = f2V[mask2V_high]
fm2V_mask = fm2V[maskm2V_high]

B0V_mask = B0V[mask0V_high]
B2V_mask = B2V[mask2V_high]
Bm2V_mask = Bm2V[maskm2V_high]

magBase0V = (f0V_mask-np.mean(f0V_zeroB)-p0V[0])/B0V_mask
magBase2V = (f2V_mask-np.mean(f2V_zeroB)-p2V[0])/B2V_mask
magBasem2V = (fm2V_mask-np.mean(fm2V_zeroB)-pm2V[0])/Bm2V_mask

plotP = np.linspace(-tol_final, tol_final)
lowM0 = plotP*p0V[1]
lowM2 = plotP*p2V[1]
lowMm2 = plotP*pm2V[1]

cf = 5E7
plt.scatter(B0V_mask, magBase0V*cf, s=5, linewidth = .5,  zorder=1, marker = '.', edgecolors = 'k', color = 'black', label = '0 V')
plt.plot(plotP, lowM0*cf, color = 'black')
plt.scatter(B2V_mask, magBase2V*cf, s=5, linewidth = .5,  zorder=1, marker = '.', edgecolors = 'red', color = 'red', label = '2 V')
plt.plot(plotP, lowM2*cf, color = 'red')
plt.scatter(Bm2V_mask, magBasem2V*cf, s=5, linewidth = .5,  zorder=1, marker = '.', edgecolors = 'green', color = 'green', label = '-2 V')
plt.plot(plotP, lowMm2*cf, color = 'green')
plt.axvline(-tol_final, alpha=.5)
plt.axvline(tol_final, alpha=.5)
plt.grid()
plt.xlabel('B (T)')
plt.ylabel(r'm ($\mu_{B})$')
#plt.ylim(-2, 2)
#plt.xlim(-.5, .5)
plt.legend(loc='best', prop={'size':10}).set_title("Gate Voltage", prop = {'size':10})
plt.savefig('/home/sam/Documents/DriftTests/SweepBs_mag.pdf', bbox_inches="tight")
plt.show()