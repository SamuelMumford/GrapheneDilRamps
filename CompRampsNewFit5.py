#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 14:24:16 2021

@author: tiffanypaul
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
plt.rcParams.update({'font.size': 18})


def prefilter(f0, b):
    keep = (f0 > 1000)
    f0 = f0[keep]
    b = b[keep]
    keep2 = np.insert((np.diff(f0) < .05), 0, True)
    f0 = f0[keep2]
    b = b[keep2]
    return f0, b

pt0V = '/Users/tiffanypaul/Desktop/KGB/Group meeting 7_2021/GrapheneDilRamps-main/BRamp0V_6_21_slow.txt'
pt2V = '/Users/tiffanypaul/Desktop/KGB/Group meeting 7_2021/GrapheneDilRamps-main/BRamp2V_6_21_slow.txt'
ptm2V = '/Users/tiffanypaul/Desktop/KGB/Group meeting 7_2021/GrapheneDilRamps-main/BRamp_min2V_6_7_21.txt'
pt2K = '/Users/tiffanypaul/Desktop/KGB/Group meeting 7_2021/GrapheneDilRamps-main/BRamp_2K_6_29_21.txt'
pt01K = '/Users/tiffanypaul/Desktop/KGB/Group meeting 7_2021/GrapheneDilRamps-main/BRamp0_1K_7_6_21.txt'


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

F2K = np.loadtxt(pt2K, delimiter="\t")
B2K = np.array(F2K[:, 0])
f2K = np.array(F2K[:, 1])
f2K, B2K = prefilter(f2K, B2K)

F01K = np.loadtxt(pt01K, delimiter="\t")
B01K = np.array(F01K[:, 0])
f01K = np.array(F01K[:, 1])
f01K, B01K = prefilter(f01K, B01K)

plt.scatter(B0V, f0V, s=10, linewidth = 1,  zorder=1, marker = '.', c = 'black', label = '1.1 K, O V', edgecolors = 'darkgreen')
plt.scatter(B2V, f2V, s=10, linewidth = 1,  zorder=1, marker = '.', c = 'black', label = '1.1 K, 2 V', edgecolors = 'coral')
plt.scatter(Bm2V, fm2V, s=10, linewidth = 1,  zorder=1, marker = '.', c = 'black', label = '1.1 K, -2 V', edgecolors = 'yellowgreen')
plt.scatter(B2K, f2K, s=10, linewidth = 1,  zorder=1, marker = '.', c = 'black', label = '2 K, 0 V', edgecolors = 'darkred')
plt.scatter(B01K, f01K, s=10, linewidth = 1,  zorder=1, marker = '.', c = 'black', label = '0.1 K, 0 V', edgecolors = 'black')
plt.grid()
plt.xlabel(r'$B$ (T)')
plt.ylabel(r'$f_0$ (Hz)')
plt.legend(prop={'size':10}, loc = 'lower left')
plt.savefig('rawf0vsB.pdf', bbox_inches="tight")
plt.show()

tol = .0005
m0V0B = (np.abs(B0V) < tol)
m2V0B = (np.abs(B2V) < tol)
mm2V0B = (np.abs(Bm2V) < tol)
f0V_zeroB = f0V[m0V0B]
f2V_zeroB = f2V[m2V0B]
fm2V_zeroB = fm2V[mm2V0B]

tolly = .1
m2k0B = (np.abs(B2K) < tolly)
f2K_zeroB = f2K[m2k0B]
m01k0B = (np.abs(B01K) < tolly)
f01K_zeroB = f01K[m01k0B]
#print('f01K_zeroB')
#print(f01K_zeroB)

#Fm2V = np.loadtxt(ptm2V, delimiter="\t")
#Bm2V = np.array(Fm2V[:, 0])
#fm2V = np.array(Fm2V[:, 1])
#fm2V, B2V = prefilter(fm2V, Bm2V)

tol_final =1.1
sm = .1

def cubicC(x, a, b):
    return a + b*x**2

def makeZeroB(Bs, f0s, tol, smWin):
    Bwin = (np.abs(Bs - tol) < smWin)
    Bfit = Bs[Bwin]
    f0fit = f0s[Bwin]
    popt, _ = curve_fit(cubicC, Bfit, f0fit, p0=[0,0])
    return popt
    
p0V = makeZeroB(B0V, f0V-np.mean(f0V_zeroB), tol_final, sm)
#print('p0V')
#print(p0V)
p2V = makeZeroB(B2V, f2V-np.mean(f2V_zeroB), tol_final, sm)
pm2V = makeZeroB(Bm2V, fm2V-np.mean(fm2V_zeroB), tol_final, sm)
p2K = makeZeroB(B2K, f2K-np.mean(f2K_zeroB), tol_final, sm)
p01K = makeZeroB(B01K, f01K-np.mean(f01K_zeroB), tol_final, sm)
#print('p01K')
#print(p01K)

plotP = np.linspace(0, tol_final)
plotf0 = cubicC(plotP, 0, p0V[1])
plotf2 = cubicC(plotP, 0, p2V[1])
plotfm2 = cubicC(plotP, 0, pm2V[1])
plotf2K = cubicC(plotP, 0, p2K[1])
plotf01K = cubicC(plotP, 0, p01K[1])

maski = (np.abs(B0V) > tol_final)
maskii = (np.abs(B2V) > tol_final)
maskiii = (np.abs(Bm2V) > tol_final)
maski2k = (np.abs(B2K) > tol_final)
maski01k = (np.abs(B01K) > tol_final)


al0 = .1


offy = 0
offy2 = 0
Bdata = np.array(B0V[maski])
fdata = np.array((f0V-np.mean(f0V_zeroB)-p0V[0])[maski] - offy)
fdata = np.array((f0V-np.mean(f0V_zeroB))[maski] - offy)

#What is p0V{0}? So what is being plotted with the fit?

def FF(x, a, b, d, e, f):
    c=1.22
    q1 = (2*a + 1)/(a)
    l1 = 1/(2*a)
    polyP = b*x*(q1/np.tanh(q1*c*a*x) - l1/(np.tanh(c*a*x*l1)))
    polyP += d*np.sign(x)*x
    polyP += e + f*np.sign(x)
    #polyP += g*np.abs(x)
    return polyP

ff = FF(Bdata,  2, 2, -2, 1.5, 0)
p0 = [2, 2, -2, 1.5, 0]
popt, pcov = curve_fit(FF, Bdata, fdata, p0=p0, maxfev=5000)
ff = FF(Bdata, popt[0],popt[1], popt[2], popt[3], popt[4])
offpart = FF(Bdata, popt[0], 0, popt[2], popt[3], popt[4])
print(popt)
print(np.sqrt(np.diag(pcov)))

sys0 = FF(B0V, popt[0],0, 0, 0, popt[4])

#plt.plot(Bdata, (ff-popt[3]), color = 'black', label = 'Fit')
plt.plot(Bdata, (ff), color = 'black', label = 'Fit')
#plt.plot(Bdata, (fdata-popt[3]), color = 'darkgreen', label = '1.1 K, 0 V Data')
plt.plot(Bdata, (fdata), color = 'darkgreen', label = '1.1 K, 0 V Data')
plt.ylabel(r'$\Delta f_{0}$ (Hz)')
plt.xlabel('B (T)')
plt.legend(prop={'size':10}, loc = 'best')
#plt.savefig('f0_0Vfit.pdf', bbox_inches="tight")
plt.show()

cf = 5E7
#plt.plot(Bdata, (ff-popt[3])/Bdata*cf, color = 'black', label = 'Fit')
plt.plot(Bdata, (ff)/Bdata*cf, color = 'black', label = 'Fit')
#plt.plot(Bdata, (fdata-popt[3])/Bdata*cf, color = 'darkgreen', label = '1.1 K, 0 V Data')
plt.plot(Bdata, (fdata)/Bdata*cf, color = 'darkgreen', label = '1.1 K, 0 V Data')
plt.ylabel(r'$m~(\mu_B)$')
plt.xlabel('B (T)')
plt.legend(prop={'size':10}, loc = 'best')
#plt.savefig('m_0Vfit.pdf', bbox_inches="tight")
plt.show()

bd0 = Bdata
fd0 = (fdata - ff)/Bdata

Bdata = np.array(Bm2V[maskiii])
fdata = np.array((fm2V-np.mean(fm2V_zeroB)-pm2V[0])[maskiii])
fdata = np.array((fm2V-np.mean(fm2V_zeroB))[maskiii])


ff = FF(Bdata,  2, 2, -2, 1.5, 0)
p0 = [2, 2, -2, 1.5, 0]
popt, pcov = curve_fit(FF, Bdata, fdata, p0=p0, maxfev=5000)
ff = FF(Bdata, popt[0],popt[1], popt[2], popt[3], popt[4])
print(popt)
print(np.sqrt(np.diag(pcov)))

sysm2 = FF(Bm2V, popt[0],0, 0, 0, popt[4])

bdm2 = Bdata
fdm2 = (fdata - ff)/Bdata

Bdata = np.array(B2V[maskii])
fdata = np.array((f2V-np.mean(f2V_zeroB)-p2V[0])[maskii])
fdata = np.array((f2V-np.mean(f2V_zeroB))[maskii])


ff = FF(Bdata,  2, 2, -2, 1.5, 0)
p0 = [2, 2, -2, 1.5, 0]
popt, pcov = curve_fit(FF, Bdata, fdata, p0=p0, maxfev=5000)
ff = FF(Bdata, popt[0],popt[1], popt[2], popt[3], popt[4])
print(popt)
print(np.sqrt(np.diag(pcov)))

sys2 = FF(B2V, popt[0],0, 0, 0, popt[4])

bd2 = Bdata
fd2 = (fdata - ff)/Bdata

#2K
Bdata = np.array(B2K[maski2k])
fdata = np.array((f2K-np.mean(f2K_zeroB)-p2K[0])[maski2k])
fdata = np.array((f2K-np.mean(f2K_zeroB))[maski2k])


ff = FF(Bdata,  2, 2, -2, 1.5, 0)
p0 = [2, 2, -2, 1.5, 0]
popt, pcov = curve_fit(FF, Bdata, fdata, p0=p0, maxfev=5000)
ff = FF(Bdata, popt[0],popt[1], popt[2], popt[3], popt[4])
print(popt)
print(np.sqrt(np.diag(pcov)))

bd2K = Bdata
fd2K = (fdata - ff)/Bdata

#sys2k = FF(B2K, popt[0],0, 0, popt[3], popt[4])
sys2k = FF(B2K, popt[0],0, 0, 0, popt[4])
#print('sys2k')
#print(sys2k)

#0.1K
def FF01(x, a, b, d, e, f):
    c=1.22
    f=0
    q1 = (2*a + 1)/(a)
    l1 = 1/(2*a)
    polyP = b*x*(q1/np.tanh(q1*c*a*x) - l1/(np.tanh(c*a*x*l1)))
    polyP += d*np.sign(x)*x
    polyP += e + f*np.sign(x)
    #polyP += g*np.abs(x)
    return polyP

Bdata = np.array(B01K[maski01k])
fdata = np.array((f01K-np.mean(f01K_zeroB)-p01K[0])[maski01k])
fdata = np.array((f01K-np.mean(f01K_zeroB))[maski01k])


ff = FF01(Bdata,  2, 2, -2, 1.5, 0)
p0 = [2, 2, -2, 1.5, 0]
popt, pcov = curve_fit(FF01, Bdata, fdata, p0=p0, maxfev=5000)
ff = FF01(Bdata, popt[0],popt[1], popt[2], popt[3], popt[4])
print(popt)
print(np.sqrt(np.diag(pcov)))

bd01K = Bdata
fd01K = (fdata - ff)/Bdata

plt.plot(Bdata, (ff), color = 'black', label = 'Fit')
#plt.plot(Bdata, (ff), color = 'black', label = 'Fit')
plt.plot(Bdata, (fdata), color = 'darkgreen', label = '0.1 K, 0 V Data')
#plt.plot(Bdata, (fdata), color = 'darkgreen', label = '1.1 K, 0 V Data')
plt.ylabel(r'$\Delta f_{0}$ (Hz)')
plt.xlabel('B (T)')
plt.legend(prop={'size':10}, loc = 'best')
#plt.savefig('f0_0Vfit.pdf', bbox_inches="tight")
plt.show()

cf = 5E7
plt.plot(Bdata, (ff)/Bdata*cf, color = 'black', label = 'Fit')
plt.plot(Bdata, (fdata)/Bdata*cf, color = 'darkgreen', label = '0.1 K, 0 V Data')
plt.ylabel(r'$m~(\mu_B)$')
plt.xlabel('B (T)')
plt.legend(prop={'size':10}, loc = 'best')
#plt.savefig('m_0Vfit.pdf', bbox_inches="tight")
plt.show()

#sys2k = FF(B2K, popt[0],0, 0, popt[3], popt[4])
sys01k = FF(B01K, popt[0],0, 0, 0, popt[4])
#print('sys01k')
#print(sys01k)

plt.scatter(np.abs(B0V), f0V-np.mean(f0V_zeroB)-sys0, s=10, linewidth = .5,  zorder=1, marker = '.', c = 'darkgreen', label = '1.1 K, O V')
plt.scatter(np.abs(B2V), f2V-np.mean(f2V_zeroB)-sys2, s=10, linewidth = .5,  zorder=1, marker = '.', c = 'coral', label = '1.1 K, 2 V')
plt.scatter(np.abs(Bm2V), fm2V-np.mean(fm2V_zeroB)-sysm2, s=10, linewidth = .5,  zorder=1, marker = '.', c = 'yellowgreen', label = '1.1 K, -2 V')
plt.scatter(np.abs(B2K), f2K-np.mean(f2K_zeroB)-sys2k, s=10, linewidth = .5,  zorder=1, marker = '.', c = 'darkred', label = '2 K, 0 V')
plt.scatter(np.abs(B01K), f01K-np.mean(f01K_zeroB)-sys01k, s=10, linewidth = .5,  zorder=1, marker = '.', c = 'black', label = '0.1 K, 0 V')
plt.axvline(tol_final, alpha=.5)
plt.xlabel('|B| (T)')
plt.ylabel(r'$\Delta f_{0}$ (Hz)')
plt.legend(prop={'size':10})
#plt.savefig('df0vsabsB.pdf', bbox_inches="tight")
plt.show()

plt.scatter(np.abs(B0V)[maski], (f0V-np.mean(f0V_zeroB)-sys0)[maski], s=10, linewidth = .5,  zorder=1, marker = '.', c = 'darkgreen', label = '1.1 K, O V')
plt.scatter(np.abs(B2V)[maskii], (f2V-np.mean(f2V_zeroB)-sys2)[maskii], s=10, linewidth = .5,  zorder=1, marker = '.', c = 'coral', label = '1.1 K, 2 V')
plt.scatter(np.abs(Bm2V)[maskiii], (fm2V-np.mean(fm2V_zeroB)-sysm2)[maskiii], s=10, linewidth = .5,  zorder=1, marker = '.', c = 'yellowgreen', label = '1.1 K, -2 V')
plt.scatter(np.abs(B2K)[maski2k], (f2K-np.mean(f2K_zeroB)-sys2k)[maski2k], s=10, linewidth = .5,  zorder=1, marker = '.', c = 'darkred', label = '2 K, 0 V')
plt.scatter(np.abs(B01K)[maski01k], (f01K-np.mean(f01K_zeroB)-sys01k)[maski01k], s=10, linewidth = .5,  zorder=1, marker = '.', c = 'black', label = '0.1 K, 0 V')
plt.axvline(tol_final, alpha=.5)
plt.xlabel('|B| (T)')
plt.ylabel(r'$\Delta f_{0}$ (Hz)')
plt.legend()
plt.show()

mag0V = ((f0V-np.mean(f0V_zeroB)-sys0)[maski])/(B0V[maski])
mag2V = ((f2V-np.mean(f2V_zeroB)-sys2)[maskii])/(B2V[maskii])
magm2V = ((fm2V-np.mean(fm2V_zeroB)-sysm2)[maskiii])/(Bm2V[maskiii])
mag2K = ((f2K-np.mean(f2K_zeroB)-sys2k)[maski2k])/(B2K[maski2k])
mag01K = ((f01K-np.mean(f01K_zeroB)-sys01k)[maski01k])/(B01K[maski01k])

cf = 5E7

plt.scatter((B0V)[maski], mag0V*cf, s=10, linewidth = .5,  zorder=1, marker = '.', c = 'darkgreen', label = '1.1 K, O V')
#plt.plot(plotP, plotf0, color = 'black')
plt.scatter((B2V)[maskii], mag2V*cf, s=10, linewidth = .5,  zorder=1, marker = '.', c = 'coral', label = '1.1 K, 2 V')
#plt.plot(plotP, plotf2, color = 'black')
plt.scatter((Bm2V)[maskiii], magm2V*cf, s=10, linewidth = .5,  zorder=1, marker = '.', c = 'yellowgreen', label = '1.1 K, -2 V')
#plt.plot(plotP, plotfm2, color = 'black')
plt.scatter((B2K)[maski2k], mag2K*cf, s=10, linewidth = .5,  zorder=1, marker = '.', c = 'darkred', label = '2 K, 0 V')
#plt.plot(plotP, plotf2K, color = 'black')
plt.scatter((B01K)[maski01k], mag01K*cf, s=10, linewidth = .5,  zorder=1, marker = '.', c = 'black', label = '0.1 K, 0 V')
plt.axvline(tol_final, alpha=.5)
plt.axvline(-tol_final, alpha=.5)
plt.xlabel('B (T)')
plt.ylabel(r'$m$ ($\mu_B$)')
plt.legend(prop={'size':10})
#plt.savefig('magcutoff.pdf', bbox_inches="tight")
plt.show()

#cf = 5E7
plt.scatter(1/bd0, fd0*cf, s=1, linewidth = .5,  zorder=1, color = 'darkgreen', label = '1.1 K, 0 V')
L = sorted(zip(1/bd0,fd0*cf))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'darkgreen', alpha = .25)
plt.scatter(1/bd2, fd2*cf, s=1, linewidth = .5,  zorder=1, color = 'coral', label = '1.1 K, 2 V')
L = sorted(zip(1/bd2,fd2*cf))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'coral', alpha = .25)
plt.scatter(1/bdm2, fdm2*cf, s=1, linewidth = .5,  zorder=1, color = 'yellowgreen', label = '1.1 K, -2 V')
L = sorted(zip(1/bdm2,fdm2*cf))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'yellowgreen', alpha = .25)
plt.scatter(1/bd2K, fd2K*cf, s=1, linewidth = .5,  zorder=1, color = 'darkred', label = '2 K, 0V')
L = sorted(zip(1/bd2K,fd2K*cf))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'darkred', alpha = .25)

plt.scatter(1/bd01K, fd01K*cf, s=1, linewidth = 2,  zorder=1, color = 'black', label = '0.1 K, 0V')
L = sorted(zip(1/bd01K,fd01K*cf))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'black', alpha = .25)

plt.ylabel(r'$m$  - Fit ($\mu_{B})$')
plt.xlabel(r'1/B (T$^{-1}$)')
plt.legend(prop={'size':10}, loc = 'best')
#plt.savefig('mresid.pdf', bbox_inches="tight")
plt.show()

cf = 5E7
plt.scatter(np.abs(bdm2), fdm2*cf*np.sign(bdm2), s=1, linewidth = .5,  zorder=1, color = 'green', label = '-2 V')
L = sorted(zip(np.abs(bdm2),fdm2*cf*np.sign(bdm2)))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'green', alpha = .25)
plt.scatter(np.abs(bd0), fd0*cf*np.sign(bd0), s=1, linewidth = .5,  zorder=1, color = 'black', label = '0 V')
L = sorted(zip(np.abs(bd0),fd0*cf*np.sign(bd0)))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'black', alpha = .25)
plt.scatter(np.abs(bd2), fd2*cf*np.sign(bd2), s=1, linewidth = .5,  zorder=1, color = 'red', label = '2 V')
L = sorted(zip(np.abs(bd2),fd2*cf*np.sign(bd2)))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.scatter(np.abs(bd2K), fd2K*cf*np.sign(bd2K), s=1, linewidth = .5,  zorder=1, color = 'blue', label = '2 K')
L = sorted(zip(np.abs(bd2K),fd2K*cf*np.sign(bd2K)))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'blue', alpha = .25)

plt.scatter(np.abs(bd01K), fd01K*cf*np.sign(bd01K), s=1, linewidth = 2,  zorder=1, color = 'black', label = '0.1 K')
L = sorted(zip(np.abs(bd01K),fd01K*cf*np.sign(bd01K)))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'black', alpha = .25)

plt.ylabel(r'$\mu$*Sign($B$) ($\mu_{B})$')
plt.xlabel(r'|B| (T)')
plt.legend(prop={'size':10}, loc = 'best')
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
plt.scatter(1/bd2K, fd2K*bd2K, s=1, linewidth = .5,  zorder=1, color = 'blue', label = '2 K')
L = sorted(zip(1/bd2K,fd2K*bd2K))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'blue', alpha = .25)

plt.scatter(1/bd01K, fd01K*bd01K, s=1, linewidth = 2,  zorder=1, color = 'black', label = '0.1 K')
L = sorted(zip(1/bd01K,fd01K*bd01K))#, key=operator.itemgetter(0))
new_x, new_y = zip(*L)
plt.plot(new_x, new_y, color = 'black', alpha = .25)

plt.ylabel(r'$\delta f_{0}$')
plt.xlabel(r'1/B (T$^{-1}$)')
plt.legend(prop={'size':10}, loc = 'best')
plt.show()

