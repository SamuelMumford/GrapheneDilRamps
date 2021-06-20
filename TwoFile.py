
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 14:09:11 2018

@author: KGB
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.optimize import curve_fit

def secSplit(vs, fs):
    vParts = []
    fParts = []
    fSigs = []
    val = vs[0]
    start = 0
    for i in range(1, len(vs)):
        if(vs[i] != val):
            fParts = fParts + [np.mean(fs[start:i])]
            fSigs = fSigs + [np.std(fs[start:i])/(1.0*np.sqrt(len(fs[start:i])))]
            vParts = vParts + [val]
            start = i
            val = vs[i]
    fParts = fParts + [np.mean(fs[start:])]
    fSigs = fSigs + [np.std(fs[start:])/(1.0*np.sqrt(len(fs[start:])))]
    vParts = vParts + [val]
    return vParts, fParts, fSigs

def upDwn(vP):
    up = [False]*len(vP)
    dw = [False]*len(vP)
    for i in range(1, len(vP)):
        if(vP[i] > vP[i-1]):
            up[i] = True
        else:
            dw[i] = True
    up[0] = up[1]
    dw[0] = dw[1]
    return up, dw

def split(vParts, fParts, fSigs, ce):
    vParts = np.array(vParts)
    fParts = np.array(fParts)
    fSigs = np.array(fSigs)
    up, down = upDwn(vParts[ce:])
    upF = fParts[ce:][up]
    upS = fSigs[ce:][up]
    upV = vParts[ce:][up]
    dwF = fParts[ce:][down]
    dwS = fSigs[ce:][down]
    dwV = vParts[ce:][down]
    stF = fParts[0:ce]
    stS = fSigs[0:ce]
    stV = vParts[0:ce]
    return np.array(upF), np.array(upS), np.array(upV), np.array(dwF), np.array(dwS), np.array(dwV), np.array(stF), np.array(stS), np.array(stV)

fname = '/home/sam/Documents/DriftTests/6_1_21_1Kcycle_long.txt'
fnamel = '/home/sam/Documents/DriftTests/6_1_21_1Kcycle_longL.txt'

flipname = '/home/sam/Documents/DriftTests/6_1_21_1Kcycle_longf.txt'
flipnamel = '/home/sam/Documents/DriftTests/6_1_21_1Kcycle_longfL.txt'

startcut = 0
ce = 8
ce2 = 24

Fdat = np.loadtxt(fname, delimiter="\t")
FL = np.loadtxt(fnamel, delimiter="\t")
endcut =-1
index = Fdat[startcut:endcut, 3]
f0 = Fdat[startcut:endcut, 2]
v = Fdat[startcut:endcut, 4]
Q = np.abs(FL[startcut:endcut, 3])
dc = Fdat[startcut:endcut, 0]

Fdat = np.loadtxt(flipname, delimiter="\t")
FL = np.loadtxt(flipnamel, delimiter="\t")
endcut =-1
indexl = Fdat[startcut:endcut, 3]
f0l = Fdat[startcut:endcut, 2]
vl = Fdat[startcut:endcut, 4]
Ql = np.abs(FL[startcut:endcut, 3])
dcl = Fdat[startcut:endcut, 0]

fa = np.mean(f0)
fb = np.mean(f0l)

vParts, fParts, fSigs = secSplit(v, f0)
upF, upS, upV, dwF, dwS, dwV, stF, stS, stV = split(vParts, fParts, fSigs, ce)

vPartsl, fPartsl, fSigsl = secSplit(vl, f0l)
upFl, upSl, upVl, dwFl, dwSl, dwVl, stFl, stSl, stVl = split(vPartsl, fPartsl, fSigsl, ce2)

plt.scatter(index, f0, s=10, linewidth = .5,  zorder=1, marker = 'o', edgecolors = 'k', color = 'black')
plt.xlabel('Run Index')
plt.ylabel('f0 (Hz)')
plt.show()

plt.scatter(index, dc, s=10, linewidth = .5,  zorder=1, marker = 'o', edgecolors = 'k', color = 'black')
plt.xlabel('Run Index')
plt.ylabel('DC Level (V)')
plt.show()

plt.scatter(index, Q, s=10, linewidth = .5,  zorder=1, marker = 'o', edgecolors = 'k', color = 'black')
plt.xlabel('Run Index')
plt.ylabel('Q')
#plt.ylim(7000, 11000)
plt.show()

plt.scatter(v, f0, s=10, linewidth = .5,  zorder=1, marker = 'o', edgecolors = 'k', color = 'black')
plt.xlabel('Gate Voltage (V)')
plt.ylabel('f0 (Hz)')
plt.show()

plt.errorbar(vParts, fParts, fSigs, marker = '.', linewidth = 0, elinewidth=1, color = 'black', capsize=3)
plt.xlabel('Gate Voltage (V)')
plt.ylabel('f0 (Hz)')
plt.show()

#plt.errorbar(upV, upF, upS, marker = '.', linewidth = 0, elinewidth=1, color = 'black', capsize=3, label = 'Up')
#plt.errorbar(dwV, dwF, dwS, marker = '.', linewidth = 0, elinewidth=1, color = 'red', capsize=3, label = 'Down')
plt.errorbar(stV, stF, stS, marker = '.', linewidth = 0, elinewidth=1, color = 'orange', capsize=3, label = 'Start')

#plt.errorbar(upVl, upFl, upSl, marker = '.', linewidth = 0, elinewidth=1, color = 'blue', capsize=3, label = 'Up')
#plt.errorbar(dwVl, dwFl, dwSl, marker = '.', linewidth = 0, elinewidth=1, color = 'green', capsize=3, label = 'Down')
plt.errorbar(stVl, stFl, stSl, marker = '.', linewidth = 0, elinewidth=1, color = 'purple', capsize=3, label = 'Start')

plt.xlabel('Gate Voltage (V)')
plt.ylabel('f0 (Hz)')
plt.legend(loc='best')
plt.show()

plt.errorbar(upV, upF-fa, upS, marker = '.', linewidth = 0, elinewidth=1, color = 'black', capsize=3, label = 'Up')
#plt.errorbar(dwV, dwF, dwS, marker = '.', linewidth = 0, elinewidth=1, color = 'red', capsize=3, label = 'Down')
#plt.errorbar(stV, stF, stS, marker = '.', linewidth = 0, elinewidth=1, color = 'orange', capsize=3, label = 'Start')

plt.errorbar(upVl, upFl-fb, upSl, marker = '.', linewidth = 0, elinewidth=1, color = 'blue', capsize=3, label = 'Up')
#plt.errorbar(dwVl, dwFl, dwSl, marker = '.', linewidth = 0, elinewidth=1, color = 'green', capsize=3, label = 'Down')
#plt.errorbar(stVl, stFl, stSl, marker = '.', linewidth = 0, elinewidth=1, color = 'purple', capsize=3, label = 'Start')

plt.xlabel('Gate Voltage (V)')
plt.ylabel('f0 (Hz)')
plt.legend(loc='best')
plt.show()

#plt.errorbar(upV, upF, upS, marker = '.', linewidth = 0, elinewidth=1, color = 'black', capsize=3, label = 'Up')
plt.errorbar(dwV, dwF-fa, dwS, marker = '.', linewidth = 0, elinewidth=1, color = 'red', capsize=3, label = 'Down')
#plt.errorbar(stV, stF, stS, marker = '.', linewidth = 0, elinewidth=1, color = 'orange', capsize=3, label = 'Start')

#plt.errorbar(upVl, upFl, upSl, marker = '.', linewidth = 0, elinewidth=1, color = 'blue', capsize=3, label = 'Up')
plt.errorbar(dwVl, dwFl-fb, dwSl, marker = '.', linewidth = 0, elinewidth=1, color = 'green', capsize=3, label = 'Down')
#plt.errorbar(stVl, stFl, stSl, marker = '.', linewidth = 0, elinewidth=1, color = 'purple', capsize=3, label = 'Start')

plt.xlabel('Gate Voltage (V)')
plt.ylabel('f0 (Hz)')
plt.legend(loc='best')
plt.show()

plt.errorbar(upV, upF-fa, upS, marker = '.', linewidth = 0, elinewidth=1, color = 'red', capsize=3, label = 'Up')
plt.errorbar(dwV, dwF-fa, dwS, marker = '.', linewidth = 0, elinewidth=1, color = 'black', capsize=3, label = 'Down')
plt.errorbar(stV, stF-fa, stS, marker = '.', linewidth = 0, elinewidth=1, color = 'orange', capsize=3, label = 'Start')

plt.errorbar(upVl, upFl-fb, upSl, marker = 'x', linewidth = 0, elinewidth=1, color = 'red', capsize=3, label = 'Up')
plt.errorbar(dwVl, dwFl-fb, dwSl, marker = 'x', linewidth = 0, elinewidth=1, color = 'black', capsize=3, label = 'Down')
plt.errorbar(stVl, stFl-fb, stSl, marker = 'x', linewidth = 0, elinewidth=1, color = 'orange', capsize=3, label = 'Start')

plt.xlabel('Gate Voltage (V)')
plt.ylabel(r'$\delta f_0$ (Hz)')
plt.legend(loc='best')
plt.savefig('/home/sam/Documents/DriftTests/V_switches.pdf')
plt.show()