
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
from scipy.signal import savgol_filter

def secSplit(vs, fs, inds):
    vParts = []
    fParts = []
    fSigs = []
    indi = []
    val = vs[0]
    start = 0
    for i in range(1, len(vs)):
        if(vs[i] != val):
            fParts = fParts + [np.mean(fs[start:i])]
            fSigs = fSigs + [np.std(fs[start:i])/(1.0*np.sqrt(len(fs[start:i])))]
            vParts = vParts + [val]
            indi = indi + [np.mean(inds[start:i])]
            start = i
            val = vs[i]
    fParts = fParts + [np.mean(fs[start:])]
    fSigs = fSigs + [np.std(fs[start:])/(1.0*np.sqrt(len(fs[start:])))]
    vParts = vParts + [val]
    indi = indi + [np.mean(inds[start:])]
    return np.array(vParts), np.array(fParts), np.array(fSigs), np.array(indi)

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

def prefilter(index, f0, v, Q, dc):
    keep = (f0 > 1000)
    index = index[keep]
    f0 = f0[keep]
    v = v[keep]
    dc = dc[keep]
    Q = Q[keep]
    return index, f0, v, Q, dc

fname = '/home/sam/Documents/DriftTests/6_9_21_SRS830Test3.txt'
fnamel = '/home/sam/Documents/DriftTests/6_9_21_SRS830Test3L.txt'

startcut = 0
endcut = -1
ce = 0

Fdat = np.loadtxt(fname, delimiter="\t")
FL = np.loadtxt(fnamel, delimiter="\t")
index = np.array(Fdat[startcut:endcut, 3])
f0 = np.array(Fdat[startcut:endcut, 2])
v = np.array(Fdat[startcut:endcut, 4])
Q = np.array(np.abs(FL[startcut:endcut, 3]))
dc = np.array(Fdat[startcut:endcut, 0])

index, f0, v, Q, dc =  prefilter(index, f0, v, Q, dc)

vParts, fParts, fSigs, indi = secSplit(v, f0, index)

colors = cm.jet(np.linspace(0, 1, len(vParts)))

plt.scatter(index, dc, s=10, linewidth = .5,  zorder=1, marker = 'o', edgecolors = 'k', color = 'black')
plt.xlabel('Run Index')
plt.ylabel('DC Level (V)')
plt.show()

plt.scatter(index, Q, s=10, linewidth = .5,  zorder=1, marker = 'o', edgecolors = 'k', color = 'black')
plt.xlabel('Run Index')
plt.ylabel('Q')
#plt.savefig('/home/sam/Documents/DriftTests/Qstab.pdf')
plt.show()

win = 30
smbool = (len(f0) > 2*win+1)

plt.scatter(index, f0-np.mean(f0), s=40, linewidth = .5,  zorder=1, marker = '.', edgecolors = 'k', color = 'black')
if(smbool):
    smf0 = savgol_filter(f0, 2*win + 1, 2)
    plt.plot(index, smf0-np.mean(f0))
plt.grid()
plt.xlabel('Run Index')
plt.ylabel('f0 - <f0> (Hz)')
plt.savefig('/home/sam/Documents/DriftTests/lowLPcurt.pdf')
plt.show()

#plt.scatter(v, f0, s=10, linewidth = .5,  zorder=1, marker = 'o', edgecolors = 'k', color = 'black')
#plt.xlabel('Gate Voltage (V)')
#plt.ylabel('f0 (Hz)')
#plt.show()

#vList = []
#for i in range(0, len(vParts)):
#    vt = vParts[i]
#    if(vt not in vList):
#        booly = (vParts == vt)
#        fv = fParts[booly]
#        fs = fSigs[booly]
#        vList = vList + [vt]
#        inds = np.arange(0, len(fv))
#        plt.errorbar(inds, fv, fs, marker = '.', linewidth = 0, elinewidth=1, capsize=3, label = str(vt))
#        plt.xlabel('Run Index')
#        plt.ylabel('DC Level (V)')
#        plt.legend()
#        plt.show()
#
#vList = []        
#for i in range(0, len(vParts)):
#    vt = vParts[i]
#    if(vt not in vList):
#        booly = (vParts == vt)
#        fv = fParts[booly]
#        fs = fSigs[booly]
#        indp = indi[booly]
#        vList = vList + [vt]
#        inds = np.arange(0, len(fv))
#        plt.errorbar(indp, fv, fs, marker = '.', linewidth = 0, elinewidth=1, capsize=3, label = str(vt))
#plt.xlabel('Run Index')
#plt.ylabel('DC Level (V)')
#plt.legend()
#plt.show()
