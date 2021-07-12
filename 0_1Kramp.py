#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 14:03:17 2021

@author: tiffanypaul
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.optimize import curve_fit
plt.style.use('default')
plt.rcParams.update({'font.size': 18})

def Bfunc(x, Bstart, Bend, indStart, indEnd):
    numer = (Bend - Bstart)*(x - indStart)/(indEnd - indStart) + Bstart
    return numer

def BfuncR(x, Bstart, indStart, rate):
    numer = (x - indStart)*rate + Bstart
    return numer

def makeBs(index, BList, Bstart, Bend, indStart, indEnd):
    loc0 = np.where(index == indStart)[0][0]
    if(indEnd + 1 < np.amax(index)):
        loc1 = np.where(index == (indEnd+1))[0][0]
    else:
        loc1 = np.where(index == (indEnd))[0][0]
    indSec = index[loc0:loc1]
    BList[loc0:loc1] = Bfunc(indSec, Bstart, Bend, indStart, indEnd)
    return BList

def cutCrash(index, BList, Bstart, indStart, rate):
    loc0 = np.where(index == indStart)[0][0]
    indSec = index[loc0:]
    BList[loc0:] = BfuncR(indSec, Bstart, indStart, rate)
    return BList

def maskOff(mask, index, BList, f0):
    return index[mask], BList[mask], f0[mask]

def prefilter(index, f0, v, Q, dc):
    keep = (f0 > 1000)
    index = index[keep]
    f0 = f0[keep]
    v = v[keep]
    dc = dc[keep]
    Q = Q[keep]
    return index, f0, v, Q, dc

startcut = 0
endcut =-1

pt4 = '/Users/tiffanypaul/Desktop/KGB/Group meeting 7_2021/GrapheneDilRamps-main/7_6_21_cooldown3.txt'
pt4l = '/Users/tiffanypaul/Desktop/KGB/Group meeting 7_2021/GrapheneDilRamps-main/7_6_21_cooldown3L.txt'

F4 = np.loadtxt(pt4, delimiter="\t")
F4L = np.loadtxt(pt4l, delimiter="\t")
index4 = np.array(F4[startcut:endcut, 3])
f04 = np.array(F4[startcut:endcut, 2])
v4 = np.array(F4[startcut:endcut, 4])
Q4 = np.array(np.abs(F4L[startcut:endcut, 3]))
dc4 = np.array(F4[startcut:endcut, 0])

#to subtract off the jump from run 619-620
#jump4 = f04[620]-f04[619]
#f04[620:-1] = f04[620:-1]-jump4

#plt.plot(index4, Q4)
#plt.show()
print(len(index4))

index4, f04, v4, Q4, dc4 =  prefilter(index4, f04, v4, Q4, dc4)
BList4 = -1000*np.ones(len(index4))

BList4 = makeBs(index4, BList4, 0, 1, 40, 90)
BList4 = makeBs(index4, BList4, 1, 4, 90, 190)
BList4 = makeBs(index4, BList4, 4, 7.4, 190, 274)

mask4 = (BList4 > -999)
index4, BList4, f04 = maskOff(mask4, index4, BList4, f04)
tol = .0002
mask0 = (np.abs(BList4) <= tol)
print(BList4[mask0])
print(f04[mask0])
f04_zeroB = f04[mask0]

plt.scatter(BList4, f04-np.mean(f04), s=10, linewidth = .5,  zorder=1, marker = '.', edgecolors = 'k', color = 'black')
plt.grid()
plt.xlabel('B (T)')
plt.ylabel('f0 (Hz)')
plt.show()

plt.scatter(np.abs(BList4), f04-np.mean(f04), s=10, linewidth = .5,  zorder=1, marker = '.', edgecolors = 'k', color = 'black')
plt.grid()
plt.xlabel('|B| (T)')
plt.ylabel('f0 (Hz)')
plt.show()

#plt.scatter(BList4, (f04)/BList4, s=10, linewidth = .5,  zorder=1, marker = '.', edgecolors = 'k', color = 'black')
#plt.grid()
#plt.xlabel('B (T)')
#plt.ylabel('f0/B (Hz/T)')
#plt.ylim(-30E5, 30E5)
#plt.xlim(-.5, .5)
##plt.savefig('G:\Tiffany\plot.pdf', bbox_inches="tight")
#plt.show()

plt.scatter(BList4, (f04-np.mean(f04_zeroB))/BList4, s=10, linewidth = .5,  zorder=1, marker = '.', edgecolors = 'k', color = 'black')
plt.grid()
plt.xlabel('B (T)')
plt.ylabel(r'$\Delta f_0$/B (Hz/T)')
plt.ylim(-2, 2)
#plt.xlim(-.5, .5)
#plt.savefig('G:\Graphene_cooldown_6_2021\plot2V.pdf', bbox_inches="tight")
plt.show()

wrA = np.vstack((BList4, f04)).T
np.savetxt('BRamp0_1K_7_6_21.txt', wrA, delimiter='\t')

