# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 14:10:12 2022

@author: lennak
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import special
from scipy.integrate import solve_ivp
from math import pi, sqrt
from scipy import integrate
import helpers
from scipy import ndimage
from scipy import fftpack

# user specifies starting parameters 1/B (picosec) 3000picoseconds
# 1picosecond,10picosecond, 13 microseconds (cannot fit)
f = .5 # FIT
B2=1 # FIT
B3=.1 # FIT
B4=1/1100000 # FIT
# for gaussian
sigma= 1 # FIT -> tau1 from XPM fit (average)

# F = lambda t, s: np.dot(np.array([[0,0,0,B4], [0,-B2,0,0], [0,B2,-B3,0],[0,0,B3,-B4]]), s)
# G = lambda t, s: np.dot(np.array([[-f*(t/(sigma**2))*np.exp(-t**2/(2*sigma**2)),0,0,B4], [f*(t/(sigma**2))*np.exp(-t**2/(2*sigma**2)),-B2,0,0], [0,B2,-B3,0],[0,0,B3,-B4]]), s)

 # xdot(1) =  B(4)*x(4);
 # xdot(2) =  - B(2)*x(2);
 # xdot(3) =  B(2)*x(2)-B(3)*x(3);
 # xdot(4) =  B(3)*x(3)-B(4)*x(4);

initial = [1-f,f,0,0]

precision = .001
t_eval = np.arange(-5, 20, precision)
t_span = [-25,25]
# sol = solve_ivp(F, t_span, initial, t_eval=t_eval)
# sol2 = solve_ivp(G, t_span, initial, t_eval=t_eval)
# area = 0
# for el in sol.y.T[:, 0]:
#     area += el*precision
# print(area)

def solve_diffeq(t, initial):
    C2 = initial[1]
    C3 = initial[2]+initial[1]*B2/(B2-B3)
    C4 = initial[3]-initial[1]*B2*B3/((B2-B3)*(B2-B4))+C3*B3/(B3-B4)
    C1 = initial[0]-C3*B4/(B3-B4)-B4*C2*B3/((B3-B2)*(B2-B4))+C4
    x1 = B4*C3*np.exp((-B3)*t)/(B3-B4) + B4*C2*B3*np.exp(-B2*t)/((B3-B2)*(B2-B4))-C4*np.exp(-B4*t) + C1
    x2 = C2*np.exp(-B2*t)
    x3 = C3*np.exp(-B3*t)-C2*B2/(B2-B3)*np.exp(-B2*t)
    x4 = C2*B2*B3*np.exp(-B2*t)/((B2-B3)*(B2-B4))-C3*B3*np.exp(-B3*t)/(B3-B4)+C4*np.exp(-B4*t)
    return x1, x2, x3, x4


x1, x2, x3, x4 = solve_diffeq(t_eval, initial)
g = helpers.gaussian(np.arange(-len(t_eval)*precision,len(t_eval)*precision,precision),sigma)
x1_blurred = ndimage.convolve(x1,g, mode='constant', cval=0.0)
x2_blurred = ndimage.convolve(x2,g, mode='constant', cval=0.0)
x3_blurred = ndimage.convolve(x3,g, mode='constant', cval=0.0)
x4_blurred = ndimage.convolve(x4,g, mode='constant', cval=0.0)

# B4=0 # versions: remove, set to constant from literature, fit parameter
# x12, x22, x33, x44 = solve_diffeq(t_eval, initial)

plt.figure(figsize = (12, 8))
plt.plot(t_eval, x1)
plt.plot(t_eval, x2)
plt.plot(t_eval, x3)
plt.plot(t_eval, x4)
plt.plot(t_eval, x1_blurred)
plt.plot(t_eval, x2_blurred)
plt.plot(t_eval, x3_blurred)
plt.plot(t_eval, x4_blurred)

# PARAMETERS
A1 = 1 # FIT
A2 = 1 # FIT
A3 = 1 # FIT
A4 = 1 # FIT

fitted = A1*x1_blurred + A2*x2_blurred + A3*x3_blurred + A4*x4_blurred

# plt.plot(t_eval, x12)
# plt.plot(t_eval, x22)
# plt.plot(t_eval, x33)
# plt.plot(t_eval, x44)

# plt.plot(t_eval, special.erf(t_eval))
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()

# wavelengths = np.arange(350,700,1)
# A1 = []
# A2 = []
# A3 = []
# A4 = []
# for w in wavelengths:
#     A1.append(0)