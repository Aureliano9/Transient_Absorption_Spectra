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

# user specifies starting parameters 1/B (picosec) 3000picoseconds
# 1picosecond,10picosecond, 13 microseconds (cannot fit)
f = .5
B2=2
B3=3
B4=4

sigma = 10 #~0.1
F = lambda t, s: np.dot(np.array([[0,0,0,B4], [0,-B2,0,0], [0,B2,-B3,0],[0,0,B3,-B4]]), s)
G = lambda t, s: np.dot(np.array([[-f*(t/(sigma**2))*np.exp(-t**2/(2*sigma**2)),0,0,B4], [f*(t/(sigma**2))*np.exp(-t**2/(2*sigma**2)),-B2,0,0], [0,B2,-B3,0],[0,0,B3,-B4]]), s)

 # xdot(1) =  B(4)*x(4);
 # xdot(2) =  - B(2)*x(2);
 # xdot(3) =  B(2)*x(2)-B(3)*x(3);
 # xdot(4) =  B(3)*x(3)-B(4)*x(4);

initial = [1-f,f,0,0]
A2 = initial[1]
A3 = initial[2]+initial[1]*B2/(B2-B3)
A4 = initial[3]-initial[1]*B2*B3/((B2-B3)*(B2-B4))+A3*B3/(B3-B4)
A1 = initial[0]-A3*B4/(B3-B4)-B4*A2*B3/((B3-B2)*(B2-B4))+A4

precision = .001
t_eval = np.arange(0, 20.01, precision)
t_span = [0,25]
sol = solve_ivp(F, t_span, initial, t_eval=t_eval)
# sol2 = solve_ivp(G, t_span, initial, t_eval=t_eval)
# area = 0
# for el in sol.y.T[:, 0]:
#     area += el*precision
# print(area)

x1 = B4*A3*np.exp((-B3)*t_eval)/(B3-B4) + B4*A2*B3*np.exp(-B2*t_eval)/((B3-B2)*(B2-B4))-A4*np.exp(-B4*t_eval) + A1
x2 = A2*np.exp(-B2*t_eval)
x3 = A3*np.exp(-B3*t_eval)-A2*B2/(B2-B3)*np.exp(-B2*t_eval)
x4 = A2*B2*B3*np.exp(-B2*t_eval)/((B2-B3)*(B2-B4))-A3*B3*np.exp(-B3*t_eval)/(B3-B4)+A4*np.exp(-B4*t_eval)

plt.figure(figsize = (12, 8))
plt.plot(t_eval, x1)
plt.plot(t_eval, x2)
plt.plot(t_eval, x3)
plt.plot(t_eval, x4)
plt.plot(t_eval, sol.y.T[:, 0])
plt.plot(t_eval, sol.y.T[:, 1])
plt.plot(t_eval, sol.y.T[:, 2])
plt.plot(t_eval, sol.y.T[:, 3])

# plt.plot(t_eval, sol2.y.T[:, 0])
# plt.plot(t_eval, sol2.y.T[:, 1])
# plt.plot(t_eval, sol2.y.T[:, 2])
# plt.plot(t_eval, sol2.y.T[:, 3])
# def gaussian(x):
#     return np.exp(-x**2/sigma**2)/(sigma*sqrt(pi))
# def error_func(x):
#     integrate.quad(lambda x: gaussian(x), -float('inf'), float('inf'))
# t_eval = np.arange(-20, 20, 0.01)
# plt.plot(t_eval, gaussian(t_eval))
# integral = []
# for t in t_eval:
#     integral.append(integrate.quad(lambda x: gaussian(x), -float('inf'), t)[0])
# integral = np.array(integral)
# plt.plot(t_eval, integral)
# result = integrate.quad(lambda x: gaussian(x), -float('inf'), float('inf'))
# print(result)
# plt.plot(t_eval, special.erf(t_eval))
plt.xlabel('x')
plt.ylabel('y')
plt.show()

wavelengths = np.arange(350,700,1)
A1 = []
A2 = []
A3 = []
A4 = []
for w in wavelengths:
    A1.append(0)