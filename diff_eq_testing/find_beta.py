# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 10:24:17 2022

@author: lennak
"""

import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import sys
import keyboard
import helpers
import os

from scipy.optimize import curve_fit

wavelengths = []
freqs = []
times = []
with open("find_beta.txt", mode ='r') as file:
   
  lines = file.readlines()
  # displaying the contents of the CSV file
  for line in lines:
      line_split = line.split(" ")
      wavelength = float(line_split[0])
      time = float(line_split[1])
      wavelengths.append(wavelength)
      times.append(time)
      speed_of_light = 2.99792458e5 # nm / ps
      freqs.append(speed_of_light/wavelength)
      
def ax_b(x,a,b):
    return a*x + b
      
wavelengths = np.array(wavelengths)
times = np.array(times)
freqs = np.array(freqs)

popt, pcov = curve_fit(ax_b, freqs, times)
print(popt)
slope = popt[0]
beta = 1/(2*slope)
print(beta, "1/ps-2")

plt.figure()
plt.plot(freqs, times)
plt.show()

# 569.802 0.9
# 579.829 0.9
# 589.868 0.9
# 599.918 0.9
# 609.979 0.9

# Beta: <=360 don't use, >=569

# XPM: >=360, <=600

# BETA 162.966315 1/ps-2