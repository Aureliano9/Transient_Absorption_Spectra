# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 09:58:12 2022

@author: lennak
"""

import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import sys
import keyboard
from scipy.optimize import curve_fit

cdict = {'red': ((0.0, 0.0, 0.0),
                 (0.1, 0.5, 0.5),
                 (0.2, 0.0, 0.0),
                 (0.4, 0.2, 0.2),
                 (0.6, 0.0, 0.0),
                 (0.8, 1.0, 1.0),
                 (1.0, 1.0, 1.0)),
        'green':((0.0, 0.0, 0.0),
                 (0.1, 0.0, 0.0),
                 (0.2, 0.0, 0.0),
                 (0.4, 1.0, 1.0),
                 (0.6, 1.0, 1.0),
                 (0.8, 1.0, 1.0),
                 (1.0, 0.0, 0.0)),
        'blue': ((0.0, 0.0, 0.0),
                 (0.1, 0.5, 0.5),
                 (0.2, 1.0, 1.0),
                 (0.4, 1.0, 1.0),
                 (0.6, 0.0, 0.0),
                 (0.8, 0.0, 0.0),
                 (1.0, 0.0, 0.0))}

my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

def plot_color(w, t, dA, current_w, current_t, w_bounds=None, t_bounds=None):
    plt.figure()
    plt.pcolor(w, t, dA, vmin=-.001, vmax=.003, cmap=my_cmap)
    plt.xlabel("Wavelength")
    plt.ylabel("Time")
    plt.colorbar()
    plt.plot([current_w,current_w],[t.min(),t.max()])
    plt.plot([w.min(),w.max()],[current_t, current_t])
    if w_bounds==None:
        w_bounds = [np.min(w),np.max(w)]
    if t_bounds==None:
        t_bounds = [np.min(t),np.max(t)]
    plt.axis([w_bounds[0], w_bounds[1], t_bounds[0], t_bounds[1]])
    plt.show()

def plot_crosssection(w,t,dA,cut_index,wavelength_flag,bounds=None):
    plt.figure()
    if wavelength_flag:
        plt.plot(w, dA[cut_index, :])
        plt.xlabel("Wavelength")
        plt.title("Time = " + str(t[cut_index]))
    else:
        plt.plot(t, dA[:, cut_index])
        plt.xlabel("Time")
        plt.title("Wavelength = " + str(w[cut_index]))
    if bounds!=None:
        plt.xlim(bounds[0],bounds[1])
    plt.ylabel("\Delta A")
    plt.show()

# return matrix from file
def readFile(filename):
    wavelengths = []
    times = []
    signal = []
     
    # opening the CSV file
    with open(filename, mode ='r')as file:
       
      # reading the CSV file
      csvFile = csv.reader(file)
     
      # displaying the contents of the CSV file
      row_counter = 0
      for lines in csvFile:
          column_counter = 0
          
          #break condition
          if len(lines)==1:
              break
          
          for line in lines:
              if row_counter == 0:
                  if column_counter>0:
                      times.append(float(line))
              else:
                  if column_counter==0:
                      wavelengths.append(float(line))
                      signal_row = []
                  else:
                      signal_row.append(float(line))
              column_counter += 1
          if row_counter>0:
              signal.append(signal_row)
          row_counter += 1
          
    wavelengths = np.array(wavelengths)
    times = np.array(times)
    signal = np.array(signal).transpose()
    
    return (wavelengths,times,signal)

def find_index(arr, value):
    return np.argmin(np.abs(arr-value))

def heaviside(x, shift, magnitude):
    return magnitude * np.heaviside(x,shift)

def fit_heaviside(x,y):
    cleaned_x = x[~np.isnan(y)]
    cleaned_y = y[~np.isnan(y)]
    
    popt, pcov = curve_fit(heaviside, cleaned_x, cleaned_y)
    best_shift = popt[0]
    best_magnitude = popt[1]
    
    return best_shift

def quadratic(x, a, b, c):
    return a*x**2 + b*x + c

def fit_quadratic(x,y):
    cleaned_x = x[~np.isnan(y)]
    cleaned_y = y[~np.isnan(y)]
    
    popt, pcov = curve_fit(quadratic, cleaned_x, cleaned_y)
    a = popt[0]
    b = popt[1]
    c = popt[2]
    
    return (a,b,c)