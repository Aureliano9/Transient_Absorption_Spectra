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

def plot_color(w, t, dA, current_w, current_t, w_bounds=None, t_bounds=None, c_bounds=None):
    plt.figure()
    if c_bounds==None or c_bounds[0]==None or c_bounds[1]==None:
        # vmin=.001, vmax=.004
        avg = np.nanmean(dA)
        std = np.nanstd(dA)
        color_width_std = .35
        c_bounds = [avg - color_width_std * std, avg + color_width_std * std]
    plt.pcolor(w, t, dA, vmin=c_bounds[0], vmax=c_bounds[1], cmap=my_cmap)
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

def plot_crosssection(w,t,dA,cut_value,wavelength_flag,bounds=None):
    plt.figure()
    if wavelength_flag:
        cut_index = find_index(t, cut_value)
        plt.plot(w, dA[cut_index, :])
        plt.xlabel("Wavelength")
        plt.title("Time = " + str(t[cut_index]))
    else:
        cut_index = find_index(w, cut_value)
        plt.plot(t, dA[:, cut_index])
        plt.xlabel("Time")
        plt.title("Wavelength = " + str(w[cut_index]))
    if bounds!=None:
        plt.xlim(bounds[0],bounds[1])
    plt.ylabel("\Delta A")
    plt.show()
    
def avg(a,b):
    return (a+b)/2.

# return matrix from file
def read_file(filename):
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
    output_index = 0 #default to zero shift
    best_loss = float('inf')
    for i in range(1,len(x)):
        # calculate loss if jump is here (assume jump to MEAN value)
        loss = np.sum((y[:i])**2) + np.sum((y[i:]-np.nanmean(y[i:]))**2)
        if loss<best_loss:
            output_index = i
            best_loss = loss
    return output_index

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

def chirp_correction(times,wavelengths,delta_A):
    shifts = []
    for i in range(len(wavelengths)):   
        shift = fit_heaviside(times,delta_A[:,i])
        shifts.append(shift)
        # if i==len(wavelengths)//2:
        #     plt.figure()
        #     plt.plot(times,delta_A[:,i])
        #     plt.title(times[shift])
        #     plt.xlim(-10,10)
        #     plt.show()
        
    shifts = np.array(shifts)
    
    # display shifts for user reference
    plt.figure()
    plt.plot(wavelengths,shifts)
    plt.show()
    
    print("Determing range to fit for chirp correction")
    chirp_min = int(input("Min for range: "))
    chirp_max = int(input("Max for range: "))
    chirp_min_index = find_index(wavelengths, chirp_min)
    chirp_max_index = find_index(wavelengths, chirp_max)
    
    a,b,c = fit_quadratic(wavelengths[chirp_min_index:chirp_max_index], shifts[chirp_min_index:chirp_max_index])
    
    fitted = np.zeros(len(wavelengths))
    fitted[chirp_min_index:chirp_max_index] = quadratic(wavelengths[chirp_min_index:chirp_max_index],a,b,c)
    
    plt.figure()
    plt.plot(wavelengths,shifts)
    plt.plot(wavelengths,fitted)
    plt.show()
    
    output = delta_A
    for w in range(len(wavelengths)):
        replace = np.zeros(len(times))
        begin_index = int(min(max(fitted[w],0),len(fitted)-1))
        replace[:len(replace)-begin_index] = delta_A[begin_index:,w]
        output[:,w] = replace
    return output

def ask_value(type, label="value", default = None, override_text = None):
    if override_text==None:
        value = input("Enter " + label + ": ")
    else:
        value = input(override_text)
    if value=="":
        return default
    else:
        return type(value)

def ask_range(type, default=(None,None)):
    min_value = input("Enter new min: ")
    if min_value=="":
        min_value = default[0]
    else:
        min_value = type(min_value)
    max_value = input("Enter new max: ")
    if max_value=="":
        max_value = default[1]
    else:
        max_value = type(max_value)
    return (min_value, max_value)

def remove_nan(delta_A):
    num_t, num_w = delta_A.shape
    for w in range(num_w):
            for t in range(num_t):
                if np.isnan(delta_A[t,w]):
                    if w-1>=0 and w+1<num_w and not np.isnan(delta_A[t,w-1]) and not np.isnan(delta_A[t,w+1]):
                        delta_A[t,w] = (delta_A[t,w-1] + delta_A[t,w+1]) / 2.
                    if t-1>=0 and t+1<num_t and not np.isnan(delta_A[t-1,w]) and not np.isnan(delta_A[t+1,w]):
                        delta_A[t,w] = (delta_A[t-1,w] + delta_A[t+1,w]) / 2.
    return delta_A

def remove_spikes(delta_A, width, factor):
    half_width = int(width/2);
    for t in range(delta_A.shape[0]):
        num_w = delta_A.shape[1]
        for w in range(half_width,num_w-half_width):
            mean = np.nanmean(delta_A[t,w-half_width:w+half_width])
            std = np.nanstd(delta_A[t,w-half_width:w+half_width])
            if (not np.isnan(mean) and abs(delta_A[t,w]-mean)>factor*std):
                delta_A[t,w] = float("nan")
    return delta_A

def gaussian(x, sigma, mu):
    return np.exp(-.5*(x-mu)**2/sigma**2)/(sigma*math.sqrt(2*math.pi))