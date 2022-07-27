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
import os
from scipy import ndimage
from scipy import fftpack


# return matrix from file
def read_file(filename):
    wavelengths = []
    times = []
    signal = []
     
    # opening the CSV file
    with open(filename, mode ='r')as file:
       
      if filename[-3:]=="dat":
          lines = file.readlines()
      else:
          # reading the CSV file
          lines = csv.reader(file)
     
      # displaying the contents of the CSV file
      row_counter = 0
      for line in lines:
          column_counter = 0
          
          if filename[-3:]=="dat":
              line = line.split("\t")
          
          #break condition
          if len(times)>0 and len(line)<len(times):
              break
            
          for el in line:
              if row_counter == 0:
                  if column_counter>0:
                      times.append(float(el))
              else:
                  if column_counter==0:
                      wavelengths.append(float(el))
                      signal_row = []
                  else:
                      signal_row.append(float(el))
              column_counter += 1
          if row_counter>0:
              signal.append(signal_row)
          row_counter += 1
          
    wavelengths = np.array(wavelengths)
    times = np.array(times)
    signal = np.array(signal).transpose()
    
    if filename[-3:]=="dat":
        return (times,wavelengths,signal.transpose())
    
    return (wavelengths,times,signal)
    
def avg(a,b):
    return (a+b)/2.

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

def ask_range(type, default=(None,None), add_text = None):
    if add_text!=None:
        print(add_text)
        
    min_value = input("Enter new min: ")
    if min_value=="":
        print("Using default value")
        min_value = default[0]
    else:
        min_value = type(min_value)
    max_value = input("Enter new max: ")
    if max_value=="":
        print("Using default value")
        max_value = default[1]
    else:
        max_value = type(max_value)
    return (min_value, max_value)

def ask_which_layer(list_of_data):
    filenames = []
    for data in list_of_data:
        filenames.append(data.get_name())
    print("Options:", filenames)
    display_index = ask_value(int, override_text="Which layer? ", default=len(filenames))
    
    if display_index>=-1 and display_index<len(filenames):
        return display_index
    else:
        print("Error: Must specify an index between 0 and " + str(len(filenames)))
        return None
    
def ask_yes_no(text, default=False):
    ans = input(text + "(y/n) ")
    if ans=="y":
        return True
    elif ans=="n":
        return False
    else:
        return default

def ask_for_indices(list_of_data):
    filenames = []
    for data in list_of_data:
        filenames.append(data.get_name())
    print("Options:", filenames)
    
    index = input("Specify index to include: ")
    output = []
    while index!="":
        index = int(index)
        if index>=0 and index<len(list_of_data):
            output.append(int(index))
        else:
            print("Error: Only specify index between 0 and " + str(len(list_of_data)-1))
        index = input("Specify index to include: ")
    return output

# def remove_nan(delta_A):
#     num_t, num_w = delta_A.shape
#     for w in range(num_w):
#             for t in range(num_t):
#                 if np.isnan(delta_A[t,w]):
#                     if w-1>=0 and w+1<num_w and not np.isnan(delta_A[t,w-1]) and not np.isnan(delta_A[t,w+1]):
#                         delta_A[t,w] = (delta_A[t,w-1] + delta_A[t,w+1]) / 2.
#                     if t-1>=0 and t+1<num_t and not np.isnan(delta_A[t-1,w]) and not np.isnan(delta_A[t+1,w]):
#                         delta_A[t,w] = (delta_A[t-1,w] + delta_A[t+1,w]) / 2.
#     return delta_A

def gaussian(x, sigma, mu=0,factor=1):
    g = np.exp(-(x-mu)**2/(2*sigma**2)) #/(sigma*math.sqrt(2*math.pi))
    g *= factor
    # g /= np.trapz(g)
    return g

def convert_to_ang_freq(wavelength):
    # convert from wavelength nm to rad/ps-1
    speed_of_light = 2.99792458e5 # nm / ps
    return 2*math.pi*speed_of_light/wavelength # rad/ps^-1


def ours(times, c1, c2, c3, tau1, t0):
    # beta = 1.7e-3 * 10**6 # chirp rate [ps^-2]
    # tau1 = 50e-3 #ps  ##???
    # target_wavelength = 400 #### ADJUST
    # center_wavelength = 460 # nm ##???
    # speed_of_light = 2.99792458e5 # nm / ps
    # omega2 = 2*math.pi*speed_of_light/target_wavelength # rad/ps^-1
    # Omega2 = 2*math.pi*speed_of_light/center_wavelength # rad/ps^-1
    # t0 = (omega2-Omega2)/(2*beta)# frequency dependent
    # c1 = t0/(2*beta)
    # c2 = t0/(2*beta)
    # c3 = -1/(4*beta)
    return np.exp(-(times-t0)**2/tau1**2)*(c1-c2*2*(times-t0)/tau1**2-c3*(2/tau1**2-4*(times-t0)**2/tau1**4))

def fifth(x, a, b, c, d, e, f):
    return a + b*x + c*x**2 + d*x**3 + e*x**4 + f*x**5

def fourth(x, a, b, c, d, e):
    return a + b*x + c*x**2 + d*x**3 + e*x**4

def third(x, a, b, c, d):
    return a + b*x + c*x**2 + d*x**3

def second(x, a, b, c):
    return a + b*x + c*x**2

def first(x, a, b):
    return a + b*x

def solve_diffeq(t, B2, B3, B4):
    f = .5 #???
    initial = [1-f,f,0,0]
    C2 = initial[1]
    C3 = initial[2]+initial[1]*B2/(B2-B3)
    C4 = initial[3]-initial[1]*B2*B3/((B2-B3)*(B2-B4))+C3*B3/(B3-B4)
    C1 = initial[0]-C3*B4/(B3-B4)-B4*C2*B3/((B3-B2)*(B2-B4))+C4
    x1 = B4*C3*np.exp((-B3)*t)/(B3-B4) + B4*C2*B3*np.exp(-B2*t)/((B3-B2)*(B2-B4))-C4*np.exp(-B4*t) + C1
    x2 = C2*np.exp(-B2*t)
    x3 = C3*np.exp(-B3*t)-C2*B2/(B2-B3)*np.exp(-B2*t)
    x4 = C2*B2*B3*np.exp(-B2*t)/((B2-B3)*(B2-B4))-C3*B3*np.exp(-B3*t)/(B3-B4)+C4*np.exp(-B4*t)
    
    # set negative time
    for i in range(len(t)):
        if t[i]<=0:
            x1[i] = 1
            x2[i] = 0
            x3[i] = 0
            x4[i] = 0
        else:
            break
    return x1, x2, x3, x4

def convolve_gaussian(t,x,sigma):
    area = np.trapz(gaussian(t,sigma))
    output = []
    for i in range(len(t)):
        value = np.sum(x*gaussian(t,sigma,mu=t[i],factor=1/area))
        output.append(value)
    return np.array(output)

def all_convolve_gaussian(t,extended_t, x1,x2,x3,x4,sigma):
    # g = gaussian(extended_t,sigma)
    # x1_blurred = ndimage.convolve(x1,g, mode='constant', cval=0.0)
    # x2_blurred = ndimage.convolve(x2,g, mode='constant', cval=0.0)
    # x3_blurred = ndimage.convolve(x3,g, mode='constant', cval=0.0)
    # x4_blurred = ndimage.convolve(x4,g, mode='constant', cval=0.0)
    x1_blurred = convolve_gaussian(extended_t, x1, sigma)
    x2_blurred = convolve_gaussian(extended_t, x2, sigma)
    x3_blurred = convolve_gaussian(extended_t, x3, sigma)
    x4_blurred = convolve_gaussian(extended_t, x4, sigma)
    min_index = find_index(extended_t, min(t))
    max_index = find_index(extended_t, max(t))
    
    # plt.figure()
    # plt.plot(extended_t, ndimage.convolve(np.ones(len(extended_t)),g, mode='constant', cval=0.0))
    # plt.plot(extended_t,x1_blurred)
    # plt.plot(extended_t,x2_blurred)
    # plt.plot(extended_t,x3_blurred)
    # plt.plot(extended_t,x4_blurred)
    # plt.show()
    
    x1_blurred = x1_blurred[min_index:max_index+1]
    x2_blurred = x2_blurred[min_index:max_index+1]
    x3_blurred = x3_blurred[min_index:max_index+1]
    x4_blurred = x4_blurred[min_index:max_index+1]
    return x1_blurred,x2_blurred,x3_blurred,x4_blurred

def rateModel(t, B2, B3, B4, A1, A2, A3, A4, sigma):
    assert(any(t<0) and any(t>0))
    # sigma related to tau1
    left_t = []
    right_t = []
    left_precision = abs(t[0]-t[1])
    right_precision = abs(t[len(t)-1]-t[len(t)-2])
    extended_left = min(t)-(len(t))*left_precision
    for i in range(len(t)):
        left_t.append(extended_left + i*left_precision)
        right_t.append(max(t) + i*right_precision)
    left_t = np.array(left_t)
    right_t = np.array(right_t)
    extended_t = np.append(np.append(left_t,t),right_t)
    
    x1, x2, x3, x4 = solve_diffeq(extended_t, B2, B3, B4)    
    
    plt.figure()
    plt.plot(extended_t,x1)
    plt.plot(extended_t,x2)
    plt.plot(extended_t,x3)
    plt.plot(extended_t,x4)
    plt.show()
    
    
    x1, x2, x3, x4 = all_convolve_gaussian(t,extended_t,x1,x2,x3,x4,sigma)
    
    plt.figure()
    plt.plot(t,x1)
    plt.plot(t,x2)
    plt.plot(t,x3)
    plt.plot(t,x4)
    plt.show()
    
    superimposed = A1*x1 + A2*x2 + A3*x3 + A4*x4
    return superimposed