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

def ask_which_layer(list_of_data, default=None):
    filenames = []
    for data in list_of_data:
        filenames.append(data.get_name())
    print("Options:", filenames)
    display_index = ask_value(int, override_text="Which layer? ")
    
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

def gaussian(x, sigma):
    g = np.exp(-x**2/(2*sigma**2))/(sigma*math.sqrt(2*math.pi))
    g /= np.trapz(g)
    return g