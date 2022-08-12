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
    '''
    Return a tuple of (wavelengths, times, signal) based on file called filename
    File should in the format where first row is time axis and first column is wavelength axis
    wavelengths and times are these axes
    signal is the 2D signal for the corresponding wavelengths and times
    '''
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
    '''
    return average of two values

    '''
    return (a+b)/2.

def find_index(arr, value):
    '''
    find the index of element in arr that is closest to the value value

    '''
    return np.argmin(np.abs(arr-value))

def heaviside(x, shift, magnitude):
    '''
    return the heaviside step function where the step occurs at shift with magnitude magnitude

    '''
    return magnitude * np.heaviside(x,shift)

def fit_heaviside(x,y):
    '''
    fit the heaviside function based on minimum sum of squared differences, provided an axis x and signal y
    return the index of where the step should be
    '''
    output_index = 0 # the output index of where the heaviside function jumps (default is 0)
    best_loss = float('inf') # best loss (minimum sum of squared differences)
    for i in range(1,len(x)):
        # i represents which index the jump could be
        # calculate loss if jump is at i (assume jump to MEAN value)
        loss = np.sum((y[:i])**2) + np.sum((y[i:]-np.nanmean(y[i:]))**2)
        
        if loss<best_loss: # update new best_loss and output_index
            output_index = i
            best_loss = loss
    return output_index

def fifth_poly(x, a, b, c, d, e, f):
    '''
    return fifth-degree polynomial
    '''
    return a + b*x + c*x**2 + d*x**3 + e*x**4 + f*x**5

def fourth_poly(x, a, b, c, d, e):
    '''
    return fourth-degree polynomial
    '''
    return a + b*x + c*x**2 + d*x**3 + e*x**4

def third_poly(x, a, b, c, d):
    '''
    return third-degree polynomial
    '''
    return a + b*x + c*x**2 + d*x**3

def second_poly(x, a, b, c):
    '''
    return second-degree polynomial
    '''
    return a + b*x + c*x**2

def first_poly(x, a, b):
    '''
    return first-degree polynomial
    '''
    return a + b*x

def fit_quadratic(x,y):
    cleaned_x = x[~np.isnan(y)]
    cleaned_y = y[~np.isnan(y)]
    
    popt, pcov = curve_fit(second_poly, cleaned_x, cleaned_y)
    a = popt[0]
    b = popt[1]
    c = popt[2]
    
    return (a,b,c)

def chirp_correction(times,wavelengths,delta_A):
    # shifts will store the t0 based on the hea
    shifts = []
    for i in range(len(wavelengths)):   
        shift = fit_heaviside(times,delta_A[:,i])
        shifts.append(shift)
        
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
    fitted[chirp_min_index:chirp_max_index] = second_poly(wavelengths[chirp_min_index:chirp_max_index],a,b,c)
    
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
    '''
    asks for a single value from user and returns it
    returns default if invalid input is provided
    '''
    
    if override_text==None:
        value = input("Enter " + label + ": ")
    else:
        value = input(override_text)
        
    if value=="":
        return default
    else:
        return type(value)

def ask_range(type, default=(None,None), add_text = None):
    '''
    asks for a range from user and returns it as a tuple (min,max)
    returns default if invalid input is provided
    '''
    
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
    '''
    asks for which item in list_of_data user would like and returns the index
    If user inputs invalid index value, returns None
    '''
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
    '''
    asks for a boolean from user
    returns default if invalid input is provided
    '''
    ans = input(text + "(y/n) ")
    if ans=="y":
        return True
    elif ans=="n":
        return False
    else:
        return default

def ask_for_indices(list_of_data):
    '''
    Asks for multiple index values in list_of_data user would like
    If user specifies invalid index, will print an error and ask for more indices
    Will continue to ask for more indices until empty string is entered
    '''
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

def gaussian(x, sigma, mu=0, factor=1):
    '''
    returns gaussian function
    multiplies by factor
    is not normalized

    '''
    g = np.exp(-(x-mu)**2/(2*sigma**2)) #/(sigma*math.sqrt(2*math.pi))
    g = g*factor
    # g /= np.trapz(g)
    return g

def convert_to_ang_freq(wavelength):
    '''
    convert wavelength nm to angular frequency rad/ps-1
    '''
    speed_of_light = 2.99792458e5 # nm / ps
    return 2*math.pi*speed_of_light/wavelength # rad/ps^-1


def ours(times, c1, c2, c3, tau1, t0):
    '''
    return XPM signal function that we wrote out based on Kovalenko's paper
    note: there is another XPM function that Kovalenko writes out explicitly but is not symmetric about center, so we primarily use our XPM signal model currently
    '''
    ### below are commented out since we take c1,c2,c3,tau1,t0 as fitting parameters, but these are theoretical values/estimates of those values
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

def solve_diffeq(t, B2, B3, B4):
    '''
    returns solution of differential equation before convolving with gaussian
    '''
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
    
    # for negative times, all population should be at ground state
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
    '''
    convolve provided signal x with axis t with the gaussian with radius sigma
    
    note: there are 2 methods for convolution - manual and using fourier transform
    as our time axis may not be uniform we cannot use FFT, so it will be easier to implement with manual convolution, which is what we do here
    '''
    area = np.trapz(gaussian(t,sigma)) # calculate area of gaussian on entire time axis for normalization
    precision = abs(t[0]-t[1]) # assuming t has uniform precision
    centered_t = np.arange(-(len(t)/2)*precision/2,(len(t)/2)*precision,precision)
    convolved = np.convolve(x,gaussian(centered_t,sigma),mode='same')/area
    return convolved
    
    # output = []
    # for i in range(len(t)):
    #     value = np.sum(x*gaussian(t,sigma,mu=t[i],factor=1/area)) # convolve manually
    #     output.append(value)
    
    # plt.figure()
    # plt.plot(t,x)
    # plt.plot(t,convolved)
    # plt.plot(t,output)
    # plt.show()
    
    # return np.array(output)

def all_convolve_gaussian(t_axis, x1,x2,x3,x4,sigma):
    '''
    Convolve x1, x2, x3, x4, which share a time axis t_axis, with a gaussian of size sigma 

    '''
    
    x1_blurred = convolve_gaussian(t_axis, x1, sigma)
    x2_blurred = convolve_gaussian(t_axis, x2, sigma)
    x3_blurred = convolve_gaussian(t_axis, x3, sigma)
    x4_blurred = convolve_gaussian(t_axis, x4, sigma)
    
    # plotting convolution result
    # plt.figure()
    # plt.plot(extended_t,x1_blurred)
    # plt.plot(extended_t,x2_blurred)
    # plt.plot(extended_t,x3_blurred)
    # plt.plot(extended_t,x4_blurred)
    # plt.show()
    
    return x1_blurred,x2_blurred,x3_blurred,x4_blurred

def rateModel(t, B2, B3, B4, A1, A2, A3, A4, sigma):
    '''
    The rate model equation that was manually solved by hand
    t is the time axis we would like values for
    B2, B3, B4 is the decay rate of the 2nd, 3rd, 4th excited state respectively
    A1, A2, A3, A4 is the magnitude of each excited state for the signal
    sigma is the size of the gaussian we will use to smear/blur the differential equation solutions
    '''
    
    # note that precision is not uniform on time axis
    # recreate a time axis that is the smallest precision or tenth of sigma
    # double length of time axis for convolution padding purposes
    precisions = np.array(t[1:] - t[:len(t)-1]) # stagger array to get precision at each time step
    smallest_precision = min(precisions) # find smallest precision
    if smallest_precision>=sigma: # if precision is smaller than sigma our precision should be smaller
        smallest_precision = sigma/10
    extended_left = min(t)-(len(t)//2+1)*smallest_precision # min value for extended time axis
    extended_right = max(t)+(len(t)//2+1)*smallest_precision # max value for extended time axis
    t_extended = np.arange(extended_left, extended_right, smallest_precision) # our new extended time axis
    
    # create solutions to differential equation
    x1_extended, x2_extended, x3_extended, x4_extended = solve_diffeq(t_extended, B2, B3, B4)    
    
    # plt.figure()
    # plt.plot(t_extended,x1_extended)
    # plt.plot(t_extended,x2_extended)
    # plt.plot(t_extended,x3_extended)
    # plt.plot(t_extended,x4_extended)
    # plt.show()
    
    # blur solution with gaussian
    x1_extended = convolve_gaussian(t_extended, x1_extended, sigma)
    x2_extended = convolve_gaussian(t_extended, x2_extended, sigma)
    x3_extended = convolve_gaussian(t_extended, x3_extended, sigma)
    x4_extended = convolve_gaussian(t_extended, x4_extended, sigma)
    
    # based on x1_extended, x2_extended, x3_extended, x4_extended, interpolate values for time values we actually want in t
    x1 = []
    x2 = []
    x3 = []
    x4 = []
    for new_t in t:
        x1.append(np.interp(new_t, t_extended, x1_extended))
        x2.append(np.interp(new_t, t_extended, x2_extended))
        x3.append(np.interp(new_t, t_extended, x3_extended))
        x4.append(np.interp(new_t, t_extended, x4_extended))
    x1 = np.array(x1)
    x2 = np.array(x2)
    x3 = np.array(x3)
    x4 = np.array(x4)
    
    # plt.figure()
    # plt.plot(t,x1)
    # plt.plot(t,x2)
    # plt.plot(t,x3)
    # plt.plot(t,x4)
    # plt.show()
    
    # super impose with provided amplitude A1, A2, A3, A4
    superimposed = A1*x1 + A2*x2 + A3*x3 + A4*x4
    return superimposed

def rateModel2(t, B2, B3, B4, A1, A2, A3, A4, sigma, t0):
    '''
    same as rateModel but the differential equation solution shifted by t0

    '''
    
    # rate model with no chirp correction
    
    # sigma related to tau1
    left_precision = abs(t[0]-t[1])
    right_precision = abs(t[len(t)-1]-t[len(t)-2])
    smallest_precision = min(left_precision, right_precision)
    if smallest_precision>=sigma:
        smallest_precision = sigma/10
    extended_left = min(t)-(len(t)//2+1)*smallest_precision
    extended_right = max(t)+(len(t)//2+1)*smallest_precision
    t_extended = np.arange(extended_left, extended_right, smallest_precision)
    
    x1_extended, x2_extended, x3_extended, x4_extended = solve_diffeq(t_extended-t0, B2, B3, B4)    # shift our differential equation solution 
    
    # plt.figure()
    # plt.plot(t_extended,x1_extended)
    # plt.plot(t_extended,x2_extended)
    # plt.plot(t_extended,x3_extended)
    # plt.plot(t_extended,x4_extended)
    # plt.show()
    
    x1_extended, x2_extended, x3_extended, x4_extended = all_convolve_gaussian(t_extended-t0, x1_extended, x2_extended, x3_extended, x4_extended, sigma)
    
    x1 = []
    x2 = []
    x3 = []
    x4 = []
    for new_t in t:
        x1.append(np.interp(new_t, t_extended, x1_extended))
        x2.append(np.interp(new_t, t_extended, x2_extended))
        x3.append(np.interp(new_t, t_extended, x3_extended))
        x4.append(np.interp(new_t, t_extended, x4_extended))
    x1 = np.array(x1)
    x2 = np.array(x2)
    x3 = np.array(x3)
    x4 = np.array(x4)
    
    # plt.figure()
    # plt.plot(t,x1)
    # plt.plot(t,x2)
    # plt.plot(t,x3)
    # plt.plot(t,x4)
    # plt.show()
    
    superimposed = A1*x1 + A2*x2 + A3*x3 + A4*x4
    return superimposed

def rateModel3(t, B2, B3, B4, A1, A2, A3, A4, sigma, t0, c1, c2, c3):
    '''
    same as rateModel2 but add XPM sigmal with parameters c1, c2, c3 on top

    '''
    # rate model with no chirp correction and without subtracting XPM
    signal = rateModel2(t,B2,B3,B4,A1,A2,A3,A4,sigma,t0)
    XPM = ours(t,c1,c2,c3,math.sqrt(2)*sigma,t0) # add XPM signal
    return signal + XPM