# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:22:43 2022

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

# PARAMETERS
# if supplying pump off/on files
pump_off_filename = None
pump_on_filename = None
# if supplying delta A file
delta_A_filenames = ["sample1/CudmpDPEphosBF4ACN_1_scan1.csv"]
subtract_surface_file = None
time_zero_correction = (-100,-.5) #time units

# READ IN DATA
if pump_off_filename!=None and pump_on_filename!=None:
    wavelengths_off, times_off, pump_off = helpers.readFile(pump_off_filename)
    wavelengths_on, times_on, pump_on = helpers.readFile(pump_on_filename)
    if not np.allclose(wavelengths_off,wavelengths_on) or not np.allclose(times_off,times_on):
        print("Pump off and pump on raw data does not have matching axes")
        sys.exit()
    delta_A = np.log(pump_off/pump_on)
    wavelengths = wavelengths_on
    times = times_on
elif len(delta_A_filenames)!=0:
    delta_As = []
    for filename in delta_A_filenames:
        wavelengths, times, delta_A = helpers.readFile(filename)
        delta_As.append(delta_A)
    delta_As = np.array(delta_As)
    delta_A = np.average(delta_As, axis=0)

# SUBTRACT SURFACE IF NEEDED
if subtract_surface_file!=None:
    wavelengths_subtract, times_subtract, subtract_surface = helpers.readFile(subtract_surface_file)
    delta_A -= subtract_surface
    


# RECORD KEEPING
wavelength_index = helpers.find_index(wavelengths, (wavelengths.min()+wavelengths.max())/2)
time_index = helpers.find_index(times, (times.min()+times.max())/2)
wavelength_bounds = [np.min(wavelengths),np.max(wavelengths)]
time_bounds = [np.min(times),np.max(times)]

# TIME ZERO CORRECTION
# time_zero_index = (helpers.find_index(times,time_zero_correction[0]),helpers.find_index(times,time_zero_correction[1]))
# print(time_zero_index)
# masked_data = np.ma.masked_array(delta_A, np.isnan(delta_A))[time_zero_index[0]:time_zero_index[1]+1,:]
# time_zero_avg = np.average(delta_A[time_zero_index[0]:time_zero_index[1]+1,:], axis=0, weights=masked_data)
# delta_A -= time_zero_avg

helpers.plot_crosssection(wavelengths,times,delta_A,wavelength_index,False)

# shifts = []
# for i in range(len(wavelengths)):
#     shifts.append(helpers.fit_heaviside(times,delta_A[:,i]))
# shifts = np.array(shifts)
# a,b,c = helpers.fit_quadratic(times, shifts)
# fitted = []
# for time in times:
#     fitted.plot(helpers.quadratic(time,a,b,c))
# fitted = np.array(fitted)
# plt.figure()
# plt.plot(times,fitted)
# plt.show()

# LIVE INTERACTION
while True:
    action = input("Action? (q=quit, wc=wavelength change, tc=time change, tp=time plot, wp=wavelength plot, cp=color plot, waxis=change wavelength axis, play=play over time)\n")
    if action=="q":
        print("quitting...")
        break
    elif action=="wc":
        print("changing wavelength value...")
        desired_wavelength = input("Enter new wavelength: ")
        if desired_wavelength!="":
            desired_wavelength = float(desired_wavelength)
            wavelength_index = helpers.find_index(wavelengths,desired_wavelength)
    elif action=="tc":
        print("changing time value...")
        desired_time = input("Enter new time: ")
        if desired_time!="":
            desired_time = float(desired_time)
            time_index = helpers.find_index(times,desired_time)
    elif action=="wp":
        print("plotting wavelength plot...")
        helpers.plot_crosssection(wavelengths,times,delta_A,time_index,True,wavelength_bounds)
    elif action=="tp":
        print("plotting time plot...")
        helpers.plot_crosssection(wavelengths,times,delta_A,wavelength_index,False)
    elif action=="cp":
        print("plotting color plot...")
        helpers.plot_color(wavelengths, times, delta_A, wavelengths[wavelength_index], times[time_index], w_bounds=wavelength_bounds, t_bounds=time_bounds)
    elif action=="waxis":
        print("changing wavelength axis...")
        wavelength_min = input("Enter new min: ")
        if wavelength_min!="":
            wavelength_bounds[0] = float(wavelength_min)
        wavelength_max = float(input("Enter new max: "))
        if wavelength_max!="":
            wavelength_bounds[1] = float(wavelength_max)
    elif action=="taxis":
        print("changing time axis...")
        time_min = input("Enter new min: ")
        if time_min!="":
            time_bounds[0] = float(time_min)
        time_max = input("Enter new max: ")
        if time_max!="":
            time_bounds[1] = float(time_max)
    elif action=="reset":
        print("reseting axis...")
        wavelength_bounds = [np.min(wavelengths),np.max(wavelengths)]
        time_bounds = [np.min(times),np.max(times)]
    elif action=="play":
        print("playing with time...")
        for i in range(len(times)):
            helpers.plot_crosssection(wavelengths,times,delta_A,i,True,wavelength_bounds)
        repeat_flag = input("Repeat? (y/n)")
        if repeat_flag!="y":
            break
    elif action=="cut w":
        print("cutting wavelength range...")
        cut_min = float(input("Enter min of cut range: "))
        cut_max = float(input("Enter max of cut range: "))
        cut_min_index = helpers.find_index(wavelengths, cut_min)
        cut_max_index = helpers.find_index(wavelengths, cut_max)
        delta_A[:,cut_min_index:cut_max_index+1] = 0
    else:
        print("error, did not recognize command")
    
# adjust wavelength range