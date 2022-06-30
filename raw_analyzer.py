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
pump_off_filename = "continuum spectrum_channel1 raw Pump off.csv"
pump_on_filename = "continuum spectrum_channel1 raw Pump on.csv"
subtract_surface_file = "continuum spectrum_channel1 raw Pump on.csv"
time_zero_correction = (0,50) #time units

# READ IN DATA
wavelengths_off, times_off, pump_off = helpers.readFile(pump_off_filename)
wavelengths_on, times_on, pump_on = helpers.readFile(pump_on_filename)
if not np.allclose(wavelengths_off,wavelengths_on) or not np.allclose(times_off,times_on):
    print("Pump off and pump on raw data does not have matching axes")
    sys.exit()
delta_A = np.log(pump_off/pump_on)
wavelengths = wavelengths_on
times = times_on

# SUBTRACT SURFACE IF NEEDED
if subtract_surface_file!=None:
    wavelengths_subtract, times_subtract, subtract_surface = helpers.readFile(subtract_surface_file)
    delta_A -= subtract_surface

# REMOVE NAN OR INVALID VALUES
delta_A = np.nan_to_num(delta_A)
    
# TIME ZERO CORRECTION
time_zero_index = (helpers.find_index(times,time_zero_correction[0]),helpers.find_index(times,time_zero_correction[1]))
time_zero_avg = np.average(delta_A[time_zero_index[0]:time_zero_index[1]+1,:], axis=0)
delta_A -= time_zero_avg

# RECORD KEEPING
wavelength_index = len(wavelengths)//2
time_index = len(times)//2
wavelength_bounds = [np.min(wavelengths),np.max(wavelengths)]
time_bounds = [np.min(times),np.max(times)]

helpers.plot_color(wavelengths, times, delta_A, wavelengths[wavelength_index], times[time_index])

# LIVE INTERACTION
while True:
    action = input("Action? (q=quit, wc=wavelength change, tc=time change, tp=time plot, wp=wavelength plot, cp=color plot, waxis=change wavelength axis, play=play over time)\n")
    if action=="q":
        print("quitting...")
        break
    elif action=="wc":
        print("changing wavelength value...")
        desired_wavelength = float(input("Enter new wavelength: "))
        wavelength_index = np.find_index(wavelengths,desired_wavelength)
    elif action=="tc":
        print("changing time value...")
        desired_time = float(input("Enter new time: "))
        time_index = np.find_index(times,desired_time)
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