# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 17:35:36 2022

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
from data_object import DataObject
from data_object import DataHandler

from scipy.optimize import curve_fit

# # READ IN DATA
# if pump_off_filename!=None and pump_on_filename!=None:
#     wavelengths_off, times_off, pump_off = helpers.read_file(pump_off_filename)
#     wavelengths_on, times_on, pump_on = helpers.read_file(pump_on_filename)
#     if not np.allclose(wavelengths_off,wavelengths_on) or not np.allclose(times_off,times_on):
#         print("Pump off and pump on raw data does not have matching axes")
#         sys.exit()
#     delta_A = np.log(pump_off/pump_on)
#     wavelengths = wavelengths_on
#     times = times_on


######################


# precision = .001
# t_eval = np.arange(-.4, .6, precision)
def kovalenko(times, beta, tau1, beta_tau2_sq, D0, t0):
    # beta = 1.7e-3 * 10**6 # chirp rate [ps^-2]
    # tau1 = 50e-3 #ps  ##???
    # # beta_tau1 = 2.2;
    # beta_tau2_sq = 42;
    # D0 = 1
    # target_wavelength = 400 ### ADJUST
    # center_wavelength = 460 # nm ##???
    # speed_of_light = 2.99792458e5 # nm / ps
    # omega2 = 2*math.pi*speed_of_light/target_wavelength # rad/ps^-1
    # Omega2 = 2*math.pi*speed_of_light/center_wavelength # rad/ps^-1
    # t0 = (omega2-Omega2)/(2*beta)# frequency dependent
    return D0*np.exp(-(times+t0)**2/tau1**2)*np.sin(1/(2*beta*tau1**2)-((times+t0)**2/(beta*tau1**4))-((times+t0)*t0/(beta_tau2_sq*tau1**2)))
# target_wavelengths = [400,450,500,550,600,650,700,750] # nm
# for target_wavelength in target_wavelengths:
#     Sk = kovalenko(t_eval, target_wavelength)
#     plt.figure()
#     plt.plot(t_eval,Sk)
#     plt.title("Kovalenko: Wavelength = "+ str(target_wavelength))
#     plt.show()

# def ours(times, c1, c2, c3, tau1, t0):
#     # beta = 1.7e-3 * 10**6 # chirp rate [ps^-2]
#     # tau1 = 50e-3 #ps  ##???
#     # target_wavelength = 400 #### ADJUST
#     # center_wavelength = 460 # nm ##???
#     # speed_of_light = 2.99792458e5 # nm / ps
#     # omega2 = 2*math.pi*speed_of_light/target_wavelength # rad/ps^-1
#     # Omega2 = 2*math.pi*speed_of_light/center_wavelength # rad/ps^-1
#     # t0 = (omega2-Omega2)/(2*beta)# frequency dependent
#     # c1 = t0/(2*beta)
#     # c2 = t0/(2*beta)
#     # c3 = -1/(4*beta)
#     return np.exp(-(times-t0)**2/tau1**2)*(c1-c2*2*(times-t0)/tau1**2-c3*(2/tau1**2-4*(times-t0)**2/tau1**4))

# for target_wavelength in target_wavelengths:
#     So = ours(t_eval, target_wavelength)
#     plt.figure()
#     plt.plot(t_eval,So)
#     plt.title("Ours: Wavelength = "+ str(target_wavelength))
#     plt.show()

# CREATE TEMP DATA
# temp_deltaA = ref_surfaces[1]

# wavelength_index = helpers.find_index(temp_deltaA.wavelengths,500) # *** 360~570 (check notes) -> 410-430 avoid
# time_min_index = helpers.find_index(temp_deltaA.times,-1)
# time_max_index = helpers.find_index(temp_deltaA.times,2)
# times = temp_deltaA.times[time_min_index:time_max_index]
# sliced_signal = temp_deltaA.signal[time_min_index:time_max_index,wavelength_index]

# FIT KOVALENKO
# popt, pcov = curve_fit(kovalenko, times, sliced_signal, [3e3, 50e-3, 42, .006, .75])
# print(popt)
# fitted = kovalenko(times, 1.7e3, 50e-3, 42, .006,.18)
# fitted = kovalenko(times, *popt)
# plt.figure()
# plt.plot(times, sliced_signal)
# plt.plot(times, fitted)
# plt.show()

# FIT OURS
# popt, pcov = curve_fit(ours, times, sliced_signal, [5.3e-5, 5.3e-5, -0.0001, .18, .0004])
# # print(popt)
# initial = ours(times, 5.3e-5, 5.3e-5, -0.0001, .05, -0.7)
# fitted = ours(times, *popt)
# plt.figure()
# plt.plot(times, sliced_signal)
# plt.plot(times, fitted)
# plt.plot(times, initial)
# plt.show()

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

# FIT OURS
# temp_deltaA = ref_surfaces[1]
# min_index = helpers.find_index(temp_deltaA.wavelengths,400)
# max_index = helpers.find_index(temp_deltaA.wavelengths,410)
# for i in range(min_index,max_index):
#     time_min_index = helpers.find_index(temp_deltaA.times,-1)
#     time_max_index = helpers.find_index(temp_deltaA.times,2)
#     times = temp_deltaA.times[time_min_index:time_max_index]
#     sliced_signal = temp_deltaA.signal[time_min_index:time_max_index,i]
    
#     t0 = temp_deltaA.find_t_peak(temp_deltaA.wavelengths[i],(-1,2))
#     # print(t0)
    
#     popt, pcov = curve_fit(ours, times, sliced_signal, [-4.79227948e-04, 3.87793025e-05, -6.49229051e-06, -t0, 1.36241557e-01])
#     print(popt)
#     fitted = ours(times, *popt)
#     plt.figure()
#     plt.plot(times, sliced_signal)
#     plt.plot(times, fitted)
#     plt.title(str(temp_deltaA.wavelengths[i]) + "," + str(t0) + "," + str(-popt[3]))
#     plt.show()


def lorenc(times):
    alpha = 1
    omega = 1
    beta = 1.7e-3 * 10**6
    tau = 50e-3
    tau_gvd = 150e-3 #fs
    return 2*np.log(1+alpha*omega/(beta*tau**2*tau_gvd)*(times*np.exp(-2*times**2/tau**2)-(times-tau_gvd)*np.exp(-2*(times-tau_gvd)**2/tau**2)))
# Sl = lorenc(t_eval)
# plt.figure()
# plt.plot(t_eval,Sl)
# plt.title("Lorenc")
# plt.show()