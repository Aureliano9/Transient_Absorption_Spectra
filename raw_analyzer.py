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
import os
from data_object import DataObject

from scipy.optimize import curve_fit



# PARAMETERS
# if supplying pump off/on files
pump_off_filename = None
pump_on_filename = None
# if supplying delta A file
delta_A_filenames = ["Data (1)/CudmpDPEphosBF4ACN_1_scan1.csv", "Data (1)/CudmpDPEphosBF4ACN_1_scan2.csv"]  # ["Four approaches for XPM treatment/CudmpDPEphosBF4ACN_1_scan1.csv"] #["sample1/CudmpDPEphosBF4ACN_1_scan2.csv","sample1/CudmpDPEphosBF4ACN_1_scan3.csv","sample1/CudmpDPEphosBF4ACN_1_scan4.csv"]
subtract_surface_filenames = ["Data (1)/Acetonitrile.csv","Data (1)/Acetonitrile2.csv"] # "Four approaches for XPM treatment/Acetonitrile2_scan1.csv", 
# time_zero_correction = (-100,-.5) #time units

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
# elif len(delta_A_filenames)!=0:
#     if len(delta_A_filenames)==1:
#         wavelengths, times, delta_A = helpers.read_file(delta_A_filenames[0])
#     else:
#         delta_As = []
#         for filename in delta_A_filenames:
#             wavelengths, times, delta_A = helpers.read_file(filename)
#             if delta_As.shape()[0]>0 and delta_As[0,:,:].shape()!=delta_A.shape():
#                 print("ERROR: all files for delta_A_filenames should have the same dimension")
#                 os.abort()
#             delta_As.append(delta_A)
#         delta_A = np.array(delta_As)
if len(delta_A_filenames)==0:
    print("Error: No data to read! Change variable delta_A_filenames")
delta_As = []
for filename in delta_A_filenames:
    delta_As.append(DataObject.CreateFromFile(filename))

# READ IN REFERENCE
ref_surfaces = [] # each element is 3 items
for filename in subtract_surface_filenames:
    ref_surfaces.append(DataObject.CreateFromFile(filename))


def apply(list_of_data, method, kwargs):
    if len(list_of_data)==1:
        return list_of_data[0].__getattribute__(method)(**kwargs)
    else:
        specify_index = helpers.ask_yes_no("Apply only to a specific layer?")
        if specify_index:
            display_index = helpers.ask_which_layer(list_of_data)
            if display_index!=None:
                return list_of_data[display_index].__getattribute__(method)(**kwargs)
        else:
            output = []
            for delta_A in list_of_data:
                output.append(delta_A.__getattribute__(method)(**kwargs))
            return output

def apply_one(list_of_data, method, kwargs):
    if len(list_of_data)==1:
        return list_of_data[0].__getattribute__(method)(**kwargs)
    else:
        display_index = helpers.ask_which_layer(list_of_data)
        if display_index!=None:
            return list_of_data[display_index].__getattribute__(method)(**kwargs)
    
# speed_of_light = 2.99792458e5
# wavelengths = speed_of_light/wavelengths


# # RECORD KEEPING
switched_data_ref = False;

precision = .001
t_eval = np.arange(-.4, .6, precision)
def kovalenko(times, beta, tau1, beta_tau2_sq, D0):
    # beta = 1.7e-3 * 10**6 # chirp rate [ps^-2]
    # tau1 = 50e-3 #ps  ##???
    # # beta_tau1 = 2.2;
    # beta_tau2_sq = 42;
    # D0 = 1
    target_wavelength = 400 ### ADJUST
    center_wavelength = 460 # nm ##???
    speed_of_light = 2.99792458e5 # nm / ps
    omega2 = 2*math.pi*speed_of_light/target_wavelength # rad/ps^-1
    Omega2 = 2*math.pi*speed_of_light/center_wavelength # rad/ps^-1
    t0 = (omega2-Omega2)/(2*beta)# frequency dependent
    return D0*np.exp(-(times+t0)**2/tau1**2)*np.sin(1/(2*beta*tau1**2)-((times+t0)**2/(beta*tau1**4))-((times+t0)*t0/(beta_tau2_sq*tau1**2)))
# target_wavelengths = [400,450,500,550,600,650,700,750] # nm
# for target_wavelength in target_wavelengths:
#     Sk = kovalenko(t_eval, target_wavelength)
#     plt.figure()
#     plt.plot(t_eval,Sk)
#     plt.title("Kovalenko: Wavelength = "+ str(target_wavelength))
#     plt.show()

def ours(times, c1, c2, c3):
    beta = 1.7e-3 * 10**6 # chirp rate [ps^-2]
    tau1 = 50e-3 #ps  ##???
    target_wavelength = 400 #### ADJUST
    center_wavelength = 460 # nm ##???
    speed_of_light = 2.99792458e5 # nm / ps
    omega2 = 2*math.pi*speed_of_light/target_wavelength # rad/ps^-1
    Omega2 = 2*math.pi*speed_of_light/center_wavelength # rad/ps^-1
    t0 = (omega2-Omega2)/(2*beta)# frequency dependent
    # c1 = t0/(2*beta)
    # c2 = t0/(2*beta)
    # c3 = -1/(4*beta)
    return np.exp(-(times+t0)**2/tau1**2)*(c1-c2*2*(times+t0)/tau1**2-c3*(2/tau1**2-4*(times+t0)**2/tau1**4))

# for target_wavelength in target_wavelengths:
#     So = ours(t_eval, target_wavelength)
#     plt.figure()
#     plt.plot(t_eval,So)
#     plt.title("Ours: Wavelength = "+ str(target_wavelength))
#     plt.show()

# popt, pcov = curve_fit(kovalenko, times, delta_A[:,helpers.find_index(wavelengths,400)], [163, 50e-3, 42, .001])
# print(popt)
# fitted = kovalenko(times, popt[0], popt[1], popt[2], popt[3])
# plt.figure()
# plt.plot(times, delta_A[:,helpers.find_index(wavelengths,400)])
# plt.plot(times, fitted)
# plt.show()

# def lorenc(times):
#     alpha = 1
#     omega = 1
#     beta = 1.7e-3 * 10**6
#     tau = 50e-3
#     tau_gvd = 150e-3 #fs
#     return 2*np.log(1+alpha*omega/(beta*tau**2*tau_gvd)*(times*np.exp(-2*times**2/tau**2)-(times-tau_gvd)*np.exp(-2*(times-tau_gvd)**2/tau**2)))
# Sl = lorenc(t_eval)
# plt.figure()
# plt.plot(t_eval,Sl)
# plt.title("Lorenc")
# plt.show()

# MENU
menu = {}
menu["q"] = "quit"
menu["wc"] = "change current wavelength"
menu["tc"] = "change current time"
menu["tp"] = "plot cross-section at current wavelength s.t. time is x-axis"
menu["wp"] = "plot cross-section at current time s.t. wavelength is x-axis"
menu["cp"] = "plot 2D color signal"
menu["waxis"] = "change wavelength axis bounds"
menu["taxis"] = "change time axis bounds"
menu["caxis"] = "change color axis bounds"
menu["reset axis"] = "reset wavelength/time/color axis bounds"
menu["p"] = "play wavelength cross-sections over all time"
menu["subtract"] = "subtract surfact from current working surface"
menu["avg"] = "take average of signal across multiple files"
menu["shift time"] = "shift time values"
menu["cut w"] = "cut out wavelength range"
menu["spikes"] = "remove spikes"
menu["nan w"] = "remove nan wavelength spectra"
menu["nan t"] = "remove nan time spectra"
menu["background"] = "perform background correction"
menu["chirp"] = "perform chirp correction"
menu["reset data"] = "reset all data changes"
menu["reset ref"] = "reset all data changes to reference surface"
menu["background ref"] = "perform background correction to reference surface"
menu["switch data ref"] = "switch data and reference for display purposes"

# LIVE INTERACTION
while True:
    
    action = input("Action? (enter 'm' to see command menu)\n")
    if action=="q":
        print("quitting...")
        break
    
    elif action=="m":
        print("displaying menu...\n")
        width  = 20
        print(f"{'Type' : <20}{'Command' : <20}")
        for key in menu:
            value = menu[key]
            print(f"{key : <20}{value : <20}")
            
    elif action=="wc":
        print("changing wavelength value...")
        new_w = helpers.ask_value(float, "wavelength")
        if new_w!=None:
            apply(delta_As,"change_current_w",{"new_w": new_w})
            
    elif action=="tc":
        print("changing time value...")
        new_t = helpers.ask_value(float, "time")
        if new_t!=None:
            apply(delta_As,"change_current_t",{"new_t": new_t})
            
    elif action=="wp":
        print("plotting wavelength plot...")
        apply_one(delta_As,"plot_wavelength_crosssection", {})
        
    elif action=="tp":
        print("plotting time plot...")
        apply_one(delta_As,"plot_time_crosssection", {})
        
    elif action=="cp":
        print("plotting color plot...")
        apply_one(delta_As,"plot_color", {})
        
    elif action=="waxis":
        print("changing wavelength axis...")
        new_range = helpers.ask_range(float)
        apply(delta_As,"change_waxis",{'new_range': new_range})
        
    elif action=="taxis":
        print("changing time axis...")
        new_range = helpers.ask_range(float)
        apply(delta_As,"change_taxis",{'new_range': new_range})
        
    elif action=="caxis":
        print("changing color axis...")
        new_range = helpers.ask_range(float)
        apply(delta_As,"change_caxis",{'new_range': new_range})
        
    elif action=="reset axis":
        print("reseting axis...")
        axis_index = helpers.ask_value(int, default=None, override_text="Which axis to reset? (0=wavelength, 1=time, 2=color, 3=all) ")
        apply(delta_As,"reset_axis",{})
        
    elif action=="play":
        print("playing with time...")
        apply_one(delta_As,"play_over_time",{})
        
    elif action=="t peak":
        print("Find peak in given range")
        peak_time = apply_one("find_t_peak",{})
        print("Peak at time: ", peak_time)
    
    elif action=="add data":
        print("Adding new data surface...")
        filename = input("Filename? ")
        delta_As.append(DataObject.CreateFromFile(filename))
        
    elif action=="add reference":
        print("Adding new reference surface...")
        filename = input("Filename? ")
        ref_surfaces.append(DataObject.CreateFromFile(filename))
        
    # MUTATING DATA
    elif action=="subtract":
        # SUBTRACT SURFACE IF NEEDED        
        subtract_index = helpers.ask_which_layer(ref_surfaces)
        if subtract_index!=None:
            surface_to_subtract = ref_surfaces[subtract_index]
            
            Er = helpers.ask_value(float, default=None, override_text="Energy for subtract surface: ") #650
            Es = helpers.ask_value(float, default=None, override_text="Energy for original surface: ") #500
            f = helpers.ask_value(float, default=None, override_text="Fraction f: ")
            
            if Es!=None and Er!=None and f!=None:
                # subtract surface 
                apply(delta_As,"subtract_surface",{"surface_to_subtract": surface_to_subtract, "Es": Es, "Er": Er, "f": f})
            else:
                print("Error: Enter valid value")
        else:
            print("Error: no file specified")
            
    elif action=="avg":
        indices = helpers.ask_for_indices(delta_As)
        if len(indices)>0:
            datas=[]
            for index in indices:
                datas.append(delta_As[index])
            delta_As.append(DataObject.average(datas))
        # if delta_A.ndim==3:
        #     delta_A = np.nanmean(delta_A, axis=0)
        # else:
        #     print("only one file")
        
    elif action=="shift time":
        uniform_time_shift = helpers.ask_value(float, default=0)
        apply(delta_As,"time_shift",{"shift_time": uniform_time_shift})
        
    elif action=="cut w":
        print("cutting wavelength range...") # ask user nan/delete
        cut_min, cut_max = helpers.ask_range(float)
        apply(delta_As,"cut_w",{"cut_min": cut_min, "cut_max": cut_max})
        
    elif action=="spikes":
        print("removing spikes...")
        width = helpers.ask_value(int,"width for spikes")
        factor = helpers.ask_value(float,"number of standard deviations allowed")
        apply(delta_As,"remove_spikes",{"width": width, "factor": factor})
        
    elif action=="nan w":
        print("removing nan wavelength spectra...")
        apply(delta_As,"remove_nan_w",{})
        
    elif action=="nan t":
        print("removing nan time spectra...")
        apply(delta_As,"remove_nan_t",{})
    
    elif action=="chirp":
        print("performing chirp correction...")
        print("ERROR: not implemented yet")
        # delta_A = helpers.chirp_correction(times,wavelengths,delta_A)
        
    # want to perform background correction for data and reference
    elif action=="background":
        print("performing background correction...")
        back_min, back_max = helpers.ask_range(float)
        apply(delta_As,"background_correction",{"back_min":back_min, "back_max":back_max})
    
    elif action=="background ref":
        print("performing background correction on reference surface...")
        back_min, back_max = helpers.ask_range(float)
        apply(ref_surfaces,"background_correction",{"back_min":back_min, "back_max":back_max})
        
    elif action=="reset data":
        print("reseting to original data...")
        apply(delta_As,"reset_data",{})
    
    elif action=="reset ref":
        print("reseting to original data...")
        apply(ref_surfaces,"reset_data",{})
        
    elif action=="switch data ref":
        print("switching data and reference for display purposes...")
        temp = delta_As
        delta_As = ref_surfaces
        ref_surfaces = temp
        switched_data_ref = ~switched_data_ref
        print("Data/Ref inverted" if switched_data_ref else "Data/Ref NOT inverted")
    else:
        print("error, did not recognize command")
    