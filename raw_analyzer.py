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
from data_object import DataHandler

from scipy.optimize import curve_fit



# PARAMETERS
# if supplying pump off/on files
pump_off_filename = None
pump_on_filename = None
# if supplying delta A file
delta_A_filenames = ["Data (1)/CudmpDPEphosBF4ACN_1_scan1.csv", "Data (1)/CudmpDPEphosBF4ACN_1_scan2.csv", "Four approaches for XPM treatment/CudmpDPEphosBF4ACN_1_scan1.csv"]  # ["Four approaches for XPM treatment/CudmpDPEphosBF4ACN_1_scan1.csv"] #["sample1/CudmpDPEphosBF4ACN_1_scan2.csv","sample1/CudmpDPEphosBF4ACN_1_scan3.csv","sample1/CudmpDPEphosBF4ACN_1_scan4.csv"]
ref_surface_filenames = ["Data (1)/Acetonitrile.csv","Data (1)/Acetonitrile2.csv"] # "Four approaches for XPM treatment/Acetonitrile2_scan1.csv", 
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

if len(delta_A_filenames)==0:
    print("Error: No data to read! Change variable delta_A_filenames")

data_handler = DataHandler(delta_A_filenames, ref_surface_filenames)


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

# HARDCODE VARIABLES
back_min_default = -100 # SET
back_max_default = -.5 # SET
ref_index_default = 1 # SET
t_range_default = (-1.5,2)
w_range_default = (346, 780) #570)
Er_default = 650
Es_default = 500
f_default = 0.8258

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
            data_handler.apply(True,"change_current_w",{"new_w": new_w})
            
    elif action=="tc":
        print("changing time value...")
        new_t = helpers.ask_value(float, "time")
        if new_t!=None:
            data_handler.apply(True,"change_current_t",{"new_t": new_t})
            
    elif action=="wp":
        print("plotting wavelength plot...")
        data_handler.apply_one(True,"plot_wavelength_crosssection", {})
        
    elif action=="tp":
        print("plotting time plot...")
        data_handler.apply_one(True,"plot_time_crosssection", {})
        
    elif action=="cp":
        print("plotting color plot...")
        data_handler.apply_one(True,"plot_color", {})
        
    elif action=="waxis":
        print("changing wavelength axis...")
        new_range = helpers.ask_range(float)
        data_handler.apply(True,"change_waxis",{'new_range': new_range})
        
    elif action=="taxis":
        print("changing time axis...")
        new_range = helpers.ask_range(float)
        data_handler.apply(True,"change_taxis",{'new_range': new_range})
        
    elif action=="caxis":
        print("changing color axis...")
        new_range = helpers.ask_range(float)
        data_handler.apply(True,"change_caxis",{'new_range': new_range})
        
    elif action=="reset axis":
        print("reseting axis...")
        which_axis = helpers.ask_value(int, default=None, override_text="Which axis to reset? (0=wavelength, 1=time, 2=color, 3=all) ")
        data_handler.apply(True,"reset_axis",{"which_axis": which_axis})
        
    elif action=="play":
        print("playing with time...")
        data_handler.apply_one(True,"play_over_time",{})
        
    elif action=="t peak":
        print("Find peak in given range")
        peak_time = data_handler.apply_one(True,"find_t_peak",{})
        print("Peak at time: ", peak_time)
    
    elif action=="add data":
        print("Adding new data surface...")
        filename = input("Filename? ")
        data_handler.add_data(DataObject.CreateFromFile(filename))
        
    elif action=="add reference":
        print("Adding new reference surface...")
        filename = input("Filename? ")
        data_handler.add_ref(DataObject.CreateFromFile(filename))
        
    # MUTATING DATA
    elif action=="subtract":
        # SUBTRACT SURFACE IF NEEDED        
        subtract_index = helpers.ask_which_layer(data_handler.reference_surfaces)
        if subtract_index!=None:
            surface_to_subtract = data_handler.reference_surfaces[subtract_index]
            
            Er = helpers.ask_value(float, default=None, override_text="Energy for subtract surface: ") #650
            Es = helpers.ask_value(float, default=None, override_text="Energy for original surface: ") #500
            f = helpers.ask_value(float, default=None, override_text="Fraction f: ")
            
            if Es!=None and Er!=None and f!=None:
                # subtract surface 
                data_handler.apply(True,"subtract_surface",{"surface_to_subtract": surface_to_subtract, "Es": Es, "Er": Er, "f": f})
            else:
                print("Error: Enter valid value")
        else:
            print("Error: no file specified")
            
    elif action=="avg":
        indices = helpers.ask_for_indices(data_handler.delta_As)
        if len(indices)>0:
            datas=[]
            for index in indices:
                datas.append(data_handler.delta_As[index])
            data_handler.add_data(DataObject.average(datas))
        
    elif action=="shift time":
        uniform_time_shift = helpers.ask_value(float, default=0)
        data_handler.apply(True,"time_shift",{"shift_time": uniform_time_shift})
        
    elif action=="cut w":
        print("cutting wavelength range...") # ask user nan/delete
        cut_min, cut_max = helpers.ask_range(float)
        data_handler.apply(True,"cut_w",{"cut_min": cut_min, "cut_max": cut_max})
        
    elif action=="spikes":
        print("removing spikes...")
        width = helpers.ask_value(int,"width for spikes")
        factor = helpers.ask_value(float,"number of standard deviations allowed")
        data_handler.apply(True,"remove_spikes",{"width": width, "factor": factor})
        
    elif action=="nan w":
        print("removing nan wavelength spectra...")
        data_handler.apply(True,"remove_nan_w",{})
        
    elif action=="nan t":
        print("removing nan time spectra...")
        data_handler.apply(True,"remove_nan_t",{})
    
    elif action=="fit ours":
        print("fitting our function to data...")
        
        ref_index = helpers.ask_which_layer(data_handler.reference_surfaces)
        if ref_index!=None:
            t_range = helpers.ask_range(int,default=t_range_default,add_text="Specify time range to fit")
            w_range = helpers.ask_range(int,default=w_range_default,add_text="Specify range of wavelength to explore") #360-570
            data_handler.fitXPM(ref_index, t_range, w_range)
    
    elif action=="chirp":
        print("performing chirp correction...")
        try:
            data_handler.apply(True,"chirp_correction",{"func": data_handler.t0_func, "popt": data_handler.t0_popt})
        except Exception as e:
            print("Error: Make sure to run fit first to define t0 parameters")    
            print(e)
    
    # want to perform background correction for data and reference
    elif action=="background":
        print("performing background correction...")
        back_min, back_max = helpers.ask_range(float, default=(-100,-.5))
        data_handler.apply(True,"background_correction",{"back_min":back_min, "back_max":back_max})
    
    elif action=="background ref":
        print("performing background correction on reference surface...")
        back_min, back_max = helpers.ask_range(float, default=(-100,-.5))
        data_handler.apply(False,"background_correction",{"back_min":back_min, "back_max":back_max})
        
    elif action=="reset data":
        print("reseting to original data...")
        data_handler.apply(True,"reset_data",{})
    
    elif action=="reset ref":
        print("reseting to original ref...")
        data_handler.apply(False,"reset_data",{})
        
    elif action=="switch data ref":
        print("switching data and reference for display purposes...")
        data_handler.switch_data_ref()
        print("Data/Ref inverted" if data_handler.switched_data_ref else "Data/Ref NOT inverted")
        
    elif action=="fit rate model":
        print("fitting rate model...")
        folder_name = "rate model fitting (chirp_corr, subtracted)/model1"
        
        w_min, w_max = helpers.ask_range(float, default=(480,600), add_text="Specify wavelength range to fit")
        interval = helpers.ask_value(float, default=10, label="interval for wavelength")
        t_min, t_max = helpers.ask_range(float, default=(-1,1.5), add_text="Specify time range to fit")
        
        if w_min!=None and w_max!=None and interval!=None and t_min!=None and t_max!=None:
            fit_params1 = data_handler.apply_one(True,"fitRateModel1",{"w_min":w_min,"w_max":w_max,"t_min":t_min,"t_max":t_max,"interval":interval, "folder_name": folder_name})
            DataHandler.write_out_params(folder_name+"/params1.txt", fit_params1)
        else:
            print("Improper inputs")
    
    elif action=="fit rate model 2":
        print("fitting rate model...")
        folder_name = "rate model fitting (chirp_corr, subtracted)/model2"
        
        w_min, w_max = helpers.ask_range(float, default=(490,600), add_text="Specify wavelength range to fit")
        interval = helpers.ask_value(float, default=10, label="interval for wavelength")
        t_min, t_max = helpers.ask_range(float, default=(-2,5), add_text="Specify time range to fit")
        
        if w_min!=None and w_max!=None and interval!=None and t_min!=None and t_max!=None:
            ref = data_handler.reference_surfaces[1]
            fit_params1 = data_handler.apply_one(True,"fitRateModel2",{"w_min":w_min,"w_max":w_max,"t_min":t_min,"t_max":t_max,"interval":interval, "ref": ref, "folder_name": folder_name})
            DataHandler.write_out_params(folder_name+"/params2.txt", fit_params1)
        else:
            print("Improper inputs")
    
    elif action=="fit rate model 3":
        print("fitting rate model...")
        folder_name = "rate model fitting (chirp_corr, subtracted)/model3"
        
        w_min, w_max = helpers.ask_range(float, default=(480,600), add_text="Specify wavelength range to fit")
        interval = helpers.ask_value(float, default=10, label="interval for wavelength")
        t_min, t_max = helpers.ask_range(float, default=(-2,5), add_text="Specify time range to fit")
        
        if w_min!=None and w_max!=None and interval!=None and t_min!=None and t_max!=None:
            ref = data_handler.reference_surfaces[1]
            fit_params1 = data_handler.apply_one(True,"fitRateModel3",{"w_min":w_min,"w_max":w_max,"t_min":t_min,"t_max":t_max,"interval":interval, "ref": ref, "folder_name": folder_name})
            DataHandler.write_out_params(folder_name+"/params3.txt", fit_params1)
        else:
            print("Improper inputs")
    
    ### SHORTCUTS
    elif action=="bfcs":
        print("This is a shortcut: Make sure you hard code proper constants or will crash")
        
        try:
            #apply background to data and reference
            data_handler.apply_one(True,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=0)
            data_handler.apply_one(False,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=1)
            
            #fit XPM
            data_handler.fitXPM(ref_index_default,t_range_default,w_range_default,skip_plot_prompt=True)
            
            #apply chirp correction to both data and reference
            data_handler.apply_one(True,"chirp_correction",{"func": data_handler.t0_func, "popt": data_handler.t0_popt},default=0)
            data_handler.switch_data_ref()
            data_handler.apply_one(True,"chirp_correction",{"func": data_handler.t0_func, "popt": data_handler.t0_popt},default=1)
            data_handler.switch_data_ref()
            
            #subtract
            surface_to_subtract = data_handler.reference_surfaces[ref_index_default]
            data_handler.apply_one(True,"subtract_surface",{"surface_to_subtract": surface_to_subtract, "Es": Es_default, "Er": Er_default, "f": f_default},default=0)
            
        except Exception as e:
            print("ERROR: Crashed. check constants defined in code for shortcuts")
            print(e)
    
    elif action=="bfc":
        print("This is a shortcut: Make sure you hard code proper constants or will crash")
        
        try:
            #apply background to data and reference
            data_handler.apply_one(True,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=0)
            data_handler.apply_one(False,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=1)
            
            #fit XPM
            data_handler.fitXPM(ref_index_default,t_range_default,w_range_default,skip_plot_prompt=True)
            
            data_handler.delta_As[0].original_signal = data_handler.delta_As[0].signal
            
            #apply chirp correction to both data and reference
            data_handler.apply_one(True,"chirp_correction",{"func": data_handler.t0_func, "popt": data_handler.t0_popt},default=0)
            data_handler.switch_data_ref()
            data_handler.apply_one(True,"chirp_correction",{"func": data_handler.t0_func, "popt": data_handler.t0_popt},default=1)
            data_handler.switch_data_ref()
            
        except Exception as e:
            print("ERROR: Crashed. check constants defined in code for shortcuts")
            print(e)
        
        
    elif action=="bfs":
        print("This is a shortcut: Make sure you hard code proper constants or will crash")
        
        try:
            #apply background to data and reference
            data_handler.apply_one(True,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=0)
            data_handler.apply_one(False,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=1)
            
            #fit XPM
            data_handler.fitXPM(ref_index_default,t_range_default,w_range_default,skip_plot_prompt=True)
            
            #subtract
            surface_to_subtract = data_handler.reference_surfaces[ref_index_default]
            data_handler.apply_one(True,"subtract_surface",{"surface_to_subtract": surface_to_subtract, "Es": Es_default, "Er": Er_default, "f": f_default},default=0)
            
        except Exception as e:
            print("ERROR: Crashed. check constants defined in code for shortcuts")
            print(e)
            
    elif action=="bf":
        print("This is a shortcut: Make sure you hard code proper constants or will crash")
        
        try:
            #apply background to data and reference
            data_handler.apply_one(True,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=0)
            data_handler.apply_one(False,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=1)
            
            #fit XPM
            data_handler.fitXPM(ref_index_default,t_range_default,w_range_default,skip_plot_prompt=True)
            
        except Exception as e:
            print("ERROR: Crashed. check constants defined in code for shortcuts")
            print(e)
        
    else:
        print("error, did not recognize command")
    