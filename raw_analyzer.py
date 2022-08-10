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

plt.ion()

# PARAMETERS
# if supplying pump off/on files
pump_off_filename = None
pump_on_filename = None
# if supplying delta A file
delta_A_filenames = ["Data (1)/CudmpDPEphosBF4ACN_1_scan1.csv", "Data (1)/CudmpDPEphosBF4ACN_1_scan2.csv", "Four approaches for XPM treatment/CudmpDPEphosBF4ACN_1_scan1.csv"]  # ["Four approaches for XPM treatment/CudmpDPEphosBF4ACN_1_scan1.csv"] #["sample1/CudmpDPEphosBF4ACN_1_scan2.csv","sample1/CudmpDPEphosBF4ACN_1_scan3.csv","sample1/CudmpDPEphosBF4ACN_1_scan4.csv"]
ref_surface_filenames = ["Data (1)/Acetonitrile.csv","Data (1)/Acetonitrile2.csv"] # "Four approaches for XPM treatment/Acetonitrile2_scan1.csv", 

if len(delta_A_filenames)==0:
    print("Error: No data to read! Change variable delta_A_filenames")

data_handler = DataHandler(delta_A_filenames, ref_surface_filenames)

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
menu["fit rate model"] = "fit rate model assuming chirp correction and XPM subtraction"
menu["fit rate model 2"] = "fit rate model assuming XPM subtraction but no chirp correction"
menu["fit rate model 3"] = "fit rate model assuming neither chirp correction or XPM subtraction"
menu["bfcs"] = "shortcut: background correction, fit ours, chirp correction, subtract XPM"
menu["bfc"] = "shortcut: background correction, fit ours, chirp correction"
menu["bfs"] = "shortcut: background correction, fit ours, subtract XPM"
menu["bf"] = "shortcut: background correction, fit ours"

plt.close('all')

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
        data_handler.apply(True,"plot_wavelength_crosssection", {})
        
    elif action=="tp":
        print("plotting time plot...")
        data_handler.apply(True,"plot_time_crosssection", {})
        
    elif action=="cp":
        print("plotting color plot...")
        data_handler.apply(True,"plot_color", {})
        
    elif action=="waxis":
        print("changing wavelength axis...")
        new_range = helpers.ask_range(float, "Specify range of wavelength to display")
        data_handler.apply(True,"change_waxis",{'new_range': new_range})
        
    elif action=="taxis":
        print("changing time axis...")
        new_range = helpers.ask_range(float, "Specify range of time to display")
        data_handler.apply(True,"change_taxis",{'new_range': new_range})
        
    elif action=="caxis":
        print("changing color axis...")
        new_range = helpers.ask_range(float, "Specify color range to display")
        data_handler.apply(True,"change_caxis",{'new_range': new_range})
        
    elif action=="reset axis":
        print("reseting axis...")
        which_axis = helpers.ask_value(int, default=None, override_text="Which axis to reset? (0=wavelength, 1=time, 2=color, 3=all) ")
        data_handler.apply(True,"reset_axis",{"which_axis": which_axis})
        
    elif action=="play":
        print("playing with time...")
        data_handler.apply(True,"play_over_time",{})
        
    elif action=="t peak":
        print("Find peak in given range")
        peak_time = data_handler.apply(True,"find_t_peak",{})
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
        cut_min, cut_max = helpers.ask_range(float, "Specify range to cut out")
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
        
    elif action=="fit rate model 1":
        print("fitting rate model...")
        folder_name = "rate model fitting (chirp_corr, subtracted)/model1"
        
        w_min, w_max = helpers.ask_range(float, default=(480,600), add_text="Specify wavelength range to fit")
        interval = helpers.ask_value(float, default=10, label="interval for wavelength")
        t_min, t_max = helpers.ask_range(float, default=(-1,1.5), add_text="Specify time range to fit")
        
        if w_min!=None and w_max!=None and interval!=None and t_min!=None and t_max!=None:
            fit_params1 = data_handler.apply(True,"fit_rate_model1",{"w_min":w_min,"w_max":w_max,"t_min":t_min,"t_max":t_max,"interval":interval, "folder_name": folder_name})
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
            fit_params1 = data_handler.apply(True,"fit_rate_model2",{"w_min":w_min,"w_max":w_max,"t_min":t_min,"t_max":t_max,"interval":interval, "ref": ref, "folder_name": folder_name})
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
            fit_params1 = data_handler.apply(True,"fit_rate_model3",{"w_min":w_min,"w_max":w_max,"t_min":t_min,"t_max":t_max,"interval":interval, "ref": ref, "folder_name": folder_name})
            DataHandler.write_out_params(folder_name+"/params3.txt", fit_params1)
            
        else:
            print("Improper inputs")
    
    ### SHORTCUTS
    elif action=="bfcs":
        print("This is a shortcut: Make sure you hard code proper constants or will crash")
        
        try:
            #apply background to data and reference
            data_handler.apply(True,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=0)
            data_handler.apply(False,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=1)
            
            #fit XPM
            data_handler.fitXPM(ref_index_default,t_range_default,w_range_default,skip_plot_prompt=True)
            
            #apply chirp correction to both data and reference
            data_handler.apply(True,"chirp_correction",{"func": data_handler.t0_func, "popt": data_handler.t0_popt},default=0)
            data_handler.switch_data_ref()
            data_handler.apply(True,"chirp_correction",{"func": data_handler.t0_func, "popt": data_handler.t0_popt},default=1)
            data_handler.switch_data_ref()
            
            #subtract
            surface_to_subtract = data_handler.reference_surfaces[ref_index_default]
            data_handler.apply(True,"subtract_surface",{"surface_to_subtract": surface_to_subtract, "Es": Es_default, "Er": Er_default, "f": f_default},default=0)
            
        except Exception as e:
            print("ERROR: Crashed. check constants defined in code for shortcuts")
            print(e)
    
    elif action=="bfc":
        print("This is a shortcut: Make sure you hard code proper constants or will crash")
        
        try:
            #apply background to data and reference
            data_handler.apply(True,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=0)
            data_handler.apply(False,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=1)
            
            #fit XPM
            data_handler.fitXPM(ref_index_default,t_range_default,w_range_default,skip_plot_prompt=True)
            
            data_handler.delta_As[0].original_signal = data_handler.delta_As[0].signal
            
            #apply chirp correction to both data and reference
            data_handler.apply(True,"chirp_correction",{"func": data_handler.t0_func, "popt": data_handler.t0_popt},default=0)
            data_handler.switch_data_ref()
            data_handler.apply(True,"chirp_correction",{"func": data_handler.t0_func, "popt": data_handler.t0_popt},default=1)
            data_handler.switch_data_ref()
            
        except Exception as e:
            print("ERROR: Crashed. check constants defined in code for shortcuts")
            print(e)
        
        
    elif action=="bfs":
        print("This is a shortcut: Make sure you hard code proper constants or will crash")
        
        try:
            #apply background to data and reference
            data_handler.apply(True,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=0)
            data_handler.apply(False,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=1)
            
            #fit XPM
            data_handler.fitXPM(ref_index_default,t_range_default,w_range_default,skip_plot_prompt=True)
            
            #subtract
            surface_to_subtract = data_handler.reference_surfaces[ref_index_default]
            data_handler.apply(True,"subtract_surface",{"surface_to_subtract": surface_to_subtract, "Es": Es_default, "Er": Er_default, "f": f_default},default=0)
            
        except Exception as e:
            print("ERROR: Crashed. check constants defined in code for shortcuts")
            print(e)
            
    elif action=="bf":
        print("This is a shortcut: Make sure you hard code proper constants or will crash")
        
        try:
            #apply background to data and reference
            data_handler.apply(True,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=0)
            data_handler.apply(False,"background_correction",{"back_min":back_min_default, "back_max":back_max_default},default=1)
            
            #fit XPM
            data_handler.fitXPM(ref_index_default,t_range_default,w_range_default,skip_plot_prompt=True)
            
        except Exception as e:
            print("ERROR: Crashed. check constants defined in code for shortcuts")
            print(e)
        
    else:
        print("error, did not recognize command")