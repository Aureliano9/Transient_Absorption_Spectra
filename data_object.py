# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 10:29:55 2022

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
import helpers


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

class DataObject:
    name = ""
    signal = np.array([])
    times = np.array([])
    wavelengths = np.array([])
    current_w = None
    current_t = None
    w_bounds = (None,None)
    t_bounds = (None,None)
    c_bounds = (None,None)
    original_signal = np.array([])
    original_wavelengths = np.array([])
    original_times = np.array([])       
    def __init__(self,wavelengths_in,times_in,signal_in,name):
        self.name = name
        self.signal = signal_in
        self.times = times_in
        self.wavelengths = wavelengths_in
        self.current_w = helpers.avg(self.wavelengths.min(),self.wavelengths.max())
        self.current_t = helpers.avg(self.times.min(),self.times.max())
        self.w_bounds = (np.min(self.wavelengths),np.max(self.wavelengths))
        self.t_bounds = (np.min(self.times),np.max(self.times))
        self.c_bounds = (None,None)
        
        # keep copy of original that will never be mutated
        self.original_signal = np.copy(self.signal)
        self.original_wavelengths = np.copy(self.wavelengths)
        self.original_times = np.copy(self.times)
        
        # check dimensions for wavelengths and times match signal
        if self.signal.shape[0]!=len(self.times) or self.signal.shape[1]!=len(self.wavelengths):
            print("The dimension of wavelength/time does not match signal")
            os.abort()
    def CreateFromFile(filename):
        wavelengths_in, times_in, signal_in = helpers.read_file(filename)
        return DataObject(wavelengths_in, times_in, signal_in,filename)
    def average(datas):
        if len(datas)==0:
            print("Error: No data")
            return
        axes_match = True
        for i in range(1,len(datas)):
            if np.any(~np.isclose(datas[0].wavelengths,datas[i].wavelengths)) or np.any(~np.isclose(datas[0].times,datas[1].times)):
                axes_match = False
        if axes_match:
            sumData = datas[0].signal
            name = "Avg ("
            for i in range(1,len(datas)):
                sumData += datas[i].get_signal()
                name += datas[i].get_name()
            name += ")"
            w = datas[0].wavelengths
            t = datas[0].times
            signal = sumData/len(datas)
            return DataObject(w, t, signal, name)
        else:
            print("Error: Axes do not match for average candidates")
    def get_name(self):
        return self.name
    def get_w(self):
        return self.wavelengths
    def get_t(self):
        return self.times
    def get_signal(self):
        return self.signal
    def find_t_peak(self):
        # return the time at which the max occurs in the range of t_bounds
        time_min_index = helpers.find_index(self.times,self.t_bounds[0])
        time_max_index = helpers.find_index(self.times,self.t_bounds[1])
        wavelength_index = helpers.find_index(self.wavelengths,self.current_w)
        peak_index = time_min_index + np.argmax(self.signal[time_min_index:time_max_index,wavelength_index])
        return self.times[peak_index]
    def plot_color(self):
        plt.figure()
        if self.c_bounds==None or self.c_bounds[0]==None or self.c_bounds[1]==None:
            # vmin=.001, vmax=.004
            avg = np.nanmean(self.signal)
            std = np.nanstd(self.signal)
            color_width_std = .35
            c_bounds = [avg - color_width_std * std, avg + color_width_std * std]
        plt.pcolor(self.wavelengths, self.times, self.signal, vmin=c_bounds[0], vmax=c_bounds[1], cmap=my_cmap)
        plt.title(self.name)
        plt.xlabel("Wavelength")
        plt.ylabel("Time")
        plt.colorbar()
        plt.plot([self.current_w,self.current_w],[self.times.min(),self.times.max()])
        plt.plot([self.wavelengths.min(),self.wavelengths.max()],[self.current_t, self.current_t])
        plt.axis([self.w_bounds[0], self.w_bounds[1], self.t_bounds[0], self.t_bounds[1]])
        plt.show()
    def plot_wavelength_crosssection(self,label=None):
        plt.figure()
        cut_index = helpers.find_index(self.times, self.current_t)
        plt.plot(self.wavelengths, self.signal[cut_index, :])
        plt.xlabel("Wavelength")
        plt.title(self.name + ": Time = " + str(self.times[cut_index]))
        plt.xlim(self.w_bounds[0],self.w_bounds[1])
        if label!=None:
            plt.text(0,0,label)
        plt.ylabel("\Delta A")
        plt.show()
    def plot_time_crosssection(self,label=None):
        plt.figure()
        cut_index = helpers.find_index(self.wavelengths, self.current_w)
        plt.plot(self.times, self.signal[:, cut_index])
        plt.xlabel("Time")
        plt.title(self.name + ": Wavelength = " + str(self.wavelengths[cut_index]))
        plt.xlim(self.t_bounds[0],self.t_bounds[1])
        if label!=None:
            plt.text(0,0,label)
        plt.ylabel("\Delta A")
        plt.show()
    def change_waxis(self, new_range):
        # accepts None, (None,#), (#,None), (#,#)
        # None will be replaced with current bound
        if new_range!=None:
            low = new_range[0]
            high = new_range[1]
            if low==None:
                low=self.w_bounds[0]
            if high==None:
                high=self.w_bounds[1]
            self.w_bounds = (low,high)
    def change_taxis(self, new_range):
        # accepts None, (None,#), (#,None), (#,#)
        # None will be replaced with current bound
        if new_range!=None:
            low = new_range[0]
            high = new_range[1]
            if low==None:
                low=self.t_bounds[0]
            if high==None:
                high=self.t_bounds[1]
            self.t_bounds = (low,high)
    def change_caxis(self, new_range):
        # accepts None, (None,#), (#,None), (#,#)
        # None will be replaced with current bound
        if new_range!=None:
            low = new_range[0]
            high = new_range[1]
            if low==None:
                low=self.c_bounds[0]
            if high==None:
                high=self.c_bounds[1]
            self.c_bounds = (low,high)
    def change_current_t(self, new_t):
        self.current_t = new_t
    def change_current_w(self, new_w):
        self.current_w = new_w
    def reset_axis(self, which_axis):
        if which_axis==0: #waxis
            self.w_bounds = (np.min(self.wavelengths),np.max(self.wavelengths))
        elif which_axis==1: #taxis
            self.t_bounds = (np.min(self.times),np.max(self.times))
        elif which_axis==2: #caxis
            self.c_bounds = (None,None)
        elif which_axis==3: #all
            self.w_bounds = (np.min(self.wavelengths),np.max(self.wavelengths))
            self.t_bounds = (np.min(self.times),np.max(self.times))
            self.c_bounds = (None,None)
        else:
            print("Error: Value must be 0,1,2,3")
    def reset_data(self):
        self.signal = self.original_signal
        self.times = self.original_times
        self.wavelengths = self.original_wavelengths
        self.current_w = helpers.avg(self.wavelengths.min(),self.wavelengths.max())
        self.current_t = helpers.avg(self.times.min(),self.times.max())
        self.w_bounds = (np.min(self.wavelengths),np.max(self.wavelengths))
        self.t_bounds = (np.min(self.times),np.max(self.times))
        self.c_bounds = (None,None)
    
    def remove_spikes(self, width, factor):
        half_width = int(width/2);
        for t in range(self.signal.shape[0]):
            num_w = self.signal.shape[1]
            for w in range(half_width,num_w-half_width):
                mean = np.nanmean(self.signal[t,w-half_width:w+half_width])
                std = np.nanstd(self.signal[t,w-half_width:w+half_width])
                if (not np.isnan(mean) and abs(self.signal[t,w]-mean)>factor*std):
                    self.signal[t,w] = float("nan")
    def shift_time(self,uniform_time_shift):
        self.times += uniform_time_shift
        self.time_bounds = (self.time_bounds[0]+uniform_time_shift,self.time_bounds[1]+uniform_time_shift)
    def subtract_surface(self,surface_to_subtract, Er, Es, f):
        if (len(surface_to_subtract.wavelengths)!=len(self.wavelengths)):
            print("Error: Reference surface must have same dimension for wavelength")
        else:
            # subtract surface
            min_time = min(self.get_t())
            max_time = max(self.get_t())
            for i in range(len(surface_to_subtract.times)):
                subtract_time = surface_to_subtract.get_t()[i]
                if subtract_time>=min_time and subtract_time<=max_time: # check if subtract surface is out of bounds for self
                    delta_A_index = helpers.find_index(self.times, subtract_time)
                    self.signal[delta_A_index,:] -= (Es*f/Er) * surface_to_subtract.signal[i,:]
    def cut_w(self,cut_min, cut_max):
        cut_min_index = helpers.find_index(self.wavelengths, cut_min)
        cut_max_index = helpers.find_index(self.wavelengths, cut_max)
        nan_flag = input("replace with nan? otherwise, will delete spectrum. (y/n)")
        if nan_flag=="y":
            self.signal[:,cut_min_index:cut_max_index] = np.NaN
        else:
            self.signal = np.concatenate((self.signal[:,:cut_min_index],self.signal[:,cut_max_index:]), axis=1)
            self.wavelengths = np.concatenate((self.wavelengths[:cut_min_index],self.wavelengths[cut_max_index:]))        
    def remove_nan_t(self):
        remove = []
        for t in range(len(self.times)):
            if np.any(np.isnan(self.signal[t,:])):
                remove.append(t)
        remove.reverse()
        for t in remove:
            if t==len(self.times)-1:
                self.signal = self.signal[:t,:]
                self.times = self.times[:t]
            else:
                self.signal = np.concatenate((self.signal[:t,:],self.signal[t+1:,:]),axis=0)
                self.times = np.concatenate((self.times[:t],self.times[t+1:]))
        print("New dimension:", self.signal.shape)
    def remove_nan_w(self):
        remove = []
        for w in range(len(self.wavelengths)):
            if np.any(np.isnan(self.signal[:,w])):
                remove.append(w)
        remove.reverse()
        for w in remove:
            if w==len(self.wavelengths)-1:
                self.signal = self.signal[:,:w]
                self.wavelengths = self.wavelengths[:w]
            else:
                self.signal = np.concatenate((self.signal[:,:w],self.signal[:,w+1:]),axis=1)
                self.wavelengths = np.concatenate((self.wavelengths[:w],self.wavelengths[w+1:]))
        print("New dimension is " + str(self.signal.shape))
    def background_correction(self, back_min, back_max):
        # background correction for data
        background_index = (helpers.find_index(self.times,back_min),helpers.find_index(self.times,back_max))
        self.signal -= np.nanmean(self.signal[background_index[0]:background_index[1],:],axis=0)
    def play_over_time(self):
        original_time = self.current_t
        for time in self.get_t():
            self.current_t = time
            self.plot_wavelength_crosssection()
        self.current_t = original_time