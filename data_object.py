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

#### NOTE ######
# data handler contains multiple data objects 
# data objects represent 2D signal for data or reference and stores independent time axes for different wavelengths
# To store these independent time axes, we will store a TimeTrace per wavelength which is composed of the time spetra and time axis


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

def my_assert(condition, msg=None):
    '''
    Aborts if condition is false
    '''
    if ~condition:
        print("Assertion false")
        if msg!=None:
            print(msg)
        os.abort()

class TimeTrace:
    # TimeTrace internal variables
    wavelength = None
    times = np.array([])
    signal = np.array([])
    
    def __init__(self,wavelength_in,times_in,signal_in):
        '''
        Constructor for TimeTrace that represents a time trace of a specific wavelength
        wavelength_in : wavelength of the time trace
        times_in : time axis (not necessarilly uniformly distributed)
        signal_in : time trace signal corresponding to times_in
        '''
        self.wavelength = wavelength_in
        self.times = times_in
        self.signal = signal_in
    
    def check_internal(self):
        '''
        asserts that internal representation is not broken
        if TimeTrace conditions are broken, system will abort

        '''
        my_assert(len(self.times)==len(self.signal), "Time trace time and signal do not have the same dimensions")
        my_assert(self.wavelength!=None and self.wavelength>=0, "Wavelength value is none or negative")
        
    def get_wavelength(self):
        '''
        return wavelength value that this time trace represents
        '''
        return self.wavelength
    def get_signal(self):
        '''
        return time trace signal as numpy array
        '''
        return self.signal
    def get_times(self):
        '''
        return time axis as numpy array
        '''
        return self.times
    def shift_time(self,shift):
        '''
        Mutator of TimeTrace
        Shift every time value by shift
        '''
        self.times += shift
    def set_times(self,new_times):
        '''
        Completely reset time axis to new_times
        '''
        self.times = new_times
        self.check_internal()
    def set_signal(self,new_signal):
        '''
        Completely reset time trace signal to new_signal
        '''
        self.signal = new_signal
        self.check_internal()

class DataObject:
    # DataObject internal variables
    name = ""
    signal = np.array([])
    times = np.array([])
    wavelengths = np.array([])
    time_traces = []
    current_w = None
    current_t = None
    w_bounds = (None,None)
    t_bounds = (None,None)
    c_bounds = (None,None)
    chirp_corrected = False
    
    def reset_time_traces(self):
        '''
        Recreate all TimeTraces in time_traces from current self.signal and self.times
        '''
        # reset time traces
        self.time_traces = []
        for w_i in range(len(self.wavelengths)):
           wavelength = self.wavelengths[w_i]
           self.time_traces.append(TimeTrace(wavelength,self.times,self.signal[:,w_i]))
    def shift_time_traces(self, shift):
        '''
        Shift all the time axes of TimeTraces by shift
        '''
        for time_trace in self.time_traces:
            time_trace.shift_time(shift)
    def time_traces_cut_wavelengths(self, min_index, max_index):
        '''
        Cut out the wavelengths from time_traces between index min_index and max_index from original self.wavelengths
        '''
        # assumes self.wavelengths already has wavelengths cut from it
        self.time_traces = self.time_traces[:min_index]+self.time_traces[max_index:]
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
        
        # create dictionary for time slices
        self.reset_time_traces()
        
        # keep copy of original that will never be mutated
        self.original_signal = np.copy(self.signal)
        self.original_wavelengths = np.copy(self.wavelengths)
        self.original_times = np.copy(self.times)
        
        # check dimensions for wavelengths and times match signal
        if self.signal.shape[0]!=len(self.times) or self.signal.shape[1]!=len(self.wavelengths):
            print("The dimension of wavelength/time does not match signal")
            os.abort()
            
    def create_from_file(filename):
        '''
        returns new DataObject based on the file filneame
        filename should have the format where first column represents the wavelength axis and first row represents the time axis
        all other value should be the corresponding delta A signal values corresponding to the respective wavelength row and time column
        '''
        wavelengths_in, times_in, signal_in = helpers.read_file(filename)
        return DataObject(wavelengths_in, times_in, signal_in,filename)
    def average(datas):
        '''
        create the average across the list of datas
        The wavelength axes and time axes must match. We will print an error if the axes values are not close enough.
        
        '''
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
        '''
        return the name of th
    
    
    
    
    data object which is its filename
        '''
        return self.name
    def get_w(self):
        '''
        return the wavelength axis of this data object
        '''
        return self.wavelengths
    def get_t(self):
        '''
        return the time axis of this data object
        '''
        return self.times
    def get_signal(self):
        '''
        return the signal matrix of this data object
        '''
        return self.signal
    def find_t_peak(self, wavelength=None,t_bounds=None):
        '''
        return the time at which the maximum signal occurs for wavelength within t_bounds
        if wavelength is None, assume current_w
        if t_bounds is None, assume t_axis_bounds
        '''
        # return the time at which the max occurs in the range of t_bounds
        if wavelength==None:
            wavelength = self.current_w
        w_i = helpers.find_index(self.wavelengths,wavelength)
        timeTrace = self.time_traces[w_i]
        
        if t_bounds==None:
            time_min_index = helpers.find_index(timeTrace.times,self.t_bounds[0])
            time_max_index = helpers.find_index(timeTrace.times,self.t_bounds[1])
        else:
            time_min_index = helpers.find_index(timeTrace.times,t_bounds[0])
            time_max_index = helpers.find_index(timeTrace.times,t_bounds[1])

        peak_index = time_min_index + np.argmax(timeTrace.signal[time_min_index:time_max_index])
        return timeTrace.times[peak_index]
    def plot_color(self):
        '''
        plot the 2D plot of the signal where time is vertical axis and wavelength is horizontal axis
        If chirp correction has been performed, the time axis across different wavelengths is not shared, so we will interpolate values for each wavelength spectra
        '''
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
        plt.show(block=True)
    def plot_wavelength_crosssection(self,label=None):
        '''
        plot wavelength crosssection at t=current_t
        If chirp correction was performed, all time axes are shifted by different amounts, so we will interpolate
        '''
        plt.figure()
        if self.chirp_corrected:
            interpolated_spectra = []
            for w_i in range(len(self.wavelengths)):
                time_trace = self.time_traces[w_i]
                interpolated_spectra.append(np.interp(self.current_t, time_trace.times, time_trace.signal))            
            plt.scatter(self.wavelengths, interpolated_spectra, s=1, marker=".")
            plt.title(self.name + ": Interpolated Time = " + str(self.current_t))
        else:
            cut_index = helpers.find_index(self.times, self.current_t)
            plt.plot(self.wavelengths, self.signal[cut_index, :])
            plt.title(self.name + ": Time = " + str(self.times[cut_index]))
        plt.xlabel("Wavelength")
        plt.xlim(self.w_bounds[0],self.w_bounds[1])
        if label!=None:
            plt.text(0,0,label)
        plt.ylabel("\Delta A")
        plt.show()
    def plot_time_crosssection(self,label=None):
        '''
        plot time crosssection at w=current_w
        If chirp correction was performed, the exact time axis will be used, no interpolation
        '''
        plt.figure()
        w_i = helpers.find_index(self.wavelengths, self.current_w)
        time_trace = self.time_traces[w_i]
        if self.chirp_corrected:
            plt.scatter(time_trace.times, time_trace.signal, s=1, marker=".")
        else:
            plt.plot(time_trace.times, time_trace.signal)
        plt.xlabel("Time")
        plt.title(self.name + ": Wavelength = " + str(self.wavelengths[w_i]))
        plt.xlim(self.t_bounds[0],self.t_bounds[1])
        if label!=None:
            plt.text(0,0,label)
        plt.ylabel("$\Delta A$")
        plt.show()
    def change_waxis(self, new_range):
        '''
        change the wavelength axis to provided new_range
        new_range should be None, (None,#), (#,None), or (#,#)
        None will be replaced with current bound values
        '''
        if new_range!=None:
            low = new_range[0]
            high = new_range[1]
            if low==None:
                low=self.w_bounds[0]
            if high==None:
                high=self.w_bounds[1]
            self.w_bounds = (low,high)
    def change_taxis(self, new_range):
        '''
        change the time axis to the provided new_range
        new_range should be None, (None,#), (#,None), or (#,#)
        None will be replaced with current bound values
        '''
        if new_range!=None:
            low = new_range[0]
            high = new_range[1]
            if low==None:
                low=self.t_bounds[0]
            if high==None:
                high=self.t_bounds[1]
            self.t_bounds = (low,high)
    def change_caxis(self, new_range):
        '''
        changes the color axis to the provided new_range
        accepts None, (None,#), (#,None), (#,#)
        None will be replaced with current bound
        '''
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
        '''
        Change current cursor location for time to new_t
        '''
        self.current_t = new_t
    def change_current_w(self, new_w):
        '''
        Change current cursor location for wavelength to new_w
        '''
        self.current_w = new_w
    def reset_axis(self, which_axis):
        '''
        Reset axis specified by which_axis
        0 -> Reset wavelength axis
        1 -> Reset time axis
        2 -> Reset color axis
        3 -> Reset all axis
        '''
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
        '''
        Reset data from original data read from file filename
        '''
        wavelengths_in, times_in, signal_in = helpers.read_file(self.filename)
        self.__init__(wavelengths_in, times_in, signal_in, self.filename)
    
    def remove_spikes(self, width, factor):
        '''
        Reset data from original data read from file filename
        '''
        assert ~self.chirp_corrected, "Chirp correction was already performed so time traces are shifted"
        
        half_width = int(width/2);
        for t in range(self.signal.shape[0]):
            num_w = self.signal.shape[1]
            for w in range(half_width,num_w-half_width):
                mean = np.nanmean(self.signal[t,w-half_width:w+half_width])
                std = np.nanstd(self.signal[t,w-half_width:w+half_width])
                if (not np.isnan(mean) and abs(self.signal[t,w]-mean)>factor*std):
                    self.signal[t,w] = float("nan")
        self.reset_time_traces()
    def shift_time(self,uniform_time_shift):
        '''
        Shift time axis by uniform_time_shift
        '''
        self.times += uniform_time_shift
        self.time_bounds = (self.time_bounds[0]+uniform_time_shift,self.time_bounds[1]+uniform_time_shift)
        self.shift_time_traces(uniform_time_shift)
    def subtract_surface(self,surface_to_subtract, Er, Es, f):
        '''
        Subtract surface_to_subtract from this surface with the factor Es*f/Er
        '''
        assert ~self.chirp_corrected, "Chirp correction was already performed so time traces are shifted"
        
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
        self.reset_time_traces()
    def cut_w(self,cut_min, cut_max):
        '''
        Remove the range cut_min and cut_max from wavelength axis (will mutate 2D signal and time_traces)
        Will prompt user whether they should be cut out or replaced Nan
        '''
        cut_min_index = helpers.find_index(self.wavelengths, cut_min)
        cut_max_index = helpers.find_index(self.wavelengths, cut_max)
        nan_flag = input("replace with nan? otherwise, will delete spectrum. (y/n)")
        if nan_flag=="y":
            self.signal[:,cut_min_index:cut_max_index] = np.NaN
        else:
            self.signal = np.concatenate((self.signal[:,:cut_min_index],self.signal[:,cut_max_index:]), axis=1)
            self.wavelengths = np.concatenate((self.wavelengths[:cut_min_index],self.wavelengths[cut_max_index:]))       
        self.time_traces_cut_wavelengths(cut_min_index,cut_max_index)
    def remove_nan_t(self):
        '''
        Remove any times with any Nan values in its spectrum
        Should be before chirp correction since time-axis won't be shared otherwise
        '''
        assert ~self.chirp_corrected, "Chirp correction was already performed so time traces are shifted"
        
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
        self.reset_time_traces()
        print("New dimension:", self.signal.shape)
    def remove_nan_w(self):
        '''
        Remove any wavelengths with any Nan values in its spectrum
        Should be before chirp correction since wavvelength-axis won't be shared otherwise
        '''
        assert ~self.chirp_corrected, "Chirp correction was already performed so time traces are shifted"
        
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
        self.reset_time_traces()
        print("New dimension is " + str(self.signal.shape))
    def background_correction(self, back_min, back_max):
        '''
        Perform background correction by taking times between back_min and back_max
        Should be performed before chirp correction (techincally code would work but not recommended)
        '''
        assert ~self.chirp_corrected, "Chirp correction was already performed so time traces are shifted"
        
        # background correction for data
        background_index = (helpers.find_index(self.times,back_min),helpers.find_index(self.times,back_max))
        self.signal -= np.nanmean(self.signal[background_index[0]:background_index[1],:],axis=0)
        self.reset_time_traces()

    def play_over_time(self):
        '''
        Display the wavelength crosssections over time
        Warning, could take long if there are many times
        '''
        original_time = self.current_t
        for time in self.get_t():
            self.current_t = time
            self.plot_wavelength_crosssection()
        self.current_t = original_time
    def chirp_correction(self, func, popt):
        '''
        Chirp correct based on the t0 function determined by fitting to reference surface
        i.e. make sure to call 'fit ours' beforehand
        '''
        new_time_min = float('inf')
        new_time_max = -float('inf')
        for w_i in range(len(self.wavelengths)):
            current_freq = helpers.convert_to_ang_freq(self.wavelengths[w_i])
            t0 = func(current_freq, *popt)
            
            current_min = np.min(self.times)-t0
            current_max = np.max(self.times)-t0
            
            if new_time_min>current_min:
                new_time_min = current_min
            if new_time_max<current_max:
                new_time_max = current_max
        
        # precision is not consistent so take the low/high ends separately and extend
        low_end_precision = abs(self.times[1]-self.times[0])
        high_end_precision = abs(self.times[-1]-self.times[-2])
        
        new_times = np.append(np.arange(new_time_min,min(self.times),low_end_precision),np.copy(self.times))
        new_times = np.append(new_times, np.arange(max(self.times),new_time_max,high_end_precision))
            
        new_signal = np.empty((len(new_times),len(self.wavelengths)))
        new_signal[:] = np.nan
        
        for w_i in range(len(self.wavelengths)):
            wavelength = self.wavelengths[w_i]
            current_freq = helpers.convert_to_ang_freq(wavelength) # rad/ps^-1
            t0 = func(current_freq, *popt)
            current_times = self.times-t0
            current_signal = self.time_traces[w_i].signal
            
            new_signal[:, w_i] = np.interp(new_times, current_times, current_signal)
            
            self.time_traces[w_i].shift_time(-t0)
        
        self.times = new_times
        self.signal = new_signal
        
        # reset all axis
        self.reset_axis(3)
        self.chirp_corrected = True
    def fit_rate_model1(self, w_min, w_max, t_min, t_max, interval, folder_name):
        '''
        fit rate model assuming current signal has been background corrected, fitted, chirp corrected, and reference-subtracted
        returns the fitted parameters and saves fitted plot if success
        prints error if fit did not succed
        '''
        outputParams = {}
        
        wavelengths_to_fit = np.arange(w_min,w_max,interval)
        
        for wavelength_to_fit in wavelengths_to_fit:
            w_i = helpers.find_index(self.wavelengths, wavelength_to_fit)
            current_time_trace = self.time_traces[w_i]
            w = current_time_trace.wavelength
            
            t_min_index = helpers.find_index(current_time_trace.times, t_min)
            t_max_index = helpers.find_index(current_time_trace.times, t_max)
    
            focus_times = current_time_trace.times[t_min_index:t_max_index]
            focus_signal = current_time_trace.signal[t_min_index:t_max_index]
            
            plt.figure()
            plt.plot(focus_times, focus_signal)
            plt.title(str(w) + " nm")
            plt.show()
            
            B2 = .8
            B3 = .15
            B4 = .05
            A1 = 0
            A2 = 2e-3
            A3 = 1.5e-3
            A4 = 6e-3
            sigma = 0.03
            try:
                popt, pcov = curve_fit(helpers.rateModel, focus_times[~np.isnan(focus_signal)], focus_signal[~np.isnan(focus_signal)], [B2, B3, B4, A1, A2, A3, A4, sigma])
                fitted = helpers.rateModel(focus_times, *popt)
                print(*popt)
                plt.figure()
                plt.plot(focus_times, focus_signal)
                plt.plot(focus_times, fitted)
                plt.title("Rate model fitted for " + str(w) + " nm")
                plt.savefig(folder_name + "/" + str(w) + ".png")
                plt.show()
            except Exception as e:
                print("Could not fit")
                print(e)
            
            outputParams[w] = popt
        
        return outputParams
    
    def fit_rate_model2(self, w_min, w_max, t_min, t_max, interval, ref, folder_name):
        '''
        fit rate model assuming current signal has been background corrected, fitted, and reference-subtracted
        assumes NOT CHIRP CORRECTED
        returns the fitted parameters and saves fitted plot if success
        prints error if fit did not succed
        '''
        
        outputParams = {}
        
        wavelengths_to_fit = np.arange(w_min,w_max,interval)
        
        for wavelength_to_fit in wavelengths_to_fit:
            w_i = helpers.find_index(self.wavelengths, wavelength_to_fit)
            current_time_trace = self.time_traces[w_i]
            w = current_time_trace.wavelength
            
            t_min_index = helpers.find_index(current_time_trace.times, t_min)
            t_max_index = helpers.find_index(current_time_trace.times, t_max)
    
            focus_times = current_time_trace.times[t_min_index:t_max_index]
            focus_signal = current_time_trace.signal[t_min_index:t_max_index]
            
            t0 = ref.find_t_peak(w,t_bounds=[t_min,t_max])
            
            # plt.figure()
            # plt.plot(focus_times, focus_signal)
            # plt.title(str(w) + " nm " + str(t0))
            # plt.show()
            
            # plt.figure()
            # plt.plot(ref.times, ref.signal[:,w_i])
            # plt.title(str(w) + " nm " + str(t0))
            # plt.show()
            
            B2 = .8
            B3 = .01
            B4 = 5.6
            A1 = -2e5
            A2 = 2e-3
            A3 = 5e-3
            A4 = 6e-3
            sigma = 4e-2
            
            # plt.figure()
            # plt.plot(focus_times, focus_signal)
            # plt.plot(focus_times, helpers.rateModel2(focus_times, B2, B3, B4, A1, A2, A3, A4, sigma, t0))
            # plt.title("Before " + str(w) + " nm. to=" + str(t0))
            # plt.show()
            # # plt.savefig(folder_name + "/" + str(w) + ".png")
            
            try:
                popt, pcov = curve_fit(helpers.rateModel2, focus_times[~np.isnan(focus_signal)], focus_signal[~np.isnan(focus_signal)], [B2, B3, B4, A1, A2, A3, A4, sigma, t0])
                fitted = helpers.rateModel2(focus_times, *popt)
                print(*popt)
                plt.figure()
                plt.plot(focus_times, focus_signal)
                plt.plot(focus_times, fitted)
                plt.title("Rate model fitted for " + str(w) + " nmm, t0 = " + str(t0))
                plt.savefig(folder_name + "/" + str(w) + ".png")
                plt.show()
            except Exception as e:
                print("Could not fit")
                # plt.figure()
                # plt.plot(focus_times, focus_signal)
                # plt.plot(focus_times, helpers.rateModel2(focus_times, B2, B3, B4, A1, A2, A3, A4, sigma, t0))
                # plt.title("Failed for " + str(w) + " nm, t0 = " + str(t0))
                # # plt.savefig(folder_name + "/" + str(w) + ".png")
                # plt.show()
                print(e)
            
            outputParams[w] = popt
        
        return outputParams
    
    def fit_rate_model3(self, w_min, w_max, t_min, t_max, interval, ref, folder_name):
        '''
        fit rate model assuming current signal has been background corrected, fitted
        assumes NOT CHIRP CORRECTED NOR REFERENCE-SUBTRACTED
        returns the fitted parameters and saves fitted plot if success
        prints error if fit did not succed
        '''
        
        outputParams = {}
        
        wavelengths_to_fit = np.arange(w_min,w_max,interval)
        
        for wavelength_to_fit in wavelengths_to_fit:
            w_i = helpers.find_index(self.wavelengths, wavelength_to_fit)
            current_time_trace = self.time_traces[w_i]
            w = current_time_trace.wavelength
            
            t_min_index = helpers.find_index(current_time_trace.times, t_min)
            t_max_index = helpers.find_index(current_time_trace.times, t_max)
    
            focus_times = current_time_trace.times[t_min_index:t_max_index]
            focus_signal = current_time_trace.signal[t_min_index:t_max_index]
            
            t0 = ref.find_t_peak(w,t_bounds=[t_min,t_max])
            
            # plt.figure()
            # plt.plot(focus_times, focus_signal)
            # plt.title(str(w) + " nm " + str(t0))
            # plt.show()
            
            # plt.figure()
            # plt.plot(ref.times, ref.signal[:,w_i])
            # plt.title(str(w) + " nm " + str(t0))
            # plt.show()
            
            B2 = .8
            B3 = .01
            B4 = 5.6
            A1 = -2e5
            A2 = 2e-3
            A3 = 5e-3
            A4 = 6e-3
            sigma = 4e-2
            
            c1 = -.0005
            c2 = 0.00005
            c3 = -1
            
            # plt.figure()
            # plt.plot(focus_times, focus_signal)
            # plt.plot(focus_times, helpers.rateModel2(focus_times, B2, B3, B4, A1, A2, A3, A4, sigma, t0))
            # plt.title("Before " + str(w) + " nm. to=" + str(t0))
            # plt.show()
            # # plt.savefig(folder_name + "/" + str(w) + ".png")
            
            try:
                popt, pcov = curve_fit(helpers.rateModel3, focus_times[~np.isnan(focus_signal)], focus_signal[~np.isnan(focus_signal)], [B2, B3, B4, A1, A2, A3, A4, sigma, t0, c1, c2, c3])
                fitted = helpers.rateModel3(focus_times, *popt)
                print(*popt)
                plt.figure()
                plt.plot(focus_times, focus_signal)
                plt.plot(focus_times, fitted)
                plt.title("Rate model fitted for " + str(w) + " nmm, t0 = " + str(t0))
                plt.savefig(folder_name + "/" + str(w) + ".png")
                plt.show()
            except Exception as e:
                print("Could not fit")
                plt.figure()
                plt.plot(focus_times, focus_signal)
                plt.plot(focus_times, helpers.rateModel2(focus_times, B2, B3, B4, A1, A2, A3, A4, sigma, t0))
                plt.title("Failed for " + str(w) + " nm, t0 = " + str(t0))
                # plt.savefig(folder_name + "/" + str(w) + ".png")
                plt.show()
                print(e)
            
            outputParams[w] = popt
        
        return outputParams

class DataHandler:
    
    # Data Handler internal variables
    delta_As = []
    reference_surfaces = []
    parameters = {}
    switched_data_ref = False
    
    # For fit
    t0_func = staticmethod(helpers.second)
    t0_popt = None
    c1_func = staticmethod(helpers.fifth)
    c1_popt = None
    c2_func = staticmethod(helpers.fifth)
    c2_popt = None
    c3_func = staticmethod(helpers.fifth)
    c3_popt = None
    tau1_func = staticmethod(helpers.fifth)
    tau1_popt = None
    
    def __init__(self, delta_A_filenames, ref_surface_filenames):
        for filename in delta_A_filenames:
            self.delta_As.append(DataObject.create_from_file(filename))
        for filename in ref_surface_filenames:
            self.reference_surfaces.append(DataObject.create_from_file(filename))
    def apply(self,is_delta_A, method, kwargs, apply_all=False):
        list_of_data = self.delta_As if is_delta_A else self.reference_surfaces
        if len(list_of_data)==1:
            return list_of_data[0].__getattribute__(method)(**kwargs)
        else:
            specify_index = False if apply_all else helpers.ask_yes_no("Apply only to a specific layer?")
            if specify_index:
                display_index = helpers.ask_which_layer(list_of_data)
                if display_index!=None:
                    return list_of_data[display_index].__getattribute__(method)(**kwargs)
            else:
                output = []
                for delta_A in list_of_data:
                    output.append(delta_A.__getattribute__(method)(**kwargs))
                return output
    def apply_one(self,is_delta_A, method, kwargs,default=None):
        list_of_data = self.delta_As if is_delta_A else self.reference_surfaces
        if len(list_of_data)==1:
            return list_of_data[0].__getattribute__(method)(**kwargs)
        elif default==None:
            display_index = helpers.ask_which_layer(list_of_data)
            if display_index!=None:
                return list_of_data[display_index].__getattribute__(method)(**kwargs)
        else:
            return list_of_data[default].__getattribute__(method)(**kwargs)
    def add_data(self,data_object):
        self.delta_As.append(data_object)
    def add_reference(self,data_object):
        self.reference_surfaces.append(data_object)
    def switch_data_ref(self):
        temp = self.delta_As
        self.delta_As = self.reference_surfaces
        self.reference_surfaces = temp
        self.switched_data_ref = ~self.switched_data_ref
    def write_out_params(name, params):
        f = open(name, "w")
        for wavelength in params:
            f.write(str(wavelength) + "," + str(params[wavelength]) + "\n")
        f.close()
    def fitXPM(self, reference_index, t_range, w_range, skip_plot_prompt=False):
        # use this function to fit XPM signal with solvent data
        
        if reference_index<0 or reference_index>=len(self.reference_surfaces):
            print("Error: reference index is not in range")
            return
        
        ref = self.reference_surfaces[reference_index]
            
        # determine time range
        min_t_index = helpers.find_index(ref.times,t_range[0])
        max_t_index = helpers.find_index(ref.times,t_range[1])
        
        # determine wavelength range
        min_w_index = helpers.find_index(ref.wavelengths,w_range[0])
        max_w_index = helpers.find_index(ref.wavelengths,w_range[1])
        min_w_index_ignore = helpers.find_index(ref.wavelengths,410)
        max_w_index_ignore = helpers.find_index(ref.wavelengths,430)
        
        if min_w_index<min_w_index_ignore and max_w_index>max_w_index_ignore:
            w_indices = list(range(min_w_index,min_w_index_ignore))
            w_indices += list(range(max_w_index_ignore,max_w_index))
        elif min_w_index<min_w_index_ignore and max_w_index<max_w_index_ignore:
            w_indices = range(min_w_index,min_w_index_ignore)
        elif min_w_index>min_w_index_ignore and max_w_index>max_w_index_ignore:
            w_indices = range(max_w_index_ignore, max_w_index)
        elif min_w_index>min_w_index_ignore and max_w_index<max_w_index_ignore:
            w_indices = range(0,0)
            print("Specify a wavelength range outside of forbidden region (410 to 430)")
        else:
            print("Error: have a bug with intervals")
        
        plot_flag = False if skip_plot_prompt else helpers.ask_yes_no("Plot all wavelengths?", default=False)
        
        for w_i in w_indices:
            wavelength = ref.wavelengths[w_i]
            sliced_times = ref.times[min_t_index:max_t_index]
            sliced_signal = ref.signal[min_t_index:max_t_index,w_i]
            
            t0 = ref.find_t_peak(wavelength,t_range)
            
            initial_c1 = -4.79227948e-04
            initial_c2 = 3.87793025e-05
            initial_c3 = -6.49229051e-06
            tau1 = 1.36241557e-01
            
            try:
                popt, pcov = curve_fit(helpers.ours, sliced_times, sliced_signal, [initial_c1, initial_c2, initial_c3, tau1, t0])
                pcov_mag = np.linalg.norm(np.array(pcov))
                if pcov_mag > .01:
                    raise RuntimeError() 
                
                # fit worked
                self.parameters[wavelength] = popt
                if plot_flag:
                    fitted = helpers.ours(sliced_times, *popt)
                    plt.figure()
                    plt.plot(sliced_times, sliced_signal)
                    plt.plot(sliced_times, fitted)
                    plt.title("$\lambda=" + str(wavelength) + " , initial\_t_0=" + str(t0) + " , final\_t_0=" + str(popt[4]) + " , cov=" + str(pcov_mag) + "$")
                    plt.show()
            except Exception:
                print("Ignoring wavelength " + str(wavelength))
    
        chosen_freqs = []
        c1s = []
        c2s = []
        c3s = []
        tau1s = []
        t0s = []
        for wavelength in self.parameters:
            omega = helpers.convert_to_ang_freq(wavelength)
            chosen_freqs.append(omega)
            c1s.append(self.parameters[wavelength][0])
            c2s.append(self.parameters[wavelength][1])
            c3s.append(self.parameters[wavelength][2])
            tau1s.append(self.parameters[wavelength][3])
            t0s.append(self.parameters[wavelength][4])
        chosen_freqs = np.array(chosen_freqs)
        c1s = np.array(c1s)
        c2s = np.array(c2s)
        c3s = np.array(c3s)
        tau1s = np.array(tau1s)
        t0s = np.array(t0s)
        
        #C1
        self.c1_popt, _ = curve_fit(self.c1_func, chosen_freqs, c1s)
        print("C1:",self.c1_popt)
        plt.figure()
        plt.plot(chosen_freqs, c1s)
        plt.plot(chosen_freqs, self.c1_func(chosen_freqs, *self.c1_popt))
        plt.title("c1")
        plt.xlabel("$rad/ps^-1$")
        plt.show()
        
        #C2
        self.c2_popt, _ = curve_fit(self.c2_func, chosen_freqs, c2s)
        print("C2:",self.c2_popt)
        plt.figure()
        plt.plot(chosen_freqs, c2s)
        plt.plot(chosen_freqs, self.c2_func(chosen_freqs, *self.c2_popt))
        plt.title("c2")
        plt.xlabel("$rad/ps^-1$")
        plt.show()
        
        #C3
        self.c3_popt, _ = curve_fit(self.c3_func, chosen_freqs, c3s)
        print("C3:",self.c3_popt)
        plt.figure()
        plt.plot(chosen_freqs, c3s)
        plt.plot(chosen_freqs, self.c3_func(chosen_freqs, *self.c3_popt))
        plt.title("c3")
        plt.xlabel("$rad/ps^-1$")
        plt.show()
        
        #tau1
        self.tau1_popt, _ = curve_fit(self.tau1_func, chosen_freqs, tau1s)
        print("tau1:",self.tau1_popt)
        plt.figure()
        plt.plot(chosen_freqs, tau1s)
        plt.plot(chosen_freqs, self.tau1_func(chosen_freqs, *self.tau1_popt))
        plt.title("tau1")
        plt.xlabel("$rad/ps^-1$")
        plt.show()
        
        #t0
        self.t0_popt, _ = curve_fit(self.t0_func, chosen_freqs, t0s)
        print("t0:",self.t0_popt)
        plt.figure()
        plt.plot(chosen_freqs, t0s)
        plt.plot(chosen_freqs, self.t0_func(chosen_freqs, *self.t0_popt))
        plt.title("t0")
        plt.xlabel("$rad/ps^-1$")
        plt.show()