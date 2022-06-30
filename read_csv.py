from turtle import color
import wave
import numpy as np
import matplotlib.pyplot as plt
import csv
import matplotlib.colors as mcolors
import matplotlib.cm as cm

def read_csv(filename):
    f = open(filename)
    csvreader = csv.reader(f)

    row_counter=0
    times = []
    wavelengths = []
    signal = []
    for row in csvreader:
        if len(row)==1:
            break
        col_counter=0
        signal_row = []
        for element in row:
            element = float(element)
            if row_counter==0 and col_counter==0:
                pass #ignore
            elif row_counter==0:
                times.append(element)
            elif col_counter==0:
                wavelengths.append(element)
            else:
                signal_row.append(element)
            col_counter+=1
        if row_counter!=0:
            signal.append(signal_row)
        row_counter+=1
    signal = np.array(signal)
    signal = signal.transpose()
    times = np.array(times)
    wavelengths = np.array(wavelengths)
    return signal, times, wavelengths

pump_off, times, wavelengths = read_csv("continuum spectrum_channel1 raw Pump off.csv")
pump_on, times, wavelengths = read_csv("continuum spectrum_channel1 raw Pump on.csv")

delta_A = np.log(pump_on/pump_off)

def plot_2D(wavelengths,times,delta_A):
    plt.figure()
    # cMap = ListedColormap(['white', 'green', 'blue','red'])
    plt.pcolor(wavelengths,times,delta_A)
    plt.colorbar()
    plt.show()

def plot_wavelength_CS(wavelengths,times,time_index,delta_A):
    plt.figure()
    plt.plot(wavelengths,delta_A[time_index,:])
    plt.title("Time = " + str(times[time_index]))
    plt.show()

def plot_time_CS(wavelengths,times,wavelength_index,delta_A):
    plt.figure()
    plt.plot(times,delta_A[:,wavelength_index])
    plt.title("Wavelength = " + str(wavelengths[wavelength_index]))
    plt.show()

while True:
    user = input("Action?\n")
    if user=="2d":
        print("plotting 2D...")
        plot_2D(wavelengths,times,delta_A)
    elif user=="cs":
        print("plotting CS...")
        wavelength_or_time = input("Wavelength (w) or Time (t) CS?\n")
        if wavelength_or_time=="w":
            time_req = float(input("Time to take CS of?\n"))
            time_index = np.argmin(np.abs(times-time_req))
            plot_wavelength_CS(wavelengths,times,time_index,delta_A)
        elif wavelength_or_time=="t":
            wavelength_req = float(input("Wavelength to take CS of?\n"))
            wavelength_index = np.argmin(np.abs(wavelengths-wavelength_req))
            plot_time_CS(wavelengths,times,wavelength_index,delta_A)
    elif user=="play":
        print("running through all wavelength CS...")
        for i in range(len(times)):
            plot_wavelength_CS(wavelengths,times,i,delta_A)
    elif user=="q":
        print("quitting...")
        break
    else:
        print("wrong command...")
        
    print("\n")

    plt.close('all')


### chirp correction
# every wavelength, plot time trace. fit time trace. heavyside step function: x shift -> time, magnitude
# time shift vs wavelength plot (2nd order polynomial f(x) = ax^2+bx+c)

# use f(x) to figure out time shift for all wavelengths,
# apply this to 2d plot for all wavelengths