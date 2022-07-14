# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 17:27:14 2022

@author: lennak
"""


import csv
import numpy as np
import matplotlib.pyplot as plt
import math

# return x,y from file
def readFile(filename, numberOfLinesToSkip):
    x = []
    y = []
     
    # opening the CSV file
    with open(filename, mode ='r')as file:
       
      # reading the CSV file
      csvFile = csv.reader(file)
     
      # displaying the contents of the CSV file
      counter = 0
      for lines in csvFile:
          if counter>=numberOfLinesToSkip:
              if lines[0]=='':
                  break
              x.append(float(lines[0]))
              y.append(float(lines[1]))
          counter += 1
    
    x = np.array(x)
    y = np.array(y)
    
    return (x,y)

# find width from data
def find_index_e_sq(x,y):
    #find index with max y
    max_index = np.argmax(y)
    
    #find diff with e squared
    e_sq_val = y[max_index]/7.3890561
    diff = np.abs(y-e_sq_val)
    
    #find first index
    first_index = np.argmin(diff)
    
    #find second index
    diff[first_index] = float('inf') # no longer need this diff value anymore
    second_index = np.argmin(diff)
    
    return (max_index, first_index, second_index)

# find theoretical width
def gaussian_width(a,b,c,d):
    return 2*math.sqrt(2*d**2*(math.log(b-a)-math.log(b/math.e**2-a)))

############################################

## ENTER INFO ##
#Fill Filename
filename = 'ProbeVSHorizontal.csv'
numberOfLinesToSkip = 5
#Gaussian Fit Values
a = 6.52944
b = 198.49146
c = 218.93036
d = 18.20977

#Set Flags
gaussian_flag = False
raw_e_sq_flag = False
subtract_noise_flag = True
noise_range = 20

# read file
x,y = readFile(filename, numberOfLinesToSkip)

if subtract_noise_flag:
    noise_constant = np.average(y[0:noise_range])
    print("Noise Constant:", noise_constant)
    y -= noise_constant
    
# plot
plt.figure()
plt.plot(x,y)

# width from raw data
max_index, first_index, second_index = find_index_e_sq(x,y)
plt.plot([x[max_index]],y[max_index],'ro')
plt.plot([x[first_index]],y[first_index],'ro')
plt.plot([x[second_index]],y[second_index],'ro')
print("Width from Raw Data:", abs(x[second_index]-x[first_index]), "microns")

# width from theoretical gaussian form
th_width = gaussian_width(a,b,c,d)
print("Theoretical Width:", th_width, "microns")
    
