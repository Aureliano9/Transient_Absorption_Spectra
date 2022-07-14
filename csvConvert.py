# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 13:27:00 2022

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

def convertDatToCsv(filename_dat):
    filename_csv = filename_dat[:-4] + ".csv"
    with open(filename_csv, mode ='w', newline='\n')as csv_file:
        csvwriter = csv.writer(csv_file, dialect='excel')
        with open(filename_dat, mode ='r') as dat_file:
           
          lines = dat_file.readlines()
          # displaying the contents of the CSV file
          for line in lines:
              row = []
              for el in line.split('\t'):
                  row.append(el)
              csvwriter.writerow(["hi","hello"])
          
convertDatToCsv("Four approaches for XPM treatment/Acetonitrile_scan1_RAW_pumped signal.dat")