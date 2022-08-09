# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 16:14:11 2022

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

# READ IN DATA
filename = "Four approaches for XPM treatment/Acetonitrile_scan1_RAW_pumped signal.dat"
wavelengths, times, delta_A = helpers.read_file(filename)
helpers.plot_color(wavelengths, times, delta_A, min(wavelengths), min(times))
