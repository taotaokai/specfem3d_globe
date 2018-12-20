#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np 

import matplotlib as mpl
mpl.use("pdf")
import matplotlib.pyplot as plt

SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18
font = {'family': 'Times New Roman', 'size': SMALL_SIZE, }
plt.rc('font', **font)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#====== read in model
FWEA18_file = "FWEA18_average_vsv.txt"
with open(FWEA18_file, 'r') as f:
  lines = [ l.split() for l in f.readlines() if '#' not in l ]
dep_FWEA18 = np.array([float(l[0]) for l in lines])
vsv_FWEA18 = np.array([float(l[1]) for l in lines])

Yang2012_file = "Yang2012_average_vsv.txt"
with open(Yang2012_file, 'r') as f:
  lines = [ l.split() for l in f.readlines() if '#' not in l ]
dep_Yang2012 = np.array([float(l[0]) for l in lines])
vsv_Yang2012 = np.array([float(l[1]) for l in lines])

Shen2016_file = "Shen2016_average_vsv.txt"
with open(Shen2016_file, 'r') as f:
  lines = [ l.split() for l in f.readlines() if '#' not in l ]
dep_Shen2016 = np.array([float(l[0]) for l in lines])
vsv_Shen2016 = np.array([float(l[1]) for l in lines])

#====== plot
#fig = plt.figure(figsize=(11,8.5))
fig = plt.figure(figsize=(8.27, 11.7))
ax = fig.add_axes([0.1, 0.1, 0.4, 0.4])

line_FWEA18, = ax.plot(vsv_FWEA18, dep_FWEA18, 'k-', linewidth=2.0)
line_Yang2012, = ax.plot(vsv_Yang2012, dep_Yang2012, 'r-', linewidth=2.0)
line_Shen2016, = ax.plot(vsv_Shen2016, dep_Shen2016, 'b-', linewidth=2.0)

ax.set_xlim([3.0, 4.8])
ax.set_ylim([0,100])
ax.grid()
ax.set_xlabel('Vs (km/s)')
ax.set_ylabel('Depth (km)')
ax.legend([line_FWEA18, line_Yang2012, line_Shen2016], ['FWEA18_avg', 'Yang2012_avg', 'Shen2016_avg'])
ax.invert_yaxis()

#plt.show()
plt.savefig("compare_average_vsv.pdf")
