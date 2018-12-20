#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot histograms of dt,cc from misfit.txt

misfit.txt: 
#station window weight CC0 CCmax dt_cc SNR AR0 ARmax
...
"""
import sys

import numpy as np
import matplotlib as mpl
mpl.rcParams['font.family'] = "times new roman"
#label_size = 14
#mpl.rcParams['xtick.labelsize'] = label_size 
#mpl.rcParams['ytick.labelsize'] = label_size 
mpl.use("pdf")

import matplotlib.pyplot as plt
SMALL_SIZE = 16
MEDIUM_SIZE = 18
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# read command line args
misfit_file = str(sys.argv[1])
window_id = str(sys.argv[2])
out_file = str(sys.argv[3])

# read in misfit_file
with open(misfit_file, 'r') as f:
  #lines = [ l.split() for l in f.readlines() if not(l.startswith('#')) and (l[1] == window_id) ]
  lines = [ l.split() for l in f.readlines() if not(l.startswith('#'))]

lines = [l for l in lines if l[1] == window_id]

#----- create figure
fig = plt.figure(figsize=(8.5, 11)) # US Letter

# a matrix of sub plot 
nrow = 2
ncol = 1
subplot_size = np.array([1.0/ncol, 1.0/nrow])

# axis position relative to the subplot region
ax_origin_subplot = np.array([0.2, 0.2])
# size: [width, height]
ax_size_subplot = np.array([0.7, 0.7])

# hist parameters
nbins = 30
max_dt_plot = 15
min_SNR = 5
max_dt_cc = 10
# cutoff after remove linear trend
num_sigma_cutoff = 3

#----- 
gcarc = np.array([ float(l[2]) for l in lines])
dt_cc = np.array([ float(l[7]) for l in lines ])
SNR = np.array([ float(l[8]) for l in lines ])

# filter SNR, dt_cc
idx = (SNR >= min_SNR) & (np.abs(dt_cc) < max_dt_cc)
gcarc = gcarc[idx]
dt_cc = dt_cc[idx]

# remove linear trend (iterate to remove large outliers)
for iter in range(10):
  n = len(dt_cc)

  polycoef = np.polyfit(gcarc, dt_cc, 1)
  p = np.poly1d(polycoef)
  dt_cc_rtrend = dt_cc - p(gcarc)
  
  # filter dt_cc_rtrend within 3-sigma
  mean_dt_cc = np.mean(dt_cc_rtrend)
  std_dt_cc = np.std(dt_cc_rtrend)
  min_dt_cc = mean_dt_cc - num_sigma_cutoff*std_dt_cc
  max_dt_cc = mean_dt_cc + num_sigma_cutoff*std_dt_cc
  idx = (dt_cc_rtrend > min_dt_cc) & (dt_cc_rtrend < max_dt_cc)
  gcarc = gcarc[idx]
  dt_cc = dt_cc[idx]
  
  print(len(dt_cc))
  if len(dt_cc) == n: break

# plot histogram of dt_cc
nrow = 1
ncol = 1

subplot_origin = np.array([ncol-1, nrow-1])*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.hist((dt_cc, dt_cc_rtrend), nbins, histtype='step', color=('b','r'))
#ax.hist(dt_cc_rtrend, nbins, histtype='step', color='r')
ax.set_xlabel('Tobs-Tsyn (s)')
ax.set_ylabel('Window number')
ax.set_xlim([-max_dt_plot, max_dt_plot])
title_str = "%s (%.3f$\pm$%.3f s)" % (window_id, np.mean(dt_cc_rtrend), np.std(dt_cc_rtrend))
ax.set_title(title_str)

# plot histogram of dt_cc
nrow = 2
ncol = 1
subplot_origin = np.array([ncol-1, nrow-1])*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.scatter(gcarc, dt_cc, marker='o', alpha=0.5, s=7, edgecolors='none', color='b')
ax.set_xlabel('Epidistance (degree)')
ax.set_ylabel('Tobs-Tsyn (s)')
#polycoef = np.polyfit(gcarc, dt_cc, 1)
#p = np.poly1d(polycoef)
xlim = [np.min(gcarc), np.max(gcarc)]
line = ax.plot(xlim, p(xlim), 'r-')
#ax.legend(line, "%5.2f s/deg"%(polycoef[0]))
#print(polycoef)
ax.set_ylim([-max_dt_plot, max_dt_plot])
title_str = "%s (%5.2f s/deg)" % (window_id, polycoef[0])
#title_str = "%s" % (window_id)
ax.set_title(title_str)

#------ save figure
plt.savefig(out_file, format='pdf')