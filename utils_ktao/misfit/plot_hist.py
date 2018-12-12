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
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

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
nbins = 50
max_dt = 10 
min_SNR = 5

#----- 
gcarc = np.array([ float(l[2]) for l in lines])
dt_cc = np.array([ float(l[7]) for l in lines ])
SNR = np.array([ float(l[8]) for l in lines ])

idx = SNR >= min_SNR
gcarc = gcarc[idx]
dt_cc = dt_cc[idx]

# filter dt_cc within 3-sigma
mean_dt_cc = np.mean(dt_cc)
std_dt_cc = np.std(dt_cc)
min_dt_cc = mean_dt_cc - 3*std_dt_cc
max_dt_cc = mean_dt_cc + 3*std_dt_cc
idx = (dt_cc > min_dt_cc) & (dt_cc < max_dt_cc)
gcarc = gcarc[idx]
dt_cc = dt_cc[idx]

# plot histogram of dt_cc
nrow = 1
ncol = 1
subplot_origin = np.array([ncol-1, nrow-1])*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.hist(dt_cc, nbins, histtype='step')
ax.set_xlabel('Tobs-Tsyn (s)')
ax.set_ylabel('Window number')
#ax.set_xlim([-max_dt, max_dt])
title_str = "%s (%.3f$\pm$%.3f s)" % (window_id, np.mean(dt_cc), np.std(dt_cc))
ax.set_title(title_str)

# plot histogram of dt_cc
nrow = 2
ncol = 1
subplot_origin = np.array([ncol-1, nrow-1])*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.plot(gcarc, dt_cc, 'o')
ax.set_xlabel('Epidistance (degree)')
ax.set_ylabel('Tobs-Tsyn (s)')
polycoef = np.polyfit(gcarc, dt_cc, 1)
p = np.poly1d(polycoef)
xlim = [np.min(gcarc), np.max(gcarc)]
line = ax.plot(xlim, p(xlim), 'r-')
#ax.legend(line, "%5.2f s/deg"%(polycoef[0]))
#print(polycoef)
#ax.set_xlim([-max_dt, max_dt])
title_str = "%s (%5.2f s/deg)" % (window_id, polycoef[0])
ax.set_title(title_str)

#------ save figure
plt.savefig(out_file, format='pdf')