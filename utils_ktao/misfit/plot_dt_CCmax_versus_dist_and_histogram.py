#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot histograms of dt,cc from misfit.txt

misfit.txt: 
#station window weight CC0 CCmax dt_cc SNR AR0 ARmax
...
"""
import sys
import scipy.stats

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
plot_title = str(sys.argv[4])

# read in misfit_file
with open(misfit_file, 'r') as f:
  #lines = [ l.split() for l in f.readlines() if not(l.startswith('#')) and (l[1] == window_id) ]
  lines = [ l.split() for l in f.readlines() if not(l.startswith('#'))]

lines = [l for l in lines if l[1] == window_id]

#----- create figure
fig = plt.figure(figsize=(8.5, 11)) # US Letter

# a matrix of sub plot 
nrow = 4
ncol = 1
subplot_size = np.array([1.0/ncol, 1.0/nrow])

# axis position relative to the subplot region
ax_origin_subplot = np.array([0.2, 0.2])
# size: [width, height]
ax_size_subplot = np.array([0.7, 0.7])

# hist parameters
#nbins = 30
min_cc_plot=0.5
max_dt_plot = 10
min_SNR = 5
max_dt_cc = 10
min_cc_max = 0.5
# cutoff after remove linear trend
num_sigma_cutoff = 3

#----- 
gcarc = np.array([ float(l[2]) for l in lines])
az = np.array([ float(l[3]) for l in lines])
baz = np.array([ float(l[4]) for l in lines])
weight = np.array([ float(l[5]) for l in lines])
cc0 = np.array([ float(l[6]) for l in lines ])
ccmax = np.array([ float(l[7]) for l in lines ])
dt_cc = np.array([ float(l[8]) for l in lines ])
SNR = np.array([ float(l[9]) for l in lines ])

# filter SNR
print(len(dt_cc))
idx = (SNR >= min_SNR) #& (np.abs(dt_cc) < max_dt_cc) #& (ccmax >= min_cc_max)
gcarc = gcarc[idx]
az = az[idx]
baz = baz[idx]
weight = weight[idx]
cc0 = cc0[idx]
ccmax = ccmax[idx]
dt_cc = dt_cc[idx]
SNR = SNR[idx]
print(len(dt_cc))

# remove linear trend (iterate to remove large outliers)
for i in range(10):
  n = len(dt_cc)

  idx = (np.abs(dt_cc) < max_dt_cc) #& (ccmax >= min_cc_max)
  print(len(dt_cc[idx]))
  polycoef = np.polyfit(gcarc[idx], dt_cc[idx], 1)
  #polycoef = np.polyfit(gcarc, dt_cc, 1)
  p = np.poly1d(polycoef)
  print(i, p)
  dt_cc_rtrend = dt_cc - p(gcarc)
  
  # filter dt_cc_rtrend within 3-sigma
  mean_dt_cc = np.mean(dt_cc_rtrend[idx])
  std_dt_cc = np.std(dt_cc_rtrend[idx])
  min_dt_cc = mean_dt_cc - num_sigma_cutoff*std_dt_cc
  max_dt_cc = mean_dt_cc + num_sigma_cutoff*std_dt_cc
  idx = (dt_cc_rtrend > min_dt_cc) & (dt_cc_rtrend < max_dt_cc)

  gcarc = gcarc[idx]
  az = az[idx]
  baz = baz[idx]
  weight = weight[idx]
  cc0 = cc0[idx]
  ccmax = ccmax[idx]
  dt_cc = dt_cc[idx]
  SNR = SNR[idx]
  dt_cc_rtrend = dt_cc_rtrend[idx]

  print(len(dt_cc))
  if len(dt_cc) == n: break

#====== plot dt_cc V.s. gcarc 
nrow = 4
ncol = 1
subplot_origin = np.array([ncol-1, nrow-1])*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.scatter(gcarc, dt_cc, marker='o', alpha=0.5, s=7, edgecolors='none', color='b')
ax.set_xlabel('Epidistance (degree)')
ax.set_ylabel('Tobs-Tsyn (s)')
xlim = [0, 16.5]
line = ax.plot(xlim, p(xlim), 'r-')
#ax.legend(line, "%5.2f s/deg"%(polycoef[0]))
ax.plot(xlim,[0,0],'k')
ax.set_xlim(xlim)
ax.set_ylim([-max_dt_plot, max_dt_plot])

title_str = "%s: %s" % (plot_title, window_id)
ax.set_title(title_str)

#====== plot CC0 against gcarc
nrow = 3
ncol = 1
subplot_origin = np.array([ncol-1, nrow-1])*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))

ax.hist([dt_cc, dt_cc_rtrend], 'auto', histtype='step', color=['b', 'r'],
    label=["%.3f$\pm$%.3f"%(np.mean(dt_cc),np.std(dt_cc)), 
           "0$\pm$%.3f"%(np.std(dt_cc_rtrend))])
ax.legend()
ax.set_xlim([-max_dt_plot, max_dt_plot])
ylim = ax.get_ylim()
ax.plot([0,0],ylim,'k')
ax.set_ylim(ylim)
ax.set_xlabel('Tobs-Tsyn (s)')
ax.set_ylabel('Window number')
ax.legend()

#====== plot CCmax against gcarc
nrow = 2
ncol = 1
subplot_origin = np.array([ncol-1, nrow-1])*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.scatter(gcarc, ccmax, marker='o', alpha=0.5, s=7, edgecolors='none', color='b')
ax.set_xlabel('Epidistance (degree)')
ax.set_ylabel('CCmax')
bin_edges = np.arange(min(gcarc), max(gcarc), 1.0)
bin_edges[-1] = max(gcarc)
nbin = len(bin_edges) - 1
cc_bin = np.zeros(nbin)
gcarc_bin = np.zeros(nbin)
cc_error = np.zeros(nbin)
for ibin in range(nbin):
  idx = (gcarc >= bin_edges[ibin]) & (gcarc < bin_edges[ibin+1])
  if np.any(idx):
    cc_bin[ibin] = np.mean(ccmax[idx])
    #gcarc_bin[ibin] = np.mean(gcarc[idx])
    gcarc_bin[ibin] = (bin_edges[ibin] + bin_edges[ibin+1])/2.0
    cc_error[ibin] = np.std(ccmax[idx])
  else:
    print("no points in this bin:", ibin)
    cc_bin[ibin] = np.nan
    gcarc_bin[ibin] = np.nan
    cc_error[ibin] = np.nan
ax.errorbar(gcarc_bin, cc_bin, yerr=cc_error, fmt='o', capsize=5)
#xlim = ax.get_xlim()
ax.set_ylim([0.5,1])
ax.set_xlim([0,16.5])

#====== plot histogram of ccmax
nrow = 1
ncol = 1

subplot_origin = np.array([ncol-1, nrow-1])*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.hist(ccmax, 'auto', histtype='step', color='b')
ax.set_xlim([min_cc_plot, 1])
ax.set_xlabel('CCmax')
ax.set_ylabel('Window number')

#------ save figure
plt.savefig(out_file, format='pdf')