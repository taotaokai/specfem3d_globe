#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Merge Xin (2018) model within FWEA18
"""
import sys
import numpy as np
from pyproj import Geod
import xarray as xr

import matplotlib as mpl
#mpl.rcParams['font.family'] = "times new roman"
mpl.use("pdf")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18
plt.rc('font', size=SMALL_SIZE, family='times new roman')          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#====== parameters
nc_file = str(sys.argv[1])
out_fig = str(sys.argv[2])

#--- read nc file 
data = xr.open_dataset(nc_file)

lons = data['lon'].values
lats = data['lat'].values

ds = data['ds'].values
path_density = data['path_density'].values
dt = data['dt'].values
dt_res = data['dt_res'].values

max_dt_cc = data.attrs['max_dt_cc']
window_id = data.attrs['window_id']

map_parallels = np.arange(0.,81,5.)
map_meridians = np.arange(0.,351,5.)

#--- plot geological blocks 
block_line_file = 'zhangpz_block.txt'
block_lines = []
with open(block_line_file, 'r') as f:
  lon = []
  lat = []
  for l in f.readlines():
    if not l.startswith('>'):
      x = l.split()
      lon.append(float(x[0]))
      lat.append(float(x[1]))
    else:
      block_lines.append([lon, lat])
      lon = []
      lat = []
#for l in block_lines:
#  #x, y = m(l[0], l[1])
#  #ax.plot(x, y, lw=0.2, color='gray')
#  fig.plot(x=l[0], y=l[1], W="thick,gray")

#--- plot STATIONS
STATIONS_file = 'STATIONS'
with open(STATIONS_file, 'r') as f:
  lines = [ l.split() for l in f.readlines() if not l.startswith('#') ] 
  stla = np.array([float(l[2]) for l in lines])
  stlo = np.array([float(l[3]) for l in lines])

#--- plot virtual sources 
STATIONS_file = 'STATIONS_select'
with open(STATIONS_file, 'r') as f:
  lines = [ l.split() for l in f.readlines() if not l.startswith('#') ] 
  evla = np.array([float(l[2]) for l in lines])
  evlo = np.array([float(l[3]) for l in lines])

#====== plot
#fig = plt.figure(figsize=(8.5, 11)) # US Letter
fig = plt.figure(figsize=(8.27, 11.7)) # A4

# a matrix of sub plot 
nrow = 3
ncol = 1
subplot_size = np.array([1.0/ncol, 1.0/nrow])
# axis position relative to the subplot region
ax_origin_subplot = np.array([0.2, 0.2])
# size: [width, height]
ax_size_subplot = np.array([0.7, 0.7])

#--- plot histogram 
nrow = 1
ncol = 1

subplot_origin = np.array([ncol-1, nrow-1])*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))

ax.hist((dt, dt_res), 30, histtype='step', color=('b','r'), label=('dt','dt_res'))
ax.legend()
ax.set_xlim([-1*max_dt_cc, max_dt_cc])
ax.set_xlabel("dt (s)")
ax.set_ylabel("number of measurements")
ax.set_title(window_id)

#--- plot spatial density of weighted path length
nrow = 2
ncol = 1

subplot_origin = np.array([ncol-1, nrow-1])*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.set_title("path density")

m = Basemap(ax=ax, projection='merc',
    resolution='l', area_thresh=30000,
    llcrnrlat=min(lats), llcrnrlon=min(lons), 
    urcrnrlat=max(lats), urcrnrlon=max(lons), lat_ts=20) 
m.drawcoastlines(linewidth=0.2)
m.drawcountries(linewidth=0.2)
m.drawparallels(map_parallels, linewidth=0.1, labels=[1,0,0,0], fontsize=12)
m.drawmeridians(map_meridians, linewidth=0.1, labels=[0,0,0,1], fontsize=12)

cmap = plt.get_cmap('rainbow')
x, y = np.meshgrid(lons,lats)
cs = m.contourf(x, y, path_density, latlon=True, cmap=cmap, extend='both')

# colorbar for contourfill
cb = m.colorbar(cs, location='right', pad="10%")
cb.ax.set_title('(km-1)', fontsize=10)
#cb.set_label('path density')

for l in block_lines:
  m.plot(l[0], l[1], latlon=True, lw=1, color='gray')
m.scatter(stlo,stla,latlon=True,s=10,marker='^', color='k')
m.scatter(evlo,evla,latlon=True,s=80,marker='*', color='r')

#--- plot slowness map
nrow = 3
ncol = 1

subplot_origin = np.array([ncol-1, nrow-1])*subplot_size
ax_origin = ax_origin_subplot*subplot_size + subplot_origin
ax_size = ax_size_subplot*subplot_size
ax = fig.add_axes(np.concatenate((ax_origin, ax_size)))
ax.set_title("slowness")

m = Basemap(ax=ax, projection='merc',
    resolution='l', area_thresh=30000,
    llcrnrlat=min(lats), llcrnrlon=min(lons), 
    urcrnrlat=max(lats), urcrnrlon=max(lons), lat_ts=20) 
m.drawcoastlines(linewidth=0.2)
m.drawcountries(linewidth=0.2)
m.drawparallels(map_parallels, linewidth=0.1, labels=[1,0,0,0], fontsize=12)
m.drawmeridians(map_meridians, linewidth=0.1, labels=[0,0,0,1], fontsize=12)

cmap = plt.get_cmap('jet')
levels = np.linspace(-0.02,0.02,101)
x, y = np.meshgrid(lons,lats)
cs = m.contourf(x, y, ds, latlon=True, cmap=cmap, levels=levels, extend='both')
cs.cmap.set_over('purple')
cs.cmap.set_under('black')

# colorbar for contourfill
ticks = np.arange(-0.02,0.02001, 0.005)
cb = m.colorbar(cs, location='right', ticks=ticks, pad="10%")
cb.ax.set_title('(s*km-1)', fontsize=10)
#cb.set_label('slowness')

for l in block_lines:
  m.plot(l[0], l[1], latlon=True, lw=1, color='gray')
m.scatter(stlo,stla,latlon=True,s=10,marker='^', color='k')
m.scatter(evlo,evla,latlon=True,s=80,marker='*', color='r')

#--- save figure
plt.savefig(out_fig, format='pdf')
