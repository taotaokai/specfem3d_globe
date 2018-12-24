#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Merge Xin (2018) model within FWEA18
"""
import sys
import numpy as np
import scipy.sparse.linalg
from pyproj import Geod
import xarray

#====== parameters

misfit_file = str(sys.argv[1])
window_id = str(sys.argv[2]) #"surface_10-25_Z"
out_file = str(sys.argv[3])

min_SNR = 5
max_dt_cc = 10

R_EARTH = 6371.0

# regular 2-D grid
# should be in strictly ascending order
lon_grid = np.arange(94,110.1,0.5)
lat_grid = np.arange(21,35.1,0.5)

#
path_integration_interval = 0.05 # degree
damping_ratio = 0.01
laplace_ratio = 2.0

# read in misfit_file
with open(misfit_file, 'r') as f:
  lines = [ l.split() for l in f.readlines() if not(l.startswith('#'))]
lines = [l for l in lines if l[1] == window_id]

# filter SNR, dt_cc
dt_cc = np.array([ float(l[7]) for l in lines ])
SNR = np.array([ float(l[8]) for l in lines ])
idx = (SNR >= min_SNR) & (np.abs(dt_cc) < max_dt_cc)

nlines = len(lines)
lines = [lines[i] for i in range(nlines) if idx[i]]

ccmax = np.array([ float(l[6]) for l in lines ])
dt_cc = np.array([ float(l[7]) for l in lines ])
SNR = np.array([ float(l[8]) for l in lines ])
stlo = np.array([ float(l[12]) for l in lines ])
stla = np.array([ float(l[13]) for l in lines ])
evlo = np.array([ float(l[14]) for l in lines ])
evla = np.array([ float(l[15]) for l in lines ])

npath = dt_cc.size 
nlon = lon_grid.size
nlat = lat_grid.size
ngrid = nlon*nlat # (nlon,nlat)

# convert ECEF to Lon,Lat,Height
g = Geod(ellps='WGS84') 

# calculate weighted path length matrix 
weight_dt = np.zeros(npath)
weight_path_matrix = np.zeros((npath,ngrid))
for ipath in range(npath):
  az, baz, dist = g.inv(evlo[ipath], evla[ipath], stlo[ipath], stla[ipath])
  dist = dist/1000.0 # km
  npts = int(dist/R_EARTH/np.deg2rad(path_integration_interval)) - 1
  path_length = dist/npts
  lonlats = g.npts(evlo[ipath], evla[ipath], stlo[ipath], stla[ipath], npts)
  lonlats.insert(0, (evlo[ipath], evla[ipath]))
  lonlats.append((stlo[ipath], stla[ipath]))
  npts = len(lonlats)

  # path weight
  #weight_path = ccmax[ipath]
  #if weight_path < 0.0:
  #  weight_path = 0.0
  weight_path = 1.0

  # weighted dt_cc
  weight_dt[ipath] = weight_path*dt_cc[ipath]

  # bilinear interpolation coeff.
  for ipts in range(npts):
    # begin/end point only affects half interval
    if ipts == 0 or ipts == npts-1:
      factor = 0.5
    else:
      factor = 1.0
    lon, lat = lonlats[ipts]
    ilon_right = np.sum(lon_grid < lon)
    ilon_left = ilon_right -1
    ilat_upper = np.sum(lat_grid < lat)
    ilat_bottom = ilat_upper -1
    dlon = lon_grid[ilon_right] - lon_grid[ilon_left]
    dlat = lat_grid[ilat_upper] - lat_grid[ilat_bottom]
    plon = (lon_grid[ilon_right] - lon)/dlon
    plat = (lat_grid[ilat_upper] - lat)/dlat
    ind_bl = np.ravel_multi_index((ilon_left,ilat_bottom),(nlon,nlat))
    weight_path_matrix[ipath,ind_bl] += weight_path*factor*path_length*plat*plon
    ind_br = np.ravel_multi_index((ilon_right,ilat_bottom),(nlon,nlat))
    weight_path_matrix[ipath,ind_br] += weight_path*factor*path_length*plat*(1.0-plon)
    ind_ul = np.ravel_multi_index((ilon_left,ilat_upper),(nlon,nlat))
    weight_path_matrix[ipath,ind_ul] += weight_path*factor*path_length*(1.0-plat)*plon
    ind_ur = np.ravel_multi_index((ilon_right,ilat_upper),(nlon,nlat))
    weight_path_matrix[ipath,ind_ur] += weight_path*factor*path_length*(1.0-plat)*(1.0-plon)

#------ damped least square
#grid_density = np.sum(weight_path_matrix,axis=0)
#median_density = np.median(grid_density)
#damping_factor = median_density * damping_ratio

#A = np.zeros((npath+ngrid, ngrid))
#A[0:npath,:] = weight_path_matrix
#np.fill_diagonal(A[npath:,:], damping_factor)
#b = np.zeros(npath+ngrid)
#b[0:npath] = weight_dt
#ds, res = np.linalg.lstsq(A, b)

# obj = 1/2*(weight*L*ds - dt)^2 + 1/2*(Damp*ds)^2 + 1/2*gamma*(Laplace*ds)^2

A = np.dot(np.transpose(weight_path_matrix), weight_path_matrix)
b = np.dot(np.transpose(weight_path_matrix), weight_dt)

# add damping matrix
diag_A = np.diag(A).copy()
median_diag_A = np.median(diag_A)
damping_factor = median_diag_A * damping_ratio
idx = diag_A < damping_factor
diag_A[idx] = damping_factor
di = np.diag_indices(ngrid)
A[di] += damping_factor*(median_diag_A/diag_A) # larger damping for smaller diag_A

# add laplacian matrix
laplace = np.zeros((ngrid,ngrid))
for igrid in range(ngrid):
  ilon, ilat = np.unravel_index(igrid, (nlon,nlat))
  ilon4 = np.array([ilon-1, ilon+1, ilon, ilon])
  ilat4 = np.array([ilat, ilat, ilat+1, ilat-1])
  idx = (ilon4 >= 0) & (ilon4 < nlon) & (ilat4 >= 0) & (ilat4 < nlat)
  #ind = np.ravel_multi_index((ilon, ilat),(nlon,nlat))
  laplace[igrid,igrid] = 1.0
  ind = np.ravel_multi_index((ilon4[idx], ilat4[idx]),(nlon,nlat))
  laplace[igrid,ind] = -1/np.sum(idx)
A += median_diag_A*laplace_ratio * np.dot(np.transpose(laplace),laplace)

# preconditioner
preconditioner = np.diag(1.0/(diag_A+median_diag_A))

# inverted slowness perturbation
ds, info = scipy.sparse.linalg.cg(A, b, x0=np.zeros(ngrid), maxiter=50, M=preconditioner)
print("number of CG iterations: ", info)

predict_weight_dt = np.dot(weight_path_matrix, ds)
residual_weight_dt = weight_dt - predict_weight_dt

# approximate area from each grid (better way: bilinear weigthed area for each grid)
grid_area = np.zeros((nlon,nlat))
for ilon in range(nlon):
  for ilat in range(nlat):
    ilon3 = np.array([ilon-1, ilon, ilon+1])
    ilat3 = np.array([ilat-1, ilat, ilat+1])
    idx = (ilon3 >= 0) & (ilon3 < nlon)
    ilon0 = np.min(ilon3[idx]) 
    ilon1 = np.max(ilon3[idx]) 
    idx = (ilat3 >= 0) & (ilat3 < nlat)
    ilat0 = np.min(ilat3[idx]) 
    ilat1 = np.max(ilat3[idx]) 
    az, baz, ew_dist = g.inv(lon_grid[ilon0], lat_grid[ilat], lon_grid[ilon1], lat_grid[ilat])
    az, baz, ns_dist = g.inv(lon_grid[ilon], lat_grid[ilat0], lon_grid[ilon], lat_grid[ilat1])
    grid_area[ilon,ilat] = ew_dist*ns_dist/1000.0**2 # km^2

#====== make xarray
# note: in order for gmt.Figure.grdimage to plot 2-D grid, the grid dimensions must be arranged in (lat,lon).
slowness_map = xarray.DataArray(np.transpose(np.reshape(ds,(nlon,nlat))), coords={'lon':lon_grid, 'lat':lat_grid}, dims=('lat', 'lon'))
slowness_map.attrs = {'long_name':"slowness perturbation", 'unit':'s*km-1'}

path_density = np.reshape(np.sum(weight_path_matrix,axis=0),(nlon,nlat))/grid_area
path_density = xarray.DataArray(np.transpose(path_density), dims=('lat', 'lon'))
path_density.attrs = {'unit':'km-1'}

weight_dt = xarray.DataArray(weight_dt, coords={'path':np.arange(npath)}, dims=('path'))
weight_dt.attrs = {'unit':'s'}

residual_weight_dt = xarray.DataArray(residual_weight_dt, dims=('path'))
residual_weight_dt.attrs = {'unit':'s'}

outdata = xarray.Dataset({
  'ds':slowness_map, 
  'path_density':path_density,
  'dt':weight_dt,
  'dt_res':residual_weight_dt,
  })

outdata.attrs = {
  'min_SNR':min_SNR,
  'max_dt_cc':max_dt_cc,
  'window_id': window_id,
  'path_integration_interval_degree':path_integration_interval, # degree
  'damping_ratio':damping_ratio,
  'laplace_ratio':laplace_ratio,
  }

outdata.to_netcdf(out_file)
