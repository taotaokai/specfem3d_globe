#!/usr/bin/env python

# ktao: conver grd file to specfem binary

import numpy as np
from netCDF4 import Dataset
from scipy.io import FortranFile

grdfile = 'ETOPO1_Ice_g_smooth_b15km_I4m.grd'
output_file = 'ETOPO1_Ice_g_smooth_b15km_I4m.bin'

fh = Dataset(grdfile,'r')

topo = np.array(fh.variables['z'][:], dtype='i2')

fid = open(output_file, 'wb')

# Add a byte-order mark
byteorder = np.array([0x1234], dtype=np.int16)
byteorder.tofile(fid)

topo = np.ravel(topo, order='C') # src/shared/model_topo_bathy.f90:read_topo_bathy_file()
topo.tofile(fid)
fid.close()