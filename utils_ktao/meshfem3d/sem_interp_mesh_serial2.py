#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" create horizontal slice of SEM model at a given depth
"""
import sys
import warnings
import time

import numpy as np
from scipy.io import FortranFile

#from mpi4py import MPI

from meshfem3d_constants import NGLLX,NGLLY,NGLLZ,GAUSSALPHA,GAUSSBETA
from gll_library import zwgljd, lagrange_poly
from meshfem3d_utils import sem_mesh_read, sem_locate_points_hex27

#====== parameters

# mesh files to interpolate from
nproc_source = int(sys.argv[1])
mesh_dir_source = str(sys.argv[2]) # <mesh_dir>/proc******_external_mesh.bin
model_dir_source = str(sys.argv[3]) # <model_dir>/proc******_<model_name>.bin

# mesh files for which to get interpolated values
nproc_target = int(sys.argv[4])
mesh_dir_target = str(sys.argv[5]) # <mesh_dir>/proc******_external_mesh.bin
#model_dir_target = str(sys.argv[6]) # <model_dir>/proc******_<model_name>.bin

model_names = str(sys.argv[6]) # comma delimited e.g. vp,vs,rho,qmu,qkappa
out_dir = str(sys.argv[7])

# model names
model_names = model_names.split(',')
nmodel = len(model_names)

#====== interpolate
#comm = MPI.COMM_WORLD
#mpi_size = comm.Get_size()
#mpi_rank = comm.Get_rank()

xigll, wx = zwgljd(NGLLX,GAUSSALPHA,GAUSSBETA)
yigll, wy = zwgljd(NGLLY,GAUSSALPHA,GAUSSBETA)
zigll, wz = zwgljd(NGLLZ,GAUSSALPHA,GAUSSBETA)

#--- loop over each slice of target SEM mesh
#for iproc_target in range(mpi_rank,nproc_target,mpi_size):
for iproc_target in range(nproc_target):

  print("====== iproc_target ", iproc_target)
  sys.stdout.flush()

  # read in target SEM mesh
  mesh_file = "%s/proc%06d_reg1_solver_data.bin"%(mesh_dir_target, iproc_target)
  mesh_data_target = sem_mesh_read(mesh_file)
  nspec_target = mesh_data_target['nspec']
  ibool_target = mesh_data_target['ibool']
  idoubling_target = mesh_data_target['idoubling']
  xyz_glob_target = mesh_data_target['xyz_glob']

  #gll_dims = mesh_data_target['gll_dims']
  #misloc_gll_target = np.zeros(gll_dims)
  #is_inside_gll_target = np.zeros(gll_dims, dtype='bool')
  #model_gll_target = np.zeros((nmodel,) + gll_dims)

  dims = (3,NGLLX*NGLLY*NGLLZ,nspec_target)
  xyz_gll_target = xyz_glob_target[:,ibool_target-1].reshape(dims)
  #idoubling_ext = np.zeros(ibool_target.shape) + idoubling_target
  #idoubling_ext = idoubling_ext.ravel()

  #npoints = xyz_target.shape[1]
  dims = (NGLLX*NGLLY*NGLLZ,nspec_target)
  misloc_gll_target = np.zeros(dims)
  misloc_gll_target[:] = np.inf
  is_inside_gll_target = np.zeros(dims,dtype='bool')
  model_gll_target = np.zeros((nmodel,)+dims)

  # loop over each slice of source SEM mesh
  #for iproc_source in range(nproc_source):
  for iproc_source in range(63,nproc_source):

    print("iproc_source ", iproc_source)
    sys.stdout.flush()

    # read in source SEM mesh
    mesh_file = "%s/proc%06d_reg1_solver_data.bin"%(mesh_dir_source, iproc_source)
    source_mesh_data = sem_mesh_read(mesh_file)

    # read in source model
    gll_dims = source_mesh_data['gll_dims']
    source_model_gll = np.zeros((nmodel,)+gll_dims)
    for imodel in range(nmodel):
      model_tag = model_names[imodel]
      model_file = "%s/proc%06d_reg1_%s.bin"%(model_dir_source, iproc_source, model_tag)
      with FortranFile(model_file, 'r') as f:
        # note: must use fortran convention when reshape to N-D array!!!
        source_model_gll[imodel,:,:,:,:] = np.reshape(f.read_ints(dtype='f4'), 
            gll_dims, order='F')

    # locate target points
    for ispec_target in range(nspec_target):

      xyz_target = xyz_gll_target[:,:,ispec_target]

      loc_data = sem_locate_points_hex27(source_mesh_data, 
          xyz_target, idoubling_target[ispec_target])
  
      # record interpolation results based on misloc and is_inside
      for ipoint in range(len(loc_data)):

        loc = loc_data[ipoint]

        if loc['is_inside'] and is_inside_gll_target[ipoint,ispec_target]:
          warnings.warn("point is located inside more than one element", 
              xyz_target[:,ipoint])

        if ( loc['misloc'] > misloc_gll_target[ipoint,ispec_target] 
            and is_inside_gll_target[ipoint,ispec_target] ):
          warnings.warn("point located inside an element but with a larger misloc")

        if ( loc['misloc'] < misloc_gll_target[ipoint,ispec_target] 
            or is_inside_gll_target[ipoint,ispec_target] ):
          misloc_gll_target[ipoint,ispec_target] = loc['misloc']
          is_inside_gll_target[ipoint,ispec_target] = loc['is_inside']
          uvw = loc['uvw']
          hlagx = lagrange_poly(xigll, uvw[0])
          hlagy = lagrange_poly(yigll, uvw[1])
          hlagz = lagrange_poly(zigll, uvw[2])
          model_gll_target[:,ipoint,ispec_target] = np.sum(source_model_gll[:,:,:,:,loc['ispec']]*hlagx[:,None,None]*hlagy[None,:,None]*hlagz[None,None,:], axis=(1,2,3))
        
  #end for loop over each slice of source SEM mesh

  # rehape results
  gll_dims = mesh_data_target['gll_dims']
  misloc_gll_target = np.reshape(misloc_gll_target, gll_dims)
  is_inside_gll_target = np.reshape(is_inside_gll_target, gll_dims)
  gll_dims = (nmodel,) + gll_dims
  model_gll_target = np.reshape(model_gll_target,gll_dims)

  # save interpolated model
  for imodel in range(nmodel):
    model_tag = model_names[imodel]
    model_file = "%s/proc%06d_reg1_%s.bin"%(out_dir, iproc_target, model_tag)
    with FortranFile(model_file, 'w') as f:
      f.write_record(np.array(model_gll_target[imodel,:,:,:,:], dtype='f4'))
