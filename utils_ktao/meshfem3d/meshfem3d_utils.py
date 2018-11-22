#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np

from meshfem3d_constants import *

#///////////////////////////////////////////////
# constants.h
#NGLLX = 5
#NGLLY = NGLLX
#NGLLZ = NGLLX
#
#MIDX = int((NGLLX-1)/2)
#MIDY = int((NGLLY-1)/2)
#MIDZ = int((NGLLZ-1)/2)
#
#GAUSSALPHA = 0
#GAUSSBETA = 0
#
## proc000*_reg1_solver_data.bin
#MESH_ARRAY_LIST = [                       
#  ('nspec','i4')                      ,
#  ('nglob','i4')                      ,
#  ('x','f4')                          ,
#  ('y','f4')                          ,
#  ('z','f4')                          ,
#  ('ibool','i4')                      ,
#  ('idoubling','i4')                  ,
#  ('ispec_is_tiso','i4')              ,
#  ('DxiDx','f4')                      ,
#  ('DxiDy','f4')                      ,
#  ('DxiDz','f4')                      ,
#  ('DetaDx','f4')                     ,
#  ('DetaDy','f4')                     ,
#  ('DetaDz','f4')                     ,
#  ('DgammaDx','f4')                   ,
#  ('DgammaDy','f4')                   ,
#  ('DgammaDz','f4')                   ,
#  ]

#//////////////////////////////////////////////
def rotmat_enu_to_ecef(lon,lat):
  """ rotation matrix from local ENU to ECEF coordinate basises
  rotmat[:,0] = Ve # column vector of the Easting direction in ECEF coordinate
  rotmat[:,1] = Vn # column vector of the Northing direction in ECEF coordinate
  rotmat[:,2] = Vu # column vector of the Up (ellipsoid height) direction in ECEF coordinate
  
  xyz_ecef = xyz0_ecef + rotmat * enu
  enu = transpose(rotmat) * (xyz_ecef - xyz0_ecef)

  , where xyz0_ecef is the reference point at (lon,lat,alt).
  """
  coslat = np.cos(np.deg2rad(lat))
  sinlat = np.sin(np.deg2rad(lat))
  coslon = np.cos(np.deg2rad(lon))
  sinlon = np.sin(np.deg2rad(lon))
  
  rotmat = np.zeros((3,3))
  rotmat[0,:] = [ -sinlon, -sinlat*coslon, coslat*coslon ]
  rotmat[1,:] = [  coslon, -sinlat*sinlon, coslat*sinlon ]
  rotmat[2,:] = [     0.0,  coslat,        sinlat        ]

  return rotmat


#///////////////////////////////////////////////////
def sem_mesh_read(mesh_file):
  """ read in SEM mesh slice
  """
  from scipy.io import FortranFile

  mesh_data = {}

  with FortranFile(mesh_file, 'r') as f:
    for field in MESH_ARRAY_LIST:
      field_name = field[0]
      data_type = field[1]
      mesh_data[field_name] = f.read_ints(dtype=data_type)
  
  mesh_data['nspec'] = mesh_data['nspec'][0]
  mesh_data['nglob'] = mesh_data['nglob'][0]

  # GLL dims
  gll_dims = (NGLLX,NGLLY,NGLLZ,mesh_data['nspec'])
  mesh_data['gll_dims'] = gll_dims

  # reshape
  for field_name in ['ibool', 'DxiDx','DxiDy','DxiDz','DetaDx','DetaDy','DetaDz','DgammaDx','DgammaDy','DgammaDz',]:
    mesh_data[field_name] = np.reshape(mesh_data[field_name], gll_dims, order='F')
    #NB: binary files are written in Fortran column-major convention !!!

  # jacobian: det( d(x,y,z)/d(xi,eta,gamma))
  jacobian = 1.0 / ( 
      mesh_data['DxiDx']*(mesh_data['DetaDy']*mesh_data['DgammaDz']-mesh_data['DetaDz']*mesh_data['DgammaDy'])
      -mesh_data['DxiDy']*(mesh_data['DetaDx']*mesh_data['DgammaDz']-mesh_data['DetaDz']*mesh_data['DgammaDx'])
      +mesh_data['DxiDz']*(mesh_data['DetaDx']*mesh_data['DgammaDy']-mesh_data['DetaDy']*mesh_data['DgammaDx'])
      )

  del mesh_data['DxiDx']
  del mesh_data['DxiDy'] 
  del mesh_data['DxiDz']
  del mesh_data['DetaDx']
  del mesh_data['DetaDy']
  del mesh_data['DetaDz']
  del mesh_data['DgammaDx']
  del mesh_data['DgammaDy']
  del mesh_data['DgammaDz']

  # use xyz_glob
  nglob = mesh_data['nglob']
  x = mesh_data['x'].reshape((1,nglob))
  y = mesh_data['y'].reshape((1,nglob))
  z = mesh_data['z'].reshape((1,nglob))
  mesh_data['xyz_glob'] = np.r_[x,y,z]

  del mesh_data['x']
  del mesh_data['y']
  del mesh_data['z']

  # add xyz_elem
  iglob_elem = mesh_data['ibool'][MIDX,MIDY,MIDZ,:] - 1
  mesh_data['xyz_elem'] = mesh_data['xyz_glob'][:,iglob_elem]

  #FIXME bad idea due to 410 undulation. Need to modify the specfem code
  ## separate mesh layers across 410-km
  ## 40: above 410, 41: below 410
  #idoubling = mesh_data['idoubling']
  #depth = (1.0 - np.sum(mesh_data['xyz_elem']**2, axis=0)*0.5) * R_EARTH_KM

  ## this is dangerous due to 410 undulation
  #ii = (idoubling == IFLAG_670_220) & (depth < 410)
  #idoubling[ii] = 10*IFLAG_670_220

  #ii = (idoubling == IFLAG_670_220) & (depth > 410)
  #idoubling[ii] = 10*IFLAG_670_220 + 1

# nspec = int(mesh_data['nspec'])
#  for ispec in range(nspec):
#    for i in range(NGLLX):
#      for j in range(NGLLY):
#        for k in range(NGLLZ):
#          iglob = mesh_data['ibool'][i,j,k,ispec] - 1
#          xyz_gll[0,i,j,k,ispec] = mesh_data['x'][iglob]
#          xyz_gll[1,i,j,k,ispec] = mesh_data['y'][iglob]
#          xyz_gll[2,i,j,k,ispec] = mesh_data['z'][iglob]
#  xyz_gll = np.zeros((3,NGLLX,NGLLY,NGLLZ,nspec))

  return mesh_data


#///////////////////////////////////////////////////
def sem_mesh_get_vol_gll(mesh_data):
  """ get xyz and volumen facotr of each gll point
  """

  from gll_library import zwgljd

  #--- quadrature weights on GLL points
  zx, wx = zwgljd(NGLLX,GAUSSALPHA,GAUSSBETA)
  zy, wy = zwgljd(NGLLY,GAUSSALPHA,GAUSSBETA)
  zz, wz = zwgljd(NGLLZ,GAUSSALPHA,GAUSSBETA)

  #wgll_cube = wx.reshape((NGLLX,1,1))*wy.reshape((1,NGLLY,1))*wx.reshape((1,1,NGLLZ))

  #--- jacobian * gll_quad_weights
  #vol_gll = mesh_data['jacobian']*wgll_cube.reshape((NGLLX,NGLLY,NGLLZ,1))
  vol_gll = mesh_data['jacobian'] * wx[NGLLX,None,None,None] * wy[None,NGLLY,None,None] * wz[None,None,NGLLZ,None]

  return vol_gll

#///////////////////////////////////////////////////
def sem_locate_points_hex27(mesh_data, xyz, idoubling=-1):
  """ locate points in the SEM mesh. 
  mesh_data: return value from sem_mesh_read()
  xyz(3,n): locations of n points
  idoubling: integer or integer array of size n. idoubling for the n points which denotes mesh regions (surface-Moho,Moho-410,410-660,etc)
    -1 means no specific region and interpolation will be done to all the elements in the mesh, otherwise interpolation is only done for those elements with the same idoubling value. 

  output:
    status(n): -1=not located,0=close to but outside the element,1=inside element
    ispec(n): element num that located
    uvw(3,n): local coordinate located
    misloc(n): location residual
    misratio(n): misloc/element_size
  """
  from scipy import spatial
  #from gll_library import zwgljd, lagrange_poly
  from jacobian_hex27 import xyz2cube_bounded_hex27, anchor_index_hex27

  npoints = xyz.shape[1]
  idoubling = np.array(idoubling, dtype='int')
  if idoubling.size == 1:
    idoubling = np.ones(npoints,dtype='int')*int(idoubling)
  elif idoubling.size != npoints:
    raise Exception("idoubling must either be an integer or an integer array of npoints ")

  nspec = mesh_data['nspec']
  ibool = mesh_data['ibool']
  source_idoubling = mesh_data['idoubling']
  xyz_glob = mesh_data['xyz_glob']
  xyz_elem = mesh_data['xyz_elem']

  #--- kdtree search nearby elements around each target point
  tree_elem = spatial.cKDTree(np.column_stack(
    (xyz_elem[0,:],xyz_elem[1,:],xyz_elem[2,:])))

  tree_xyz = spatial.cKDTree(np.column_stack(
    (xyz[0,:],xyz[1,:],xyz[2,:])))
  
  # determine element size (approximately)
  element_size = np.zeros(nspec)
  for ispec in range(nspec):
    # distance between gll points and the central gll point 
    iglob1 = ibool[:,:,:,ispec].ravel() - 1
    dist = np.sum((xyz_elem[:,ispec:ispec+1] - xyz_glob[:,iglob1])**2, axis=0)**0.5
    element_size[ispec] = np.max(dist)

  # get neighbouring elements around each target location xyz
  neighbor_lists = tree_xyz.query_ball_tree(tree_elem, 1.2*max(element_size))

  #--- loop over each point, get the location info 
  iax, iay, iaz = anchor_index_hex27(NGLLX,NGLLY,NGLLZ)

  #xigll, wx = zwgljd(NGLLX,GAUSSALPHA,GAUSSBETA)
  #yigll, wy = zwgljd(NGLLY,GAUSSALPHA,GAUSSBETA)
  #zigll, wz = zwgljd(NGLLZ,GAUSSALPHA,GAUSSBETA)

  status_all = np.zeros(npoints, dtype='int')
  status_all[:] = -1
  ispec_all = np.zeros(npoints,dtype='int')
  uvw_all = np.zeros((3,npoints))
  misloc_all = np.zeros(npoints)
  misloc_all[:] = np.inf
  misratio_all = np.zeros(npoints) 

  ipoint_select = [ ipoint for ipoint in range(npoints) if neighbor_lists[ipoint] ]
  #for ipoint in range(npoints):
  for ipoint in ipoint_select:
    #if not neighbor_lists[ipoint]: continue
    # get neibouring elements
    ispec_list = np.array(neighbor_lists[ipoint]) # convert list to numpy array to have index slicing
    # get ratio between distance to center and element size 
    dist_ratio = np.sum((xyz_elem[:,ispec_list] - xyz[:,ipoint:ipoint+1])**2, axis=0)**0.5 / element_size[ispec_list]
    # remove elements too far away from target point
    idx = dist_ratio < 1.1
    # skip elements that does NOT have the same idoubling as xyz
    if idoubling[ipoint] != -1:
      idx = idx & (source_idoubling[ispec_list] == idoubling[ipoint])
    ispec_list = ispec_list[idx]
    dist_ratio = dist_ratio[idx]
    # loop each element, start from the closest element
    for ispec in ispec_list[np.argsort(dist_ratio)]:
      #if (idoubling[ipoint] != -1 and
      #    idoubling[ipoint] != source_idoubling[ispec]):
      #  continue
      iglob = ibool[iax,iay,iaz,ispec] - 1
      xyz_anchor = xyz_glob[:,iglob]
      uvw, misloc, is_inside = xyz2cube_bounded_hex27(xyz_anchor, xyz[:,ipoint])
      ##DEBUG
      #if is_inside and status_all[ipoint]==1:
      #  warnings.warn("point is located inside more than one element", 
      #      xyz[:,ipoint], xyz_anchor)
      if misloc > misloc_all[ipoint] and is_inside:
        warnings.warn("point located inside an element but with a larger misloc(loc/previous)", misloc, misloc_all[ipoint])
      if misloc < misloc_all[ipoint] or is_inside:
        status_all[ipoint] = 0
        ispec_all[ipoint] = ispec
        uvw_all[:,ipoint] = uvw
        misloc_all[ipoint] = misloc
        misratio_all[ipoint] = misloc/element_size[ispec]
      # skip the rest elements since points already located inside an element
      # this means if multiple elements overlap (should not occur) we only take the first found element where the point locates inside
      if is_inside: 
        status_all[ipoint] = 1
        break
    #if 'uvw' in loc_data[ipoint]:
    #  hlagx = lagrange_poly(xigll, uvw[0])
    #  hlagy = lagrange_poly(yigll, uvw[1])
    #  hlagz = lagrange_poly(zigll, uvw[2])
    #  loc_data[ipoint]['lagrange'] = hlagx[:,None,None]*hlagy[None,:,None]*hlagz[None,None,:]

  return status_all, ispec_all, uvw_all, misloc_all, misratio_all
