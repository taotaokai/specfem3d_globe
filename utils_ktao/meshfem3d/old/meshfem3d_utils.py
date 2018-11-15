#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np

#///////////////////////////////////////////////
# constants.h
NGLLX = 5
NGLLY = NGLLX
NGLLZ = NGLLX

MIDX = int((NGLLX-1)/2)
MIDY = int((NGLLY-1)/2)
MIDZ = int((NGLLZ-1)/2)

GAUSSALPHA = 0
GAUSSBETA = 0

field_list = [                       
  ('nspec','i4')                      ,
  ('nglob','i4')                      ,
  ('x','f4')                          ,
  ('y','f4')                          ,
  ('z','f4')                          ,
  ('ibool','i4')                      ,
  ('idoubling','i4')                  ,
  ('ispec_is_tiso','i4')              ,
  ('DxiDx','f4')                      ,
  ('DxiDy','f4')                      ,
  ('DxiDz','f4')                      ,
  ('DetaDx','f4')                     ,
  ('DetaDy','f4')                     ,
  ('DetaDz','f4')                     ,
  ('DgammaDx','f4')                   ,
  ('DgammaDy','f4')                   ,
  ('DgammaDz','f4')                   ,
  ]

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
    for field in field_list:
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
    #NB: binary files are written in Fortran !!!

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

  wgll_cube = wx.reshape((NGLLX,1,1))*wy.reshape((1,NGLLY,1))*wx.reshape((1,1,NGLLZ))

  #--- jacobian * gll_quad_weights
  jacobian = 1.0/mesh_data[D]
  vol_gll = mesh_data['jacobian']*wgll_cube.reshape((NGLLX,NGLLY,NGLLZ,1))

  return vol_gll

#///////////////////////////////////////////////////
def sem_locate_points_hex8(mesh_data, xyz):
  """ locate points in the SEM mesh. 
  The anchor points are located on the 8 corners of each element.
  xyz(3,n)
  """
  from scipy import spatial
  from gll_library import zwgljd, lagrange_poly
  from jacobian_hex8 import xyz2cube_bounded_hex8, anchor_index_hex8

  nspec = mesh_data['nspec']
  ibool = mesh_data['ibool']
  xyz_glob = mesh_data['xyz_glob']
  xyz_elem = mesh_data['xyz_elem']

  #--- kdtree search nearby elements around each target point
  tree_elem = spatial.cKDTree(np.column_stack(
    (xyz_elem[0,:],xyz_elem[1,:],xyz_elem[2,:])))

  tree_xyz = spatial.cKDTree(np.column_stack(
    (xyz[0,:],xyz[1,:],xyz[2,:])))
  
  # determine maximum search radius
  max_element_size = 0.0 
  for ispec in range(nspec):
    iglob1 = ibool[0,0,0,ispec]-1
    iglob2 = ibool[-1,-1,-1,ispec]-1
    xyz1 = xyz_glob[:,iglob1]
    xyz2 = xyz_glob[:,iglob2]
    max_element_size = max(max_element_size, sum((xyz1-xyz2)**2)**0.5)

  neighbor_lists = tree_xyz.query_ball_tree(tree_elem, 1.5*max_element_size)

  #--- loop over each point, get the location info 
  iax, iay, iaz = anchor_index_hex8(NGLLX,NGLLY,NGLLZ)
  xigll, wx = zwgljd(NGLLX,GAUSSALPHA,GAUSSBETA)
  yigll, wy = zwgljd(NGLLY,GAUSSALPHA,GAUSSBETA)
  zigll, wz = zwgljd(NGLLZ,GAUSSALPHA,GAUSSBETA)

  npoints = xyz.shape[1]
  loc_data = [None] * npoints
  for ipoint in range(npoints):
    loc_data[ipoint] = {}
    loc_data[ipoint]['misloc'] = np.inf
    loc_data[ipoint]['is_inside'] = False
    if not neighbor_lists[ipoint]: continue
    # sort distance, start from the neareast element
    ispec_list = np.array(neighbor_lists[ipoint]) # covnert list to numpy array to have index slicing
    dist2 = np.sum((xyz_elem[:,ispec_list] - xyz[:,ipoint].reshape((3,1)))**2, axis=0)
    idx = np.argsort(dist2)
    for ispec in ispec_list[idx]:
      iglob = ibool[iax,iay,iaz,ispec] - 1
      xyz_anchor = xyz_glob[:,iglob]
      uvw, misloc, is_inside = xyz2cube_bounded_hex8(xyz_anchor, xyz[:,ipoint])
      if misloc < loc_data[ipoint]['misloc']:
        loc_data[ipoint]['uvw'] = uvw
        loc_data[ipoint]['misloc'] = misloc
        loc_data[ipoint]['ispec'] = ispec
        loc_data[ipoint]['is_inside'] = is_inside
      if is_inside: break
    if 'uvw' in loc_data[ipoint]:
      hlagx = lagrange_poly(xigll, uvw[0])
      hlagy = lagrange_poly(yigll, uvw[1])
      hlagz = lagrange_poly(zigll, uvw[2])
      loc_data[ipoint]['lagrange'] = hlagx[:,None,None]*hlagy[None,:,None]*hlagz[None,None,:]

  return loc_data
