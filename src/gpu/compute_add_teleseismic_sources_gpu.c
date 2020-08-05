/*
!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
*/

#include "mesh_constants_gpu.h"


extern EXTERN_LANG
void FC_FUNC_ (compute_add_teleseismic_sources_gpu,
               compute_add_teleseismic_sources_gpu) (long *Mesh_pointer_f,
                                                     realw *field_teleseismic,
                                                     int *itype) {

  TRACE ("compute_add_teleseismic_sources_gpu");

  int num_teleseismic_boundary_faces = 0;

  gpu_int_mem d_teleseismic_boundary_ispec;
  gpu_realw_mem d_field_teleseismic;

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // teleseismic boundary type
  int interface_type = *itype;
  switch (interface_type) {
  case 0:
    // xmin
    num_teleseismic_boundary_faces = mp->nspec2D_teleseismic_xmin;
    d_teleseismic_boundary_ispec = mp->d_ibelm_teleseismic_xmin;
    d_field_teleseismic = mp->d_field_teleseismic_xmin;
    break;
  case 1:
    // xmax
    num_teleseismic_boundary_faces = mp->nspec2D_teleseismic_xmax;
    d_teleseismic_boundary_ispec = mp->d_ibelm_teleseismic_xmax;
    d_field_teleseismic = mp->d_field_teleseismic_xmax;
    break;
  case 2:
    // ymin
    num_teleseismic_boundary_faces = mp->nspec2D_teleseismic_ymin;
    d_teleseismic_boundary_ispec = mp->d_ibelm_teleseismic_ymin;
    d_field_teleseismic = mp->d_field_teleseismic_ymin;
    break;
  case 3:
    // ymax
    num_teleseismic_boundary_faces = mp->nspec2D_teleseismic_ymax;
    d_teleseismic_boundary_ispec = mp->d_ibelm_teleseismic_ymax;
    d_field_teleseismic = mp->d_field_teleseismic_ymax;
    break;
  case 4:
    // zmin
    num_teleseismic_boundary_faces = mp->nspec2D_teleseismic_zmin;
    d_teleseismic_boundary_ispec = mp->d_ibelm_teleseismic_zmin;
    d_field_teleseismic = mp->d_field_teleseismic_zmin;
    break;
  default:
    exit_on_error ("compute_add_teleseismic_sources_gpu: unknown interface type");
  }

  // checks if anything to do
  if (num_teleseismic_boundary_faces == 0) return;

  // way 1
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //int blocksize = 32;

  // way 2: seems slightly faster
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (num_teleseismic_boundary_faces, &num_blocks_x, &num_blocks_y);

  // copies teleseismic source array to GPU
  gpuCopy_todevice_realw (&d_field_teleseismic, field_teleseismic, NDIM * NGLL2 * num_teleseismic_boundary_faces);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    // absorbing boundary contributions
    clCheck (clSetKernelArg (mocl.kernels.compute_add_teleseismic_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_accel_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_teleseismic_sources_kernel, idx++, sizeof (cl_mem), (void *) &d_field_teleseismic.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_teleseismic_sources_kernel, idx++, sizeof (int), (void *) &interface_type));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_teleseismic_sources_kernel, idx++, sizeof (int), (void *) &num_teleseismic_boundary_faces));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_teleseismic_sources_kernel, idx++, sizeof (cl_mem), (void *) &d_teleseismic_boundary_ispec.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_teleseismic_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_add_teleseismic_sources_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // absorbing boundary contributions
    compute_add_teleseismic_sources_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel_crust_mantle.cuda,
                                                                                  d_field_teleseismic.cuda,
                                                                                  interface_type,
                                                                                  num_teleseismic_boundary_faces,
                                                                                  d_teleseismic_boundary_ispec.cuda,
                                                                                  mp->d_ibool_crust_mantle.cuda);
  }
#endif

  GPU_ERROR_CHECKING ("compute_add_teleseismic_sources_gpu");
}


/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_add_teleseismic_sources_backward_gpu,
               compute_add_teleseismic_sources_backward_gpu) (long *Mesh_pointer_f,
                                                              realw *field_teleseismic,
                                                              int *itype) {

  TRACE ("compute_add_teleseismic_sources_backward_gpu");

  int num_teleseismic_boundary_faces = 0;

  gpu_int_mem d_teleseismic_boundary_ispec;
  gpu_realw_mem d_b_field_teleseismic;

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // teleseismic boundary type
  int interface_type = *itype;
  switch (interface_type) {
  case 0:
    // xmin
    num_teleseismic_boundary_faces = mp->nspec2D_teleseismic_xmin;
    d_teleseismic_boundary_ispec = mp->d_ibelm_teleseismic_xmin;
    d_b_field_teleseismic = mp->d_field_teleseismic_xmin;
    break;
  case 1:
    // xmax
    num_teleseismic_boundary_faces = mp->nspec2D_teleseismic_xmax;
    d_teleseismic_boundary_ispec = mp->d_ibelm_teleseismic_xmax;
    d_b_field_teleseismic = mp->d_field_teleseismic_xmax;
    break;
  case 2:
    // ymin
    num_teleseismic_boundary_faces = mp->nspec2D_teleseismic_ymin;
    d_teleseismic_boundary_ispec = mp->d_ibelm_teleseismic_ymin;
    d_b_field_teleseismic = mp->d_field_teleseismic_ymin;
    break;
  case 3:
    // ymax
    num_teleseismic_boundary_faces = mp->nspec2D_teleseismic_ymax;
    d_teleseismic_boundary_ispec = mp->d_ibelm_teleseismic_ymax;
    d_b_field_teleseismic = mp->d_field_teleseismic_ymax;
    break;
  case 4:
    // zmin
    num_teleseismic_boundary_faces = mp->nspec2D_teleseismic_zmin;
    d_teleseismic_boundary_ispec = mp->d_ibelm_teleseismic_zmin;
    d_b_field_teleseismic = mp->d_field_teleseismic_zmin;
    break;
  default:
    exit_on_error ("compute_add_teleseismic_sources_backward_gpu: unknown interface type");
  }

  // checks if anything to do
  if (num_teleseismic_boundary_faces == 0) return;

  // way 1
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //int blocksize = 32;

  // way 2: seems slightly faster
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (num_teleseismic_boundary_faces, &num_blocks_x, &num_blocks_y);

  // copies teleseismic source array to GPU
  gpuCopy_todevice_realw (&d_b_field_teleseismic, field_teleseismic, NDIM * NGLL2 * num_teleseismic_boundary_faces);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    // absorbing boundary contributions
    clCheck (clSetKernelArg (mocl.kernels.compute_add_teleseismic_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_accel_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_teleseismic_sources_kernel, idx++, sizeof (cl_mem), (void *) &d_b_field_teleseismic.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_teleseismic_sources_kernel, idx++, sizeof (int), (void *) &interface_type));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_teleseismic_sources_kernel, idx++, sizeof (int), (void *) &num_teleseismic_boundary_faces));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_teleseismic_sources_kernel, idx++, sizeof (cl_mem), (void *) &d_teleseismic_boundary_ispec.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_teleseismic_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_add_teleseismic_sources_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // absorbing boundary contributions
    compute_add_teleseismic_sources_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel_crust_mantle.cuda,
                                                                                  d_b_field_teleseismic.cuda,
                                                                                  interface_type,
                                                                                  num_teleseismic_boundary_faces,
                                                                                  d_teleseismic_boundary_ispec.cuda,
                                                                                  mp->d_ibool_crust_mantle.cuda);
  }
#endif

  GPU_ERROR_CHECKING ("compute_add_teleseismic_sources_backward_gpu");
}


