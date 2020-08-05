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

  subroutine compute_add_teleseismic_sources()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: i,j,k,iglob,ispec,ispec2D
  integer :: it_tmp

  ! checks
  if (SIMULATION_TYPE /= 1) return

  ! sets iteration step
  if (USE_LDDRK) then
    !it_tmp = it -1 + C_LDDRK(istage) ! C_LDDRK is fractional step lengths from 0 to 1
    !TODO interpolate source between it-1 and it at it_tmp for LDDRK
    call exit_MPI(myrank,'LDDRK is not yet supported in compute_add_teleseismic_source')
  else
    it_tmp = it
  endif

  !   xmin
  ! if two chunks exclude this face for one of them
  if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then
    if (nspec2D_teleseismic_xmin > 0) then
      call read_teleseismic(0,field_teleseismic_xmin,reclen_teleseismic_xmin,it_tmp)
    endif
    if (.not. GPU_MODE) then
      ! on CPU
      do ispec2D = 1,nspec2D_teleseismic_xmin
        ispec = ibelm_teleseismic_xmin(ispec2D)
        i = 1
        do k = 1,NGLLZ
          do j = 1,NGLLY
            iglob = ibool_crust_mantle(i,j,k,ispec)
            accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &
              + field_teleseismic_xmin(:,j,k,ispec2D)
          enddo
        enddo
      enddo
    else
      ! on GPU
      if (nspec2D_teleseismic_xmin > 0 ) then
        call compute_add_teleseismic_sources_gpu(Mesh_pointer, field_teleseismic_xmin, 0)
      endif
    endif
  endif ! NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC

  !   xmax
  ! if two chunks exclude this face for one of them
  if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then
    if (nspec2D_teleseismic_xmax > 0) then
      call read_teleseismic(0,field_teleseismic_xmax,reclen_teleseismic_xmax,it_tmp)
    endif
    if (.not. GPU_MODE) then ! on CPU
      do ispec2D = 1,nspec2D_teleseismic_xmax
        ispec = ibelm_teleseismic_xmax(ispec2D)
        i = NGLLX
        do k = 1,NGLLZ
          do j = 1,NGLLY
            iglob = ibool_crust_mantle(i,j,k,ispec)
            accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &
              + field_teleseismic_xmax(:,j,k,ispec2D)
          enddo
        enddo
      enddo
    else ! on GPU
      if (nspec2D_teleseismic_xmax > 0 ) then
        call compute_add_teleseismic_sources_gpu(Mesh_pointer, field_teleseismic_xmax, 1)
      endif
    endif
  endif ! NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB

  ! ymin
  if (nspec2D_teleseismic_ymin > 0) then
    call read_teleseismic(0,field_teleseismic_ymin,reclen_teleseismic_ymin,it_tmp)
  endif
  if (.not. GPU_MODE) then
    ! on CPU
    do ispec2D = 1,nspec2D_teleseismic_ymin
      ispec = ibelm_teleseismic_ymin(ispec2D)
      j = 1
      do k = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &
            + field_teleseismic_ymin(:,j,k,ispec2D)
        enddo
      enddo
    enddo
  else
    ! on GPU
    if (nspec2D_teleseismic_ymin > 0 ) then
      call compute_add_teleseismic_sources_gpu(Mesh_pointer, field_teleseismic_ymin, 2)
    endif
  endif

  ! ymax
  if (nspec2D_teleseismic_ymax > 0) then
    call read_teleseismic(0,field_teleseismic_ymax,reclen_teleseismic_ymin,it_tmp)
  endif
  if (.not. GPU_MODE) then
    ! on CPU
    do ispec2D = 1,nspec2D_teleseismic_ymax
      ispec = ibelm_teleseismic_ymax(ispec2D)
      j = NGLLY
      do k = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &
            + field_teleseismic_ymax(:,j,k,ispec2D)
        enddo
      enddo
    enddo
  else
    ! on GPU
    if (nspec2D_teleseismic_ymax > 0 ) then
      call compute_add_teleseismic_sources_gpu(Mesh_pointer, field_teleseismic_ymax, 3)
    endif
  endif

  ! zmin
  if (nspec2D_teleseismic_zmin > 0) then
    call read_teleseismic(0,field_teleseismic_zmin,reclen_teleseismic_ymin,it_tmp)
  endif
  if (.not. GPU_MODE) then
    ! on CPU
    do ispec2D = 1,nspec2D_teleseismic_zmin
      ispec = ibelm_teleseismic_zmin(ispec2D)
      k = 1
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &
            + field_teleseismic_zmin(:,j,k,ispec2D)
        enddo
      enddo
    enddo
  else
    ! on GPU
    if (nspec2D_teleseismic_zmin > 0 ) then
      call compute_add_teleseismic_sources_gpu(Mesh_pointer, field_teleseismic_zmin, 4)
    endif
  endif

!  if (.not. GPU_MODE) then
!    ! on CPU
!    !xmin
!    do ispec2D = 1, nspec2D_teleseismic_xmin
!      ispec = ibelm_teleseismic_xmin(ispec2D)
!      i = 1
!      do k = 1,NGLLZ
!        do j = 1,NGLLY
!          iglob = ibool_crust_mantle(i,j,k,ispec)
!          accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &
!            + teleseismic_source_xmin(:,j,k,ispec2D,it_tmp)
!        enddo
!      enddo
!    enddo
!    !xmax
!    do ispec2D = 1, nspec2D_teleseismic_xmax
!      ispec = ibelm_teleseismic_xmax(ispec2D)
!      i = NGLLX
!      do k = 1,NGLLZ
!        do j = 1,NGLLY
!          iglob = ibool_crust_mantle(i,j,k,ispec)
!          accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &
!            + teleseismic_source_xmax(:,j,k,ispec2D,it_tmp)
!        enddo
!      enddo
!    enddo
!    !ymin
!    do ispec2D = 1, nspec2D_teleseismic_ymin
!      ispec = ibelm_teleseismic_ymin(ispec2D)
!      j = 1
!      do k = 1,NGLLZ
!        do i = 1,NGLLX
!          iglob = ibool_crust_mantle(i,j,k,ispec)
!          accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &
!            + teleseismic_source_ymin(:,i,k,ispec2D,it_tmp)
!        enddo
!      enddo
!    enddo
!    !ymax
!    do ispec2D = 1, nspec2D_teleseismic_ymax
!      ispec = ibelm_teleseismic_ymax(ispec2D)
!      j = 1
!      do k = 1,NGLLZ
!        do i = 1,NGLLX
!          iglob = ibool_crust_mantle(i,j,k,ispec)
!          accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &
!            + teleseismic_source_ymax(:,i,k,ispec2D,it_tmp)
!        enddo
!      enddo
!    enddo
!    !zmin
!    do ispec2D = 1, nspec2D_teleseismic_zmin
!      ispec = ibelm_teleseismic_zmin(ispec2D)
!      k = 1
!      do j = 1,NGLLY
!        do i = 1,NGLLX
!          iglob = ibool_crust_mantle(i,j,k,ispec)
!          accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &
!            + teleseismic_source_zmin(:,i,j,ispec2D,it_tmp)
!        enddo
!      enddo
!    enddo
!
!  else
!    ! on GPU
!    !call compute_add_sources_gpu(Mesh_pointer, NSOURCES,stf_pre_compute)
!    call exit_MPI(myrank,'GPU is not yet supported in compute_add_teleseismic_source')
!  endif

  end subroutine compute_add_teleseismic_sources

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_add_teleseismic_sources_backward()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: i,j,k,iglob,ispec,ispec2D
  integer :: it_tmp

  ! checks
  if (SIMULATION_TYPE /= 3) return

  ! get the iteration step corresponding to forward simulation
  if (UNDO_ATTENUATION) then
    ! forward reconstruction subsets of time steps starts from the last snapshot 
    ! (iteration_on_subset = 1) to the first/begining snapshot
    ! (iteration_on_subset = NSUBSET_ITERATIONS)
    it_tmp = (NSUBSET_ITERATIONS - iteration_on_subset)*NT_DUMP_ATTENUATION + it_of_this_subset
  else
    ! backward reconstruction
    it_tmp = NSTEP - it + 1
  endif

  ! sets time/index for force term
  if (USE_LDDRK) then
  !  if (UNDO_ATTENUATION) then
  !    !stepping forward from snapshot position
  !    !current_time = (it_tmp-1)*DT-t0. 
  !    !LDDRK uses force terms before current_time, from it-1 to it.
  !    time_t = -t0 + dble(it_tmp-1)*DT - DT + dble(C_LDDRK(istage))*DT
  !    it_tmp = it_tmp - 1 + C_LDDRK(istage)
  !  else ! iterate_time.F90
  !    ! stepping backward
  !    ! use force terms in the future of current time, from it+1 to it
  !    !FIXME need confirmation of this formula
  !    time_t = -t0 + dble(it_tmp-1)*DT + DT - dble(C_LDDRK(istage))*DT
  !    it_tmp = it_tmp + 1 - C_LDDRK(istage)
  !  endif
  ! !NOTE it_tmp is not an integer and needs interpolation.
    call exit_MPI(myrank, &
    'LDDRK is not yet supported in compute_add_teleseismic_sources_backward')
  !else
  !  ! Newmark only explicitly uses force term at it
  !  ! forces at it-1 or it are implicitly contained in the acceleration term.
  !  time_t = dble(it_tmp-1)*DT - t0
  !  it_tmp = it_tmp
  endif

  !   xmin
  ! if two chunks exclude this face for one of them
  if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then
    if (nspec2D_teleseismic_xmin > 0) then
      call read_teleseismic(0,field_teleseismic_xmin,reclen_teleseismic_xmin,it_tmp)
    endif
    if (.not. GPU_MODE) then
      ! on CPU
      do ispec2D = 1,nspec2D_teleseismic_xmin
        ispec = ibelm_teleseismic_xmin(ispec2D)
        i = 1
        do k = 1,NGLLZ
          do j = 1,NGLLY
            iglob = ibool_crust_mantle(i,j,k,ispec)
            b_accel_crust_mantle(:,iglob) = b_accel_crust_mantle(:,iglob) &
              + field_teleseismic_xmin(:,j,k,ispec2D)
          enddo
        enddo
      enddo
    else
      ! on GPU
      if (nspec2D_teleseismic_xmin > 0 ) then
        call compute_add_teleseismic_sources_backward_gpu(Mesh_pointer, field_teleseismic_xmin, 0)
      endif
    endif
  endif ! NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC

  !   xmax
  ! if two chunks exclude this face for one of them
  if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then
    if (nspec2D_teleseismic_xmax > 0) then
      call read_teleseismic(0,field_teleseismic_xmax,reclen_teleseismic_xmax,it_tmp)
    endif
    if (.not. GPU_MODE) then ! on CPU
      do ispec2D = 1,nspec2D_teleseismic_xmax
        ispec = ibelm_teleseismic_xmax(ispec2D)
        i = NGLLX
        do k = 1,NGLLZ
          do j = 1,NGLLY
            iglob = ibool_crust_mantle(i,j,k,ispec)
            b_accel_crust_mantle(:,iglob) = b_accel_crust_mantle(:,iglob) &
              + field_teleseismic_xmax(:,j,k,ispec2D)
          enddo
        enddo
      enddo
    else ! on GPU
      if (nspec2D_teleseismic_xmax > 0 ) then
        call compute_add_teleseismic_sources_backward_gpu(Mesh_pointer, field_teleseismic_xmax, 1)
      endif
    endif
  endif ! NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB

  ! ymin
  if (nspec2D_teleseismic_ymin > 0) then
    call read_teleseismic(0,field_teleseismic_ymin,reclen_teleseismic_ymin,it_tmp)
  endif
  if (.not. GPU_MODE) then
    ! on CPU
    do ispec2D = 1,nspec2D_teleseismic_ymin
      ispec = ibelm_teleseismic_ymin(ispec2D)
      j = 1
      do k = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          b_accel_crust_mantle(:,iglob) = b_accel_crust_mantle(:,iglob) &
            + field_teleseismic_ymin(:,j,k,ispec2D)
        enddo
      enddo
    enddo
  else
    ! on GPU
    if (nspec2D_teleseismic_ymin > 0 ) then
      call compute_add_teleseismic_sources_backward_gpu(Mesh_pointer, field_teleseismic_ymin, 2)
    endif
  endif

  ! ymax
  if (nspec2D_teleseismic_ymax > 0) then
    call read_teleseismic(0,field_teleseismic_ymax,reclen_teleseismic_ymin,it_tmp)
  endif
  if (.not. GPU_MODE) then
    ! on CPU
    do ispec2D = 1,nspec2D_teleseismic_ymax
      ispec = ibelm_teleseismic_ymax(ispec2D)
      j = NGLLY
      do k = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          b_accel_crust_mantle(:,iglob) = b_accel_crust_mantle(:,iglob) &
            + field_teleseismic_ymax(:,j,k,ispec2D)
        enddo
      enddo
    enddo
  else
    ! on GPU
    if (nspec2D_teleseismic_ymax > 0 ) then
      call compute_add_teleseismic_sources_backward_gpu(Mesh_pointer, field_teleseismic_ymax, 3)
    endif
  endif

  ! zmin
  if (nspec2D_teleseismic_zmin > 0) then
    call read_teleseismic(0,field_teleseismic_zmin,reclen_teleseismic_ymin,it_tmp)
  endif
  if (.not. GPU_MODE) then
    ! on CPU
    do ispec2D = 1,nspec2D_teleseismic_zmin
      ispec = ibelm_teleseismic_zmin(ispec2D)
      k = 1
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          b_accel_crust_mantle(:,iglob) = b_accel_crust_mantle(:,iglob) &
            + field_teleseismic_zmin(:,j,k,ispec2D)
        enddo
      enddo
    enddo
  else
    ! on GPU
    if (nspec2D_teleseismic_zmin > 0 ) then
      call compute_add_teleseismic_sources_backward_gpu(Mesh_pointer, field_teleseismic_zmin, 4)
    endif
  endif

  end subroutine compute_add_teleseismic_sources_backward
