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


  subroutine prepare_teleseismic_incidence()

! sets up arrays for teleseismic incidence

  use specfem_par
  !use shared_parameters, only: NSTEP,TELESEISMIC_INCIDENCE
  use specfem_par_crustmantle

  implicit none
  ! local parameters
  integer :: ier
  integer :: filesize

  ! checks if anything to do
  if (.not. TELESEISMIC_INCIDENCE) return

  ! sets up absorbing boundary buffer arrays
  if (myrank == 0) then
    write(IMAIN,*) "preparing teleseismic incidence"
    call flush_IMAIN()
  endif

  ! crust_mantle
  ! create name of database
  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

  ! allocates buffers
  allocate(field_teleseismic_xmin(NDIM,NGLLY,NGLLZ,nspec2D_teleseismic_xmin), &
           field_teleseismic_xmax(NDIM,NGLLY,NGLLZ,nspec2D_teleseismic_xmax), &
           field_teleseismic_ymin(NDIM,NGLLX,NGLLZ,nspec2D_teleseismic_ymin), &
           field_teleseismic_ymax(NDIM,NGLLX,NGLLZ,nspec2D_teleseismic_ymax), &
           field_teleseismic_zmin(NDIM,NGLLX,NGLLY,nspec2D_teleseismic_zmin), &
           stat=ier)
  if (ier /= 0 ) then
    call exit_MPI(myrank,'Error allocating field_teleseismic_xmin... arrays')
  endif

  ! file I/O for teleseismic wavefields
  if (nspec2D_teleseismic_xmin > 0) then
    ! size of single record
    reclen_teleseismic_xmin = CUSTOM_REAL*NDIM*NGLLY*NGLLZ*nspec2D_teleseismic_xmin
    ! total file size
    filesize = reclen_teleseismic_xmin*NSTEP
    if (SIMULATION_TYPE == 2) then
      call open_file_teleseismic_w(0, &
        trim(prname)//'field_teleseismic_xmin.bin', &
        len_trim(trim(prname)//'field_teleseismic_xmin.bin'), filesize)
    else
      call open_file_teleseismic_r(0, &
        trim(prname)//'field_teleseismic_xmin.bin', &
        len_trim(trim(prname)//'field_teleseismic_xmin.bin'), filesize, myrank)
    endif
  endif

  if (nspec2D_teleseismic_xmax > 0) then
    ! size of single record
    reclen_teleseismic_xmax = CUSTOM_REAL*NDIM*NGLLY*NGLLZ*nspec2D_teleseismic_xmax
    ! total file size
    filesize = reclen_teleseismic_xmax*NSTEP
    if (SIMULATION_TYPE == 2) then
      call open_file_teleseismic_w(1, &
        trim(prname)//'field_teleseismic_xmax.bin', &
        len_trim(trim(prname)//'field_teleseismic_xmax.bin'), filesize)
    else
      call open_file_teleseismic_r(1, &
        trim(prname)//'field_teleseismic_xmax.bin', &
        len_trim(trim(prname)//'field_teleseismic_xmax.bin'), filesize, myrank)
    endif
  endif

  if (nspec2D_teleseismic_ymin > 0) then
    ! size of single record
    reclen_teleseismic_ymin = CUSTOM_REAL*NDIM*NGLLX*NGLLZ*nspec2D_teleseismic_ymin
    ! total file size
    filesize = reclen_teleseismic_ymin*NSTEP
    if (SIMULATION_TYPE == 2) then
      call open_file_teleseismic_w(2, &
        trim(prname)//'field_teleseismic_ymin.bin', &
        len_trim(trim(prname)//'field_teleseismic_ymin.bin'), filesize)
    else
      call open_file_teleseismic_r(2, &
        trim(prname)//'field_teleseismic_ymin.bin', &
        len_trim(trim(prname)//'field_teleseismic_ymin.bin'), filesize, myrank)
    endif
  endif

  if (nspec2D_teleseismic_ymax > 0) then
    ! size of single record
    reclen_teleseismic_ymax = CUSTOM_REAL*NDIM*NGLLX*NGLLZ*nspec2D_teleseismic_ymax
    ! total file size
    filesize = reclen_teleseismic_ymax*NSTEP
    if (SIMULATION_TYPE == 2) then
      call open_file_teleseismic_w(3, &
        trim(prname)//'field_teleseismic_ymax.bin', &
        len_trim(trim(prname)//'field_teleseismic_ymax.bin'), filesize)
    else
      call open_file_teleseismic_r(3, &
        trim(prname)//'field_teleseismic_ymax.bin', &
        len_trim(trim(prname)//'field_teleseismic_ymax.bin'), filesize, myrank)
    endif
  endif

  if (nspec2D_teleseismic_zmin > 0) then
    ! size of single record
    reclen_teleseismic_zmin = CUSTOM_REAL*NDIM*NGLLX*NGLLZ*nspec2D_teleseismic_zmin
    ! total file size
    filesize = reclen_teleseismic_zmin*NSTEP
    if (SIMULATION_TYPE == 2) then
      call open_file_teleseismic_w(4, &
        trim(prname)//'field_teleseismic_zmin.bin', &
        len_trim(trim(prname)//'field_teleseismic_zmin.bin'), filesize)
    else
      call open_file_teleseismic_r(4, &
        trim(prname)//'field_teleseismic_zmin.bin', &
        len_trim(trim(prname)//'field_teleseismic_zmin.bin'), filesize, myrank)
    endif
  endif

!  ! read teleseismic source
!  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
!    ! user output
!    if (myrank == 0) then
!      write(IMAIN,*) "read teleseismic_source"
!      call flush_IMAIN()
!    endif
!
!    ! read in teleseismic source files
!    call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)
!    open(unit=IIN,file=prname(1:len_trim(prname))//'teleseismic_source.bin', &
!          status='old',form='unformatted',action='read',iostat=ier)
!    if (ier /= 0 ) call exit_mpi(myrank,'Error opening teleseismic_source.bin file')
!
!    read(IOUT) teleseismic_source_xmin
!    read(IOUT) teleseismic_source_xmax
!    read(IOUT) teleseismic_source_ymin
!    read(IOUT) teleseismic_source_ymax
!    read(IOUT) teleseismic_source_zmin
!  endif
!
!
!  ! allocate arrays for teleseismic gradient
!  if (SIMULATION_TYPE == 2) then
!    ! user output
!    if (myrank == 0) then
!      write(IMAIN,*) "allocate arrays for teleseismic gradient"
!      call flush_IMAIN()
!    endif
!
!    ! allocates arrays for teleseismic gradient
!    allocate(teleseismic_gradient_xmin(NDIM,NGLLY,NGLLZ,nspec2D_teleseismic_xmin,NSTEP), &
!             teleseismic_gradient_xmax(NDIM,NGLLY,NGLLZ,nspec2D_teleseismic_xmax,NSTEP), &
!             teleseismic_gradient_ymin(NDIM,NGLLX,NGLLZ,nspec2D_teleseismic_ymin,NSTEP), &
!             teleseismic_gradient_ymax(NDIM,NGLLX,NGLLZ,nspec2D_teleseismic_ymax,NSTEP), &
!             teleseismic_gradient_zmin(NDIM,NGLLX,NGLLY,nspec2D_teleseismic_zmin,NSTEP), stat=ier)
!    if (ier /= 0 ) then
!      call exit_MPI(myrank,'Error allocating teleseismic_gradient_xmin... arrays')
!    endif
!  endif

  call synchronize_all()

  end subroutine prepare_teleseismic_incidence

