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


  subroutine compute_teleseismic_gradient()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: ispec2D,ispec,i,j,k,iglob
  integer :: it_adj

  ! checks
  if (SIMULATION_TYPE /= 2 ) return

  if (.not. GPU_MODE) then ! on CPU
    !xmin
    do ispec2D = 1, nspec2D_teleseismic_xmin
      ispec = ibelm_teleseismic_xmin(ispec2D)
      ! loop on all the points inside the element
      i = 1
      do k = 1,NGLLZ
        do j = 1,NGLLY
          iglob = ibool_crust_mantle(i,j,k,ispec)
          ! stores teleseismic gradient at each gll point (adjoint displ * area of the gll point)
          field_teleseismic_xmin(:,j,k,ispec2D) = displ_crust_mantle(:,iglob) &
            * scale_displ * area_teleseismic_xmin(j,k,ispec2D)
        enddo
      enddo
    enddo
    !xmax
    do ispec2D = 1, nspec2D_teleseismic_xmax
      ispec = ibelm_teleseismic_xmax(ispec2D)
      ! loop on all the points inside the element
      i = NGLLX
      do k = 1,NGLLZ
        do j = 1,NGLLY
          iglob = ibool_crust_mantle(i,j,k,ispec)
          ! stores displacement
          field_teleseismic_xmax(:,j,k,ispec2D) = displ_crust_mantle(:,iglob) &
            * scale_displ * area_teleseismic_xmax(j,k,ispec2D)
        enddo
      enddo
    enddo
    !ymin
    do ispec2D = 1, nspec2D_teleseismic_ymin
      ispec = ibelm_teleseismic_ymin(ispec2D)
      ! loop on all the points inside the element
      j = 1
      do k = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          ! stores displacement
          field_teleseismic_ymin(:,i,k,ispec2D) = displ_crust_mantle(:,iglob) &
            * scale_displ * area_teleseismic_ymin(i,k,ispec2D)
        enddo
      enddo
    enddo
    !ymax
    do ispec2D = 1, nspec2D_teleseismic_ymax
      ispec = ibelm_teleseismic_ymax(ispec2D)
      ! loop on all the points inside the element
      j = NGLLY
      do k = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          ! stores displacement
          field_teleseismic_ymax(:,i,k,ispec2D) = displ_crust_mantle(:,iglob) &
            * scale_displ * area_teleseismic_ymax(i,k,ispec2D)
        enddo
      enddo
    enddo
    !zmin
    do ispec2D = 1, nspec2D_teleseismic_zmin
      ispec = ibelm_teleseismic_zmin(ispec2D)
      ! loop on all the points inside the element
      k = 1
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          ! stores displacement
          field_teleseismic_zmin(:,i,j,ispec2D) = displ_crust_mantle(:,iglob) &
            * scale_displ * area_teleseismic_zmin(i,j,ispec2D)
        enddo
      enddo
    enddo

  else ! on GPU
    if (nspec2D_teleseismic_xmin > 0) then
      call compute_teleseismic_gradient_gpu( &
        Mesh_pointer,field_teleseismic_xmin,scale_displ,0)
    endif
    if (nspec2D_teleseismic_xmax > 0) then
      call compute_teleseismic_gradient_gpu( &
        Mesh_pointer,field_teleseismic_xmax,scale_displ,1)
    endif
    if (nspec2D_teleseismic_ymin > 0) then
      call compute_teleseismic_gradient_gpu( &
        Mesh_pointer,field_teleseismic_ymin,scale_displ,2)
    endif
    if (nspec2D_teleseismic_ymax > 0) then
      call compute_teleseismic_gradient_gpu( &
        Mesh_pointer,field_teleseismic_ymax,scale_displ,3)
    endif
    if (nspec2D_teleseismic_zmin > 0) then
      call compute_teleseismic_gradient_gpu( &
        Mesh_pointer,field_teleseismic_zmin,scale_displ,4)
    endif
  endif

  ! writes field_teleseismic values

  ! since adjoint simulation is in reversed time order, writing from back
  it_adj = NSTEP - it + 1
  ! xmin
  if (nspec2D_teleseismic_xmin > 0) then
    call write_teleseismic(0,field_teleseismic_xmin,reclen_teleseismic_xmin,it_adj)
  endif
  ! xmax
  if (nspec2D_teleseismic_xmax > 0) then
    call write_teleseismic(1,field_teleseismic_xmax,reclen_teleseismic_xmax,it_adj)
  endif
  ! ymin
  if (nspec2D_teleseismic_ymin > 0) then
    call write_teleseismic(2,field_teleseismic_ymin,reclen_teleseismic_ymin,it_adj)
  endif
  ! ymax
  if (nspec2D_teleseismic_ymax > 0) then
    call write_teleseismic(3,field_teleseismic_ymax,reclen_teleseismic_ymax,it_adj)
  endif
  ! zmin
  if (nspec2D_teleseismic_zmin > 0) then
    call write_teleseismic(4,field_teleseismic_zmin,reclen_teleseismic_zmin,it_adj)
  endif

  !! gets resulting array values onto CPU
  !if (GPU_MODE) then
  !  call transfer_displ_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,displ_crust_mantle,Mesh_pointer)
  !endif

  !! save files
  !if (it == NSTEP) then
  !  ! create name of database
  !  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

  !  open(unit=IOUT,file=prname(1:len_trim(prname))//'teleseismic_gradient.bin', &
  !        status='unknown',form='unformatted',action='write',iostat=ier)
  !  if (ier /= 0 ) call exit_mpi(myrank,'Error opening teleseismic_gradient.bin file')

  !  write(IOUT) teleseismic_gradient_xmin
  !  write(IOUT) teleseismic_gradient_xmax
  !  write(IOUT) teleseismic_gradient_ymin
  !  write(IOUT) teleseismic_gradient_ymax
  !  write(IOUT) teleseismic_gradient_zmin

  !  close(IOUT)

  !endif

  end subroutine compute_teleseismic_gradient
