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

  subroutine get_jacobian_teleseismic_zmin(ispec,r1,r2,r3,r4, &
      xstore,ystore,zstore,dershape2D_bottom)

  use constants

  !>>>KTAO 
  use shared_parameters, only: TELESEISMIC_INCIDENCE
  use regions_mesh_par, only: wgllwgll_xy
  use regions_mesh_par2, only: &
    r_teleseismic_zmin, &
    ispec2D_teleseismic_zmin,nspec2D_teleseismic_zmin, &
    ibelm_teleseismic_zmin, &
    area_teleseismic_zmin,normal_teleseismic_zmin
  !<<<

  implicit none

  ! input
  integer :: ispec
  double precision :: r1, r2, r3, r4
  double precision :: xstore(NGLLX,NGLLY,NGLLZ)
  double precision :: ystore(NGLLX,NGLLY,NGLLZ)
  double precision :: zstore(NGLLX,NGLLY,NGLLZ)
  double precision :: dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY)

  ! local variables
  double precision, dimension(NGNOD2D) :: xelm2, yelm2, zelm2

  ! ======================

  ! check
  if (.not. TELESEISMIC_INCIDENCE) return

  ! get coordinates of 9 nodes for the bottom surface element to compute the Jacobian
  xelm2(1)=xstore(1,1,1)
  yelm2(1)=ystore(1,1,1)
  zelm2(1)=zstore(1,1,1)
  xelm2(2)=xstore(NGLLX,1,1)
  yelm2(2)=ystore(NGLLX,1,1)
  zelm2(2)=zstore(NGLLX,1,1)
  xelm2(3)=xstore(NGLLX,NGLLY,1)
  yelm2(3)=ystore(NGLLX,NGLLY,1)
  zelm2(3)=zstore(NGLLX,NGLLY,1)
  xelm2(4)=xstore(1,NGLLY,1)
  yelm2(4)=ystore(1,NGLLY,1)
  zelm2(4)=zstore(1,NGLLY,1)
  xelm2(5)=xstore((NGLLX+1)/2,1,1)
  yelm2(5)=ystore((NGLLX+1)/2,1,1)
  zelm2(5)=zstore((NGLLX+1)/2,1,1)
  xelm2(6)=xstore(NGLLX,(NGLLY+1)/2,1)
  yelm2(6)=ystore(NGLLX,(NGLLY+1)/2,1)
  zelm2(6)=zstore(NGLLX,(NGLLY+1)/2,1)
  xelm2(7)=xstore((NGLLX+1)/2,NGLLY,1)
  yelm2(7)=ystore((NGLLX+1)/2,NGLLY,1)
  zelm2(7)=zstore((NGLLX+1)/2,NGLLY,1)
  xelm2(8)=xstore(1,(NGLLY+1)/2,1)
  yelm2(8)=ystore(1,(NGLLY+1)/2,1)
  zelm2(8)=zstore(1,(NGLLY+1)/2,1)
  xelm2(9)=xstore((NGLLX+1)/2,(NGLLY+1)/2,1)
  yelm2(9)=ystore((NGLLX+1)/2,(NGLLY+1)/2,1)
  zelm2(9)=zstore((NGLLX+1)/2,(NGLLY+1)/2,1)

  !only record elements whose bottom lie on teleseismic bottom boundary
  if (    abs(r1-r_teleseismic_zmin)/r_teleseismic_zmin < SMALLVAL &
    .and. abs(r2-r_teleseismic_zmin)/r_teleseismic_zmin < SMALLVAL &
    .and. abs(r3-r_teleseismic_zmin)/r_teleseismic_zmin < SMALLVAL &
    .and. abs(r4-r_teleseismic_zmin)/r_teleseismic_zmin < SMALLVAL ) then

    ispec2D_teleseismic_zmin = ispec2D_teleseismic_zmin + 1
    ibelm_teleseismic_zmin(ispec2D_teleseismic_zmin) = ispec

    call compute_jacobian_2D(ispec2D_teleseismic_zmin, &
      xelm2,yelm2,zelm2,dershape2D_bottom, & 
      area_teleseismic_zmin,normal_teleseismic_zmin, &
      NGLLX,NGLLY,nspec2D_teleseismic_zmin)

    area_teleseismic_zmin(:,:,ispec2D_teleseismic_zmin) = &
      area_teleseismic_zmin(:,:,ispec2D_teleseismic_zmin) * wgllwgll_xy
  endif

end subroutine get_jacobian_teleseismic_zmin
