! The model, entitled mHM-Nitrate model v2.0, is a fully distributed nitrate 
! transport and removal model. Please refer to "README.md" for more instructions   	

! Copyright (C) 2020,
! Xiaoqiang Yang and Michael Rode
! Helmholtz Centre for Environmental Research - UFZ.
! All rights reserved.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

!> \file mo_wqm_soiltexture.f90

!> \brief calculate coefficient of shading effect in riparian zone.

!> \details This module calculated the shading effect of riparizan zone according to different land use types. 
!>          the coefficient was generated from the global radiation (daily timeseries for the whole catchment) 
!>          and LAI data (land use dependent monthly data from loop up table).\n


!> \authors Xiaoqiang Yang
!> \date Nov 2016

MODULE mo_wqm_soiltexture


  USE mo_kind, ONLY: i4, sp, dp

  
  IMPLICIT NONE

  PUBLIC :: L1_soiltexture               ! 


  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         L1_soiltexture

  !     PURPOSE
  !>        \brief .

  !>        \details 

  !     CALLING SEQUENCE
  !         

  !     INTENT(IN)
  !>        None

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None


  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Nov 2016

  subroutine L1_soiltexture(fsoil)

  use mo_julian,    only: julday
  use mo_global_variables, only: soilDB
  
  implicit none

  real(dp), dimension(:,:),             intent(in)       :: fsoil  
  !local
  real(dp), dimension(:), allocatable   :: fsdsoil    !sandy share
  real(dp), dimension(:), allocatable   :: fclsoil    !clay share
  real(dp), dimension(:), allocatable   :: fbdsoil    !bulk density  

 !local
  integer(i4)   :: ncells_1, nsoilty    
  real(dp)      :: mm, nn,kk
  integer(i4)   :: i,j 
  
  ncells_1 = size(fsoil,1)
  nsoilty = size(fsoil,2)
  allocate(fsdsoil(ncells_1))
  allocate(fclsoil(ncells_1))
  allocate(fbdsoil(ncells_1))

  fsdsoil = 0.0_dp
  fclsoil = 0.0_dp
  fbdsoil = 0.0_dp
  
  do i=1,ncells_1 !each grid at L1
     do j = 1, nsoilty  !each soil type 
         
         mm = sum(soilDB%sand(j,:) * (soilDB%LD(j,:)-soilDB%UD(j,:))) &
              /maxval(soilDB%LD(j,:))
         nn = sum(soilDB%clay(j,:) * (soilDB%LD(j,:)-soilDB%UD(j,:))) &
              /maxval(soilDB%LD(j,:))
         kk = sum(soilDB%DbM(j,:) * (soilDB%LD(j,:)-soilDB%UD(j,:))) &
              /maxval(soilDB%LD(j,:))
  
         fsdsoil(i) = fsdsoil(i) + fsoil(i,j) * mm
         fclsoil(i) = fclsoil(i) + fsoil(i,j) * nn
         fbdsoil(i) = fbdsoil(i) + fsoil(i,j) * kk   

     end do
  end do


  
  end subroutine L1_soiltexture
  
  
  

END MODULE mo_wqm_soiltexture