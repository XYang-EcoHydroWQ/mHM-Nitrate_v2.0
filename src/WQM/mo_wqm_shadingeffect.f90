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

!> \file mo_wqm_shadingeffect.f90

!> \brief calculate coefficient of shading effect in riparian zone.

!> \details This module calculated the shading effect of riparizan zone according to different land use types. 
!>          the coefficient was generated from the global radiation (daily timeseries for the whole catchment) 
!>          and LAI data (land use dependent monthly data from loop up table).\n


!> \authors Xiaoqiang Yang
!> \date Nov 2016

MODULE mo_wqm_shadingeffect


  USE mo_kind, ONLY: i4, sp, dp

  
  IMPLICIT NONE

  PUBLIC :: rz_shading_coeff               ! 
  PUBLIC :: read_shading_data

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         read_shading_data

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

  subroutine read_shading_data(nbasins, global_radi, LAI, lai_daily, gr_nor)

  use mo_julian,    only: julday
  
  implicit none
  integer(i4),                          intent(in)       :: nbasins       !number of basins in total
  real(dp), dimension(:,:),             intent(in)       :: global_radi   !global radiation
  real(dp), dimension(:,:),             intent(in)       :: LAI           !monthly LAI value
  real(dp), dimension(:,:), allocatable, intent(inout)   :: lai_daily  !daily normalised lai after interpolating
  real(dp), dimension(:,:), allocatable, intent(inout)   :: gr_nor    ! normalised daily global radiation data
!  real(dp), dimension(:,:,:), allocatable, intent(inout) :: rz_coeff      !riparian zone coefficient

  !local
  integer(i4), dimension(12)            :: mondays     !for LAI interpolation, number of days in the middle of each monthly
  real(dp), dimension(:,:), allocatable       :: lai_nor   ! normalised monthly LAI values
  integer(i4)                           :: i, nn, ii, mm
  
  allocate(lai_nor(size(LAI,1), size(LAI,2)))
  allocate(lai_daily(size(LAI,1), 366))
  allocate(gr_nor(size(global_radi,1), size(global_radi,2)))
  
  !normalise LAI value where it is forest
  lai_nor = (LAI -minval(LAI(1:4,:)))/(maxval(LAI(1:4,:) )-minval(LAI(1:4,:)))
  !interpolate monthy LAI to daily  
  do i=1, 12
     mondays(i) = julday(15, i, 2017) - julday(1,1,2017) !"2017" can be changed
  end do

  do i=1,366             !day loop
  do nn =1, size(LAI,1)  !nLAIclass loop
  if (nn .lt. 5) then    !if land use is forest (id 1:4)
     if (i <= mondays(1) ) then
        lai_daily(nn, i) = lai_nor(nn, 1)
     elseif (i > mondays(12)) then
        lai_daily(nn, i) = lai_nor(nn, 12)
     else
     do mm =2, 12
        if ((i > mondays(mm-1)) .and. (i <= mondays(mm) )) then
           lai_daily(nn,i) = lai_nor(nn,mm-1) + (lai_nor(nn,mm)-lai_nor(nn,mm-1)) * (i-mondays(mm-1))/(mondays(mm)-mondays(mm-1)) 
           exit 
        end if
     end do
     end if
  else                   !if land use are not forest, no shading in stream
     lai_daily(nn,i) = 0.0_dp
  end if
  end do !nLAIclass loop
  end do !day loop
  !normalise global radiation  
  do ii=1, nbasins
     gr_nor(ii,:) = (global_radi(ii,:) - minval(global_radi(ii,:)))/(maxval(global_radi(ii,:)) - minval(global_radi(ii,:))) 
  end do
!open(223, file = "gr_data.txt")
!write(223,*) gr_nor(1,:)
!close(223)
  end subroutine read_shading_data
!--------------------------------------------------------  

  ! ------------------------------------------------------------------

  !     NAME
  !         rz_shading_coeff

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

  subroutine rz_shading_coeff(i, noday, nor_gr, rz_coeff_i)

  use mo_wqm_global_variables,   only: norlai_daily, L11_fLAI
  implicit none
  integer(i4),            intent(in)    :: i             !reach id 
  integer(i4),            intent(in)    :: noday
  real(dp),               intent(in)    :: nor_gr 
  real(dp),               intent(inout) :: rz_coeff_i    !shading coefficient in the 'i'th reach

  !local
  real(dp)               :: lai_coe   !total lai coefficient, considering weight of land use area
  integer(i4)            :: nn
  
  !lai_coe = sum(norlai_daily(:, noday)* L11_fLAI(i,:), )
  lai_coe = 0.0_dp
  do nn =1, size(norlai_daily, 1)
     lai_coe = lai_coe + norlai_daily(nn,noday)* L11_fLAI(i, nn)
  end do
  !
!  if (sum(L11_fLAI(i,:) ) > 1 ) then
!  write(23,*)  i, sum(L11_fLAI(i,:) )
!  end if
  rz_coeff_i = nor_gr * (1- lai_coe)
!  if (i==11) then
!  write(40,*) lai_coe, nor_gr, rz_coeff_i
!  write(23,*) L11_fLAI(i,:)
!  end if

  
  end subroutine rz_shading_coeff
  

END MODULE mo_wqm_shadingeffect