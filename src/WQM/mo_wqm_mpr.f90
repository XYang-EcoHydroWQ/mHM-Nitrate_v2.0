!> \file mo_wqm_mpr.f90

!> \brief Upscale land-use dependent parameters to model levels

!> \details 

!> \author Xiaoqiang Yang
!> \date Jul 2017 

MODULE mo_wqm_mpr

 
  USE mo_kind, ONLY: i4, sp, dp


  IMPLICIT NONE

 
  PUBLIC :: wqm_mpr  ! parameters of water quality model


CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         wqm_mpr

  !     PURPOSE
  !>        \brief Water quality parameters 

  !>        \details Water quality parameters are previously read-in from "mhm_parameter.nml" and stored in "global_parameter".
  !>        This module assigns those values to parameter variables in each L1 cell. \n
  !>        Notes: WQ parameter regionalisation has not been developed yet. Same value is assigned to all L1 cells. 
  !>               Currently, we take the five soil phase parameters as land-use dependent parameter, and the land use 
  !>               differences are read-in as input initial values (see mo_wqm_read::wqm_readdata and "initial_value.txt").
  !     CALLING SEQUENCE
  !         

  !     INTENT(IN)
  !>        \param[in] "integer(i4), dimension(:,:)        :: procMat"     process setting dim1=11 for water quality
  !>        \param[in] "real(dp), dimension(:)             :: param"       global parameters from parameter namelist        

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !>        \param[out] "real(dp), dimension(:)             :: degradN_rate"   parameter 
  !>        \param[out] "real(dp), dimension(:)             :: mineraN_rate"   parameter   
  !>        \param[out] "real(dp), dimension(:)             :: dissolN_rate"   parameter   
  !>        \param[out] "real(dp), dimension(:)             :: sdenitr_rate"   parameter   
  !>        \param[out] "real(dp), dimension(:)             :: gwresid_time"   parameter  
  !>        \param[out] "real(dp), dimension(:)             :: adenitr_rate"   parameter   
  !>        \param[out] "real(dp), dimension(:)             :: priprod_rate"   parameter   

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>       

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         
  !         

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Jun 2017

  subroutine wqm_mpr( nNodes, L1id_on_L11, L11id_on_L1, map_flag, procMat,param, fLAI, fLAI_11, degradN_rate, &
                     mineraN_rate, dissolN_rate, sdenitr_rate, adenitr_rate, priprod_rate, nLink_from, &
                     areaCell1, areaCell11 )

  use mo_message,   only: message
  implicit none
  !integer(i4),                       intent(in)     :: ncells1      ! number of cells at model level (L1)
  integer(i4),                       intent(in)     :: nNodes       ! number of cells at routing level (L11)
  integer(i4), dimension(:),         intent(in)     :: L1id_on_L11  ! L1 cell id on L11 level
  integer(i4), dimension(:),         intent(in)     :: L11id_on_L1  ! L11 cell id on L1 level
  logical,                           intent(in)     :: map_flag     ! L11 is bigger than L1 or not
  integer(i4), dimension(:,:),       intent(in)     :: procMat      ! process setting dim1=11 for water quality
  real(dp), dimension(:),            intent(in)     :: param        ! global parameters from parameter namelist
  real(dp), dimension(:,:),          intent(in)     :: fLAI         ! area fraction of each landuse type
  real(dp), dimension(:,:),          intent(inout)  :: fLAI_11      ! area fraction of each landuse type at L11
  real(dp), dimension(:),            intent(out)    :: degradN_rate ! parameter 
  real(dp), dimension(:),            intent(out)    :: mineraN_rate ! parameter   
  real(dp), dimension(:),            intent(out)    :: dissolN_rate ! parameter   
  real(dp), dimension(:),            intent(out)    :: sdenitr_rate ! parameter    
  real(dp), dimension(:),            intent(out)    :: adenitr_rate ! parameter   
  real(dp), dimension(:),            intent(out)    :: priprod_rate ! parameter
  integer(i4), dimension(:),         intent(in)     :: nLink_from   ! id of Node where current reach connects from
  real(dp), dimension(:),            intent(in)     :: areaCell1    ! cell area at L1
  real(dp), dimension(:),            intent(in)     :: areaCell11   ! cell area at 	L11

  !local
  integer(i4)    :: istart, iend
  integer(i4)    :: k, nn, mm, ii
  !integer(i4), dimension(size(nLink_from))  :: locating  
  !real(dp), dimension(nNodes, size(fLAI,2))  :: tmp_fLAI_11   ! area fraction of different land use at L11, temporal variable

  istart = 1_i4
  iend   = procMat(11,3) - procMat(10,3)

  if (iend .ne. (istart +10)) then
     call message()
     call message('***ERROR: the number of nitrogen parameters in mhm_parameter.nml does not match the model structure!')
	 
  stop
  end if

  !general parameter
  adenitr_rate(:) = param(istart)

  !considering the difference of arable denitri_rate and other denitri_rate
  !instream paramters, at L11 level
  !priprod_rate(:) = param(istart + 1) 
  do k =1, nNodes
     !get the reach id "mm" from node id "k"   where(nLink_from .eq. k)
     !locating(:) = k
     !mm = minloc(abs(nLink_from - k),1)	
     do ii = 1, nNodes -1
     if (nLink_from(ii) == k ) then
     mm = ii  
	
     if (map_flag) then   !L11 >= L1
      do nn=1, size(fLAI,2)
      fLAI_11(mm,nn) = sum(fLAI(:,nn) * areaCell1(:) / areaCell11(k) , L11id_on_L1(:) .eq. k ) !* nNodes / real(ncells1, dp)
      if ((nn == 9) .or. (nn == 7)) then
         priprod_rate(mm) = priprod_rate(mm) + fLAI_11(mm,nn) * param(istart + 2)
      else
         priprod_rate(mm) = priprod_rate(mm) + fLAI_11(mm,nn) * param(istart + 1) 
     end if
      end do

     else         ! L11 < L1
      do nn=1, size(fLAI,2)
      fLAI_11(mm,nn) = fLAI(L1id_on_L11(k),nn)
      if ((nn == 9) .or. (nn == 7)) then
         priprod_rate(mm) = priprod_rate(mm) + fLAI_11(mm,nn) * param(istart + 2)
      else
         priprod_rate(mm) = priprod_rate(mm) + fLAI_11(mm,nn) * param(istart + 1)
      end if
      end do
     end if

     end if
     end do

  !
  end do
  
  ! soil phase, at L1 level
  do k = 1, size(fLAI,1) 
  do nn=1, size(fLAI,2) 
  if ((nn==9) .or. (nn == 7)) then
     degradN_rate(k) = degradN_rate(k) + fLAI(k,nn) * param(istart + 4) 
     mineraN_rate(k) = mineraN_rate(k) + fLAI(k,nn) * param(istart + 6) 
     dissolN_rate(k) = dissolN_rate(k) + fLAI(k,nn) * param(istart + 8) 
     sdenitr_rate(k) = sdenitr_rate(k) + fLAI(k,nn) * param(istart + 10)
  else
     degradN_rate(k) = degradN_rate(k) + fLAI(k,nn) * param(istart + 3)
     mineraN_rate(k) = mineraN_rate(k) + fLAI(k,nn) * param(istart + 5) 
     dissolN_rate(k) = dissolN_rate(k) + fLAI(k,nn) * param(istart + 7) 
     sdenitr_rate(k) = sdenitr_rate(k) + fLAI(k,nn) * param(istart + 9)
  end if
  end do    
  end do
  
  end subroutine wqm_mpr
END MODULE mo_wqm_mpr
