!> \file mo_mrm_mpr.f90

!> \brief Perform Multiscale Parameter Regionalization on Routing Parameters

!> \details This module contains the subroutine for calculating the regionalized
!> routing parameters (beta-parameters) given the five global routing parameters
!> (gamma) at the level 0 scale.

!> \author Luis Samaniego, Stephan Thober
!> \date Aug 2015
module mo_mrm_mpr
  use mo_kind, only: dp
  implicit none
  public :: reg_rout
  private
contains
  
  ! ----------------------------------------------------------------------------

  !      NAME
  !         reg_rout

  !>        \brief Regionalized routing

  !>        \details sets up the Regionalized Routing parameters\n
  !>                 Global parameters needed (see mhm_parameter.nml):\n
  !>                    - param(1) = muskingumTravelTime_constant    \n
  !>                    - param(2) = muskingumTravelTime_riverLength \n
  !>                    - param(3) = muskingumTravelTime_riverSlope  \n
  !>                    - param(4) = muskingumTravelTime_impervious  \n
  !>                    - param(5) = muskingumAttenuation_riverSlope \n

  !      INTENT(IN)
  !>        \param[in] "real(dp) :: param(5)"  - five input parameters
  !>        \param[in] "real(dp) :: length(:)" - [m] total length
  !>        \param[in] "real(dp) :: slope(:)"  - average slope
  !>        \param[in] "real(dp) :: fFPimp(:)" - fraction of the flood plain with
  !>                                             impervious layer
  !>        \param[in] "real(dp) :: TS"        - [h] time step in

  !      INTENT(INOUT)
  !          None
  
  !      INTENT(OUT)
  !>        \param[out] "real(dp) :: C1(:)"    - routing parameter C1 (Chow, 25-41)
  !>        \param[out] "real(dp) :: C2(:)"    - routing parameter C2 (")

  !      INTENT(IN), OPTIONAL
  !          None

  !      INTENT(INOUT), OPTIONAL
  !          None

  !      INTENT(OUT), OPTIONAL
  !          None

  !      RETURN
  !          None

  !      RESTRICTIONS
  !          None

  !      EXAMPLE
  !          None

  !      LITERATURE
  !          None

  !      HISTORY
  !>        \author Stephan Thober, Rohini Kumar
  !>        \date Dec 2012
  !         Written  Stephan Thober, Dec 2012
  !         Modified Xiaoqiang Yang, Aug 2016  !deleted K constrains and changed the time unit from hourly to daily. 
                                               !temporally for daily nitrate sumilation.

  subroutine reg_rout( param, length, slope, fFPimp, TS, &
       C1, C2 )

    implicit none

    ! Input
    real(dp), dimension(5), intent(in)  :: param  ! input parameter
    real(dp), dimension(:), intent(in)  :: length ! [m] total length
    real(dp), dimension(:), intent(in)  :: slope  ! average slope
    real(dp), dimension(:), intent(in)  :: fFPimp ! fraction of the flood plain with
    !                                                ! impervious layer
    real(dp),               intent(in)  :: TS     ! [h] time step in

    ! Output
    real(dp), dimension(:), intent(out) :: C1     ! routing parameter C1 (Chow, 25-41)
    real(dp), dimension(:), intent(out) :: C2     ! routing parameter C2 (")

    ! local variables
    real(dp)                            :: ssMax  ! stream slope max
    real(dp), dimension(size(fFPimp,1)) :: K      ! [d] Muskingum travel time parameter
    real(dp), dimension(size(fFPimp,1)) :: xi     ! [1] Muskingum diffusion parameter (attenuation)
	!for daily simulation-- added by yangx 2018-08-05
    real(dp)                            :: tds    ![d] timestep transform unit
    tds = TS/24.0_dp
	
    ! ! normalize stream bed slope
    ! ssMax = maxval( slope(:) )

    ! ! New regional relationship; K = f(length, slope, & fFPimp)
    ! K = param(1) + param(2) * (length * 0.001_dp) &
         ! + param(3) * slope &
         ! + param(4) * fFPimp

    ! ! Xi = f(slope)
    ! xi = param(5)*(1.0_dp + slope / ssMax)

    ! ! constraints on Xi
    ! xi = merge( 0.5_dp, xi, xi > 0.5_dp )
    ! xi = merge( 0.005_dp, xi, xi < 0.005_dp )

    ! constrains on Ki
    !K = merge( 0.5_dp * TS / xi,            K, K > 0.5_dp * TS / xi )
    !K = merge( 0.5_dp * TS / (1.0_dp - xi), K, K < 0.5_dp * TS / (1.0_dp - xi))

    ! Muskingum parameters
    !C1 = tds / ( K * (1.0_dp - xi) + 0.5_dp * tds )
    !C2 = 1.0_dp - C1 * K / tds

    !A simple linear- reservoir, delay based method is used for "daily time step + 1km routing"
	C1 = param(1) !the coefficient Q2 = C1*I2+ (1-C1)*Q1
    C2 = 0




  end subroutine reg_rout
end module mo_mrm_mpr
