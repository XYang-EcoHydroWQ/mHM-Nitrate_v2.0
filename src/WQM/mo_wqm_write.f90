!> \file mo_wqm_write.f90

!> \brief write out water quality model information

!> \details 

!> \author Xiaoqiang Yang
!> \date  Sep 2017

MODULE mo_wqm_write


  USE mo_kind, ONLY: i4, sp, dp
  USE mo_wqm_write_fluxes_states   !, ONLY: OutputDatasetwqm

  IMPLICIT NONE


  PUBLIC :: wqm_write   ! Constant Pi in single precision
  PUBLIC :: wqm_write_output_fluxes   !State variables and Fluxes


CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         wqm_write_output_fluxes

  !     PURPOSE
  !>        \brief write fluxes to netcdf output files
  !
  !>        \details This subroutine creates a netcdf data set
  !>        for writing water quality model related state variables and fluxes for different time averages.
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !       None
  !
  !     LITERATURE
  !       None

  !     HISTORY
  !>        \author Stephan Thober
  !>        \date Aug 2015
  !         Modified, Xiaoqiang Yang, Jul 2017 - major modified to water qulaity model

  subroutine wqm_write_output_fluxes( &
       iBasin, &
       s1,e1,s11,e11, & !grid ids
       timeStep_model_outputs, & ! timestep of model outputs
       warmingDays, & ! number of warming days
       newTime, & ! julian date of next time step
       nTimeSteps, & ! number of total timesteps
       nTStepDay, & ! number of timesteps per day
       tt, & ! current model timestep
       day, & ! current day of the year
       month, & ! current month of the year
       year, & ! current year
       timestep, & ! current model time resolution
       mask1, & ! mask at level 1
       mask11, & ! mask at level 11
       ncn)
    use mo_kind, only: i4, dp
    use mo_julian, only: caldat
    use mo_wqm_global_variables, only: L1_csoilMoist, L1_cfastRunoff, L1_cslowRunoff, &
                                       L1_cbaseflow, L1_ctotal_runoff, L1_soilUptakeN, &
                                       L1_soilDenitri, L1_soilMineralN,L1_soilINfrtmanapp,   &
                                       L11_concMod, L11_aquaticDenitri,L11_aquaticAssimil
		
    implicit none
    ! input variables
    integer(i4), intent(in) :: iBasin
    integer(i4), intent(in) :: s1,e1,s11,e11
    integer(i4), intent(in) :: timeStep_model_outputs
    integer(i4), intent(in) :: warmingDays
    real(dp), intent(in) :: newTime
    integer(i4), intent(in) :: nTimeSteps
    integer(i4), intent(in) :: nTStepDay
    integer(i4), intent(in) :: tt
    integer(i4), intent(in) :: day
    integer(i4), intent(in) :: month
    integer(i4), intent(in) :: year
    integer(i4), intent(in) :: timestep
    logical, intent(in) :: mask1(:,:)
    logical, intent(in) :: mask11(:,:)
    type(OutputDatasetwqm), intent(inout)  :: ncn !for L1 , L11        

    ! local variables
    integer(i4)          :: tIndex_out
    logical              :: writeout
    integer(i4)          :: new_year ! year of next timestep (newTime)
    integer(i4)          :: new_month ! month of next timestep (newTime)
    integer(i4)          :: new_day ! day of next timestep (newTime)

    !
    ! if ( tt .EQ. 1 ) then
       ! day_counter   = day
       ! month_counter = month
       ! year_counter  = year
       ! average_counter = 0_i4
    ! end if
    
    ! update the counters
    ! if (day_counter   .NE. day  ) day_counter   = day
    ! if (month_counter .NE. month) month_counter = month
    ! if (year_counter  .NE. year)  year_counter  = year
    call caldat(int(newTime), yy=new_year, mm=new_month, dd=new_day)

    ! output only for evaluation period
    tIndex_out = (tt-warmingDays*NTSTEPDAY) ! tt if write out of warming period

    if ((tIndex_out .gt. 0_i4)) then
       ! average_counter = average_counter + 1

       ! create output dataset	   
       if ( tIndex_out .EQ. 1 ) then
       	   
#ifdef pgiFortran154
            ncn = newOutputDatasetwqm(iBasin, mask1, mask11)
#else
            ncn = OutputDatasetwqm(iBasin, mask1, mask11)
#endif
       !nc11 = OutputDatasetwqm(iBasin, mask11)
       end if
       ! update Dataset
       call ncn%updateDatasetwqm( &
            s1,e1, s11, e11 , &
            L1_csoilMoist   , &
            L1_cfastRunoff  , &
            L1_cslowRunoff  , &
            L1_cbaseflow    , &
            L1_ctotal_runoff, &
            L1_soilUptakeN  , &
            L1_soilDenitri  , &
            L1_soilMineralN , &
            L1_soilINfrtmanapp, &
            L11_concMod         , &
            L11_aquaticDenitri  , &
            L11_aquaticAssimil    &
            )
                
       ! determine write flag
       writeout = .false.
       if (timeStep_model_outputs .gt. 0) then
          if ((mod(tIndex_out, timeStep_model_outputs) .eq. 0) .or. (tt .eq. nTimeSteps)) writeout = .true.
       else
          select case(timeStep_model_outputs)
          case(0) ! only at last time step
             if (tt .eq. nTimeSteps) writeout = .true.
          case(-1) ! daily
             if (((tIndex_out .gt. 1) .and. (day .ne. new_day)) .or. (tt .eq. nTimeSteps))     writeout = .true.
          case(-2) ! monthly
             if (((tIndex_out .gt. 1) .and. (month .ne. new_month)) .or. (tt .eq. nTimeSteps)) writeout = .true.
          case(-3) ! yearly
             if (((tIndex_out .gt. 1) .and. (year .ne. new_year)) .or. (tt .eq. nTimeSteps))   writeout = .true.
          case default ! no output at all
             continue
          end select
       endif

       ! write data
       if (writeout) then
         call ncn%writeTimestepwqm(tIndex_out*timestep-1)
         !call nc11%writeTimestep(tIndex_out*timestep-1)
       end if
       ! close dataset
       if (tt .eq. nTimeSteps) then
         call ncn%close()
         !call nc11%close()
       end if
    end if
    
  end subroutine wqm_write_output_fluxes
  
  
  
  ! ------------------------------------------------------------------

  !     NAME
  !         wqm_write

  !     PURPOSE
  !>        \brief 

  !>        \details 

  !     CALLING SEQUENCE
  !         None

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
  !         

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang - Modified from the corresponding nHM routines
  !>        \date Sep 2016 

  subroutine wqm_write()

    use mo_mrm_global_variables,    only: &
          evalPer, nGaugesTotal, simPer, nBasins, NTSTEPDAY, &
          warmingDays_mrm, basin_mrm

    use mo_wqm_global_variables,    only: &
          WQM_nutrient, maxcols, basin_wqm


    implicit none
    integer(i4)     :: ii,gg
    integer(i4)     :: iDay, tt, iS, iE,intertt
    integer(i4)     :: nTimeSteps
    real(dp), dimension(:,:,:), allocatable   :: d_ConcMod



    ii = maxval( evalPer(1:nBasins)%julEnd - evalPer(1:nBasins)%julStart + 1 )
    allocate( d_ConcMod(ii, nGaugesTotal, maxcols) ) 
    d_ConcMod = 0.0_dp

    ! loop over basins
    do ii = 1, nBasins
       nTimeSteps = ( simPer(ii)%julEnd - simPer(ii)%julStart + 1 ) * NTSTEPDAY
       iDay = 0

       ! loop over timesteps
       do tt = warmingDays_mrm(ii)*NTSTEPDAY+1, nTimeSteps, NTSTEPDAY
          iS = tt
          iE = tt + NTSTEPDAY - 1
          iDay = iDay + 1
          ! over gauges
          do gg = 1, basin_mrm%nGauges(ii)
             do intertt=iS, iE
             d_ConcMod(iDay, basin_mrm%gaugeIndexList(ii,gg),:) = d_ConcMod(iDay, basin_mrm%gaugeIndexList(ii,gg),:) &
                  + WQM_nutrient(intertt, basin_mrm%gaugeIndexList(ii,gg),:)
             end do
             d_ConcMod(iDay, basin_mrm%gaugeIndexList(ii,gg),:)=d_ConcMod(iDay, basin_mrm%gaugeIndexList(ii,gg),:) &
                  / real(NTSTEPDAY,dp)
          end do
          !
       end do

    end do

    ! write in an ASCII file          ! OBS[nModeling_days X nGauges_total] , SIM[nModeling_days X nGauges_total] 
    call write_daily_obs_sim_conc( basin_wqm%GaugeConc(:,:,1), d_ConcMod(:,:,1) )

    ! free space
    deallocate(d_ConcMod)   
	
  end subroutine wqm_write

  ! ------------------------------------------------------------------

  !     NAME
  !         write_daily_obs_sim_conc

  !     PURPOSE
  !>        \brief 

  !>        \details 

  !     CALLING SEQUENCE
  !         None

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
  !         

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang - Modified from the corresponding nHM routines
  !>        \date Sep 2016 
 
  subroutine write_daily_obs_sim_conc( Concobs, Concsim )

  use mo_errormeasures,          only: kge, nse, sse, rmse, bias
  use mo_utils,                  only: ge
  use mo_mrm_global_variables,   only: &
       nBasins, basin_mrm, dirOut, evalPer, gauge

  use mo_wqm_global_variables,   only: &
       file_daily_conc  
  use mo_julian,                 only: dec2date  
  use mo_message,                only: message
  use mo_string_utils,           only: num2str
  implicit none
  real(dp), dimension(:,:), intent(in)  :: Concobs
  real(dp), dimension(:,:), intent(in)  :: Concsim
  !local
    character(256) :: fName, formHeader, formData, dummy
    integer(i4) :: bb, gg, tt, err
    integer(i4) :: igauge_start, igauge_end
    integer(i4) :: day, month, year
    real(dp) :: newTime
    ! initalize igauge_start
    igauge_start = 1

    ! basin loop
    do bb = 1, nBasins
       if( basin_mrm%nGauges(bb) .lt. 1 ) cycle

       ! estimate igauge_end
       igauge_end = igauge_start + basin_mrm%nGauges(bb) - 1

       ! check the existance of file
       fName = trim(adjustl(dirOut(bb))) // trim(adjustl(file_daily_conc))
       open(89, file=trim(fName), status='unknown', action='write', iostat=err)
       if( err .ne. 0 ) then
          call message ('  IOError while opening ',trim(fName))
          call message ('  Error-Code ', num2str(err))
          stop
       end if

       ! header
       write(formHeader, *) '( 4a8, ' , basin_mrm%nGauges(bb),'(2X, a8, i7.7, 2X, a8, i7.7) )' 
       write(89, formHeader) 'No', 'Day', 'Mon', 'Year', &
            ( 'INConcobs_', gauge%gaugeId(gg), &
            'INConcsim_', gauge%gaugeId(gg), gg=igauge_start, igauge_end )

       ! form data
       write(formData, *) '( 4I8, ' , basin_mrm%nGauges(bb),'(2X,   f15.7 , 2X,  f15.7  ) )' 

       ! write data
       newTime  = real(evalPer(bb)%julStart,dp) - 0.5_dp

       do tt = 1, (evalPer(bb)%julEnd - evalPer(bb)%julStart + 1)          
          call dec2date(newTime, yy=year, mm=month, dd=day)
          !currently, only inorganic nitrogen(nitrate) is write out to **.out file 
          write(89, formData) tt, day, month, year, ( Concobs(tt,gg), Concsim(tt,gg) , gg=igauge_start, igauge_end )
          newTime = newTime + 1.0_dp
       end do

       ! close file
       close(89)
       ! ======================================================================
       ! screen output
       ! ======================================================================
       call message()
       write(dummy,'(I3)') bb
       call message('  OUTPUT: saved daily IN concentration file for basin ', trim(adjustl(dummy)))
       call message('    to ',trim(fname))
! statistical calculation is different from discharge, regarding to those missing value(nodata_dp)
! Here needs future work!!
       do gg=igauge_start, igauge_end
          if (count(ge(Concobs(:,gg), 0.0_dp)) > 1 )  then
             call message('    RMSE of daily nitrate (gauge #',trim(adjustl(num2str(gg))),'): ', &
                  trim(adjustl(num2str(rmse(Concobs(:,gg), Concsim(:,gg), mask=(ge(Concobs(:,gg), 0.0_dp)))))) )
             call message('    NSE of daily nitrate (gauge #',trim(adjustl(num2str(gg))),'): ', &
                  trim(adjustl(num2str(nse(Concobs(:,gg), Concsim(:,gg), mask=(ge(Concobs(:,gg), 0.0_dp)))))) )
             call message('    PBIAS(%) of daily nitrate (gauge #',trim(adjustl(num2str(gg))),'): ', &
                  trim(adjustl(num2str(bias(Concobs(:,gg), Concsim(:,gg), mask=(ge(Concobs(:,gg), 0.0_dp)))))) )
          end if
       end do

       ! update igauge_start
       igauge_start = igauge_end + 1  

    end do

  end subroutine write_daily_obs_sim_conc 
  
  
END MODULE mo_wqm_write
