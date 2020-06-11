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

!> \file mo_wqm_readtimeseries.f90

!> \brief read evaluation data.

!> \details 

!> \authors Xiaoqiang Yang
!> \date Jun 2017

MODULE mo_wqm_readtimeseries


  USE mo_kind, ONLY: i4, sp, dp

  IMPLICIT NONE

  PUBLIC :: read_timeseries_conc  ! 



  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         read_timeseries_conc

  !     PURPOSE
  !>        \brief read_timeseries_conc

  !>        \details 
  !>                 

  !     INTENT(IN)
  !>        ...

  !     INTENT(INOUT)
  !         ...

  !     INTENT(OUT)
  !         ...

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        None

  !     RESTRICTIONS
  !         

  !     EXAMPLE
  !         None
  !         

  !     LITERATURE
  !         

  !     HISTORY
  !>        \authors Xiaoqiang Yang
  !>        \date Jun 2017

  subroutine read_timeseries_conc(filename, fileunit, &
               periodStart, periodEnd, optimize, opti_function, maxcols, numScol,&
               head_str, data_conc, mask, nCmeasPerday)


    use mo_julian,        only: julday
    use mo_message,       only: message
	
    implicit none

    character(len=*),                                 intent(in)  :: filename         ! name of input file
    integer(i4),                                      intent(in)  :: fileunit         ! file unit
    integer(i4), dimension(3),                        intent(in)  :: periodStart      ! format (/YYYY, MM, DD/)
    integer(i4), dimension(3),                        intent(in)  :: periodEnd        ! format (/YYYY, MM, DD/)
    logical,                                      intent(in)  :: optimize         ! optimization on or off (.TRUE. or .FALSE.)
    ! number of optimization function in mHM determining the type of data used
    integer(i4),                                      intent(in)  :: opti_function
    ! type of measured data in gauge file
    ! currently for nitrogen submodel,"maxcols" = 3 
    integer(i4),                                      intent(in)  :: maxcols          ! maximum wq data columns in file
    integer(i4),                                      intent(out) :: numScol          ! number of wq data columns   
    character(256), dimension(:), allocatable,        intent(out) :: head_str         ! heading(name) of the wq data

    real(dp),    dimension(:,:), allocatable,  intent(out) :: data_conc      ! time series output (dim1= periodStart:periodEnd)
                                                                       !dim2 = maxcols, 1 for IN, 2 for ON, 3 for TN   	 
    logical,     dimension(:,:), allocatable, optional, intent(out) :: mask           ! indicating valid data (false
    !                                                                                 ! at no data value points)
    integer(i4),                            optional, intent(out) :: nCmeasPerday      ! Number of data points read in
    !                                                                                 ! per day, e.g. hourly = 24

    ! local vaiables
    real(dp)                                             :: nodata_file      ! no data value of data
    integer(i4)                                          :: timestep_file    ! [h] time resolution of input data
    integer(i4),  dimension(3)                           :: periodStart_file ! starting date of timeseries data within file
    integer(i4),  dimension(3)                           :: periodEnd_file   ! ending date of timeseries data within file
    integer(i4),  dimension(5)                           :: time_file        ! year, month, day, hour, minute
    ! for heading of data file 
    character(256)                                       :: line             ! data heading line.


    integer(i4)                                                   :: i, j, k
    integer(i4)                                                   :: idx_st_period    ! index to put data from file to data
    integer(i4)                                                   :: idx_en_period    ! index to put data from file to data
    integer(i4)                                                   :: idx_st_file      ! index to put data from file to data
    integer(i4)                                                   :: idx_en_file      ! index to put data from file to data
    integer(i4)                                                   :: startJul_file    ! start julian day of available data
    integer(i4)                                                   :: endJul_file      ! end   julian day of available data
    integer(i4)                                                   :: startJul_period  ! start julian day of needed data
    integer(i4)                                                   :: endJul_period    ! end   julian day of needed data
    integer(i4)                                                   :: length_file      ! number of days in file
    integer(i4)                                                   :: length_period    ! number of days in period
    real(dp),    dimension(:,:), allocatable                      :: data_file        ! time series output (fileStart:fileEnd)
    !                                                                                 ! points --> 0.25 [d-1]
    character(256)                                                :: dummy            ! dummy for char read in

    open(unit=fileunit, file=filename, action='read', status='old')
      ! read header
      read(fileunit,'(a256)') dummy
      read(fileunit,*)        dummy, nodata_file
      read(fileunit,*)        dummy, timestep_file
      read(fileunit,*)        dummy, (periodStart_file(i), i = 1, 3)
      read(fileunit,*)        dummy, (periodEnd_file(i),   i = 1, 3)
      read(fileunit,'(a256)') line     
      dummy = dummy//''   ! only to avoid warning
      if ((timestep_file .lt. 1_i4) .or. (timestep_file .gt. 1440_i4)) then
         call message('***ERROR: Number of measurements per day has to be between 1 (daily) and 1440 (every minute)! ',&
                       trim(filename))
         stop
      end if


      ! checking if period is covered by data in file
      startJul_period = julday(periodStart(3),      periodStart(2),      periodStart(1))
      endJul_period   = julday(periodEnd(3)  ,      periodEnd(2)  ,      periodEnd(1)  )
      startJul_file   = julday(periodStart_file(3), periodStart_file(2), periodStart_file(1) )
      endJul_file     = julday(periodEnd_file(3),   periodEnd_file(2),   periodEnd_file(1) )

      if (((startJul_period < startJul_file) .OR. (endJul_period > endJul_file )) &
         .AND. optimize .and. (opti_function .le. 9_i4 .or. opti_function .eq. 14_i4)) then
         ! adjust this whenever a new opti function on discharge is added to mhm!
         call message('***ERROR: Simulation period is not covered by observations! ', trim(filename))
         stop
      end if

      !identify which type of data provided in the file	  
      call get_heading_column(line, numScol, head_str)

      ! allocation of arrays
      allocate( data_conc((endJul_period - startJul_period + 1_i4) * timestep_file, maxcols))
      data_conc      = nodata_file
      allocate( data_file((endJul_file   - startJul_file + 1_i4)* timestep_file, numScol))
      data_file = nodata_file

      if (present(mask)) then
         allocate( mask((endJul_period - startJul_period + 1_i4) * timestep_file, maxcols))
         mask = .true.
      end if

      ! read data from file to temporal array
      do i=1, (endJul_file - startJul_file + 1_i4) * timestep_file
         ! read date and data
         read(fileunit, *) (time_file(j),j=1,5), (data_file(i, k), k=1,numScol)

      end do
      time_file(1) = time_file(2) + 1   ! only to avoid warning

      length_file   = (endJul_file   - startJul_file   + 1 )
      length_period = (endJul_period - startJul_period + 1 )

      !       |---------------------------------|    FILE
      !                  |--------------|            PERIOD
      if (( startJul_period .ge. startJul_file ) .and. ( endJul_period .le. endJul_file )) then
         idx_st_period = 1
         idx_en_period = length_period * timestep_file
         idx_st_file   = (startJul_period - startJul_file + 1 ) * timestep_file
         idx_en_file   = (idx_st_file + length_period     - 1 ) * timestep_file
      end if

      !                  |--------------|            FILE
      !       |---------------------------------|    PERIOD
      if (( startJul_period .lt. startJul_file ) .and. ( endJul_period .gt. endJul_file )) then
         idx_st_period = (startJul_file - startJul_period + 1 ) * timestep_file
         idx_en_period = (idx_st_period + length_file     - 1 ) * timestep_file
         idx_st_file   = 1
         idx_en_file   = length_file                            * timestep_file
      end if

      !  |--------------|                            FILE
      !       |---------------------------------|    PERIOD
      if (( startJul_period .ge. startJul_file ) .and. ( endJul_period .gt. endJul_file )) then
         idx_st_period = 1
         idx_en_period = ( endJul_file     - startJul_period + 1 ) * timestep_file
         idx_st_file   = ( startJul_period - startJul_file   + 1 ) * timestep_file
         idx_en_file   = length_file                               * timestep_file
      end if

      !                          |--------------|    FILE
      !  |---------------------------------|         PERIOD
      if (( startJul_period .lt. startJul_file ) .and. ( endJul_period .le. endJul_file )) then
         idx_st_period = ( startJul_file - startJul_period + 1 ) * timestep_file
         idx_en_period = ( length_period                       ) * timestep_file
         idx_st_file   = 1
         idx_en_file   = ( endJul_period - startJul_file   + 1 ) * timestep_file
      end if

      do i=1, numScol
         select case (head_str(i))
         case ('IN', 'in', 'In', 'iN')   !inorganic nitrogen(case free)
            data_conc(idx_st_period:idx_en_period,1) = data_file(idx_st_file:idx_en_file,i)
         case ('ON', 'on', 'On', 'oN')!organic nitrogen
            data_conc(idx_st_period:idx_en_period,2) = data_file(idx_st_file:idx_en_file,i)
         case ('TN', 'tn', 'Tn', 'tN')   !total nitrogen
            data_conc(idx_st_period:idx_en_period,3) = data_file(idx_st_file:idx_en_file,i)
         end select
      end do
      !data_conc(idx_st_period:idx_en_period,:) = data_file(idx_st_file:idx_en_file,:)

      if (present(mask)) then
         where ( abs(data_conc-nodata_file) .lt. tiny(1.0_dp) )
            mask = .false.
         end where
      end if

      if (present(nCmeasPerDay)) then
         nCmeasPerDay = timestep_file
      end if

    close(fileunit)  
 
  end subroutine read_timeseries_conc
  !---------------------------------------------------------

  subroutine get_heading_column(line, nsubcol,head_str)
  
  implicit none
  character(256),              intent(in)   :: line      !string variable of the file heading 
  integer(i4),                 intent(out)  :: nsubcol   !column number for wq data (minus date/time column)
  character(256),dimension(:), allocatable,  intent(out)  :: head_str  !string variable for wq data name
  
  !local
  character(256), dimension(15)    :: str  ! tmp variables for file heading, maximum 10 columns for wq data
  integer(i4)                      :: i, j, is, ie
  logical                          :: jchar, prev_jchar

  i=1_i4  
  is = 1_i4
  str =''


  do j=2, len(line)
      jchar=.NOT.(line(j:j)==char(32).OR.line(j:j)==char(9).OR. &
            line(j:j)==char(13).OR.line(j:j)==char(10))
      prev_jchar=.NOT.(line(j-1:j-1)==char(32).OR.line(j-1:j-1)==char(9).OR.  &
            line(j-1:j-1)==char(13).OR.line(j-1:j-1)==char(10))

      
	  
      if(jchar .AND. (.NOT.(prev_jchar)))then
        is=j
      elseif(.NOT.(jchar) .AND. prev_jchar) then
        ie=j-1
        str(i) =line(is:ie)
        i=i+1
      end if 

  end do
  i=i-1
  nsubcol = i-5_i4
  allocate(head_str(nsubcol) )
  head_str(1:nsubcol) = str(6:i)
 
  end subroutine get_heading_column
  
  
  
  !-----------------------------------------------------------

  
END MODULE mo_wqm_readtimeseries
