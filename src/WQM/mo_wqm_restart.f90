!> \file mo_wqm_restart.f90

!> \brief reading and writing states, fluxes and configuration for restart of water quality model.

!> \details routines are seperated for reading and writing variables for:\n
!>          - states and fluxes, and \n
!>          - configuration.\n


!> \authors Xiaoqiang Yang, Modified from original mHM code
!> \date  Jul 2017 

MODULE mo_wqm_restart



  IMPLICIT NONE

  PUBLIC :: wqm_read_restart_states     ! 

  PUBLIC :: wqm_write_restart_files     ! 

  PRIVATE

CONTAINS
  ! ------------------------------------------------------------------
  
  !      NAME
  !         wqm_write_restart_files

  !     PURPOSE
  !>        \brief write restart files for each basin

  !>        \details write restart files for each basin. For each basin
  !>        xxx_wqm_states.nc(xxx being the three digit basin index) is written.
  !>        All variables in mo_wqm_read.f90: wqm_variable_default_init() should be
  !>        included here.
  !>        If a variable is added here, it should also be added
  !>        in the read restart routines below.

  !     INTENT(IN)
  !>        \param[in] "character(256), dimension(:) :: OutPath"     Output Path for each basin

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

  !     RESTRICTIONS 
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         Modified from the original mHM and mRM

  !     HISTORY
  !>        \author   Xiaoqiang Yang
  !>        \date     Jul 2017
  !         Modified  

  ! ------------------------------------------------------------------ 
  subroutine wqm_write_restart_files(OutPath )

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_init_states,      only: get_basin_info
    use mo_mrm_tools,        only: get_basin_info_mrm
    use mo_string_utils,     only: num2str
    use mo_netcdf,           only: NcDataset, NcDimension, NcVariable
    use mo_mhm_constants,    only: nodata_dp
    use mo_wqm_global_variables, only: &
         L1_humusN, L1_fastN, L1_dissolvedIN,              & ! INOUT NS Four different Nitrate pools in each soil layer
         L1_dissolvedON, L1_csoilMoist,    & ! L1_initcsm,INOUT NS conc. in each soil layer, dim2=soillayers, dim3=substances  
         L1_csealSTW, L1_cunsatSTW, L1_csatSTW,            & ! INOUT NS conc. in each water storage, dim2=substances
         L1_soiltemp, L1_baseflow_avg,  & ! INOUT NS soil temperature, variables for calculating baseflow conc.
         L1_cinfilSoil, L1_cpreEffect,  & ! INOUT NX dim2=soillayers, dim3=substances/dim2=substances  
         L1_crain, L1_cpercolate, L1_crunoffSeal,          & ! INOUT NX dim2=substances
         L1_cfastRunoff, L1_cslowRunoff, L1_cbaseflow,     & ! INOUT NX dim2=substances 
         L1_ctotal_runoff,                                                   & ! INOUT NX dim2=substances 
         L1_rdegradN, L1_rmineralN, L1_rdissolN,                   & ! INOUT NE1 nitrate submodel parameters
         L1_rdeniSoil,                       &     ! INOUT NE1 nitrate submodel parameters  
         prevstep_sealedStorage,             & ! sealed storage in last step
         prevstep_unsatStorage,              & ! unsaturated storage in last step
         prevstep_satStorage,                & ! saturated storage in last step
         prevstep_soilMoisture,                 & ! soil moisture in last step		 
         prevstep_percol,              & !
         prevstep_baseflow,            & !
         L1_reachtemp,   & !
         L11_riverbox, L11_criverbox, L11_yravg_q,       & !inout 
         L11_concOUT, L11_interload, L11_concTIN,    & !inout dim2=substances                  
         L11_concMod,                            & !inout 
         L11_rdeniAqtc,L11_rpprodN,L11_rivert_avg10, &
         L11_rivert_avg20, L11_rzcoeff, L11_flight         
    implicit none

    character(256)                           :: Fname
    character(256),dimension(:), intent(in)  :: OutPath ! list of Output paths per Basin
    integer(i4)                              :: iBasin
    integer(i4)                              :: ii, jj
    integer(i4)                              :: s0       ! start index at level 0
    integer(i4)                              :: e0       ! end index at level 0
    integer(i4)                              :: ncols0   ! number of colums at level 0
    integer(i4)                              :: nrows0   ! number of rows at level 0
    logical, dimension(:,:), allocatable     :: mask0    ! mask at level 0
    integer(i4)                              :: s1       ! start index at level 1
    integer(i4)                              :: e1       ! end index at level 1
    integer(i4)                              :: ncols1   ! number of colums at level 1
    integer(i4)                              :: nrows1   ! number of rows at level 1
    logical, dimension(:,:), allocatable     :: mask1    ! mask at level 1
    real(dp), dimension(:,:,:), allocatable  :: dummy_d3 ! dummy variable
    real(dp), dimension(:,:,:,:), allocatable:: dummy_d4
    integer(i4)                              :: s11       ! start index at level 11
    integer(i4)                              :: e11       ! end index at level 11
    integer(i4)                              :: ncols11   ! number of colums at level 11
    integer(i4)                              :: nrows11   ! number of rows at level 11
    logical, dimension(:,:), allocatable     :: mask11    ! mask at level 11

    type(NcDataset)                          :: nc
    type(NcDimension)                        :: rows1, cols1, rows11,cols11, soil1,nsubst
    type(NcVariable)                         :: var
    
    basin_loop: do iBasin = 1, size(OutPath)

       ! get Level0 information about the basin
       call get_basin_info( iBasin, 0, nrows0, ncols0, iStart=s0, iEnd=e0, mask=mask0 )

       ! get Level1 information about the basin
       call get_basin_info( iBasin, 1, nrows1, ncols1, iStart=s1, iEnd=e1, mask=mask1 )
       ! level 11 info
       call get_basin_info_mrm( iBasin, 11, nrows11, ncols11, iStart=s11, iEnd=e11, mask=mask11 )

       ! write restart file for iBasin
       Fname = trim(OutPath(iBasin)) // "WQM_restart_" // trim(num2str(iBasin, "(i3.3)")) // ".nc"
       ! print a message
       call message("    Writing Restart-file: ", trim(adjustl(Fname))," ...")

       nc     = NcDataset(fname, "w")
       rows1  = nc%setDimension("nrows1 ", nrows1)
       cols1  = nc%setDimension("ncols1 ", ncols1)
       rows11  = nc%setDimension("nrows11", nrows11)
       cols11  = nc%setDimension("ncols11", ncols11)
       soil1  = nc%setDimension("soilhorizons", size( L1_csoilMoist, 2))
       nsubst  = nc%setDimension("substances", size( L1_csoilMoist, 3))

       !for three dimensions(row, col, nsoillyer(soil1))
       allocate( dummy_d3( nrows1, ncols1, size( L1_csoilMoist, 2) ) )
       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_humusN(s1:e1,ii), mask1, nodata_dp )
       end do
       var = nc%setVariable("L1_humusN","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","soil humusN storage at level 1")
	   
       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_fastN(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_fastN","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","soil fastN storage at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_dissolvedIN(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_dissolvedIN","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","soil dissolved IN storage at level 1")
      
       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_dissolvedON(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_dissolvedON","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","soil dissolved ON storage at level 1")
      
       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( prevstep_soilMoisture(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("prevstep_soilMoisture","f64",(/rows1,cols1,soil1/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","soil moisture of previous step at level 1")

       !for four dimensions(row, col, nsoillyer(soil1),nsubst)
       allocate( dummy_d4( nrows1,ncols1, size(L1_csoilMoist,2), size(L1_csoilMoist,3)))
       do ii = 1, size(dummy_d4,4)
           do jj =1, size(dummy_d4,3)
               dummy_d4(:,:,jj,ii) = unpack( L1_csoilMoist(s1:e1,jj,ii), mask1, nodata_dp )
           end do
       end do
       var = nc%setVariable("L1_csoilMoist","f64",(/rows1,cols1, soil1, nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d4)
       call var%setAttribute("long_name","soil moisture concentration at level 1")
	  
       do ii = 1, size(dummy_d4,4)
           do jj =1, size(dummy_d4,3)
               dummy_d4(:,:,jj,ii) = unpack( L1_cinfilSoil(s1:e1,jj,ii), mask1, nodata_dp )
           end do
       end do
       var = nc%setVariable("L1_cinfilSoil","f64",(/rows1,cols1, soil1, nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d4)
       call var%setAttribute("long_name","concentration of infiltrated water at level 1")

       !for three dimensions (row, col, nsubst)
       deallocate(dummy_d3)
       allocate(dummy_d3( nrows1, ncols1, size( L1_csoilMoist, 3) ))

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_csealSTW(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_csealSTW","f64",(/rows1,cols1,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","concentration in sealed storage at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_cunsatSTW(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_cunsatSTW","f64",(/rows1,cols1,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","concentration in unsaturated storage at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_csatSTW(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_csatSTW","f64",(/rows1,cols1,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","concentration in saturated storage at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_cpreEffect(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_cpreEffect","f64",(/rows1,cols1,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","effective precipitation concentration at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_crain(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_crain","f64",(/rows1,cols1,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","precipitation concentration wet deposition at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_cpercolate(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_cpercolate","f64",(/rows1,cols1,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","percolated water concentration at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_crunoffSeal(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_crunoffSeal","f64",(/rows1,cols1,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","sealed runoff concentration at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_cfastRunoff(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_cfastRunoff","f64",(/rows1,cols1,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","fast interflow runoff concentration at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_cslowRunoff(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_cslowRunoff","f64",(/rows1,cols1,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","slow interflow runoff concentration at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_cbaseflow(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_cbaseflow","f64",(/rows1,cols1,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","base flow concentration at level 1")

       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L1_ctotal_runoff(s1:e1,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L1_ctotal_runoff","f64",(/rows1,cols1,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","total runoff concentration at level 1")
	   
       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L11_criverbox(s11:e11,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L11_criverbox","f64",(/rows11,cols11,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","concentration in stored volume at level 11")
	   
       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L11_concOUT(s11:e11,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L11_concOUT","f64",(/rows11,cols11,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","concentration in local terrestrial output at level 11")
	   
       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L11_interload(s11:e11,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L11_interload","f64",(/rows11,cols11,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","load in stored volume at level 11")
       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L11_concTIN(s11:e11,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L11_concTIN","f64",(/rows11,cols11,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","concentration of total reach input at level 11")
       do ii = 1, size( dummy_d3, 3 )
          dummy_d3(:,:,ii) = unpack( L11_concMod(s11:e11,ii), mask1, nodata_dp )
       end do      
       var = nc%setVariable("L11_concMod","f64",(/rows11,cols11,nsubst/))
       call var%setFillValue(nodata_dp)
       call var%setData(dummy_d3)
       call var%setAttribute("long_name","concentration of total reach output at level 11")  

       !for two dimensions (row, col)
       var = nc%setVariable("L1_soiltemp","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_soiltemp(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","soil temperature at level 1")  

       var = nc%setVariable("L1_baseflow_avg","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_baseflow_avg(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","30 days averaged baseflow at level 1")

       var = nc%setVariable("L1_rdegradN","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_rdegradN(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","degradation rate at level 1")

       var = nc%setVariable("L1_rmineralN","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_rmineralN(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","mineralization rate at level 1")

       var = nc%setVariable("L1_rdissolN","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_rdissolN(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","dissolution rate at level 1")

       var = nc%setVariable("L1_rdeniSoil","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_rdeniSoil(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","soil denitrification rate at level 1")
       var = nc%setVariable("prevstep_sealedStorage","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(prevstep_sealedStorage(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","sealed storage of previous step at level 1")
       var = nc%setVariable("prevstep_unsatStorage","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(prevstep_unsatStorage(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","unsaturated storage of previous step at level 1")
       var = nc%setVariable("prevstep_satStorage","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(prevstep_satStorage(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","saturated storage of previous step at level 1")
       var = nc%setVariable("prevstep_percol","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(prevstep_percol(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","percolated volume of previous step at level 1")
       var = nc%setVariable("prevstep_baseflow","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(prevstep_baseflow(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","baseflow of previous step at level 1")
       var = nc%setVariable("L1_reachtemp","f64",(/rows1,cols1/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L1_reachtemp(s1:e1), mask1, nodata_dp))
       call var%setAttribute("long_name","surface water temperature at level 1")

       var = nc%setVariable("L11_riverbox","f64",(/rows11,cols11/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L11_riverbox(s11:e11), mask1, nodata_dp))
       call var%setAttribute("long_name","water volume stored in reach at level 11")
       var = nc%setVariable("L11_yravg_q","f64",(/rows11,cols11/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L11_yravg_q(s11:e11), mask1, nodata_dp))
       call var%setAttribute("long_name","moving averaged discharge at level 11")
       var = nc%setVariable("L11_rdeniAqtc","f64",(/rows11,cols11/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L11_rdeniAqtc(s11:e11), mask1, nodata_dp))
       call var%setAttribute("long_name","aquatic denitrification rate at level 11")
       var = nc%setVariable("L11_rpprodN","f64",(/rows11,cols11/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L11_rpprodN(s11:e11), mask1, nodata_dp))
       call var%setAttribute("long_name","aquatic primary production rate at level 11")
       var = nc%setVariable("L11_rzcoeff","f64",(/rows11,cols11/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L11_rzcoeff(s11:e11), mask1, nodata_dp))
       call var%setAttribute("long_name","shading coefficient of riparian zone at level 11")
       var = nc%setVariable("L11_flight","f64",(/rows11,cols11/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L11_flight(s11:e11), mask1, nodata_dp))
       call var%setAttribute("long_name","overall near surface light coefficient at level 11")
       var = nc%setVariable("L11_rivert_avg10","f64",(/rows11,cols11/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L11_rivert_avg10(s11:e11), mask1, nodata_dp))
       call var%setAttribute("long_name","10day moving averaged water temperature at level 11")
       var = nc%setVariable("L11_rivert_avg20","f64",(/rows11,cols11/))
       call var%setFillValue(nodata_dp)
       call var%setData(unpack(L11_rivert_avg20(s11:e11), mask1, nodata_dp))
       call var%setAttribute("long_name","20day moving averaged water temperature at level 11")

       call nc%close()
       
    end do basin_loop
    
  end subroutine wqm_write_restart_files


  ! ------------------------------------------------------------------

  !      NAME
  !         wqm_read_restart_states

  !     PURPOSE
  !>        \brief reads fluxes and state variables from file

  !>        \details read fluxes and state variables from given 
  !>        restart directory and initialises all state variables
  !>        that are initialized in the subroutine wqm_initialise,
  !>        contained in module mo_water_quality.

  !     INTENT(IN)
  !>        \param[in] "integer(i4)    :: iBasin"        number of basin


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

  !     RESTRICTIONS 
  !>        \note Restart Files must have the format, as if
  !>        it would have been written by subroutine write_restart_files

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang  -- Modified from the original mHM and mRM
  !>        \date Jun 2017

  subroutine wqm_read_restart_states( iBasin, InPath )

    use mo_kind,             only: i4, dp
    ! use mo_message,          only: message
    use mo_netcdf,           only: NcDataset, NcVariable
    use mo_string_utils,     only: num2str
    use mo_init_states,      only: get_basin_info
    use mo_mrm_tools,        only: get_basin_info_mrm
    !use mo_mhm_constants,    only: YearMonths_i4
    use mo_global_variables, only:  &
         nSoilHorizons_mHM
    use mo_wqm_global_variables, only: &
         L1_humusN, L1_fastN, L1_dissolvedIN,              & ! INOUT NS Four different Nitrate pools in each soil layer
         L1_dissolvedON, L1_csoilMoist,    & ! L1_initcsm,INOUT NS conc. in each soil layer, dim2=soillayers, dim3=substances  
         L1_csealSTW, L1_cunsatSTW, L1_csatSTW,            & ! INOUT NS conc. in each water storage, dim2=substances
         L1_soiltemp, L1_baseflow_avg,  & ! INOUT NS soil temperature, variables for calculating baseflow conc.
         L1_cinfilSoil, L1_cpreEffect,  & ! INOUT NX dim2=soillayers, dim3=substances/dim2=substances  
         L1_crain, L1_cpercolate, L1_crunoffSeal,          & ! INOUT NX dim2=substances
         L1_cfastRunoff, L1_cslowRunoff, L1_cbaseflow,     & ! INOUT NX dim2=substances 
         L1_ctotal_runoff,                                                   & ! INOUT NX dim2=substances 
         L1_rdegradN, L1_rmineralN, L1_rdissolN,                   & ! INOUT NE1 nitrate submodel parameters
         L1_rdeniSoil,                       &     ! INOUT NE1 nitrate submodel parameters  
         prevstep_sealedStorage,             & ! sealed storage in last step
         prevstep_unsatStorage,              & ! unsaturated storage in last step
         prevstep_satStorage,                & ! saturated storage in last step
         prevstep_soilMoisture,                 & ! soil moisture in last step		 
         prevstep_percol,              & !
         prevstep_baseflow,            & !
         L1_reachtemp,   & !
         L11_riverbox, L11_criverbox, L11_yravg_q,       & !inout 
         L11_concOUT, L11_interload, L11_concTIN,    & !inout dim2=substances                  
         L11_concMod,                            & !inout 
         L11_rdeniAqtc,L11_rpprodN,L11_rivert_avg10, &
         L11_rivert_avg20, L11_rzcoeff, L11_flight, &
         nsubstances

    implicit none

    integer(i4),    intent(in) :: iBasin
    character(256), intent(in) :: InPath

    character(256)                                    :: Fname
    integer(i4)                                       :: ii,jj       ! loop index
    integer(i4)                                       :: s1,s11       ! start index at level 1/11
    integer(i4)                                       :: e1,e11       ! end index at level 1/11
    integer(i4)                                       :: ncols1,ncols11   ! number of colums at level 1/11
    integer(i4)                                       :: nrows1,nrows11   ! number of rows at level 1/11
    integer(i4)                                       :: ncells1  ! number of cells at level 1
    logical, dimension(:,:), allocatable              :: mask1,mask11    ! mask at level 1/11

    real(dp), dimension(:,:),   allocatable           :: dummyD2  ! dummy, 2 dimension
    real(dp), dimension(:,:,:), allocatable           :: dummyD3  ! dummy, 3 dimension
    real(dp), dimension(:,:,:,:), allocatable         :: dummyD4

    type(NcDataset) :: nc
    type(NcVariable) :: var
    
    Fname = trim(InPath) // 'WQM_restart_' // trim(num2str(iBasin, '(i3.3)')) // '.nc'
    ! call message('    Reading states from ', trim(adjustl(Fname)),' ...')

    
    ! get basin information at level 1
    call get_basin_info( iBasin, 1, nrows1, ncols1, ncells=ncells1, &
         iStart=s1, iEnd=e1, mask=mask1 )
    ! level 11 info
    call get_basin_info_mrm( iBasin, 11, nrows11, ncols11, iStart=s11, iEnd=e11, mask=mask11 )

	
    nc = NcDataset(fname,"r")

    ! three dimensions(row,col, nsoillayer)
    var = nc%getVariable("L1_humusN")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_humusN(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_fastN")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_fastN(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_dissolvedIN")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_dissolvedIN(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_dissolvedON")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       L1_dissolvedON(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("prevstep_soilMoisture")
    call var%getData(dummyD3)
    do ii = 1, nSoilHorizons_mHM
       prevstep_soilMoisture(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
	
    !four dimensions
    var = nc%getVariable("L1_csoilMoist")
    call var%getData(dummyD4)
    do ii = 1, nsubstances
       do jj = 1, nSoilHorizons_mHM
       L1_csoilMoist(s1:e1, jj,ii) = pack( dummyD4( :,:,jj, ii), mask1)
       end do
    end do
    var = nc%getVariable("L1_cinfilSoil")
    call var%getData(dummyD4)
    do ii = 1, nsubstances
       do jj = 1, nSoilHorizons_mHM
       L1_cinfilSoil(s1:e1, jj,ii) = pack( dummyD4( :,:,jj, ii), mask1)
       end do
    end do

    !three dimensions (row, col, nsubstances)
    ! level 1
    deallocate(dummyD3)
    var = nc%getVariable("L1_csealSTW")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L1_csealSTW(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_cunsatSTW")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L1_cunsatSTW(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_csatSTW")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L1_csatSTW(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_cpreEffect")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L1_cpreEffect(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_crain")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L1_crain(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_cpercolate")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L1_cpercolate(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_crunoffSeal")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L1_crunoffSeal(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_cfastRunoff")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L1_cfastRunoff(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_cslowRunoff")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L1_cslowRunoff(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_cbaseflow")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L1_cbaseflow(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    var = nc%getVariable("L1_ctotal_runoff")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L1_ctotal_runoff(s1:e1, ii) = pack( dummyD3( :,:,ii), mask1)
    end do
    !level 11
    deallocate(dummyD3)
    var = nc%getVariable("L11_criverbox")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L11_criverbox(s11:e11, ii) = pack( dummyD3( :,:,ii), mask11)
    end do
    var = nc%getVariable("L11_concOUT")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L11_concOUT(s11:e11, ii) = pack( dummyD3( :,:,ii), mask11)
    end do
    var = nc%getVariable("L11_interload")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L11_interload(s11:e11, ii) = pack( dummyD3( :,:,ii), mask11)
    end do
    var = nc%getVariable("L11_concTIN")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L11_concTIN(s11:e11, ii) = pack( dummyD3( :,:,ii), mask11)
    end do
    var = nc%getVariable("L11_concMod")
    call var%getData(dummyD3)
    do ii = 1, nsubstances
       L11_concMod(s11:e11, ii) = pack( dummyD3( :,:,ii), mask11)
    end do

    !two dimensions
    !level 1
    var = nc%getVariable("L1_soiltemp")
    call var%getData(dummyD2)
    L1_soiltemp(s1:e1) = pack( dummyD2, mask1 ) 
    var = nc%getVariable("L1_baseflow_avg")
    call var%getData(dummyD2)
    L1_baseflow_avg(s1:e1) = pack( dummyD2, mask1 ) 
    var = nc%getVariable("L1_rdegradN")
    call var%getData(dummyD2)
    L1_rdegradN(s1:e1) = pack( dummyD2, mask1 ) 
    var = nc%getVariable("L1_rmineralN")
    call var%getData(dummyD2)
    L1_rmineralN(s1:e1) = pack( dummyD2, mask1 ) 
    var = nc%getVariable("L1_rdissolN")
    call var%getData(dummyD2)
    L1_rdissolN(s1:e1) = pack( dummyD2, mask1 ) 
    var = nc%getVariable("L1_rdeniSoil")
    call var%getData(dummyD2)
    L1_rdeniSoil(s1:e1) = pack( dummyD2, mask1 ) 
    var = nc%getVariable("prevstep_sealedStorage")
    call var%getData(dummyD2)
    prevstep_sealedStorage(s1:e1) = pack( dummyD2, mask1 ) 
    var = nc%getVariable("prevstep_unsatStorage")
    call var%getData(dummyD2)
    prevstep_unsatStorage(s1:e1) = pack( dummyD2, mask1 ) 
    var = nc%getVariable("prevstep_satStorage")
    call var%getData(dummyD2)
    prevstep_satStorage(s1:e1) = pack( dummyD2, mask1 ) 
    var = nc%getVariable("prevstep_percol")
    call var%getData(dummyD2)
    prevstep_percol(s1:e1) = pack( dummyD2, mask1 ) 
    var = nc%getVariable("prevstep_baseflow")
    call var%getData(dummyD2)
    prevstep_baseflow(s1:e1) = pack( dummyD2, mask1 ) 
    var = nc%getVariable("L1_reachtemp")
    call var%getData(dummyD2)
    L1_reachtemp(s1:e1) = pack( dummyD2, mask1 ) 
    !level 11
    deallocate(dummyD2)
    var = nc%getVariable("L11_riverbox")
    call var%getData(dummyD2)
    L11_riverbox(s11:e11) = pack( dummyD2, mask11 ) 
    var = nc%getVariable("L11_rdeniAqtc")
    call var%getData(dummyD2)
    L11_rdeniAqtc(s11:e11) = pack( dummyD2, mask11 ) 
    var = nc%getVariable("L11_rpprodN")
    call var%getData(dummyD2)
    L11_rpprodN(s11:e11) = pack( dummyD2, mask11 ) 
    var = nc%getVariable("L11_rzcoeff")
    call var%getData(dummyD2)
    L11_rzcoeff(s11:e11) = pack( dummyD2, mask11 ) 
    var = nc%getVariable("L11_flight")
    call var%getData(dummyD2)
    L11_flight(s11:e11) = pack( dummyD2, mask11 )    
    var = nc%getVariable("L11_yravg_q")
    call var%getData(dummyD2)
    L11_yravg_q(s11:e11) = pack( dummyD2, mask11 ) 
    var = nc%getVariable("L11_rivert_avg10")
    call var%getData(dummyD2)
    L11_rivert_avg10(s11:e11) = pack( dummyD2, mask11 ) 
    var = nc%getVariable("L11_rivert_avg20")
    call var%getData(dummyD2)
    L11_rivert_avg20(s11:e11) = pack( dummyD2, mask11 ) 

    call nc%close()

  end subroutine wqm_read_restart_states

END MODULE mo_wqm_restart
