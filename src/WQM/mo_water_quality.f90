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

!> \file mo_water_quality.f90

!> \brief All main processes of water quality model.

!> \details This module calls all processes of water quality model if the "water qulity process" is turned on
!>          in the configure file. The configuration is specified in the namelist mhm.nml.
!>          Currently, only Nitrogen sub-model is implemented.
!>          Notes: Routing processes for water quality modeling are not seperated, thus the routing 
!>          process (8) should be turned on beforehand.

!> \authors Xiaoqiang Yang
!> \date Apr 2018

MODULE mo_water_quality


  USE mo_kind, ONLY: i4, sp, dp

  IMPLICIT NONE


  PUBLIC :: wqm  ! water quality model 

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         wqm

  !     PURPOSE
  !>        \brief Water quality modelling 

  !>        \details Calculation of all routines involved in water quality modeling
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
  !         Routing process must be turned on beforehand

  !     EXAMPLE
  !         None
  !         

  !     LITERATURE
  !         HYPE model and INCA model

  !     HISTORY
  !>        \authors Xiaoqiang Yang
  !>        \date Jun 2017
  !>        \modified          Xiaoqiang Yang  Nov. 2018 - implemented new regionalization approach of in-stream uptake

  subroutine wqm(&
      !configuration
      processMatrix,     & ! process setting
      global_parameter,  & ! global parameters for nitrate, read in from namelist: mhm_parameter.nml: nitrateparameter
      perform_mpr,       & !
      read_restart,      & !
      nCells1,           & ! number of cells in a given basin at level L1
      nHorizons_mHM,     & ! number of soil layers defined in config file
      HorizonDepth,      & ! depth of each soil layer of mhm
      tt,                & ! simulated time step
      time,              & ! current decimal Julian day (shifted by half a day)
      timeStep,          & ! simulation time step in hour[h]
      iBasin,            & ! id for the current basin
      mask0,             & ! mask at level L0
      cellId0,           & !
      nTCells0_inL1,     & ! total number of valid L0 cells in a given L1 cell
      L0upBound_inL1,    & ! upper row of L0 block within L1 cell
      L0downBound_inL1,  & ! lower row of L0 block within L1 cell
      L0leftBound_inL1,  & ! left column of L0 block within L1 cell
      L0rightBound_inL1, & ! right column of L0 block within L1 cell
      !wet deposition
      rain,              & ! precipitation in form of rainfall 
      prec_effect,       & ! effective precipitation
      temp,              & ! air temperature(***should be adjusted to day-night corrected value***)
      !soil moisture
      fsealed,           & ! fraction of sealed area
      sealedStorage,     & ! sealed water storage
      runoff_sealed,     & ! direct runoff
      infiltration,      & ! infiltration for each soil layer
      soilMoisture,      & ! soil moisture in each layer
      wilting_point,     & ! soil water content: wilting point 
      soilmoist_sat,     & ! soil water content: saturation soil moisture
      aet_soil,          & ! actual evaptranspiration in each soil layer
      aet_sealed,        & ! actual evaptranspiration in sealed area
      !runoff generation
      satStorage,        & ! saturated storage
      unsatStorage,      & ! unsaturated storage
      snowpack,          & ! snowpack depth
      slowRunoff,        & ! slow runoff(interflow)
      fastRunoff,        & ! fast runoff(surfaceflow)  
      baseflow,          & ! baseflow
      percol,            & ! percolation
      total_runoff,      & ! total runoff
      !previous step state value
      prevstep_sealedStorage,  & ! 
      prevstep_unsatStorage,   & !
      prevstep_satStorage,     & !
      prevstep_soilMoisture,   & !
      prevstep_percol,         & !
      prevstep_baseflow,       & !
      !rotation and land use cover
      L0_cover_rotation, & ! rotaion cover from rotation map
      L0_LCover_LAI,     & ! land use cover from LAI map
      L0_soilId,         & ! soil type id
      L0_geoUnit,        & ! geounit id
      fgeounit,          & ! aera fraction of  geounit	  
      frotation,         & ! aera fraction of rotation type at L1 level
      fLAI,              & ! aera fraction of LAI type at L1 level
      fLAI_11,           & ! area fraction of LAI type at L11 level
      fsoil,             & ! area fraction of soil type at L1 level
      !nitrate in soil phase
      nsubstances,       & ! number of substances 
      nCroptation,       & ! number of crop rotation class in current basin (ii) 
      !rotation,          & ! rotation information in current basin(ii)
      humusN,            & ! humus Nitrate pool in each soil layer(organic form)
      fastN,             & ! fast Nitrate pool in each soil layer
      dissIN,            & ! dissolved inorganic pool in each soil layer
      dissON,            & ! dissolved organic pool in each soil layer
      csoilMoist,        & ! conc. in soil moisture
      csealSTW,          & ! conc. in sealed water storage
      cunsatSTW,         & ! conc. in unsaturated water storage
      csatSTW,           & ! conc. in saturated water storage
      soiltemp,          & ! soil temperature (calculated from air temperature)
      basef_avg,         & ! mean baseflow for calculating baseflow conc.(eq. introduced from INCA model) 
      cinfiltration,     & ! conc. in infiltrated soil water (between each soillayer)
      cprec_effect,      & ! conc. in effective precipitation
      crain,             & ! conc. in rainfall (wet atmospheric deposition)
      cpercol,           & ! conc. in percolated water from unsat storage to satstorage
      crunoff_sealed,    & ! conc. in direct runoff
      cfastrunoff,       & ! conc. in fast runoff (HERE mHM named fast_interflow)
      cslowrunoff,       & ! conc. in slow runoff (HERE mHM named slow_interflow)
      cbaseflow,         & ! conc. in baseflow
      ctotal_runoff,     & ! conc. in total runoff (mixture from all runoff component)
      soil_uptakeN,      & ! N uptake amount in soil phase [mg/m^2/timestep]
      soil_denitri,      & ! N denitrified amount in soil phase
      soil_mineralN,     & ! N mineralised amount in soil phase
      infrtmanapp,       & ! IN applied from fertilizer and manure
      degradN_rate,      & ! nitrate parameter degradation rate from humusN pool to fastN pool
      mineraN_rate,      & ! nitrate parameter mineralisation rate from fastN pool to dissolved inorganic pool
      dissolN_rate,      & ! nitrate parameter dissolution rate from fastN pool to dissolved organic pool
      sdenitr_rate,      & ! nitrate parameter IN denitrification rate in soil water
      !routing
      !      
      nNodes,            & ! number of nodes of river connection
      areaCell1,         & ! [km2]cell area at level L1
      areaCell11,        & ! [km2]cell area at level L11
      L11id_on_L1,       & ! L11id on L1 map
      L1id_on_L11,       & ! L1 di on L11 map
      map_flag,          & ! L11 cell size is bigger than L1 or not
      nLink_from,        & ! Link from node 
      nLink_to,          & ! Link to node
      nLink_length,      & ! Link length
      !additional inflow gauging information
      nInflowGauges,     & ! number of inflow gauges in basin ii
      InflowHeadwater,   & ! flag for headwater cell of inflow gauge
      InflowNodeList,    & ! guage node list at L11
      InflowIndexList,   & ! index list of inflow gauges
      Qinflow,           & ! discharge at inflow station at current day??
      cqInflow,          & ! concentration at inflow station at current day??
      globalradi_nor,    & ! normalised global radiation(daily)
      reachtemp1,        & ! regional reach temperature in each L1 cell (only for calcualting of rivertemp)
      nLink_rivertemp,   & ! river reach temperature of every river (link)
      nNode_qOUT,        & ! outflow from each L11 cell[m3/s]
      nNode_qTR,         & ! outflow from each reach (link)
      nPerm,             & ! L11 routing order
      nLink_riverbox,    & ! amount of water stored in each river reach
      nLink_criverbox,   & ! conc. of each riverbox (above)
      nLink_yravg_q,     & ! yearly average discharge of each link
      nNode_concOUT,     & ! conc. of runoff from L11 cell, corresponding to "L11_qOUT" 
      nNode_interload,   & ! total load of inflow from upstream reach(es)
      nNode_concTIN,     & ! conc. of total inflow in each reach node, corresponding to "L11_qTIN"
      nNode_concMod,     & ! modelled conc. in each node(ouput result), corresponding to "L11_qMod"
      aquatic_denitri,   & ! N denitrified in stream water
      aquatic_assimil,   & ! N uptake (assimilatory) amount in stream water
      adenitr_rate,      & ! parameter: denitrification rate in aquatic system (in-stream) 
      atruptk_rate,     & ! parameter: autotrophic uptake rate in-stream mgNm-2d-1 
      priprod_rate,      & ! parameter: assimilatory(primary production) rate in aquatic system	  
      nLink_rzcoeff,     & ! riparian zone shading coefficient 
      nLink_flight)      ! moving averaged shading coefficient
		  
    use mo_global_variables,     only: nLAIclass, nSoilTypes,nGeoUnits, &  ! number of LAI types (land use types)/soil type/Geounits
                                       GeoUnitList
    use mo_wqm_global_variables, only: rotation,year_start, day_preyr   ! simulation starting year		 
    use mo_julian,               only: ndays, dec2date
    use mo_upscaling_operators,  only: L0_fractionalCover_in_Lx
    use mo_wqm_mpr,              only: wqm_mpr
    !use mo_wqm_restart,          only: wqm_read_restart_states
    use mo_wqm_soiltexture,    only: L1_soiltexture

    implicit none
    !intent
    !IN
    integer(i4), dimension(:,:), intent(in) :: processMatrix
    real(dp), dimension(:),      intent(in) :: global_parameter
    logical,                     intent(in) :: perform_mpr
    logical,                     intent(in) :: read_restart
    integer(i4),                 intent(in) :: nCells1
    integer(i4),                 intent(in) :: nHorizons_mHM
    real(dp), dimension(:),      intent(in) :: HorizonDepth
    integer(i4),                 intent(in) :: tt
    real(dp),                    intent(in) :: time
    integer(i4),                 intent(in) :: timeStep

    ! Physiographic L1
    integer(i4),                   intent(in) :: iBasin
    logical, dimension(:,:),       intent(in) :: mask0
    integer(i4), dimension(:),     intent(in) :: cellId0
    integer(i4), dimension(:),     intent(in) :: nTCells0_inL1
    integer(i4), dimension(:),     intent(in) :: L0upBound_inL1
    integer(i4), dimension(:),     intent(in) :: L0downBound_inL1
    integer(i4), dimension(:),     intent(in) :: L0leftBound_inL1
    integer(i4), dimension(:),     intent(in) :: L0rightBound_inL1
	!hydrological variables
    real(dp), dimension(:),      intent(in) :: rain
    real(dp), dimension(:),      intent(in) :: prec_effect
    real(dp), dimension(:),      intent(in) :: temp
    real(dp), dimension(:),      intent(in) :: fsealed
    real(dp), dimension(:),      intent(in) :: sealedStorage
    real(dp), dimension(:),      intent(in) :: runoff_sealed
    real(dp), dimension(:,:),    intent(in) :: infiltration
    real(dp), dimension(:,:),    intent(in) :: soilMoisture
    real(dp), dimension(:,:),    intent(in) :: wilting_point
    real(dp), dimension(:,:),    intent(in) :: soilmoist_sat
    real(dp), dimension(:,:),    intent(in) :: aet_soil
    real(dp), dimension(:),      intent(in) :: aet_sealed
    real(dp), dimension(:),      intent(in) :: satStorage
    real(dp), dimension(:),      intent(in) :: unsatStorage
    real(dp), dimension(:),      intent(in) :: snowpack
    real(dp), dimension(:),      intent(in) :: slowRunoff
    real(dp), dimension(:),      intent(in) :: fastRunoff
    real(dp), dimension(:),      intent(in) :: baseflow
    real(dp), dimension(:),      intent(in) :: percol
    real(dp), dimension(:),      intent(in) :: total_runoff
    !values of pervious time step
    real(dp), dimension(:),      intent(inout) :: prevstep_sealedStorage    ! sealed storage in last step    
    real(dp), dimension(:),      intent(inout) :: prevstep_unsatStorage    ! usat storage in last step    
    real(dp), dimension(:),      intent(inout) :: prevstep_satStorage    ! sat storage in last step    
    real(dp), dimension(:,:),    intent(inout) :: prevstep_soilMoisture    ! soil moisture in last step
    real(dp), dimension(:),      intent(inout) :: prevstep_percol
    real(dp), dimension(:),      intent(inout) :: prevstep_baseflow
    !land use(LAI) and rotation type 
    integer(i4), dimension(:),   intent(in)  :: L0_cover_rotation
    integer(i4), dimension(:),   intent(in)  :: L0_LCover_LAI
    integer(i4), dimension(:),   intent(in)  :: L0_soilId !dim1=nCells1, dim2 fixed as 1 for the case flag_DB= 0
    integer(i4), dimension(:),   intent(in)  :: L0_geoUnit
    real(dp), dimension(:,:),  intent(inout) :: fgeounit
    real(dp), dimension(:,:),  intent(inout) :: frotation
    real(dp), dimension(:,:),  intent(inout) :: fLAI
    real(dp), dimension(:,:),  intent(inout) :: fLAI_11
    real(dp), dimension(:,:),  intent(inout) :: fsoil
    !INOUT
    ! Global variables for Nitrate submodel
    integer(i4),                   intent(in)    :: nsubstances
    integer(i4),                   intent(in)    :: nCroptation
    !type(rotationType),            intent(in)    :: rotation
    real(dp), dimension(:,:),      intent(inout) :: humusN
    real(dp), dimension(:,:),      intent(inout) :: fastN
    real(dp), dimension(:,:),      intent(inout) :: dissIN
    real(dp), dimension(:,:),      intent(inout) :: dissON
    real(dp), dimension(:,:,:),    intent(inout) :: csoilMoist
    real(dp), dimension(:,:),      intent(inout) :: csealSTW
    real(dp), dimension(:,:),      intent(inout) :: cunsatSTW
    real(dp), dimension(:,:),      intent(inout) :: csatSTW
    real(dp), dimension(:),        intent(inout) :: soiltemp
    real(dp), dimension(:),        intent(inout) :: basef_avg
    real(dp), dimension(:,:,:),    intent(inout) :: cinfiltration 
    real(dp), dimension(:,:),      intent(inout) :: cprec_effect
    real(dp), dimension(:,:),      intent(inout) :: crain
    real(dp), dimension(:,:),      intent(inout) :: cpercol
    real(dp), dimension(:,:),      intent(inout) :: crunoff_sealed
    real(dp), dimension(:,:),      intent(inout) :: cfastrunoff
    real(dp), dimension(:,:),      intent(inout) :: cslowrunoff
    real(dp), dimension(:,:),      intent(inout) :: cbaseflow
    real(dp), dimension(:,:),      intent(inout) :: ctotal_runoff
    real(dp), dimension(:),        intent(inout) :: soil_uptakeN
    real(dp), dimension(:),        intent(inout) :: soil_denitri
    real(dp), dimension(:),        intent(inout) :: soil_mineralN
    real(dp), dimension(:),        intent(inout) :: infrtmanapp
    real(dp), dimension(:),        intent(inout) :: degradN_rate
    real(dp), dimension(:),        intent(inout) :: mineraN_rate
    real(dp), dimension(:),        intent(inout) :: dissolN_rate
    real(dp), dimension(:),        intent(inout) :: sdenitr_rate
    !routing
    integer(i4),                   intent(in)   :: nNodes
    real(dp), dimension(:),        intent(in)   :: areaCell1
    real(dp), dimension(:),        intent(in)   :: areaCell11
    integer(i4), dimension(:),     intent(in)   :: L11id_on_L1
    integer(i4), dimension(:),     intent(in)   :: L1id_on_L11
    logical,                       intent(in)   :: map_flag
    integer(i4), dimension(:),     intent(in)   :: nLink_from
    integer(i4), dimension(:),     intent(in)   :: nLink_to
    real(dp), dimension(:),        intent(in)   :: nLink_length
    !inflow gauging station
    integer(i4),                   intent(in)   :: nInflowGauges
    logical, dimension(:),         intent(in)   :: InflowHeadwater
    integer(i4), dimension(:),     intent(in)   :: InflowNodeList
    integer(i4), dimension(:),     intent(in)   :: InflowIndexList
    real(dp), dimension(:),        intent(in)   :: Qinflow
    real(dp), dimension(:,:),      intent(inout):: cqInflow

    real(dp),                      intent(in)   :: globalradi_nor
    real(dp), dimension(:),        intent(inout):: reachtemp1
    real(dp), dimension(:),        intent(inout):: nLink_rivertemp
    real(dp), dimension(:),        intent(in)   :: nNode_qOUT
    real(dp), dimension(:),        intent(in)   :: nNode_qTR
    integer(i4), dimension(:),     intent(in)   :: nPerm
    real(dp), dimension(:),        intent(inout):: nLink_riverbox
    real(dp), dimension(:,:),      intent(inout):: nLink_criverbox 
    real(dp), dimension(:),        intent(inout):: nLink_yravg_q 
    real(dp), dimension(:,:),      intent(inout):: nNode_concOUT
    real(dp), dimension(:,:),      intent(inout):: nNode_interload
    real(dp), dimension(:,:),      intent(inout):: nNode_concTIN
    real(dp), dimension(:,:),      intent(inout):: nNode_concMod
    real(dp), dimension(:),        intent(inout):: aquatic_denitri
    real(dp), dimension(:),        intent(inout):: aquatic_assimil
    real(dp), dimension(:),        intent(inout):: adenitr_rate
    real(dp), dimension(:),        intent(inout):: atruptk_rate
    real(dp), dimension(:),        intent(inout):: priprod_rate
    real(dp), dimension(:),        intent(inout):: nLink_rzcoeff
    real(dp), dimension(:),        intent(inout):: nLink_flight

    !local
    real(dp), dimension(2)     :: potential_uptake !nitrate potential plant uptake amount in soil phase,
                                                     !  dim1=soillyrs, dim1=1(N)
    integer(i4)        :: k,i !L1 cells loop index
    !real(dp)           :: temp_reach   !river temperature at L1 
    integer(i4)        :: year,month,day,hour !date of current timestep
    logical            :: isday    !indicate day or night	

    integer(i4)        :: no_day, no_year ! day number in a year in current timestep,

	
    if (tt == 1)  then
      !area fraction of LAI(land-use) and rotation type
      !currently, regionalisation of nitrate parameters has not developed yet,
      !thus fraction of LAI type must be calculated for those land-use dependent parameters
       do i=1, nLAIclass
          fLAI(:,i) = L0_fractionalCover_in_Lx( int(L0_LCover_LAI), i, mask0, &
                                               L0upBound_inL1,    &
                                               L0downBound_inL1,  &
                                               L0leftBound_inL1,  &
                                               L0rightBound_inL1, &
                                               nTCells0_inL1)
       end do
      !new implementation for crop rotation
      !here, rotation type map is modified from land use map due to technical simplification
       do i=1, nCroptation
          frotation(:,i) = L0_fractionalCover_in_Lx( int(L0_cover_rotation), rotation(iBasin)%id(i), mask0, &
                                               L0upBound_inL1,    &
                                               L0downBound_inL1,  &
                                               L0leftBound_inL1,  &
                                               L0rightBound_inL1, &
                                               nTCells0_inL1)
       end do
      !new implementation for geounit

       do i=1, nGeoUnits
          fgeounit(:,i) = L0_fractionalCover_in_Lx( int(L0_geoUnit), GeoUnitList(i), mask0, &
                                               L0upBound_inL1,    &
                                               L0downBound_inL1,  &
                                               L0leftBound_inL1,  &
                                               L0rightBound_inL1, &
                                               nTCells0_inL1)
       end do
       !#############################################
       !New implementation for area fraction of each soil type (i.e., soil texture)
	   !the overall share of each soil texture in each grid cell of Level 1
       !added by yangx, 2019-03-13
	   !areal fraction of each soil type
       do i=1, nSoilTypes
          fsoil(:,i) = L0_fractionalCover_in_Lx( int(L0_soilId), i, mask0, &
                                               L0upBound_inL1,    &
                                               L0downBound_inL1,  &
                                               L0leftBound_inL1,  &
                                               L0rightBound_inL1, &
                                               nTCells0_inL1)
       end do
       !to get the share of each soil texture
	   !the overall percentage of sand, clay and the bulk density of each grid at Level 1
       !this information is currently coded as local variable, can be added as global variable if needed
       call L1_soiltexture(fsoil)
       !#############################################
	   
       !nitrate parameters upscaling
       if (perform_mpr) then
          call wqm_mpr( nNodes, L1id_on_L11, L11id_on_L1, map_flag, processMatrix,global_parameter, fLAI, fLAI_11, &
               degradN_rate, mineraN_rate, dissolN_rate,sdenitr_rate, adenitr_rate,atruptk_rate, priprod_rate, nLink_from, &
               areaCell1, areaCell11)
       end if

       !initial values for water quality model
       if (.not. read_restart) then
          call wqm_initialise(nCells1, nHorizons_mHM, HorizonDepth, nLAIclass, fLAI, humusN, fastN, dissIN, dissON, &
          soilMoisture, csoilMoist,cbaseflow, cunsatSTW, cpercol, fgeounit)    
       !!restart not yet implemented...
       !else
       !   call wqm_read_restart_states(iBasin)
       end if
       !special initial values needed
       call dec2date(time, yy=year_start)
       !at the first step, initialise vars from previous step(not exist)
       prevstep_sealedStorage = sealedStorage             
       prevstep_unsatStorage= unsatStorage         
       prevstep_satStorage= satStorage     
       prevstep_soilMoisture= soilMoisture
       !reachtemp1 = 0.0_dp
       !initial setting for groundwater conc. calculating
       prevstep_percol = percol
       prevstep_baseflow = baseflow
       basef_avg(:) = baseflow(:)

    end if
	!
 	
    !get date of global calendar from Julian day
    call dec2date(time, yy=year, mm=month, dd=day, hh=hour)
    isday = ( hour .gt. 6 ) .AND. ( hour .le. 18 )
	
    !get the number of days in current year and number of years
    no_day = ndays(day,month,year) - ndays(1,1,year) + 1
    no_year = year - year_start + 1
    !initial value for the days of previous year
    if (no_year == 1)  day_preyr =365_i4
    !update at the end of year
    if ((day == 31) .and. (month == 12)) then
       day_preyr = no_day   
    end if
 
    !*******************************
    !Soil phase
    !*******************************
    do k = 1, nCells1

       !preci. from snow melting is assumed as pure water
       if (prec_effect(k) > 0.0_dp) then 
       cprec_effect(k,1)=rain(k)* crain(k,1)/prec_effect(k)
       cprec_effect(k,2)=rain(k)* crain(k,2)/prec_effect(k)
       end if   

       !conc. update in sealed storage and soil moisture
       call soil_moisture_concentration(nsubstances,fsealed(k), prec_effect(k), sealedStorage(k), &
         infiltration(k,:), soilMoisture(k,:), aet_soil(k,:), aet_sealed(k), prevstep_sealedStorage(k), &
         prevstep_soilMoisture(k,:), csealSTW(k,:), cinfiltration(k,:,:), csoilMoist(k,:,:),&
         cprec_effect(k,:),crunoff_sealed(k,:))

	   !--temporally settings---
       cfastrunoff(k,:) = csoilMoist(k,1,:)
	   
       !conc. update in unsat- and saturated storage
       call runoff_concentration(infiltration(k,nHorizons_mHM), cinfiltration(k,nHorizons_mHM,:),  &       
         prevstep_unsatStorage(k), prevstep_percol(k), prevstep_baseflow(k), basef_avg(k), cunsatSTW(k,:), & 
         cpercol(k,:),cslowrunoff(k,:), cbaseflow(k,:),satStorage(k))
       
       !conc. in total runoff that comes out from L1 cell
       call totalrunoff_concentration(fsealed(k),runoff_sealed(k), fastRunoff(k), slowRunoff(k), baseflow(k),  &
         total_runoff(k), crunoff_sealed(k,:), cfastrunoff(k,:), cslowrunoff(k,:), cbaseflow(k,:),   &
         ctotal_runoff(k,:))

       !agricultural management--nutrient input from fertilizer and manure application (crop data)
       !AND potential plant uptake (potential_uptake)
      
       call agri_management(timeStep, iBasin, no_day, no_year, day_preyr,frotation(k,:), soilMoisture(k,:), &
         csoilMoist(k,:,:), fastN(k,:), humusN(k,:), temp(k), potential_uptake, nCroptation,infrtmanapp(k))       
	   
       !update dissolved ON and IN pools		  
       dissIN(k,:) = soilMoisture(k,:) * csoilMoist(k,:,1)
       dissON(k,:) = soilMoisture(k,:) * csoilMoist(k,:,2)
	   
	   !transformation between differnt pools AND DENITRIFICATION
       call soil_nutrient_transformation(timeStep, nHorizons_mHM, HorizonDepth(:), temp(k), soiltemp(k), snowpack(k), &
          wilting_point(k,:),soilmoist_sat(k,:), soilMoisture(k,:), csoilMoist(k,:,:), fastN(k,:), humusN(k,:), &
          dissIN(k,:), dissON(k,:), degradN_rate(k), mineraN_rate(k), dissolN_rate(k), soil_mineralN(k))
       
       !plant uptake and denitrification in soil phase 
       call soil_plant_uptake(wilting_point(k,:), soilMoisture(k,:), csoilMoist(k,:,:), &
          potential_uptake, soil_uptakeN(k))
       dissIN(k,:) = soilMoisture(k,:) * csoilMoist(k,:,1)

       call soil_denitrification(timeStep, nHorizons_mHM, soilmoist_sat(k,:), soiltemp(k), soilMoisture(k,:), csoilMoist(k,:,:), &
          dissIN(k,:), sdenitr_rate(k), soil_denitri(k))

       !water temperature of regional stream in L1 cells, 20-day's average of L1 cell air temperature
    reachtemp1(k) = reachtemp1(k) + (temp(k) - reachtemp1(k)) / 20.0_dp


    end do   
    !**************************
    !In-stream phase (routing)
    !**************************
    ! including the point sources from sewage plants, using the structure for inflow gauging stations
    call conc_routing_accmix(nsubstances, nNodes, timeStep, total_runoff, &
        ctotal_runoff,areaCell1,L11id_on_L1, L1id_on_L11, map_flag, nLink_from, reachtemp1,nInflowGauges,  &
       InflowHeadwater, InflowNodeList,InflowIndexList,  Qinflow, cqInflow, nLink_rivertemp, nNode_concOUT )

    !in-stream processes
    call conc_routing_process(nNodes, nNodes-1 ,no_day, timeStep, nPerm, nLink_from, nLink_to, nLink_length, &
       nNode_qOUT,nNode_qTR, nNode_concOUT,nNode_interload,nNode_concTIN,nLink_riverbox, nLink_criverbox,nLink_yravg_q,&
       nLink_rivertemp,nNode_concMod, aquatic_denitri, aquatic_assimil, adenitr_rate, atruptk_rate,priprod_rate, &
       globalradi_nor, nLink_rzcoeff, nLink_flight)


  end subroutine wqm

  !--------------------------------------------------------------------
  ! ------------------------------------------------------------------

  !     NAME
  !         conc_routing_process

  !     PURPOSE
  !>        \brief water qulity processes within routing. 

  !>        \details calculates river water concentration including in-stream processes including mixture in each node and 
  !>        denitrification and assimilatory in each reach(link).         

  !     INTENT(IN)
  !>        \param[in] "integer(i4)                     :: nNodes"
  !>        \param[in] "integer(i4)                     :: nLinks"
  !>        \param[in] "integer(i4)                     :: TS"                simulation time step [h] 
  !>        \param[in] "integer(i4)                     :: nsubstances"
  !river networks (links)
  !>        \param[in] "integer(i4), dimension(:)       :: nPerm"             basin routing order
  !>        \param[in] "integer(i4), dimension(:)       :: nLink_from"        from node
  !>        \param[in] "integer(i4), dimension(:)       :: nLink_to"          to node
  !>        \param[in] "real(dp), dimension(:)          :: nLink_length"      length of each reach(link)
  !variables from hydrological model
  !>        \param[in] "real(dp), dimension(:)          :: nNode_qOUT"        total outflow from L11 cells
  !>        \param[in] "real(dp), dimension(:)          :: nNode_qTR"      outflow from each reach, only IT=2 was called by WQM
  !variables for in-stream WQ processes
  !>        \param[in] "real(dp), dimension(:,:)        :: nNode_concOUT"     conc. of total outflow from L11 cells  
  !>        \param[in] "real(dp), dimension(:)          :: rivertemp11"       reach temperature
  !>        \param[in] "real(dp), dimension(:)          :: pardeniratew"      denitrification rate in aquatic 
  !>        \param[in] "real(dp), dimension(:),         :: parautouptkrate"   potiential autotrophic uptake mgNm-2d-1
  !>        \param[in] "real(dp), dimension(:)          :: parprodrate"       assimilatory rate in aquatic   

  !     INTENT(INOUT)
  !>        \param[inout] "real(dp), dimension(:,:)        :: nNode_interload"    variables to store input load at each node
  !>        \param[inout] "real(dp), dimension(:,:)        :: nNode_concTIN"      conc. of total inflow of each nodes
  !>        \param[inout] "real(dp), dimension(:)          :: nLink_riverbox"     water volumn of each reach (link)
  !>        \param[inout] "real(dp), dimension(:,:)        :: nLink_criverbox"    conc. in the water volumn
  !>        \param[inout] "real(dp), dimension(:)          :: nLink_yravg_q"      yearly average discharge of each reach
  !results
  !>        \param[inout] "real(dp), dimension(:,:)   :: nNode_concMod"      conc. at each node (ouput variable!)  
  !>        \param[inout] "real(dp), dimension(:)     :: aquatic_denitri"    amount of denitrification in aquatic system (instream)
  !>        \param[inout] "real(dp), dimension(:)     :: aquatic_assimil"    amount of assimilatory(primary production)

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
  !>        \author Xiaoqiang Yang
  !>        \date Sep 2016
  subroutine conc_routing_process(nNodes, nLinks, no_day, TS, nPerm, nLink_from, nLink_to, nLink_length,  &
       nNode_qOUT, nNode_qTR, nNode_concOUT, nNode_interload, nNode_concTIN, nLink_riverbox, nLink_criverbox, nLink_yravg_q, &
       rivertemp11, nNode_concMod, aquatic_denitri, aquatic_assimil,pardeniratew, parautouptkrate, parprodrate, &
       nor_gr,rz_coeff, f_light)
  
  use mo_wqm_global_variables,     only: L11_rivert_avg10, L11_rivert_avg20  !, L11_flightglobal state variables for moving average
  
  implicit none

  integer(i4),                     intent(in)    :: nNodes
  integer(i4),                     intent(in)    :: nLinks
  integer(i4),                     intent(in)    :: no_day     !number of day accounter in current year
  integer(i4),                     intent(in)    :: TS         !simulation time step [h] 

  !river networks (links)
  integer(i4), dimension(:),       intent(in)    :: nPerm      ! basin routing order
  integer(i4), dimension(:),       intent(in)    :: nLink_from ! from node
  integer(i4), dimension(:),       intent(in)    :: nLink_to   ! to node
  real(dp), dimension(:),          intent(in)    :: nLink_length ! length of each reach(link)
  !variables from hydrological model
  real(dp), dimension(:),          intent(in)    :: nNode_qOUT ! total outflow from L11 cells
  real(dp), dimension(:),          intent(in)    :: nNode_qTR  ! outflow from each reach, for dim2 only IT was called by WQM
  !variables for in-stream WQ processes
  real(dp), dimension(:,:),        intent(in)    :: nNode_concOUT     ! conc. of total outflow from L11 cells  
  real(dp), dimension(:,:),        intent(inout) :: nNode_interload   ! variables to store input load at each node
  real(dp), dimension(:,:),        intent(inout) :: nNode_concTIN     ! conc. of total inflow of each nodes
  real(dp), dimension(:),          intent(inout) :: nLink_riverbox    ! water volumn of each reach (link)
  real(dp), dimension(:,:),        intent(inout) :: nLink_criverbox   ! conc. in the water volumn
  real(dp), dimension(:),          intent(inout) :: nLink_yravg_q     ! yearly average discharge of each reach
  real(dp), dimension(:),          intent(in)    :: rivertemp11       ! reach temperature
  !results
  real(dp), dimension(:,:),        intent(inout) :: nNode_concMod     ! conc. at each node (ouput variable!)  
  real(dp), dimension(:),          intent(inout) :: aquatic_denitri   ! amount of denitrification in aquatic system (instream)
  real(dp), dimension(:),          intent(inout) :: aquatic_assimil   ! amount of assimilatory(primary production)
  !WQ parameters
  real(dp), dimension(:),          intent(in)    :: pardeniratew      !denitrification rate in aquatic 
  real(dp), dimension(:),          intent(in)    :: parautouptkrate   !potiential autotrophic uptake mgNm-2d-1
  real(dp), dimension(:),          intent(in)    :: parprodrate       !assimilatory rate in aquatic  
  real(dp),                        intent(in)    :: nor_gr            !normalised global radiation (daily)  
  real(dp), dimension(:),          intent(inout)    :: rz_coeff          !rz shading coefficient dim = number of reaches
  real(dp), dimension(:),          intent(inout)    :: f_light          !rz shading coefficient dim = number of reaches
  
  !local
  real(dp), dimension(nNodes)  :: temp_qTIN! store total inflow of each node, represents "netNode_qTIN" in hydrological routing
  integer(i4)             :: sec_TS, DT      ! seconds per timestep, number of steps in one day
  integer(i4)             :: k,i
  real(dp)                :: width, depth, benthic_area   ![m] reach(link) width, depth and beneath area
  real(dp)                :: newbox
  integer(i4)             :: iNode, tNode

  
  temp_qTIN = 0.0_dp
  nNode_interload = 0.0_dp
  !constant parameter  
  sec_TS = TS* 3600_i4    ! [s]
  DT = 24_i4 / TS

  do k=1, nLinks  
     i = nPerm(k)
     iNode = nLink_from(i)
     tNode = nLink_to(i)
     ! yearly average channel discharge (for empirial channel width and depth calculation)
     nLink_yravg_q(i) = nLink_yravg_q(i) + (nNode_qTR(iNode) - nLink_yravg_q(i)) / (90.0_dp )  
     ! empirial equations for channel width and depth
     width = 5.40_dp *nLink_yravg_q(i) **(0.50_dp)  ! from M. Rode (2016), EST
     depth = 0.27_dp *nLink_yravg_q(i) **(0.39_dp)  ! from J.A. Moody(2002), Earth Surface Processes and Landforms
     benthic_area = width * nLink_length(i)

     !total input at "iNode"
     temp_qTIN(iNode) = temp_qTIN(iNode) + nNode_qOUT(iNode)
     !conc. of total infow of node(at iNode)
     nNode_concTIN(iNode,:) = (nNode_interload(iNode,:) + &
                nNode_qOUT(iNode) * nNode_concOUT(iNode,:) * sec_TS) / (temp_qTIN(iNode) * sec_TS)
     newbox = nLink_riverbox(i) + temp_qTIN(iNode) * sec_TS
     

     !update conc. and water volumn of reach i
     nLink_criverbox(i,:) = (nLink_criverbox(i,:) * nLink_riverbox(i) + nNode_interload(iNode,:) + &
                nNode_qOUT(iNode) * nNode_concOUT(iNode,:) * sec_TS) / newbox          
     nLink_riverbox(i) = newbox
     !instream denitrification and primary production
     !10- and 20-day moving mean temperature of river water
     L11_rivert_avg10(i) = L11_rivert_avg10(i) + (rivertemp11(i) - L11_rivert_avg10(i)) / (10.0_dp )
     L11_rivert_avg20(i) = L11_rivert_avg20(i) + (rivertemp11(i) - L11_rivert_avg20(i)) / (20.0_dp )
     call instream_nutrient_processes(TS,no_day,i, rivertemp11(i), L11_rivert_avg10(i), L11_rivert_avg20(i), &
         benthic_area, depth, nLink_riverbox(i), nLink_criverbox(i,:),aquatic_denitri(i), aquatic_assimil(i), &
        pardeniratew(i),parautouptkrate(i), parprodrate(i), nor_gr, rz_coeff(i), f_light(i))

     !update water volumn of reach(link) 
     nLink_riverbox(i) = nLink_riverbox(i) - nNode_qTR(iNode) * sec_TS
     !calculate the load from each upstream reach and store as input of the to_node (tNode) 
     nNode_interload(tNode,:) = nNode_interload(tNode,:) + nLink_criverbox(i,:) * nNode_qTR(iNode) * sec_TS
     !accumulate flow to the 'to node' of the current link, as the upstream inflow of 'to node'
     temp_qTIN(tNode) = temp_qTIN(tNode) + nNode_qTR(iNode)

  end do
  !*************************************
  !accumulate at the outlet of catchment
  !*************************************
  iNode = nLink_from(nPerm(nLinks))
  tNode = nLink_to(nPerm(nLinks))
  temp_qTIN(tNode) = temp_qTIN(tNode) + nNode_qOUT(tNode)
  nNode_concTIN(tNode,:) = (nNode_interload(tNode,:) + nNode_qOUT(tNode)* nNode_concOUT(tNode,:) * sec_TS) / &
                           (temp_qTIN(tNode) * sec_TS)

  !variable for final concentration output
  do i=1, nNodes
  nNode_concMod(i,:) = nNode_concTIN(i,:)  
  end do
  
  end subroutine conc_routing_process
  !--------------------------------------------------------------------
  ! ------------------------------------------------------------------

  !     NAME
  !         instream_nutrient_processes

  !     PURPOSE
  !>        \brief In-stream processes in a reach(link). 

  !>        \details Calculates in-stream processes in a sepcific reach(link), including denitrification and assimilatory.
  !>                 

  !     INTENT(IN)
  !>        \param[in] "integer(i4)                    :: TS"            time step         
  !>        \param[in] "integer(i4)                    :: i"         from node id of this sepcific reach
  !>        \param[in] "real(dp)                       :: rivertemp11"   river water temperature
  !>        \param[in] "real(dp)                       :: rt_avg10"      10-day's average water temperature
  !>        \param[in] "real(dp)                       :: rt_avg20"      20-day'S average water temperature
  !>        \param[in] "real(dp)                       :: benthic_area"  [m2]reach benthic area
  !>        \param[in] "real(dp)                       :: depth"         [m]average reach depth 
  !>        \param[in] "real(dp)                       :: pardeniratew"  parameter
  !>        \param[in] "real(dp)                       :: parautouptkrate"  parameter
  !>        \param[in] "real(dp)                       :: parprodrate"   parameter
  !     INTENT(INOUT)
  !>        \param[inout] "real(dp)                       :: riverbox"    river water volumn in each reach  
  !>        \param[inout] "real(dp), dimension(:)         :: criverbox"   conc. in water wolumn 


  !     INTENT(OUT)
  !>        \param[out] "real(dp)                       :: aqdenitri"   amount of N denitrified in river water 
  !>        \param[out] "real(dp)                       :: aqassimil"   amount of N assimilated in river water

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
  !         HYPE model

  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Sep 2016
  !>    Modified
  !>        X. Yang Jun 2017 enabled different time-step (hourly)

  subroutine instream_nutrient_processes(TS,no_day,i, rivertemp11, rt_avg10, rt_avg20, benthic_area, depth, &
        riverbox, criverbox, aqdenitri, aqassimil, pardeniratew,parautouptkrate, parprodrate,nor_gr, rz_coeff, f_light)

  use mo_wqm_shadingeffect,      only: rz_shading_coeff
  use mo_wqm_global_variables,   only: GR_file_exist
  implicit none
  integer(i4),                    intent(in)     :: TS
  integer(i4),                    intent(in)     :: no_day     ! for light calculate -- shading effect
  integer(i4),                    intent(in)     :: i
  real(dp),                       intent(in)     :: rivertemp11
  real(dp),                       intent(in)     :: rt_avg10
  real(dp),                       intent(in)     :: rt_avg20
  real(dp),                       intent(in)     :: benthic_area  
  real(dp),                       intent(in)     :: depth    
  real(dp),                       intent(inout)  :: riverbox  
  real(dp), dimension(:),         intent(inout)  :: criverbox    
  real(dp),                       intent(out)  :: aqdenitri    
  real(dp),                       intent(out)  :: aqassimil    
  real(dp),                       intent(in)     :: pardeniratew
  real(dp),                       intent(in)     :: parautouptkrate
  real(dp),                       intent(in)     :: parprodrate  
  real(dp),                       intent(in)     :: nor_gr   ! for shading effect
  real(dp),                       intent(inout)  :: rz_coeff ! coefficient for each reach
  real(dp),                       intent(inout)  :: f_light  ! 5 days' moving average of rz_coeff


  !local
  real(dp)                 :: INpool, ONpool, tp_ave,aqassimil0
  real(dp)                 :: f_temp, f_conc, f_tp, f_temp1, f_temp2
  real(dp)              :: DT
  
  !HYPE parameter
  real(dp), parameter      :: tpmean = 0.005_dp
  ! constant variables
  real(dp), parameter      :: halfsatINwater = 1.5  ! mg/l
  real(dp), parameter      :: maxdenitriwater = 0.5_dp !
  real(dp), parameter      :: maxprodwater = 0.5_dp  ! 
  !real(dp), parameter      :: maxdegradwater = 0.5_dp  
  real(dp), parameter      :: halfsatIPwater = 0.05_dp  
  real(dp), parameter      :: activedepth = 1.0_dp  
  
  DT = 24.0_dp / TS

 
  INpool = riverbox * criverbox(1) / 1000.0_dp   !kg
  f_temp = tempfactor(rivertemp11)
  f_conc = criverbox(1) / (criverbox(1) + halfsatINwater )

  !denitrification amount
  aqdenitri = pardeniratew * f_temp * f_conc * benthic_area / DT
  aqdenitri = min(maxdenitriwater*INpool, aqdenitri) 

  !update pool and conc.
  INpool = INpool - aqdenitri
  if (riverbox >0.0_dp) then
    criverbox(1) = INpool / riverbox * 1000.0_dp   !mg/l
  else
    criverbox(1) = 0.0_dp
  end if
  !primary production and mineralisation (inverse processes)
  if (GR_file_exist) then
    !##############################################
    !##new method considering light availability ##
    !##############################################
    call rz_shading_coeff(i, no_day, nor_gr, rz_coeff)
    f_light = f_light + (rz_coeff-f_light) / 5.0_dp

    aqassimil =0.0_dp
    aqassimil0 = 0.0_dp
    if (rivertemp11 >= 0.0_dp) then
       aqassimil0 =  parautouptkrate* f_light * activedepth * benthic_area * 10E-6 /DT   ![kg]
	   ! here the parautouptkrate (potential N autotrophic uptake = 300 mg N m-2 d-1) 
       aqassimil0 = min(0.9_dp *INpool*activedepth, aqassimil0 )  
       !mineralization (respiration) assumed part of autotrophic assimilation
       aqassimil = parprodrate*aqassimil0
       !aqassimil = min(maxprodwater*INpool*activedepth, aqassimil )
    end if
  else
    !##method from HYPE##    
    tp_ave = tpmean   ! since phorsphose is not implemented at the monment..
    f_tp = tp_ave / (tp_ave + halfsatIPwater)
    ONpool = riverbox * criverbox(2) /1000.0_dp
    if (rivertemp11 > 0.0_dp) then
       f_temp1 = rivertemp11 /20.0_dp
    else
       f_temp1 = 0.0_dp
    end if  
    f_temp2 = (rt_avg10 - rt_avg20) / 5.0_dp  
    f_temp = f_temp1 * f_temp2
    aqassimil =0.0_dp
    if (rivertemp11 >= 0.0_dp) then
       aqassimil = parprodrate * f_temp * f_tp * activedepth * depth * benthic_area/ DT !
       if (aqassimil > 0.0_dp) then
          aqassimil = min(maxprodwater*INpool*activedepth, aqassimil ) 
       else
          aqassimil = max(-maxprodwater*ONpool*activedepth, aqassimil) 
       end if
    end if
  end if

  
  
  INpool = INpool - aqassimil
  ONpool = ONpool + aqassimil
  if (riverbox >0.0_dp) then
    criverbox(1) = INpool / riverbox * 1000.0_dp   !mg/l
    criverbox(2) = ONpool / riverbox * 1000.0_dp   !mg/l
  else
    criverbox(1) = 0.0_dp
    criverbox(2) = 0.0_dp
  end if 
  
  !!for writing out pure assimiloray uptake
  !aqassimil = aqassimil0
   


  end subroutine instream_nutrient_processes 
  !--------------------------------------------------------------------
  ! ------------------------------------------------------------------

  !     NAME
  !         conc_routing_accmix

  !     PURPOSE
  !>        \brief Accumulation and mixture of concentration in cells at L11 level.   

  !>        \details Based on the resolution of L1 and L11, concentration of each node at L11 is calculated from conc. of
  !>         totalrunoff and mixed with the weight of total runoff of each cell at L1.          

  !     INTENT(IN)
  !>        \param[in] "integer(i4)                    :: nsubstances"  muber of substances involved
  !>        \param[in] "integer(i4)                    :: nNodes"       number of cells at L11
  !>        \param[in] "real(dp), dimension(:)         :: qAll"         total runoff from L1 cell
  !>        \param[in] "integer(i4)                    :: TS"           timestep
  !>        \param[in] "real(dp), dimension(:,:)       :: cqAll"        conc. of total runoff (above)
  !>        \param[in] "real(dp), dimension(:)         :: areaCell1"    [km2] cell area at level L1
  !>        \param[in] "real(dp), dimension(:)         :: areaCell11"   [km2] cell area at level L11  
  !>        \param[in] "integer(i4), dimension(:)      :: L11id_on_L1"  L11 id on L1 map
  !>        \param[in] "integer(i4), dimension(:)      :: L1id_on_L11"  L1 id on L11 map
  !>        \param[in] "logical                        :: map_flag"     L11 IS bigger than L1 OR NOT
  !>        \param[in] "integer(i4), dimension(:)      :: nLink_from"   link from node 
  !>        \param[in] "real(dp), dimension(:)         :: reachtemp1"   temperature of regional reach in L1 cells        

  !     INTENT(INOUT)
  !>        \param[inout] "real(dp), dimension(:)         :: rivertemp11   wate temperature in river
  !>        \param[inout] "real(dp), dimension(:,:)       :: cqOUT         conc. of outflow of L11 cell after mixing

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
  !>        \author Xiaoqiang Yang
  !>        \date Sep 2016

  subroutine conc_routing_accmix(nsubstances, nNodes, TS, qAll, cqAll,areaCell1, L11id_on_L1,&
      L1id_on_L11,map_flag,nLink_from, reachtemp1, nInflowGauges, InflowHeadwater, InflowNodeList, InflowIndexList,  &
      Qinflow, cqInflow, rivertemp11, cqOUT)

  !use mo_mrm_constants,           only: nodata_dp
	  
  implicit none 
  integer(i4),                    intent(in)     :: nsubstances  ! muber of substances involved
  integer(i4),                    intent(in)     :: nNodes       ! number of cells at L11
  integer(i4),                    intent(in)     :: TS           ! [h]timestep
  real(dp), dimension(:),         intent(in)     :: qAll         ! total runoff from L1 cell
  real(dp), dimension(:,:),       intent(in)     :: cqAll        ! conc. of total runoff (above)
  real(dp), dimension(:),         intent(in)     :: areaCell1    ! [km2] cell area at level L1
!  real(dp), dimension(:),         intent(in)     :: areaCell11   ! [km2] cell area at level L11  
  integer(i4), dimension(:),      intent(in)     :: L11id_on_L1  ! L11 id on L1 map
  integer(i4), dimension(:),      intent(in)     :: L1id_on_L11  ! L1 id on L11 map
  logical,                        intent(in)     :: map_flag     ! L11 IS bigger than L1 OR NOT
  integer(i4), dimension(:),      intent(in)     :: nLink_from   ! link from node 
  real(dp), dimension(:),         intent(in)     :: reachtemp1   ! temperature of regional reach in L1 cells 
  !variables from routing model
  integer(i4),                    intent(in)     :: nInflowGauges  !  
  logical, dimension(:),          intent(in)     :: InflowHeadwater!
  integer(i4), dimension(:),      intent(in)     :: InflowNodeList  
  integer(i4), dimension(:),      intent(in)     :: InflowIndexList
  real(dp), dimension(:),         intent(in)     :: Qinflow
  real(dp), dimension(:,:),       intent(inout)  :: cqInflow


  real(dp), dimension(:),         intent(inout)  :: rivertemp11  ! wate temperature in river
  real(dp), dimension(:,:),       intent(inout)  :: cqOUT        ! conc. of outflow of L11 cell after mixing
  !local
  real(dp), dimension(nNodes)     :: sumq, temperature11       ! temporal variables, instead of discharge/temperature   
  integer(i4)                     :: k, i, mm, nn 
  real(dp)                        :: jj
  real(dp)                        :: tmpvar, tst       !temporal variables, transfer total runoff(mm) to discharge(m^3/s)
                                                       !timestep[s]  
  sumq = 0.0_dp
  temperature11 = 0.0_dp
  cqOUT =0.0_dp
  tst = real(TS, dp) * 3600_dp   ! transfer from Hour to Second
  
  if (map_flag) then    ! if L11 >= L1
  do k =1, nNodes  

     jj = 0.0_dp
     do i = 1, size(L11id_on_L1)
     if (L11id_on_L1(i) .eq. k ) then
        jj= jj + 1.0_dp
        tmpvar = qAll(i) * areaCell1(i) * 1000.0_dp / tst
        cqOUT(k,:) = (cqOUT(k,:) * sumq(k) + cqAll(i,:) * tmpvar) /(sumq(k) + tmpvar )
        sumq(k) = sumq(k) + tmpvar
        temperature11(k) = temperature11(k) + reachtemp1(i)
     

     end if
     end do

     temperature11(k)= temperature11(k) /jj   !average reach temperature at L11
     !
     !identify reach id
     mm = minloc(abs(nLink_from - k),1)
     !assign aggragated temperature at L11 to the linked reach
     rivertemp11(mm) = temperature11(k)

  end do

  else                     ! if L11 < L1
    do k = 1, size(qAll)
     do i =1, nsubstances
        cqOUT(:,i) = merge(cqAll(k,i), 0.0_dp, L1id_on_L11(:) .eq. k ) 
     end do
     temperature11(:) = merge(reachtemp1(k), 0.0_dp, L1id_on_L11(:) .eq. k)
    end do
	!
    do k = 1, nNodes
     mm = minloc(abs(nLink_from - k),1)
     rivertemp11(mm) = temperature11(k)
    end do     
  end if
  
  !ADDITIONAL INFLOW STATION (point source and upstream inflow guages)....
  if (nInflowGauges .gt. 0_i4 ) then
     where (cqInflow .lt. 0.0_dp) cqInflow = 0.0_dp
     do nn=1, nInflowGauges
        if (InflowHeadwater(nn)) then
           cqOUT(InflowNodeList(nn),:) = cqInflow(InflowIndexList(nn), 1:2) 
        else
           cqOUT(InflowNodeList(nn),:) = (cqOUT(InflowNodeList(nn),:)* sumq(InflowNodeList(nn))+ &
              cqInflow(InflowIndexList(nn),1:2)*Qinflow(InflowIndexList(nn)) )/ &
              (sumq(InflowNodeList(nn))+Qinflow(InflowIndexList(nn)))   
        end if
     end do
    
  end if

  end subroutine conc_routing_accmix
  !--------------------------------------------------------------------
  ! ------------------------------------------------------------------
  !     NAME
  !         wqm_initialise

  !     PURPOSE
  !>        \brief initialisation of water quality model. 

  !>        \details initialises WQM variables at the beginning of each basin. 
  !>        calculates initial values of four different nutrient(Nitrogen, at the moment) pools, which is land-use dependent        
  !>         and is read in from input file "initial_value.txt"
  !     INTENT(IN)
  !>        \param[in] "integer(i4)                            :: nCells1"
  !>        \param[in] "integer(i4)                            :: nHorizons_mHM"
  !>        \param[in] "real(dp), dimension(:)                 :: horizon_depth"
  !>        \param[in] "integer(i4)                            :: nLAIclass"
  !>        \param[in] "real(dp), dimension(:,:)               :: fLAI"
  !>        \param[in] "real(dp), dimension(:,:)               :: soilMoisture"
  !>        \param[in] "real(dp), dimension(:)                 :: sealedStorage" 
  !>        \param[in] "real(dp), dimension(:)                 :: unsatStorage"
  !>        \param[in] "real(dp), dimension(:)                 :: satStorage"

 !     INTENT(INOUT)
  !>        \param[inout] "real(dp), dimension(:,:)           :: humus_N"
  !>        \param[inout] "real(dp), dimension(:,:)           :: fast_N"
  !>        \param[inout] "real(dp), dimension(:,:)           :: diss_IN"
  !>        \param[inout] "real(dp), dimension(:,:)           :: diss_ON"
  !>        \param[inout] "real(dp), dimension(:,:,:)         :: csoilMoist"
  !>        \param[inout] "real(dp), dimension(:,:)           :: cbaseflow"
  !>        \param[inout] "real(dp), dimension(:,:)           :: cunsatSTW"
  !>        \param[inout] "real(dp), dimension(:,:)           :: cpercol"
  !>        \param[inout] "real(dp), dimension(:)             :: reachtemp1"
  !>        \param[inout] "real(dp), dimension(:)             :: prevstep_sealedStorage"
  !>        \param[inout] "real(dp), dimension(:)             :: prevstep_unsatStorage"
  !>        \param[inout] "real(dp), dimension(:)             :: prevstep_satStorage" 
  !>        \param[inout] "real(dp), dimension(:)             :: prevstep_soilMoisture"

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
  !         HYPE model

  !     HISTORY
  !>        \authors Xiaoqiang Yang 
  !>        \date Jun 2016

  subroutine wqm_initialise(nCells1, nHorizons_mHM, horizon_depth, nLAIclass, fLAI, humus_N, fast_N, diss_IN, diss_ON, &
      soilMoisture, csoilMoist,cbaseflow, cunsatSTW, cpercol, fgeounit)

  use mo_wqm_global_variables,   only: init_humusN, init_fastN, init_concIN, init_concON, hnhalf,Geoform  
                                     !read in from file: initial_value.txt
                                     !readin init_humusN&init_fastN unit[mg/m^3] init_concIN&init_concON [mg/l]
  implicit none

  integer(i4),                        intent(in)    :: nCells1
  integer(i4),                        intent(in)    :: nHorizons_mHM
  real(dp), dimension(:),             intent(in)    :: horizon_depth   ![mm] Horizon depth from surface
  integer(i4),                        intent(in)    :: nLAIclass
  real(dp), dimension(:,:),           intent(in)    :: fLAI
  real(dp), dimension(:,:),           intent(inout) :: humus_N         ![mg/m^2]
  real(dp), dimension(:,:),           intent(inout) :: fast_N          ![mg/m^2]
  real(dp), dimension(:,:),           intent(inout) :: diss_IN         ![mg/m^2]
  real(dp), dimension(:,:),           intent(inout) :: diss_ON         ![mg/m^2]
  real(dp), dimension(:,:),           intent(in)    :: soilMoisture    ![mm]
  real(dp), dimension(:,:,:),         intent(inout) :: csoilMoist      ![mg/l]
  real(dp), dimension(:,:),           intent(inout) :: cbaseflow       ![mg/l]
  real(dp), dimension(:,:),           intent(inout) :: cunsatSTW       ![mg/l]
  real(dp), dimension(:,:),           intent(inout) :: cpercol         ![mg/l]
  real(dp), dimension(:,:),           intent(in)    :: fgeounit

  !local
  integer(i4)     :: j, k,i
  real(dp)        :: init_humusN0, init_fastN0, init_concIN0, init_concON0  
  real(dp)        :: hnhalf0,half_N,f_imper
  real(dp), dimension(size(horizon_depth))  :: horizon_depth1   !temporal variable for soil layer depth
  !         

  horizon_depth1 = horizon_depth/1000.0_dp        !mm -> m
  do k=1, nCells1
     init_humusN0 = 0.0_dp              !kg/km^2 (mg/m^2)
     init_fastN0 = 0.0_dp               !kg/km^2
     init_concIN0 = 0.0_dp
     init_concON0 = 0.0_dp
     hnhalf0 = 0.0_dp
     do j=1, nLAIclass
     if (fLAI(k,j) > 0.0_dp) then
        init_humusN0 = init_humusN0 + fLAI(k, j) * init_humusN(j) * horizon_depth1(1)
        init_fastN0  = init_fastN0 + fLAI(k, j) * init_fastN(j) * horizon_depth1(1)
        init_concIN0 = init_concIN0 + fLAI(k, j) * init_concIN(j)
        init_concON0 = init_concON0 + fLAI(k, j) * init_concON(j)
        hnhalf0 = hnhalf0 + fLAI(k,j) * hnhalf(j)
     end if
     end do        
    fast_N(k,:) = init_fastN0
    humus_N(k,1) = init_humusN0
      half_N = log(2.0_dp)/ hnhalf0
    if (nHorizons_mHM > 1_i4) then  
       humus_N(k,2) = init_humusN0 * exp(-half_N * horizon_depth1(2)/2.0_dp)
    end if
    if (nHorizons_mHM > 2_i4) then
       humus_N(k,3) = init_humusN0 * exp(-half_N * (horizon_depth1(3) + horizon_depth1(2) - horizon_depth1(1))/2.0_dp)
    end if   
    csoilMoist(k,3,1) = init_concIN0
    diss_IN(k,3) = csoilMoist(k,3,1) * soilMoisture(k,3)
    !cbaseflow(k,1) = init_concIN0
    cunsatSTW(k,1) = init_concIN0
    csoilMoist(k,3,2) = init_concON0
    diss_ON(k,3) = csoilMoist(k,3,2) * soilMoisture(k,3)
    cbaseflow(k,2) = init_concON0
    cunsatSTW(k,2) = init_concON0
	
    cpercol(k,:) = cunsatSTW(k,:)
  
    ! !Special setting for test site (Selke catchment)
    ! !Because the arable land in the upper part(mountain area) do not observe high conc. in ground water
    ! if (k >= 320) then
       ! cbaseflow(k,1) =0.5_dp

    ! end if
    !####a more general implementation#########
    !considering the geological condition of each grid cell using the geounit dependent classification
    !if Geoform = 1 then the same as original
    !if Geoform = 2 then a initial value 0.5 is assigned
    f_imper = 0.0_dp
    do i = 1, size(Geoform)
       if (Geoform(i) == 2) then
          f_imper = f_imper + fgeounit(k,i)
       end if
    end do
    cbaseflow(k,1) =  0.5_dp* f_imper + init_concIN0 * (1- f_imper)
	!here 0.5 mg/l is given as a default Nitrate-N concentration for deep groundwater in forests


  end do

  end subroutine wqm_initialise
  !-------------------------------------------------------------------
  ! ------------------------------------------------------------------

  !     NAME
  !         soil_denitrification

  !     PURPOSE
  !>        \brief Nitrate denitrification in soil water 

  !>        \details calculates nitrate denitrification in soil water, considering temperature, concentration, 
  !>        and soil moisture factors.         

  !     INTENT(IN)
  !>        \param[in] "integer(i4)                    :: nHorizons_mHM"
  !>        \param[in] "real(dp), dimension(:)         :: sat"
  !>        \param[in] "real(dp)                       :: soiltemp"
  !>        \param[in] "real(dp), dimension(:)         :: soilmoist"        
  !>        \param[in] "real(dp)                       :: pardenirate"
  !     INTENT(INOUT)
  !>        \param[inout] "real(dp), dimension(:,:)    :: concsoil"
  !>        \param[inout] "real(dp), dimension(:)      :: diss_IN"

  !     INTENT(OUT)
  !>        \param[out] "real(dp)                      :: soildenitri"

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
  !         HYPE model

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Jun 2016
  !>    Modified
  !>        X. Yang Jun 2017 enabled different time-step (hourly)

  subroutine soil_denitrification(TS, nHorizons_mHM,sat, soiltemp, soilmoist, concsoil, diss_IN, pardenirate, soildenitri)

  implicit none
  integer(i4),                    intent(in)     :: nHorizons_mHM
  integer(i4),                    intent(in)     :: TS
  real(dp), dimension(:),         intent(in)     :: sat
  real(dp),                       intent(in)     :: soiltemp
  real(dp), dimension(:),         intent(in)     :: soilmoist
  real(dp), dimension(:,:),       intent(inout)  :: concsoil
  real(dp), dimension(:),         intent(inout)  :: diss_IN
  real(dp),                       intent(in)     :: pardenirate
  real(dp),                       intent(out)    :: soildenitri
  !local
  real(dp), dimension(nHorizons_mHM)  :: lyr_deni
  real(dp)    :: fct_temp, fct_sm, fct_conc  !terms for different factors
  integer(i4) :: j
  real(dp)    :: DT
  
  
  DT = 24.0_dp / TS
  ! factors for denitrification
  fct_temp = tempfactor(soiltemp)
  do j=1, nHorizons_mHM
     if (soilmoist(j) > 0.0_dp) then
     fct_sm = deni_moistfactor(soilmoist(j),sat(j) )
     fct_conc = concfactor(concsoil(j,1))
     !denitrification
     lyr_deni(j) = pardenirate * diss_IN(j) * fct_temp * fct_sm * fct_conc / DT
	 !update dissolved inorgainc pool
     diss_IN(j) = diss_IN(j) - lyr_deni(j)
     concsoil(j,1) = diss_IN(j) / soilmoist(j)
     else
     lyr_deni(j) = 0.0_dp
     concsoil(j,1) = 0.0_dp
     diss_IN(j) = 0.0_dp
     end if


  end do
  
  soildenitri = sum(lyr_deni(:) )

    
  end subroutine soil_denitrification
  !--------------------------------------------------------------------
  ! ------------------------------------------------------------------

  !     NAME
  !         soil_plant_uptake

  !     PURPOSE
  !>        \brief plant uptake of Nitrate in soil phase. 

  !>        \details calculates nitrate uptake amount in soil water, considering potential uptake(crop dependent),
  !>         concentration, and soil moisture factors.         

  !     INTENT(IN)
  !>        \param[in] "real(dp), dimension(:)        :: wp"
  !>        \param[in] "real(dp), dimension(:)        :: horizon_depth"
  !>        \param[in] "real(dp), dimension(:)        :: soilmoist"
  !>        \param[in] "real(dp), dimension(:,:)     :: potential_uptake"
  !     INTENT(INOUT)
  !>        \param[inout] "real(dp), dimension(:,:)    :: concsoil"
  !>        \param[inout] "real(dp), dimension(:)      :: diss_IN"

  !     INTENT(OUT)
  !>        \param[out] "real(dp)                      :: soiluptkN"

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
  !         HYPE model

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Jun 2016
  
  subroutine soil_plant_uptake(wp, soilmoist, concsoil, potential_uptake, soiluptkN)
  
  implicit none

  real(dp), dimension(:),         intent(in)     :: wp
  real(dp), dimension(:),         intent(in)     :: soilmoist
  real(dp), dimension(:,:),       intent(inout)  :: concsoil
  real(dp), dimension(:),       intent(in)     :: potential_uptake
  real(dp),                       intent(out)    :: soiluptkN
  !local
  integer(i4)         :: j
  real(dp), dimension(2)    :: max_uptk, lyr_uptk
  real(dp), dimension(2)    :: diss_INt
  !only happens in upper two layers

  do j = 1, 2
     if (soilmoist(j) > 0.0001_dp) then
     diss_INt(j) = soilmoist(j) * concsoil(j,1)
     !in case max_uptk less than zero
     if (soilmoist(j) > wp(j)) then
     max_uptk(j) = (soilmoist(j)- wp(j))/ soilmoist(j)
     else
     max_uptk(j) = 0.0_dp
     end if
     lyr_uptk(j) = min(potential_uptake(j), max_uptk(j)*diss_INt(j) )
    
	 !update dissolved inorganic pool and conc. in soil water
     diss_INt(j) = diss_INt(j) - lyr_uptk(j)
     concsoil(j,1) = diss_INt(j) / soilmoist(j)
     else
     lyr_uptk(j) =0.0_dp
     concsoil(j,1) = 0.0_dp
     diss_INt(j) = 0.0_dp
     end if

  end do  
	 
  soiluptkN = sum(lyr_uptk(:))

  end subroutine soil_plant_uptake

  !--------------------------------------------------------------------
  ! ------------------------------------------------------------------

  !     NAME
  !         agri_management

  !     PURPOSE
  !>        \brief agricultural management. 

  !>        \details calculates nutrient(Nitrogen at the moment) input from agricultural management (e.g. fertilizer 
  !>        and manure application) and from residuals by havesting or ploughing. Also, calcualtes potential nutrient 
  !>        uptake during the growing season (from plant date to havest date).
  !>        Notes: catch crops has not been implemented yet....   

  !     INTENT(IN)
  !>        \param[in] "integer(i4)               :: TS"
  !>        \param[in] "integer(i4)               :: noday" 
  !>        \param[in] "integer(i4)               :: noyear"
  !>        \param[in] "integer(i4)               :: days_prev"
  !>        \param[in] "real(dp), dimension(:)    :: frac_rotation"
  !>        \param[in] "real(dp), dimension(:)    :: soilmoist"
  !>        \param[in] "real(dp)                  :: temp"    

  !     INTENT(INOUT)
  !>        \param[inout] "real(dp), dimension(:,:)   :: concsoil"
  !>        \param[inout] "real(dp), dimension(:)     :: fast_N"  
  !>        \param[inout] "real(dp), dimension(:)     :: humus_N"
  !>        \param[inout] "real(dp), dimension(:,:)   :: potential_uptake"  

  !     INTENT(OUT)
  !>        None

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
  !         HYPE model

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Jun 2016
  !>    Modified
  !>        X. Yang Jun 2017 enabled different time-step (hourly)
  !>        X. Yang Aug 2019 modified potential crop uptake calculation
  
  subroutine agri_management(TS, iBasin, noday, noyear, days_prev, frac_rotation, soilmoist, concsoil, &
       fast_N, humus_N, temp, potential_uptake, nCroptation,infrtmanapp)
  
  use mo_wqm_global_variables,   only: &
       cropdata, rotation              !&
	   
  implicit none
  integer(i4),               intent(in)    :: TS  !timestep [h]
  integer(i4),               intent(in)    :: iBasin
  integer(i4),               intent(in)    :: noday, noyear, days_prev
  real(dp), dimension(:),    intent(in)    :: frac_rotation
  real(dp), dimension(:),    intent(in)    :: soilmoist
  real(dp),                  intent(in)    :: temp                            !temperature
  real(dp), dimension(:,:),  intent(inout) :: concsoil
  real(dp), dimension(:),    intent(inout) :: fast_N, humus_N                 ![mg/m^2]
  real(dp), dimension(:),    intent(inout) :: potential_uptake  
  integer(i4),               intent(in)    :: nCroptation
  real(dp),                  intent(inout) :: infrtmanapp                     ![mg/m^2]


  !local variables
  real(dp),dimension(2,2)   :: frtman_nadd, res_nadd   ! amount of sources added to nitrogen pools [kg/km^2 or mg/m^2]
  integer(i4)       :: j, num_crp
  real(dp)          :: DT  !number of steps in one day
  integer(i4)       :: ir,jc, jc_prev, jc_next    !rotation index/crop index in a specific rotation
  integer(i4)       :: remidr, remidr_prev, remidr_next    !reminder after mod function to identify cropid in a rotation
  real(dp)          :: uptk_help, uptake_N  
  real(dp)          :: f_temp              ! temperature factor of soil uptake

  
  DT = 24.0_dp / TS
  frtman_nadd = 0.0_dp 
  res_nadd    = 0.0_dp
  potential_uptake = 0.0_dp

  do ir = 1, nCroptation
     if (frac_rotation(ir) > 0.0_dp ) then 
        ! abbreviation of crop rotation variables
        num_crp = rotation(iBasin)%ncrops(ir)
        ! identify correct crop in specific year
        remidr = mod(noyear, num_crp)
        remidr_prev = mod(noyear-1, num_crp)
        remidr_next = mod(noyear+1, num_crp)  ! the year after, **modified 2019-07-17**
        if (remidr == 0_i4)   remidr = num_crp
        if (remidr_prev == 0_i4)   remidr_prev = num_crp
        if (remidr_next == 0_i4)   remidr_next = num_crp
        jc = rotation(iBasin)%crop(ir,remidr) 
        jc_prev = rotation(iBasin)%crop(ir,remidr_prev)
        jc_next = rotation(iBasin)%crop(ir,remidr_next)  !**modified 2019-07-17**
        !***get crop_id 'jc' from rotation_info.txt and then get corresponding crop management information from cropdata.txt***
        ! fertiliser application
        if ((noday >= cropdata(jc)%frtday1 .and. noday < cropdata(jc)%frtday1 + cropdata(jc)%frtperiod) .or. &
            (noday < (cropdata(jc_prev)%frtday1 + cropdata(jc_prev)%frtperiod - days_prev))) then
        frtman_nadd(1,1) = frtman_nadd(1,1) + frac_rotation(ir) * cropdata(jc)%frtn1 * (1.0 - cropdata(jc)%frtdown1) &
                                            / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,1) = frtman_nadd(2,1) + frac_rotation(ir) * cropdata(jc)%frtn1 * cropdata(jc)%frtdown1  &
                                            / (cropdata(jc)%frtperiod * DT)
        end if
        !in case fert applied twice a year
        if (cropdata(jc)%frtday2 .ne. 0_i4) then
        if (noday >= cropdata(jc)%frtday2 .and. noday < (cropdata(jc)%frtday2 + cropdata(jc)%frtperiod)) then
        frtman_nadd(1,1) = frtman_nadd(1,1) + frac_rotation(ir) * cropdata(jc)%frtn2 * (1.0 - cropdata(jc)%frtdown2) &
                                            / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,1) = frtman_nadd(2,1) + frac_rotation(ir) * cropdata(jc)%frtn2 * cropdata(jc)%frtdown2  &
                                            / (cropdata(jc)%frtperiod * DT)
        end if
        end if
        !in case fert applied twice in previous year and the second application period across the calendar year
        if (cropdata(jc_prev)%frtday2 .ne. 0_i4) then
        if (noday < (cropdata(jc_prev)%frtday2 + cropdata(jc_prev)%frtperiod - days_prev)) then
        frtman_nadd(1,1) = frtman_nadd(1,1) + frac_rotation(ir) * cropdata(jc)%frtn2 * (1.0 - cropdata(jc)%frtdown2) &
                                            / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,1) = frtman_nadd(2,1) + frac_rotation(ir) * cropdata(jc)%frtn2 * cropdata(jc)%frtdown2  &
                                            / (cropdata(jc)%frtperiod * DT)
        end if
        end if
		!manure application
        if ((noday >= cropdata(jc)%manday1 .and. noday < cropdata(jc)%manday1 + cropdata(jc)%frtperiod) .or. &
            (noday < (cropdata(jc_prev)%manday1 + cropdata(jc_prev)%frtperiod - days_prev))) then
        frtman_nadd(1,1) = frtman_nadd(1,1) + frac_rotation(ir) * cropdata(jc)%mann1 * (1.0 - cropdata(jc)%mandown1) &
                                            * cropdata(jc)%manfIN / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,1) = frtman_nadd(2,1) + frac_rotation(ir) * cropdata(jc)%mann1 * cropdata(jc)%mandown1  &
                                            * cropdata(jc)%manfIN / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(1,2) = frtman_nadd(1,2) + frac_rotation(ir) * cropdata(jc)%mann1 * (1.0 - cropdata(jc)%mandown1) &
                                            * (1.0 - cropdata(jc)%manfIN) / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,2) = frtman_nadd(2,2) + frac_rotation(ir) * cropdata(jc)%mann1 * cropdata(jc)%mandown1  &
                                            * (1.0 - cropdata(jc)%manfIN) / (cropdata(jc)%frtperiod * DT)
        end if
		!in case manure applied two times a year
        if (cropdata(jc)%manday2 .ne. 0_i4) then
        if (noday >= cropdata(jc)%manday2 .and. noday < (cropdata(jc)%manday2 + cropdata(jc)%frtperiod)) then
        frtman_nadd(1,1) = frtman_nadd(1,1) + frac_rotation(ir) * cropdata(jc)%mann2 * (1.0 - cropdata(jc)%mandown2) &
                                            * cropdata(jc)%manfIN / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,1) = frtman_nadd(2,1) + frac_rotation(ir) * cropdata(jc)%mann2 * cropdata(jc)%mandown2  &
                                            * cropdata(jc)%manfIN / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(1,2) = frtman_nadd(1,2) + frac_rotation(ir) * cropdata(jc)%mann2 * (1.0 - cropdata(jc)%mandown2) &
                                            * (1.0 - cropdata(jc)%manfIN) / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,2) = frtman_nadd(2,2) + frac_rotation(ir) * cropdata(jc)%mann2 * cropdata(jc)%mandown2  &
                                            * (1.0 - cropdata(jc)%manfIN) / (cropdata(jc)%frtperiod * DT)
        end if
        end if
        !in case manure applied twice in previous year and the second application period across the calendar year	
        if (cropdata(jc_prev)%manday2 .ne. 0_i4) then
        if (noday < (cropdata(jc_prev)%manday2 + cropdata(jc_prev)%frtperiod - days_prev))  then
        frtman_nadd(1,1) = frtman_nadd(1,1) + frac_rotation(ir) * cropdata(jc)%mann2 * (1.0 - cropdata(jc)%mandown2) &
                                            * cropdata(jc)%manfIN / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,1) = frtman_nadd(2,1) + frac_rotation(ir) * cropdata(jc)%mann2 * cropdata(jc)%mandown2  &
                                            * cropdata(jc)%manfIN / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(1,2) = frtman_nadd(1,2) + frac_rotation(ir) * cropdata(jc)%mann2 * (1.0 - cropdata(jc)%mandown2) &
                                            * (1.0 - cropdata(jc)%manfIN) / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,2) = frtman_nadd(2,2) + frac_rotation(ir) * cropdata(jc)%mann2 * cropdata(jc)%mandown2  &
                                            * (1.0 - cropdata(jc)%manfIN) / (cropdata(jc)%frtperiod * DT)
        end if
        end if
        !residuals
        if (cropdata(jc)%resday == 0_i4) cropdata(jc)%resperiod = 365_i4
        if ((noday >= cropdata(jc)%resday .and. noday < cropdata(jc)%resday + cropdata(jc)%resperiod) .or. &
            (noday < (cropdata(jc_prev)%resday + cropdata(jc_prev)%resperiod - days_prev))) then
        res_nadd(1,1) = res_nadd(1,1) + frac_rotation(ir) * cropdata(jc)%resn * (1.0 - cropdata(jc)%resdown) &
                                      * cropdata(jc)%resfast / (cropdata(jc)%resperiod * DT)
        res_nadd(1,2) = res_nadd(1,2) + frac_rotation(ir) * cropdata(jc)%resn * (1.0 - cropdata(jc)%resdown) &
                                      * (1.0 - cropdata(jc)%resfast) / (cropdata(jc)%resperiod * DT)
        res_nadd(2,1) = res_nadd(2,1) + frac_rotation(ir) * cropdata(jc)%resn *  cropdata(jc)%resdown &
                                      * cropdata(jc)%resfast / (cropdata(jc)%resperiod * DT)
        res_nadd(2,2) = res_nadd(2,2) + frac_rotation(ir) * cropdata(jc)%resn *  cropdata(jc)%resdown &
                                      * (1.0 - cropdata(jc)%resfast) / (cropdata(jc)%resperiod * DT)
        end if
		
		

		!*************	 
        !calculate Nitrate potential uptake because of plant/crop growth 
        !*************  
        uptake_N =0.0_dp
        if (temp < 5.0_dp) then
           f_temp = 0.0_dp
        else
           f_temp = min(1.0_dp, (temp-5.0)/20.0)
        end if
        !!normally plant in Spring and harvest in Autumn
        if (cropdata(jc)%plantd <= cropdata(jc)%havestd) then    
           if (noday >= cropdata(jc)%plantd .and. noday < cropdata(jc)%havestd) then
           uptk_help = (cropdata(jc)%up1 - cropdata(jc)%up2) * exp(-cropdata(jc)%up3*(noday - cropdata(jc)%plantd))
           uptake_N = cropdata(jc)%up1 *cropdata(jc)%up2*cropdata(jc)%up3 * uptk_help /(1.0+uptk_help)/(1.0+ uptk_help)
           end if
        end if
        !Winter crop, plant in Autumn of previous year and harvest in the current year, in other words, plant date is latter than harvest date
        !in "cropdata.txt" we still use the day number of a year
        ! in this case, the emerging date (emergd) should be specified in cropdata.txt
        if (cropdata(jc)%plantd > cropdata(jc)%havestd) then
           if (noday >= cropdata(jc)%emergd .and. noday < cropdata(jc)%havestd ) then
              uptk_help = (cropdata(jc)%up1 - cropdata(jc)%up2) * exp(-cropdata(jc)%up3*(noday - cropdata(jc)%emergd))   
              uptake_N = cropdata(jc)%up1 *cropdata(jc)%up2*cropdata(jc)%up3 * uptk_help/(1.0+uptk_help)/(1.0+ uptk_help)           
           elseif (noday < cropdata(jc)%emergd) then
              uptk_help = (cropdata(jc)%up1 - cropdata(jc)%up2) * exp(-cropdata(jc)%up3*(noday + 365 - cropdata(jc)%plantd))
              uptake_N = f_temp* cropdata(jc)%up1 *cropdata(jc)%up2*cropdata(jc)%up3 * uptk_help/(1.0+uptk_help)/(1.0+ uptk_help)
           else
              uptake_N = 0.0_dp  
           end if
        end if
        ! if the next year cropping winter crops, means the next year's crop will be planted in current years' winter
        if (cropdata(jc_next)%plantd > cropdata(jc_next)%havestd) then  
           if (noday >= cropdata(jc_next)%plantd ) then
           uptk_help = (cropdata(jc_next)%up1 - cropdata(jc_next)%up2) &
                       * exp(-cropdata(jc_next)%up3*(noday - cropdata(jc_next)%plantd)-25) !the "25" is directly introduced from HYPE
           uptake_N = f_temp*cropdata(jc_next)%up1 *cropdata(jc_next)%up2*cropdata(jc_next)%up3* &
                         uptk_help /(1.0+uptk_help)/(1.0+ uptk_help)
           end if        
        end if

        potential_uptake(1) =potential_uptake(1) + frac_rotation(ir) * uptake_N * cropdata(jc)%uppsoil 
        potential_uptake(2) =potential_uptake(2) + frac_rotation(ir) * uptake_N * (1.0 - cropdata(jc)%uppsoil) 
    
     end if  !end of (if frac_rotation>0)
  end do     !end of rotation loop

  !add to corresponding pools
  do j= 1, 2  !upper two layers
  !fertilizer and manure: orgN goes to fastN pool, inorgN goes to soil water(if > 0) or fastN pool 
  if (frtman_nadd(j,1) > 0.0_dp) then
     if (soilmoist(j) > 0.0_dp) then
        call add_source_water(soilmoist(j), concsoil(j,1), frtman_nadd(j,1))
     else
        fast_N(j) = fast_N(j) + frtman_nadd(j,1)
     end if 
  end if
  if (frtman_nadd(j,2) > 0.0_dp) then
     fast_N(j) = fast_N(j) + frtman_nadd(j,2)
  end if
  !residuals: go to humusN and fastN pools
  if (sum(res_nadd(j,:)) > 0.0_dp ) then
     fast_N(j) = fast_N(j) + res_nadd(j,1)
     humus_N(j)= humus_N(j) + res_nadd(j,2)
  end if
  !potential uptake
  if (potential_uptake(j) > 0.0_dp) then
     potential_uptake(j) =1000.0_dp * potential_uptake(j) / DT    !g/m^2/d --> mg/m^2/timestep 
  end if  

  end do 
  infrtmanapp =  sum(frtman_nadd(:,1))
 
  end subroutine agri_management   
  !-------------------------------------------------------------------
  ! ------------------------------------------------------------------

  !     NAME
  !         soil_nutrient_transformation

  !     PURPOSE
  !>        \brief Nutrient transformation between different pools in soil phase. 

  !>        \details calculates the transformation between four different pools in every soil layer, including 
  !>        degradation, minearalisation and dissolution.         

  !     INTENT(IN)
  !>        \param[in] "integer(i4)                   :: TS"
  !>        \param[in] "integer(i4)                   :: nhorizons"     !
  !>        \param[in] "real(dp), dimension(:)        :: horizon_depth" !
  !>        \param[in] "real(dp)                      :: airtemp"       !
  !>        \param[in] "real(dp)                      :: snowdp"        !
  !>        \param[in] "real(dp), dimension(:)        :: wp"            !
  !>        \param[in] "real(dp), dimension(:)        :: sat"           !
  !>        \param[in] "real(dp), dimension(:)        :: soilmoist"     !
  !>        \param[in] "real(dp)                      :: pardegradN"    !
  !>        \param[in] "real(dp)                      :: parmineraN"    !
  !>        \param[in] "real(dp)                      :: pardissolN"    !

  !     INTENT(INOUT)
  !>        \param[inout] "real(dp)                      :: soiltemp"      !
  !>        \param[inout] "real(dp), dimension(:,:)      :: concsoil"      !
  !>        \param[inout] "real(dp), dimension(:)        :: fast_N"        !
  !>        \param[inout] "real(dp), dimension(:)        :: humus_N"       !
  !>        \param[inout] "real(dp), dimension(:)        :: diss_IN"       !
  !>        \param[inout] "real(dp), dimension(:)        :: diss_ON"       !

  !     INTENT(OUT)
  !>        None

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
  !         HYPE model

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Jun 2016
  !>    Modified
  !>        X. Yang Jun 2017 enabled different time-step (hourly)
  
  subroutine soil_nutrient_transformation(TS, nhorizons, horizon_depth, airtemp, soiltemp, snowdp, wp, sat,  soilmoist,  &
        concsoil, fast_N, humus_N, diss_IN, diss_ON, pardegradN, parmineraN, pardissolN, mineralN_soil)


    implicit none
    integer(i4),                   intent(in)   :: TS
    integer(i4),                   intent(in)   :: nhorizons     !
    real(dp), dimension(:),        intent(in)   :: horizon_depth ![mm]
    real(dp),                      intent(in)   :: airtemp       !
    real(dp),                      intent(inout):: soiltemp      !
    real(dp),                      intent(in)   :: snowdp        ![mm]
    real(dp), dimension(:),        intent(in)   :: wp            ![mm]
    real(dp), dimension(:),        intent(in)   :: sat           ![mm]
    real(dp), dimension(:),        intent(in)   :: soilmoist     ![mm]
    real(dp), dimension(:,:),      intent(inout):: concsoil      ![mg/l]
    real(dp), dimension(:),        intent(inout):: fast_N        ![mg/m^2]
    real(dp), dimension(:),        intent(inout):: humus_N       ![mg/m^2]
    real(dp), dimension(:),        intent(inout):: diss_IN       ![mg/m^2]
    real(dp), dimension(:),        intent(inout):: diss_ON       ![mg/m^2]
    real(dp),                      intent(in)   :: pardegradN    ![/d]
    real(dp),                      intent(in)   :: parmineraN    ![/d]
    real(dp),                      intent(in)   :: pardissolN    ![/d]
    real(dp),                      intent(out)  :: mineralN_soil ![mg/m^2]
    !local variables
    real(dp)    :: DT
    real(dp)    :: fct_temp     ! temperature factor
    real(dp)    :: fct_sm       ! soil moisture factor
    real(dp)    :: degradN, mineraN, mineraN2,newconc !transformated amount in each process
    integer(i4) :: j
    real(dp), dimension(nHorizons)   :: mineralN_lyr
	
	
    DT = 24.0_dp / TS
    !calcualte soil temperature
    call calculate_soil_temperature(airtemp, soiltemp, snowdp)
    !temperature factor
    fct_temp = tempfactor(soiltemp) 
    do j =1, nhorizons
       degradN = 0.0_dp
       mineraN = 0.0_dp 
       mineraN2 = 0.0_dp
       newconc = 0.0_dp

    !soil moisture factor
       if (j == 1) then
          fct_sm = moistfactor(soilmoist(j), wp(j), sat(j), horizon_depth(j)) 
       else
          fct_sm = moistfactor(soilmoist(j), wp(j), sat(j), horizon_depth(j)- horizon_depth(j-1))
       end if
    !degradation: from humusN pool to fastN pool
       degradN =  pardegradN * fct_temp * fct_sm * humus_N(j) / DT 
       humus_N(j) = humus_N(j) - degradN
       fast_N(j) = fast_N(j) + degradN    
    !mineralisation: from fastN pool to dissolved IN
       mineraN =  parmineraN * fct_temp * fct_sm * fast_N(j) / DT
       fast_N(j)  = fast_N(j) - mineraN
       diss_IN(j) = diss_IN(j) + mineraN
	!update concentration of IN and ON
       concsoil(j,1) = diss_IN(j) / soilmoist(j)
       concsoil(j,2) = diss_ON(j) / soilmoist(j)
    !mineralisation: from dissolved ON to dissolved IN
       mineraN2 =  parmineraN * fct_temp * fct_sm * concsoil(j,2) / DT 
       concsoil(j,1) = concsoil(j,1) + mineraN2
       concsoil(j,2) = concsoil(j,2) - mineraN2
    !total mineralistaion amount
       mineralN_lyr(j) = mineraN + mineraN2*soilmoist(j) 
	!update dissolved IN and ON pools
       diss_IN(j) = concsoil(j,1) * soilmoist(j)
       diss_ON(j) = concsoil(j,2) * soilmoist(j)
    !dissolution: from fastN pool to dissolved ON
       fast_N(j) = fast_N(j) + diss_ON(j)
       if (j == 1) then
          newconc = fast_N(j) / (soilmoist(j) + pardissolN * horizon_depth(j)/1000.0_dp)
       else
          newconc = fast_N(j) / (soilmoist(j) + pardissolN * (horizon_depth(j)-horizon_depth(j-1))/1000.0_dp)
       end if
       concsoil(j,2) = concsoil(j,2) + (newconc - concsoil(j,2))*(1.0_dp - exp(-0.01_dp / DT))  !1% perday
    !update fastN and dissolved ON pools
       diss_ON(j) = concsoil(j,2) * soilmoist(j)
       fast_N(j) = fast_N(j) - diss_ON(j)
       
    end do
       mineralN_soil = sum(mineralN_lyr(:))

  end subroutine soil_nutrient_transformation

  !-------------------------------------------------------------------
  ! ------------------------------------------------------------------

  !     NAME
  !         totalrunoff_concentration

  !     PURPOSE
  !>        \brief totalrunoff_concentration. 

  !>        \details   

  !     INTENT(IN)
  !>        \param[in] 

  !     INTENT(INOUT)
  !>        \param[inout] 

  !     INTENT(OUT)
  !>        None

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
  !         HYPE model

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Jun 2016
  
  subroutine totalrunoff_concentration(frac_sealed,runoff_sealed, fastRunoff, slowRunoff, baseflow,  &
         total_runoff, crunoff_sealed, cfastrunoff, cslowrunoff, cbaseflow,   &
         ctotal_runoff)

    implicit none
    real(dp),                         intent(in) :: frac_sealed
    real(dp),                         intent(in) :: runoff_sealed
    real(dp),                         intent(in) :: fastRunoff
    real(dp),                         intent(in) :: slowRunoff
    real(dp),                         intent(in) :: baseflow
    real(dp),                         intent(in) :: total_runoff
    real(dp), dimension(:),           intent(in) :: crunoff_sealed
    real(dp), dimension(:),           intent(in) :: cfastrunoff    
    real(dp), dimension(:),           intent(in) :: cslowrunoff
    real(dp), dimension(:),           intent(in) :: cbaseflow
    real(dp), dimension(:),        intent(inout) :: ctotal_runoff
    !local
	
    ctotal_runoff(:) = (runoff_sealed * crunoff_sealed(:) + (fastRunoff * cfastrunoff(:) + &
                        slowRunoff * cslowrunoff(:) + baseflow * cbaseflow(:)) ) / total_runoff        
!    ctotal_runoff(:) = (frac_sealed * runoff_sealed * crunoff_sealed(:) + (1-frac_sealed)*(fastRunoff * cfastrunoff(:) + &
!                        slowRunoff * cslowrunoff(:) + baseflow * cbaseflow(:)) ) / total_runoff        

  end subroutine totalrunoff_concentration    
  !-------------------------------------------------------------------
  ! ------------------------------------------------------------------

  !     NAME
  !         runoff_concentration

  !     PURPOSE
  !>        \brief runoff_concentration. 

  !>        \details   

  !     INTENT(IN)
  !>        \param[in] 

  !     INTENT(INOUT)
  !>        \param[inout] 

  !     INTENT(OUT)
  !>        None

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
  !         HYPE model
  !         INCA model

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Jun 2016
 
  subroutine runoff_concentration(infil_third, cinfil_third,  &       
         prevstep_unsatStorage, pre_percol, pre_baseflow, basef_avg, cunsatSTW, cpercol, &
         cslowrunoff, cbaseflow,satStorage)

    implicit none
	

    real(dp),                         intent(in) :: infil_third
    real(dp), dimension(:),           intent(in) :: cinfil_third
    real(dp),                         intent(in) :: prevstep_unsatStorage
    real(dp),                         intent(in) :: satStorage
    real(dp),                         intent(in) :: pre_percol
    real(dp),                         intent(in) :: pre_baseflow
    real(dp),                         intent(inout) :: basef_avg
    real(dp), dimension(:),           intent(inout) :: cunsatSTW
    real(dp), dimension(:),           intent(inout) :: cpercol
    real(dp), dimension(:),           intent(inout) :: cslowrunoff
!    real(dp), dimension(:),           intent(inout) :: cfastrunoff
    real(dp), dimension(:),           intent(inout) :: cbaseflow
    !real(dp),                         intent(in) :: gwresidt_par    
                                               ! parameter:ground water resident time for gwconc. 
    !local
    !variables for the numerical solving method
    real(dp), dimension(size(cbaseflow))    :: k1, k2, k3, k4
    real(dp), dimension(size(cbaseflow))    :: y1, y2, y3
    										   
	!unsaturated zone
    if (infil_third > 0.0_dp) then
       call mix_conc(prevstep_unsatStorage, cunsatSTW(:), infil_third, cinfil_third(:))
    end if
    cslowrunoff(:) = cunsatSTW(:)

	!saturated zone
	!equation introduced from INCA model
    !Also using the fourth-order Runge-Kutta technique to solve the differential equation
    !700 times of mHM saturated storage is assigned as the size of retention storage
    basef_avg = basef_avg + (pre_baseflow -basef_avg) / 30.0_dp
    k1(:) = (cpercol(:)* pre_percol - pre_baseflow * cbaseflow(:)) / (satStorage*700.0_dp)!(1000.0_dp * gwresidt_par * basef_avg)
    y1(:) = cbaseflow(:) + k1(:) * 0.5_dp
    k2(:) = (cpercol(:)* pre_percol - pre_baseflow * y1(:)) / (satStorage*700.0_dp)!(1000.0_dp * gwresidt_par * basef_avg)	
    y2(:) = cbaseflow(:) + k2(:) * 0.5_dp
    k3(:) = (cpercol(:)* pre_percol - pre_baseflow * y2(:)) / (satStorage*700.0_dp)!(1000.0_dp * gwresidt_par * basef_avg)
    y3(:) = cbaseflow(:) + k2(:)
    k4(:) = (cpercol(:)* pre_percol - pre_baseflow * y3(:)) / (satStorage*700.0_dp)!(1000.0_dp * gwresidt_par * basef_avg)
    cbaseflow(:) = cbaseflow(:) + (k1(:)+2*k2(:)+2*k3(:)+k4(:))/6.0_dp
    cpercol(:) = cunsatSTW(:)
  end subroutine runoff_concentration

  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

  !     NAME
  !         soil_moisture_concentration

  !     PURPOSE
  !>        \brief soil_moisture_concentration. 

  !>        \details   

  !     INTENT(IN)
  !>        \param[in] 

  !     INTENT(INOUT)
  !>        \param[inout] 

  !     INTENT(OUT)
  !>        None

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
  !         HYPE model
  !         INCA model

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Jun 2016
  
  subroutine soil_moisture_concentration(nsubstances, frac_sealed, prec_effect, &
      sealedStorage, infiltration, soilMoisture, aet_soil, aet_sealed, prevstep_sealedStorage, &
      prevstep_soilMoisture, csealSTW, cinfiltration, concsoil, &
      cprec_effect, crunoff_sealed)

    implicit none
    integer(i4),                      intent(in) :: nsubstances
    real(dp),                         intent(in) :: frac_sealed
    real(dp),                         intent(in) :: prec_effect
    real(dp),                         intent(in) :: sealedStorage
    real(dp), dimension(:),           intent(in) :: infiltration
    real(dp), dimension(:),           intent(in) :: soilMoisture
    real(dp), dimension(:),           intent(in) :: aet_soil
    real(dp),                         intent(in) :: aet_sealed
    real(dp),                         intent(in) :: prevstep_sealedStorage
    real(dp), dimension(:),           intent(in) :: prevstep_soilMoisture
    real(dp), dimension(:),           intent(inout) :: csealSTW
    real(dp), dimension(:,:),         intent(inout) :: cinfiltration
    real(dp), dimension(:,:),         intent(inout) :: concsoil
    real(dp), dimension(:),           intent(inout) :: cprec_effect
    real(dp), dimension(:),           intent(inout) :: crunoff_sealed
 
    !local
    real(dp)                                      ::tmp    !temporal var for water
    real(dp), dimension(nsubstances)              ::ctmp   !temporal var for conc.

    integer(i4)                                   ::hh

    !----    
    !---conc. in sealed storage ***SHOULD BE IMPROVED***
    if (frac_sealed > 0.0_dp) then
       tmp = prevstep_sealedStorage + prec_effect * frac_sealed
       if ( tmp > 0.0_dp) then
       ctmp(:) = (csealSTW(:) * prevstep_sealedStorage + prec_effect * frac_sealed * cprec_effect(:)) / tmp
       end if
       crunoff_sealed(:) = 1.0_dp !ctmp(:) !***SHOULD BE IMPROVED LATER***
       if (sealedStorage <= 0.1_dp) then
          csealSTW(:) = 1.0_dp
       else
          csealSTW(:) = ctmp(:) * (sealedStorage + aet_sealed) / sealedStorage 
       end if
    end if

    !every step initialise conc. in infiled water of lower two layers
    !first layer infiltration is always updated by effect rainfall 
    cinfiltration(:,:) = 0.0_dp
    !---conc. in soil layers
    do hh =1, size(concsoil, 1)
     
       concsoil(hh,:) = concsoil(hh,:) * (soilMoisture(hh) + aet_soil(hh) + infiltration(hh)) / (soilMoisture(hh)+ infiltration(hh))  

       if (hh == 1) then
          call mix_conc(prevstep_soilMoisture(hh), concsoil(hh,:), prec_effect*(1-frac_sealed), cprec_effect(:))
       else
          call mix_conc(prevstep_soilMoisture(hh), concsoil(hh,:), infiltration(hh-1), cinfiltration(hh-1,:))
       end if
      


       !---conc. after soil ET

       !to avoid numerical error
       if (soilMoisture(hh) < 0.0001_dp) then
          concsoil(hh,:) = 0.0_dp
       end if
       cinfiltration(hh,:) = concsoil(hh,:)

    end do

  
  end subroutine soil_moisture_concentration

  !--------------------------------------------------------------------------  
  subroutine mix_conc(sm,csm,infil,cinfil)

    implicit none
    real(dp),               intent(in)      ::sm    ! soil moisture
    real(dp), dimension(:), intent(inout)   ::csm   ! conc. in the soil layer
    real(dp),               intent(in)      ::infil ! input water from the upper soil layer
    real(dp), dimension(:), intent(in)     ::cinfil ! conc. in the input soil water from upper layer

    csm(:) = (sm * csm(:) + infil * cinfil(:)) / (sm + infil)

  end subroutine mix_conc
  !----------------------------------------------------------------------------
  subroutine add_source_water(vol, conc, source)
     implicit none
     real(dp),              intent(in)     :: vol          ![mm]water volumn
     real(dp),              intent(inout)  :: conc         ![mg/l]concentration in water
     real(dp),              intent(in)     :: source       ![mg/m^2]amount of source  
     !local
     !vol*conc-->[mm]*[mg/l]=[mg/m^2]
     !to aviod numerical error
     if (vol > 0.0001_dp) then
        conc = (conc * vol + source) / vol
     end if
     
  end subroutine add_source_water
  !-----------------------------------------------------------------------------
  subroutine production_pool(pool, source)
  
    implicit none
    real(dp),              intent(inout)  :: pool
    real(dp),              intent(in)     :: source

    pool = pool + source
  
  end subroutine production_pool  
  !------------------------------------------------------------------------------
  subroutine calculate_soil_temperature(airt, soilt, snowdp)  
    implicit none
    real(dp),              intent(in)   ::airt     ! air temperature
    real(dp),              intent(inout):: soilt   ! soil temperature
    real(dp),              intent(in)   :: snowdp  ! snow depth [mm]
    !local
    !parameter variables, constant
    real(dp), parameter  :: temp_deep = 5.0_dp     ! temperature of deep soil
    real(dp), parameter  :: soil_mem = 30.0_dp     ! soil memory [days]
    real(dp), parameter  :: sp_frost = 10.0_dp     ! -----
    real(dp), parameter  :: w_deep = 0.001_dp  ! weight of deep soil temperature

    real(dp)  :: w_air  ! weight of air temperature

    w_air = 1.0_dp / (soil_mem + sp_frost * snowdp)
    soilt = soilt * (1.0_dp-w_air-w_deep) + airt * w_air + temp_deep * w_deep
	
  end subroutine calculate_soil_temperature
  !-------------------------------------------------------------------------------
  real(kind=dp) function tempfactor(c_temp)
  
  implicit none
  real(dp),               intent(in)  :: c_temp     ! temperature
  real(dp)                            :: f_temp     ! function value
  
  f_temp = 2**((c_temp - 20.0_dp) / 10.0_dp)
  if (c_temp < 5.0_dp) f_temp = f_temp * (c_temp / 5.0_dp)
  if (c_temp < 0.0_dp) f_temp = 0.0_dp
  tempfactor = f_temp
    
 
  end function tempfactor
  !-------------------------------------------------------------------------------
  real(kind=dp) function moistfactor(sm, wp, sat, thickness)
  implicit none
  !values of a sepcific soil layer
  real(dp),       intent(in) :: sm   !soil moisture
  real(dp),       intent(in) :: wp          !wilting point
  real(dp),       intent(in) :: sat         !saturated water content
  real(dp),       intent(in) :: thickness       ! soil thickness
  !local
  ! parameter variables constant
  real(dp), parameter        :: smf_satact = 0.6_dp  !
  real(dp), parameter        :: smf_upp    = 0.12_dp  !  
  real(dp), parameter        :: smf_low    = 0.08_dp  !  
  real(dp), parameter        :: smf_pow    = 1.0_dp  !

  real(dp)    :: f_sm  
  if (sm >= sat) then
     f_sm = smf_satact  
  elseif (sm <= wp) then
     f_sm = 0.0_dp
  else
     f_sm = min(1.0_dp, (1.0_dp-smf_satact)*((sat-sm)/ ((smf_upp/100.0_dp)*thickness))**smf_pow + smf_satact)
     f_sm = min(f_sm,((sm-wp)/((smf_low/100.0_dp)*thickness))**smf_pow)
  end if

  moistfactor = f_sm

  end function moistfactor
  !-----------------------------------------------------------------------------------
  real(kind =dp) function deni_moistfactor(sm, sat)
  implicit none 
  real(dp),          intent(in)  :: sm
  real(dp),          intent(in)  :: sat
  !constant parameter  
  real(dp), parameter         :: fsm_denilimit = 0.3_dp
  real(dp), parameter         :: fsm_denipow   = 2.5_dp
  !local
  real(dp)  :: f_denism
  
  f_denism = 0.0_dp
  if (sm/sat > fsm_denilimit) then
     f_denism = (((sm/sat) - fsm_denilimit)/ (1.0-fsm_denilimit)) ** fsm_denipow  
  end if
  deni_moistfactor = f_denism

  end function deni_moistfactor
  
  !-----------------------------------------------------------------------------------
  real(kind=dp) function concfactor(conc)
  implicit none
  real(dp),      intent(in) :: conc
  !constant parameter
  real(dp), parameter         :: halfsatINsoil = 1.0_dp
  
  concfactor = conc / (conc + halfsatINsoil)

  end function concfactor  

END MODULE mo_water_quality
