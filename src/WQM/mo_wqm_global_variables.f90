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

!> \file mo_wqm_global_variables.f90

!> \brief Global variables used in reading, writing and startup for water quality modelling.

!> \details

!> \authors Xiaoqiang Yang
!> \date Jun 2016
!>  Modified   X Yang, May 2018   added global variables for riparian zone shading effect


MODULE mo_wqm_global_variables

  
  USE mo_kind,             ONLY: i4, i8, dp
  use mo_common_variables, ONLY: period
  USE mo_mhm_constants,    ONLY: nOutFlxState, YearMonths, maxNoBasins, maxNLCovers
!  use mo_global_variables, only: nLAIclass, nBasins

  IMPLICIT NONE
  
  public :: rotationType
  public :: cropType
  public :: wqm_basinInfo
  


  
  !read-in input files for WQM
    character(256), dimension(:), allocatable, public :: dirWaterquality     ! Directory where WQ input files are located
    ! number of substances involved, currently only IN and ON
    integer(i4), parameter, public                    :: nsubstances = 2     
    integer(i4), dimension(:), allocatable, public    :: nCroptation         ! Number of crop rotation types  "nBains"
    integer(i4), public                               :: num_crops           ! number of total crop types
    !crop data and land-use fraction
    integer(i4), dimension(:), allocatable, public    :: L0_cover_rotation   !crop rotation id at level L0
    real(dp), dimension(:,:), allocatable, public     :: L1_frotation        ! area fraction at Level L1
    real(dp), dimension(:,:), allocatable, public     :: L1_fLAI             ! area fraction at Level L1
    real(dp), dimension(:,:), allocatable, public     :: L1_fsoil            ! area fraction at Level L1
    real(dp), dimension(:,:), allocatable, public     :: L1_fgeounit          ! area fraction at Level L1
    !(for crop rotation calculation)
    integer(i4), public                               :: year_start          ! simulation starting year
    integer(i4), public                               :: day_preyr           ! number of days in the previous year.

	!values of land-use denpent parameters and initial values of different pools
	!dim=nLAIclass
    real(dp), dimension(:), allocatable,  public   :: degradN_LAI    ! degradation rate (land-use factor)
    real(dp), dimension(:), allocatable,  public   :: mineraN_LAI    ! mineralisation rate (land use factor)
    real(dp), dimension(:), allocatable,  public   :: dissolN_LAI    ! dissolution rate (land use factor)
    real(dp), dimension(:), allocatable,  public   :: init_concIN    ! initial IN concentration in soil water
    real(dp), dimension(:), allocatable,  public   :: init_concON    ! initial ON concentration in soil water
    real(dp), dimension(:), allocatable,  public   :: init_humusN    ! amount of pool
    real(dp), dimension(:), allocatable,  public   :: init_fastN     ! amount of pool
    !real(dp), dimension(:), allocatable,  public   :: upsoil_part    ! fraction of plant uptake in first layer 
    real(dp), dimension(:), allocatable,  public   :: hnhalf         ! half depthe for humusN pool
    integer(i4), dimension(:), allocatable, public    :: Geoform     ! Formation 1 - permeable sedimentary materials
                                                                     !           2 - impermeable bedrock
    !variables as global arguements
    !state variables	
    real(dp), dimension(:,:),   allocatable, public :: L1_humusN         ! humus Nitrate pool in each soil layer(organic form)
    real(dp), dimension(:,:),   allocatable, public :: L1_fastN          ! fast Nitrate pool in each soil layer(orgainc form)
    real(dp), dimension(:,:),   allocatable, public :: L1_dissolvedIN    ! dissolved inorganic pool in each soil layer
    real(dp), dimension(:,:),   allocatable, public :: L1_dissolvedON    ! dissolved organic pool in each soil layer
    real(dp), dimension(:,:,:), allocatable, public :: L1_csoilMoist     ! conc. in soil moisture
    real(dp), dimension(:,:),   allocatable, public :: L1_csealSTW       ! conc. in sealed water storage
    real(dp), dimension(:,:),   allocatable, public :: L1_cunsatSTW      ! conc. in unsaturated water storage
    real(dp), dimension(:,:),   allocatable, public :: L1_csatSTW        ! conc. in saturated water storage
    real(dp), dimension(:),     allocatable, public :: L1_soiltemp       ! soil temperature (calculated from air temperature)
    real(dp), dimension(:),     allocatable, public :: L1_baseflow_avg   ! mean baseflow for calculating baseflow conc.
                                                                         !(eq. introduced from INCA model) 
    real(dp), dimension(:,:),   allocatable, public :: L1_cbaseflow_delta! changes of baseflow conc. in each time step
    !matter fluxes
    real(dp), dimension(:,:,:), allocatable, public :: L1_cinfilSoil     ! conc. in infiltrated soil water
                                                                         !  (between each soillayer)
    real(dp), dimension(:,:),   allocatable, public :: L1_cpreEffect     ! conc. in effective precipitation
    real(dp), dimension(:,:),   allocatable, public :: L1_crain          ! conc. in rainfall (wet atmospheric deposition)(dim1=ncells, dim2=nsubstances)
    real(dp), dimension(:,:),   allocatable, public :: L1_cpercolate     ! conc. in percolated water from unsat storage 
                                                                         !  to satstorage
    real(dp), dimension(:,:),   allocatable, public :: L1_crunoffSeal    ! conc. in direct runoff
    real(dp), dimension(:,:),   allocatable, public :: L1_cfastRunoff    ! conc. in fast runoff (HERE mHM named fast_interflow)
    real(dp), dimension(:,:),   allocatable, public :: L1_cslowRunoff    ! conc. in slow runoff (HERE mHM named slow_interflow)
    real(dp), dimension(:,:),   allocatable, public :: L1_cbaseflow      ! conc. in baseflow
    real(dp), dimension(:,:),   allocatable, public :: L1_ctotal_runoff  ! conc. in total runoff 
                                                                         !   (mixture from all runoff component)
    real(dp), dimension(:),     allocatable, public :: L1_soilUptakeN    ! soil N uptake amount
    real(dp), dimension(:),     allocatable, public :: L1_soilDenitri    ! soil denitrification amount
    real(dp), dimension(:),     allocatable, public :: L1_soilMineralN    ! soil N mineralisation amount
    real(dp), dimension(:),     allocatable, public :: L1_soilINfrtmanapp ! soil IN surpplied amount
    !nitrate parameters
    real(dp), dimension(:),     allocatable, public :: L1_rdegradN       ! degradation rate from humusN pool to fastN pool
    real(dp), dimension(:),     allocatable, public :: L1_rmineralN      ! mineralisation rate from fastN pool 
                                                                         ! to dissolved inorganic pool
    real(dp), dimension(:),     allocatable, public :: L1_rdissolN       ! dissolution rate from fastN pool 
                                                                         ! to dissolved organic pool
    real(dp), dimension(:),     allocatable, public :: L1_rdeniSoil      ! IN denitrification rate in soil water
    real(dp), dimension(:),     allocatable, public :: L1_gwresidT       ! groundwater resident time baseflow conc(INCA EQ.)
   
  !state variables to store last step values (previous step) for soil water conc. 
    real(dp), dimension(:),     allocatable, public :: prevstep_sealedStorage    ! sealed storage in previous step    
    real(dp), dimension(:),     allocatable, public :: prevstep_unsatStorage    ! unsaturated storage in previous step    
    real(dp), dimension(:),     allocatable, public :: prevstep_satStorage    ! saturated storage in previous step    
    real(dp), dimension(:,:),   allocatable, public :: prevstep_soilMoisture    ! soil moisture in pervious step  
    real(dp), dimension(:),     allocatable, public :: prevstep_percol          ! percolated water in previous step	
    real(dp), dimension(:),     allocatable, public :: prevstep_baseflow          ! baseflow in previous step

  type rotationType
     ! input data
     integer(i4), dimension(:), allocatable        :: id        !         rotation Id
     integer(i4), dimension(:), allocatable        :: ncrops    !Number of crops in specific rotation type, maximum 10 types
     integer(i4), dimension(:,:), allocatable      :: crop      ! crop_ids in specific rotation type
     !                                                                !            
  end type rotationType
  type(rotationType), dimension(:),allocatable, public :: rotation  ! rotation infos

  type cropType
     ! input data
     character(256)        :: cropname  ! name decription of crop
     integer(i4)           :: cropid    ! crop id
     ! data of input sources     
     real(dp)              :: frtn1       ! amount of nitrogen in first fertilizer
     integer(i4)           :: frtday1     ! number of days in a year that firstly apply fertilizer
     real(dp)              :: frtdown1    ! fraction of N that goes down to the second layers
     real(dp)              :: frtn2       ! amount of nitrogen in second fertilizer
     integer(i4)           :: frtday2     ! number of days in a year that secondly apply fertilizer
     real(dp)              :: frtdown2    ! fraction of N that goes down to the second layers
     real(dp)              :: mann1       ! amount of nitrogen in first manure
     integer(i4)           :: manday1     ! number of days in a year that firstly apply manure
     real(dp)              :: mandown1    ! fraction of N that goes down to the second layers
     real(dp)              :: mann2       ! amount of nitrogen in second manure
     integer(i4)           :: manday2     ! number of days in a year that secondly apply manure
     real(dp)              :: mandown2    ! fraction of N that goes down to the second layers
     real(dp)              :: manfIN       ! fraction of inorganic N in applied mannure
     integer(i4)           :: frtperiod   ! number of days to spread fert and man
     real(dp)              :: resn      ! amount of nitrogen in residuals
     integer(i4)           :: resday    ! number of days in a year that residuals appear
     real(dp)              :: resdown   ! fraction of N in residuals that goes down to the second layers
     real(dp)              :: resfast   ! fraction of N in residuals that goes into fast N pool
     integer(i4)           :: resperiod ! number of days to spread residuals
     ! uptake parameters
     real(dp)              :: up1
     real(dp)              :: up2
     real(dp)              :: up3
     real(dp)              :: uppsoil
     ! farming date

     integer(i4)           :: plantd    ! planting date (sowning date)
     integer(i4)           :: emergd   ! emerging date for winter crops
     integer(i4)           :: havestd   ! havesting date
     ! catch crop
     integer(i4)           :: ccrop     ! 1= including catch crop; 0= non catch crop
     integer(i4)           :: ccplantd  ! planting date of catch crop
     integer(i4)           :: cchavestd ! havesting date of catch crop (can be more than 365 if it's winter catch crop 
	                                    ! and havest in the spring of next year)
  end type cropType
  type(cropType), dimension(:), allocatable, public      :: cropdata    ! crop managements 
	! -------------------------------------------------------------------

  !VARIABLES FOR IN-STREAM PROCESSES
  real(dp), dimension(:),       allocatable, public      :: L1_reachtemp     ! regional reach water temperature at L1 level
                                                                             !(20day's average of air temperature)
  real(dp), dimension(:),       allocatable, public      :: L11_rivertemp    ! water temperature of each reach(link)
  real(dp), dimension(:),       allocatable, public      :: L11_riverbox     ! amount of water stored in each reach
  real(dp), dimension(:,:),     allocatable, public      :: L11_criverbox    ! conc. of each riverbox
  real(dp), dimension(:),       allocatable, public      :: L11_yravg_q      ! yearly average discharge of each reach
  real(dp), dimension(:,:),     allocatable, public      :: L11_concOUT      ! conc.of runoff from L11 cells,
                                                                             !   corresponding to "L11_qOUT"
  real(dp), dimension(:,:),     allocatable, public      :: L11_interload    ! total load of inflow from upstream reach(es)
  real(dp), dimension(:,:),     allocatable, public      :: L11_concTIN      ! conc. of total inflow in each reach node,
                                                                             !   corresponding to "L11_qTIN"
  real(dp), dimension(:,:),     allocatable, public      :: L11_concMod     ! modelled conc. of each node,
                                                                             !   corresponding to "L11_qMod"
  !model parameters in parameter namelist
  real(dp), dimension(:),       allocatable, public      :: L11_rdeniAqtc    ! denitrification rate in aquatic system(river)
  real(dp), dimension(:),       allocatable, public      :: L11_ratuptkN     ! autotropic N uptake mgNm-2d-1
  real(dp), dimension(:),       allocatable, public      :: L11_rpprodN      ! assimilatory rate in aquatic system
  ! state variables for uptake amount
  real(dp), dimension(:),       allocatable, public      :: L11_aquaticDenitri  ! N denitrified in stream water 
  real(dp), dimension(:),       allocatable, public      :: L11_aquaticAssimil  !N uptake (assimilatory) amount in stream water
  ! global radiation for in-stream Primarily production  
  logical, public                                        :: GR_file_exist = .FALSE.  ! if the "global_radiation.txt" presents or not
  real(dp), dimension(:,:),     allocatable, public      :: global_radiation 
                                                           ! measured global radiation data, dim1=nbasin,dim2=nosimulatedays
  real(dp), dimension(:,:),     allocatable, public      :: nor_globalradi   ! normalised global radiation
  real(dp), dimension(:,:),     allocatable, public      :: norlai_daily     ! normalised daily lai, dim1=nLAIclass, dim2=days_yr
 
  real(dp), dimension(:,:),     allocatable, public      :: L11_fLAI         ! area fraction of each land use type at L11 level

  real(dp), dimension(:),   allocatable, public      :: L11_rzcoeff  !riparian zone shading coefficient  
                                                           ! dim = number of reaches 
  real(dp), dimension(:),   allocatable, public      :: L11_flight   ! moving average of L11_recoeff, for GPP calculation

  !state variables for calculating
  real(dp), dimension(:),  allocatable, public  :: L11_rivert_avg10 ! 10-day's average temperature of river water in each reach 
  real(dp), dimension(:),  allocatable, public  :: L11_rivert_avg20 ! 20-day's average temperature of river water in each reach    

  
  !VARIABLES FOR CONFIGURE INFORMATION READ-IN
  !dim = number of basins
  ! directory of input files for water quality model
  character(256), dimension(:),  allocatable, public     :: dirInputWQM      
  ! directory of measured water quality data (evaluation, upstream and point source inflow) 
  character(256), dimension(:),  allocatable, public     :: dirGaugesWQM     
  ! 
  type wqm_basinInfo
     !dim = number of total stations
     character(256), dimension(:), allocatable  :: fname        ! file name of measured wq data at evaluation gauging station
     ! file name of additional inflow wq data(upstream and point source)   
     character(256), dimension(:), allocatable  :: Inflowfname  
     !dim1= number of observation, dim2=number of total gauges, dim3=number of measured data type
     real(dp), dimension(:,:,:), allocatable    :: GaugeConc        ! measured wq data at evaluation gauging station
     real(dp), dimension(:,:,:), allocatable    :: InflowGaugeConc  ! measured wq data at additional inflow station
  end type wqm_basinInfo  
  type(wqm_basinInfo), public                            :: basin_wqm        ! basin info for wqm
  
  integer(i4),dimension(:), allocatable, public  :: nEvalCmeasPerday     ! number of measured conc. data per day in files
  integer(i4),dimension(:), allocatable, public  :: nAddinCmeasPerday     ! number of measured conc. data per day in files
  integer(i4),dimension(:), allocatable, public  :: numEvalCol       ! number of wq data columns in evaluation data file  
  integer(i4),dimension(:), allocatable, public  :: numAddinCol      ! number of wq data columns in additional inflow data file

  !constant: maximum columns of wq data in gauging station
  ! currently, for nitrogen submodel only IN, ON, TN are available
  integer(i4), parameter, public                         :: maxcols = 3_i4

  character(256), dimension(:,:), allocatable, public  :: evalHead_str     ! name of data measured at evaluation station
  character(256), dimension(:,:), allocatable, public  :: inflowHead_str   ! name of data measured at inflow station   
  
  !for model write out
  real(dp), dimension(:,:,:), allocatable                :: WQM_nutrient
  ! file name to store the final daily concentration results.
  character(len=*), parameter :: file_daily_conc  = 'daily_concentration.out' 
  
  !State Fluxes variable for output
  integer(i4)                      :: timeStep_model_outputs_wqm ! timestep for writing model outputs
  integer(i4),      public    :: nOUTstate_wqm       ! total number of states that can be written out
  logical, dimension(:), allocatable, public :: outputFlxState_wqm         ! Define model outputs see "wqm_outputs.nml"
  character(256), parameter, public      :: Version_WQM = '1.0'          ! model version

  
  
  
  
  
  
END MODULE mo_wqm_global_variables