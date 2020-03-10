!> \file mo_wqm_write_fluxes_states.f90

!> \brief Creates NetCDF output for different fluxes and state variables of mHM.

!> \details NetCDF is first initialized and later on variables are put to the NetCDF.

!  HISTORY
!>     \authors Matthias Zink
!>     \date Apr 2013
!      Modified:
!          David Schaefer, Aug 2015 - major rewrite
!          Stephan Thober, Oct 2015 - adapted to mRM
!          Xiaoqiang Yang, Jul 2017 - Major modified to WQM
!       
module mo_wqm_write_fluxes_states

  use mo_kind                 , only: i4, dp
  use mo_string_utils         , only: num2str
  use mo_mrm_global_variables     , only: gridGeoRef
  use mo_mhm_constants        , only: nodata_dp
  use mo_netcdf               , only: NcDataset, NcDimension, NcVariable

  implicit none

  type OutputVariablewqm
     type(NcVariable)      :: ncn                 !> NcDataset which contains the variable
     logical               :: avg = .false.      !> average data before writing
     logical,  pointer     :: mask(:,:)          !> mask to reconstruct data
     real(dp), allocatable :: data(:)            !> store the data between writes
     integer(i4)           :: counter = 0        !> count the number of updateVariablewqm calls
     
   contains
     procedure, public :: updateVariablewqm
     procedure, public :: writeVariableTimestepwqm

  end type OutputVariablewqm

  ! constructor interface
  interface OutputVariablewqm
     procedure newOutputVariablewqm
  end interface OutputVariablewqm

  type OutputDatasetwqm
     integer(i4)                       :: ibasin      !> basin id
     type(NcDataset)                   :: ncn          !> NcDataset to write
     type(OutputVariablewqm), allocatable :: vars(:)     !> store all created (dynamic) variables
     integer(i4)                       :: counter = 0 !> count written time steps

   contains
     procedure, public :: updateDatasetwqm
     procedure, public :: writeTimestepwqm
     procedure, public :: close

  end type OutputDatasetwqm

  ! constructor interface
  interface OutputDatasetwqm
     procedure newOutputDatasetwqm
  end interface OutputDatasetwqm

  
  private
  
  public :: OutputDatasetwqm
#ifdef pgiFortran154
  public :: newOutputDatasetwqm
#endif

contains
  
  !------------------------------------------------------------------
  !     NAME
  !         newOutputVariablewqm
  !
  !     PURPOSE
  !>        \brief Initialize OutputVariablewqm
  !
  !     CALLING SEQUENCE
  !         var = OutputVariablewqm(ncn, ncells, mask, avg)
  !
  !     INTENT(IN)
  !>        \param[in] "type(NcDataset)   :: ncn"        -> NcDataset which contains the variable
  !>        \param[in] "integer(i4)       :: ncells"    -> number of cells in basin
  !>        \param[in] "logical, target   :: mask(:,:)" -> mask to reconstruct data
  !>        \param[in] "logical, optional :: avg"       -> average the data before writing
  !
  !     INTENT(INOUT)
  !         None
  !
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
  !         \return type(OutputVariablewqm)
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author David Schaefer
  !>        \date June 2015
  function newOutputVariablewqm(ncn, ncells, mask, avg) result(out)
    type(NcVariable), intent(in) :: ncn
    integer(i4), intent(in)      :: ncells
    logical, target              :: mask(:,:)
    logical, optional            :: avg
    type(OutputVariablewqm)         :: out

    allocate(out%data(ncells))
    out%ncn   =  ncn
    out%mask => mask
    out%data =  0
    if (present(avg)) out%avg = avg
  end function newOutputVariablewqm

  !------------------------------------------------------------------
  !     NAME
  !         updateVariablewqm
  !
  !     PURPOSE
  !>        \brief Update OutputVariablewqm
  !>        \details Add the array given as actual argument 
  !>                 to the derived type's component 'data'
  !
  !     CALLING SEQUENCE
  !         -> with ncn of type(OutputVariablewqm):
  !         call var%updateVariablewqm(data)              
  !
  !     INTENT(IN)
  !>        \param[in] "type(NcDataset)   :: ncn"        -> NcDataset which contains the variable
  !>        \param[in] "integer(i4)       :: ncells"    -> number of cells in basin
  !>        \param[in] "logical, target   :: mask(:,:)" -> mask to reconstruct data
  !>        \param[in] "logical, optional :: avg"       -> average the data before writing
  !
  !     INTENT(INOUT)
  !         None
  !
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
  !         \return type(OutputVariablewqm)
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author David Schaefer
  !>        \date June 2015
  subroutine updateVariablewqm(self, data)
    class(OutputVariablewqm), intent(inout) :: self
    real(dp)             , intent(in)    :: data(:)

    self%data    =  self%data + data
    self%counter = self%counter + 1

  end subroutine updateVariablewqm

  !------------------------------------------------------------------
  !     NAME
  !         writeVariableTimestepwqm
  !
  !     PURPOSE
  !>        \brief Write timestep to file
  !>        \details Write the content of the derived types's component
  !>                 'data' to file, average if necessary
  !
  !     CALLING SEQUENCE
  !         -> with var of type(OutputVariablewqm):
  !         call var%updateVariablewqm(data)              
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: timestep" 
  !>            -> index along the time dimension of the netcdf variable
  !
  !     INTENT(INOUT)
  !         None
  !
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
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author David Schafer
  !>        \date June 2015
  !         Modified:
  !             David Schaefer, Sep. 2015 - bugfix
  subroutine writeVariableTimestepwqm(self, timestep)
    class(OutputVariablewqm), intent(inout) :: self
    integer(i4)          , intent(in)    :: timestep

    if (self%avg) then
       self%data = self%data / real(self%counter, dp)
    end if

    call self%ncn%setData(unpack(self%data, self%mask, nodata_dp), &
         (/1,1,timestep/))
    self%data = 0
    self%counter = 0
    
  end subroutine writeVariableTimestepwqm
  
  !------------------------------------------------------------------
  !     NAME
  !         newOutputDatasetwqm
  !
  !     PURPOSE
  !>        \brief Initialize OutputDatasetwqm
  !>        \details Create and initialize the output file. If new a new output
  !>                 variable needs to be written, this is the first of two
  !>                 procedures to change (second: updateDatasetwqm)
  !
  !     CALLING SEQUENCE
  !         ncn = OutputDatasetwqm(ibasin, mask1)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: ibasin" -> basin id 
  !>        \param[in] "logical     :: mask1"  -> L1 mask to reconstruct the data
  !
  !     INTENT(INOUT)
  !         None
  !
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
  !         \return type(OutputDatasetwqm)
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Apr 2013
  !         Modified:
  !             R. Kumar & S. Thober, Aug. 2013 - code change to incorporate output timestep
  !                                               during writing of the netcdf file
  !             Matthias Zink       , Feb. 2014 - added aditional output: pet
  !             V. Prykhodk, J. Mai , Nov. 2014 - adding new variable infilSoil - case 16
  !             David Schaefer      , Jun. 2015 - major rewrite
  !             Stephan Thober      , Oct  2015 - adapted to mRM
  !             Xiaoqiang Yang      , Jul. 2017 - major modified for WQM 
  
  function newOutputDatasetwqm(ibasin, mask1, mask11) result(out)

    use mo_wqm_global_variables,  only : outputFlxState_wqm
    use mo_mrm_tools,         only : get_basin_info_mrm
    use mo_init_states,       only: get_basin_info
    use mo_global_variables,  only: nSoilHorizons_mHM

    integer(i4), intent(in) :: ibasin
    logical, target         :: mask1(:,:)
    logical, target         :: mask11(:,:)
    type(OutputDatasetwqm)     :: out 
    ! local
    integer(i4)             :: ii, nn, ncols, nrows, ncells
    integer(i4)             :: ncols11, nrows11, ncells11
    character(3)            :: dtype
    character(16)           :: dims1(3)
    character(16)           :: dims11(3)
    type(NcDataset)         :: ncn
    type(OutputVariablewqm)    :: tmpvars(size(outputFlxState_wqm)* nSoilHorizons_mHM)
    type(NcVariable)        :: xx

    call get_basin_info_mrm(ibasin, 11, ncols11, nrows11, ncells=ncells11)
    call get_basin_info (ibasin, 1, ncols, nrows, ncells=ncells)
	
    dtype = "f64"
    dims1 = (/"easting ", "northing","time    "/)
    dims11 = (/"easting11 ", "northing11","time      "/)
    ncn    = createOutputFile(ibasin)

    ii = 0
    !level 1
    if (outputFlxState_wqm(1)) then
       do nn = 1, nSoilHorizons_mHM
       ii = ii + 1
       tmpvars(ii) = OutputVariablewqm( &
            ncn%setVariable("SM_conc"//trim(num2str(nn, '(i2.2)')), dtype, dims1), &
            ncells, mask1, .true.)  !last argument define if taking the average value over setted step (in wqm_outputs.nml) or not .false. by default
       call writeVariableAttributes(tmpvars(ii), "soil moisture concentration of layer"//trim(num2str(nn, '(i2.2)')), "mg l-1")
       end do
    end if

    if (outputFlxState_wqm(2)) then
       ii = ii + 1
        
       tmpvars(ii) = OutputVariablewqm( &
            ncn%setVariable("cfastrunoff", dtype, dims1), &
            ncells, mask1, .true.)
       call writeVariableAttributes(tmpvars(ii), "fast runoff concentration", "mg l-1")
    end if

    if (outputFlxState_wqm(3)) then
       ii = ii + 1
       xx=ncn%setVariable("cslowrunoff", dtype, dims1)
       tmpvars(ii) = newOutputVariablewqm( &
            xx, &
            ncells, mask1, .true.)
       call writeVariableAttributes(tmpvars(ii), "slow runoff concentration", "mg l-1")
    end if

    if (outputFlxState_wqm(4)) then
       ii = ii + 1
       tmpvars(ii) = OutputVariablewqm( &
            ncn%setVariable("cbasefow", dtype, dims1), &
            ncells, mask1, .true.)
       call writeVariableAttributes(tmpvars(ii), "baseflow concentration", "mg l-1")
    end if

    if (outputFlxState_wqm(5)) then
       ii = ii + 1
       tmpvars(ii) = OutputVariablewqm( &
            ncn%setVariable("ctotalrunoff", dtype, dims1), &
            ncells, mask1, .true.)
       call writeVariableAttributes(tmpvars(ii), "total runoff concentration", "mg l-1")
    end if
	
    if (outputFlxState_wqm(6)) then
       ii = ii + 1
       tmpvars(ii) = OutputVariablewqm( &
            ncn%setVariable("soiluptakeN", dtype, dims1), &
            ncells, mask1)
       call writeVariableAttributes(tmpvars(ii), "total uptake amount in terrestrial phase", "mg m-2")
    end if
	
    if (outputFlxState_wqm(7)) then
       ii = ii + 1
       tmpvars(ii) = OutputVariablewqm( &
            ncn%setVariable("soildenitri", dtype, dims1), &
            ncells, mask1)
       call writeVariableAttributes(tmpvars(ii), "total denitrification amount", "mg m-2")
    end if

    if (outputFlxState_wqm(8)) then
       ii = ii + 1
       tmpvars(ii) = OutputVariablewqm( &
            ncn%setVariable("soilmineralN", dtype, dims1), &
            ncells, mask1)
       call writeVariableAttributes(tmpvars(ii), "total mineralization amount", "mg m-2")
    end if

    if (outputFlxState_wqm(9)) then
       ii = ii + 1
       tmpvars(ii) = OutputVariablewqm( &
            ncn%setVariable("INfrtmanapp", dtype, dims1), &
            ncells, mask1)
       call writeVariableAttributes(tmpvars(ii), "total IN application from fertilizer and manure ", "mg m-2")
    end if

    ! level 11
    if (outputFlxState_wqm(10)) then
       ii = ii + 1
       tmpvars(ii) = OutputVariablewqm( &
            ncn%setVariable("concmod", dtype, dims11), &
            ncells11, mask11, .true.)
       call writeVariableAttributes(tmpvars(ii), "output concentration of each reach", "mg l-1")
    end if
    if (outputFlxState_wqm(11)) then
       ii = ii + 1
       tmpvars(ii) = OutputVariablewqm( &
            ncn%setVariable("waterdenitri", dtype, dims11), &
            ncells11, mask11)
       call writeVariableAttributes(tmpvars(ii), "instream denitrification amount", "mg m-2")
    end if
    if (outputFlxState_wqm(12)) then
       ii = ii + 1
       tmpvars(ii) = OutputVariablewqm( &
            ncn%setVariable("waterassimi", dtype, dims11), &
            ncells11, mask11)
       call writeVariableAttributes(tmpvars(ii), "instream assimilatory uptake amount", "mg m-2")
    end if
	

    out = OutputDatasetwqm(ibasin, ncn, tmpvars(1:ii))
	

  end function newOutputDatasetwqm
  
  !------------------------------------------------------------------
  !     NAME
  !         updateDatasetwqm
  !
  !     PURPOSE
  !>        \brief Update all variables .
  !>        \details Call the type bound procedure updateVariablewqm for
  !>                 all output variables. If a new output
  !>                 variable needs to be written, this is the second
  !>                 of two procedures to change (first: newOutputDatasetwqm)
  !
  !     CALLING SEQUENCE
  !        with ncn of type(OutputDatasetwqm):
  !        call ncn%updateDatasetwqm(&
  !             self         , sidx         , eidx,           ,    &
  !             L11_qMod )
  !             
  !
  !     INTENT(IN)
  !>             \param[in] "sidx"        -> start index of the basin related data in L1_* arguments
  !>             \param[in] "eidx"        -> end index of the basin related data in L1_* arguments
  !>             \param[in] "L11_qMod"
  !
  !     INTENT(INOUT)
  !         None
  !
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
  !         None
  !
  !     LITERATURE
  !         None
  !
  !
  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Apr 2013
  !         Modified:
  !             R. Kumar & S. Thober, Aug. 2013 - code change to incorporate output timestep
  !                                               during writing of the netcdf file
  !             L. Samaniego et al.,  Dec  2013 - nullify pointer
  !             Matthias Zink,        Feb. 2014 - added aditional output: pet
  !             V. Prykhodk, J. Mai,  Nov. 2014 - adding new variable infilSoil - case 16
  !             David Schaefer      , Jun. 2015 - major rewrite
  !             Stephan Thober      , Oct  2015 - adapted to mRM
  !             Xiaoqiang Yang      , Jul  2017 - Major modified for WQM
  subroutine updateDatasetwqm(self, sidx1, eidx1, sidx11, eidx11, &
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
            L11_aquaticAssimil    )
    
    use mo_wqm_global_variables,  only : outputFlxState_wqm
    use mo_global_variables,      only: nSoilHorizons_mHM

    class(OutputDatasetwqm), intent(inout), target :: self
    integer(i4),          intent(in)            :: sidx1, eidx1
    integer(i4),          intent(in)            :: sidx11, eidx11
    ! fluxes,

    real(dp),             intent(in)            :: L1_csoilMoist(:,:,:)
    real(dp),             intent(in)            :: L1_cfastRunoff(:,:)
    real(dp),             intent(in)            :: L1_cslowRunoff(:,:)
    real(dp),             intent(in)            :: L1_cbaseflow(:,:)
    real(dp),             intent(in)            :: L1_ctotal_runoff(:,:)
    real(dp),             intent(in)            :: L1_soilUptakeN(:)
    real(dp),             intent(in)            :: L1_soilDenitri(:)
    real(dp),             intent(in)            :: L1_soilMineralN(:)
    real(dp),             intent(in)            :: L1_soilINfrtmanapp(:)
    real(dp),             intent(in)            :: L11_concMod(:,:)
    real(dp),             intent(in)            :: L11_aquaticDenitri(:)
    real(dp),             intent(in)            :: L11_aquaticAssimil(:)
    ! local
    type(OutputVariablewqm), pointer :: vars(:)
    integer(i4)       :: ii, nn

    ii = 0
    vars  => self%vars

    if (outputFlxState_wqm(1)) then
       do nn = 1, nSoilHorizons_mHM
       ii = ii + 1
#ifdef pgiFortran
       call updateVariablewqm(vars(ii), L1_csoilMoist(sidx1:eidx1,nn,1))
#else
       call vars(ii)%updateVariablewqm(L1_csoilMoist(sidx1:eidx1,nn,1))
#endif
       end do
    end if
    
    if (outputFlxState_wqm(2)) then
       ii = ii + 1
#ifdef pgiFortran
       call updateVariablewqm(vars(ii), L1_cfastRunoff(sidx1:eidx1,1))
#else
       call vars(ii)%updateVariablewqm(L1_cfastRunoff(sidx1:eidx1,1))
#endif
    end if

    if (outputFlxState_wqm(3)) then
       ii = ii + 1
#ifdef pgiFortran
       call updateVariablewqm(vars(ii), L1_cslowRunoff(sidx1:eidx1,1))
#else
       call vars(ii)%updateVariablewqm(L1_cslowRunoff(sidx1:eidx1,1))
#endif
    end if
    if (outputFlxState_wqm(4)) then
       ii = ii + 1
#ifdef pgiFortran
       call updateVariablewqm(vars(ii), L1_cbaseflow(sidx1:eidx1,1))
#else
       call vars(ii)%updateVariablewqm(L1_cbaseflow(sidx1:eidx1,1))
#endif
    end if
    if (outputFlxState_wqm(5)) then
       ii = ii + 1
#ifdef pgiFortran
       call updateVariablewqm(vars(ii), L1_ctotal_runoff(sidx1:eidx1,1))
#else
       call vars(ii)%updateVariablewqm(L1_ctotal_runoff(sidx1:eidx1,1))
#endif
    end if
    if (outputFlxState_wqm(6)) then
       ii = ii + 1
#ifdef pgiFortran
       call updateVariablewqm(vars(ii), L1_soilUptakeN(sidx1:eidx1))
#else
       call vars(ii)%updateVariablewqm(L1_soilUptakeN(sidx1:eidx1))
#endif
    end if
    if (outputFlxState_wqm(7)) then
       ii = ii + 1
#ifdef pgiFortran
       call updateVariablewqm(vars(ii), L1_soilDenitri(sidx1:eidx1))
#else
       call vars(ii)%updateVariablewqm(L1_soilDenitri(sidx1:eidx1))
#endif
    end if
    if (outputFlxState_wqm(8)) then
       ii = ii + 1
#ifdef pgiFortran
       call updateVariablewqm(vars(ii), L1_soilMineralN(sidx1:eidx1))
#else
       call vars(ii)%updateVariablewqm(L1_soilMineralN(sidx1:eidx1))
#endif
    end if
    if (outputFlxState_wqm(9)) then
       ii = ii + 1
#ifdef pgiFortran
       call updateVariablewqm(vars(ii), L1_soilINfrtmanapp(sidx1:eidx1))
#else
       call vars(ii)%updateVariablewqm(L1_soilINfrtmanapp(sidx1:eidx1))
#endif
    end if
 
    !Level 11 
    if (outputFlxState_wqm(10)) then
       ii = ii + 1
#ifdef pgiFortran
       call updateVariablewqm(vars(ii), L11_concMod(sidx11:eidx11, 1))
#else
       call vars(ii)%updateVariablewqm(L11_concMod(sidx11:eidx11, 1))
#endif
    end if
    if (outputFlxState_wqm(11)) then
       ii = ii + 1
#ifdef pgiFortran
       call updateVariablewqm(vars(ii), L11_aquaticDenitri(sidx11:eidx11))
#else
       call vars(ii)%updateVariablewqm(L11_aquaticDenitri(sidx11:eidx11))
#endif
    end if
    if (outputFlxState_wqm(12)) then
       ii = ii + 1
#ifdef pgiFortran
       call updateVariablewqm(vars(ii), L11_aquaticAssimil(sidx11:eidx11))
#else
       call vars(ii)%updateVariablewqm(L11_aquaticAssimil(sidx11:eidx11))
#endif
    end if

 
  end subroutine updateDatasetwqm

  !------------------------------------------------------------------
  !     NAME
  !         writeTimestepwqm
  !
  !     PURPOSE
  !>        \brief Write all accumulated data.
  !>        \details Write all accumulated and potentially averaged
  !>                 data to disk.
  !
  !     CALLING SEQUENCE
  !         -> with ncn of type(OutputDatasetwqm)
  !         call ncn%writeTimestepwqm(timestep)     
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: timestep" The model timestep to write
  !
  !     INTENT(INOUT)
  !         None
  !
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
  !         \return type(OutputVariablewqm)
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author David Schaefer
  !>        \date June 2015
  subroutine writeTimestepwqm(self, timestep)
    class(OutputDatasetwqm), intent(inout), target :: self
    integer(i4),          intent(in)    :: timestep
    integer(i4)                         :: ii
    type(NcVariable)                    :: tvar

    self%counter = self%counter + 1

    ! add to time variable
    tvar = self%ncn%getVariable("time")
    call tvar%setData(timestep, (/self%counter/))

    do ii = 1, size(self%vars)
       call self%vars(ii)%writeVariableTimestepwqm(self%counter)
    end do
    
  end subroutine writeTimestepwqm

  !------------------------------------------------------------------
  !     NAME
  !         close
  !
  !     PURPOSE
  !>        \brief Close the file
  !>        \details Close the file associated with variable of
  !>                 type(OutputDatasetwqm)
  !
  !     CALLING SEQUENCE
  !         -> with ncn of type(OutputDatasetwqm):
  !         call ncn%close()     
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !         None
  !
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
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Rohini Kumar & Stephan Thober
  !>        \date August 2013
  !         Modified:
  !             David Schaefer, June 2015 - adapted to new structure
  !             Stephan Thober, Oct  2015 - adapted to mRM
  subroutine close(self)

    use mo_String_utils,         only: num2str
    use mo_message,              only: message
    use mo_mrm_global_variables, only: dirOut 

    class(OutputDatasetwqm) :: self
    call self%ncn%close()
    call message('  OUTPUT: saved netCDF file for basin', trim(num2str(self%ibasin)))
    call message('    to ', trim(dirOut(self%ibasin)))

  end subroutine close

  !------------------------------------------------------------------
  !     NAME
  !         createOutputFile
  !
  !     PURPOSE
  !>        \brief Create and initialize output file
  !>        \details Create output file, write all non-dynamic variables
  !>                 and global attributes for the given basin.
  !>
  !
  !     CALLING SEQUENCE
  !         ncn = createOutputFile(ibasin)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4)     :: ibasin"      -> basin id
  !>        \param[in] "logical, target :: mask1(:,:)"  -> level11 mask
  !
  !     INTENT(INOUT)
  !         None
  !
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
  !         \return type(NcDataset)
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author David Schaefer
  !>        \date June 2015
  !         modified,
  !             Stephan Thober, Oct  2015 - adapted to mRM
  !             Xiaoqiang Yang, Jul  2017 - Major modified for WQM
  function createOutputFile(ibasin) result(ncn)

  
    !use mo_global_variables,      only: level1
    use mo_mrm_global_variables,  only: dirOut, evalPer, level1, level11
    use mo_wqm_global_variables,  only: Version_WQM
    use mo_julian,                only: dec2date

    integer(i4), intent(in) :: ibasin
    type(NcDataset)         :: ncn
    type(NcDimension)       :: dimids1(5)
    type(NcVariable)        :: var
    integer(i4)             :: day, month, year
    character(128)          :: fname, unit, date, time, datetime
    real(dp), allocatable   :: northing(:), easting(:), lat(:,:), lon(:,:)
    real(dp), allocatable   :: northing11(:), easting11(:), lat11(:,:), lon11(:,:)

    fname = trim(dirOut(ibasin)) // 'WQM_Fluxes_States.nc'
    call mapCoordinates(ibasin, level1, northing, easting)
    call geoCoordinates(ibasin, level1, lat, lon)
    call mapCoordinates(ibasin, level11, northing11, easting11)
    call geoCoordinates(ibasin, level11, lat11, lon11)

    ncn = NcDataset(trim(fname),"w")
    dimids1 = (/&
         ncn%setDimension("easting",  size(easting)), &
         ncn%setDimension("northing", size(northing)), &
         ncn%setDimension("easting11",  size(easting11)), &
         ncn%setDimension("northing11", size(northing11)), &
         ncn%setDimension("time",     0) &
         /)

    ! time units
    call dec2date(real(evalPer(ibasin)%julStart, dp), dd=day, mm=month, yy=year)
    write(unit,"('hours since ', i4, '-' ,i2.2, '-', i2.2, 1x, '00:00:00')") year, month, day

    ! time
    var = ncn%setVariable("time", "i32", (/ dimids1(5) /))
    call var%setAttribute("units", unit)
    call var%setAttribute("long_name","time")

    ! northing
    var = ncn%setVariable("northing", "f64", (/ dimids1(2) /))
    call var%setData(northing)
    call var%setAttribute("units","m")
    call var%setAttribute("long_name","y-coordinate in cartesian coordinates GK4")

    ! easting
    var = ncn%setVariable("easting", "f64", (/ dimids1(1) /))
    call var%setData(easting)
    call var%setAttribute("units","m")
    call var%setAttribute("long_name","x-coordinate in cartesian coordinates GK4")
    ! northing at L11
    var = ncn%setVariable("northing11", "f64", (/ dimids1(4) /))
    call var%setData(northing11)
    call var%setAttribute("units","m")
    call var%setAttribute("long_name","y-coordinate in cartesian coordinates GK4")

    ! easting at L11
    var = ncn%setVariable("easting11", "f64", (/ dimids1(3) /))
    call var%setData(easting11)
    call var%setAttribute("units","m")
    call var%setAttribute("long_name","x-coordinate in cartesian coordinates GK4")

    ! lon
    var = ncn%setVariable("lon","f64",dimids1(1:2))
    call var%setData(lon)
    call var%setAttribute("units","degerees_east")
    call var%setAttribute("long_name","longitude")
    call var%setAttribute("missing_value","-9999.")

    ! lat
    var = ncn%setVariable("lat","f64",dimids1(1:2))
    call var%setData(lat)
    call var%setAttribute("units","degerees_north")
    call var%setAttribute("long_name","latitude")
    call var%setAttribute("missing_value","-9999.")
    ! lon at L11
    var = ncn%setVariable("lon11","f64",dimids1(3:4))
    call var%setData(lon11)
    call var%setAttribute("units","degerees_east")
    call var%setAttribute("long_name","longitude")
    call var%setAttribute("missing_value","-9999.")

    ! lat at L11
    var = ncn%setVariable("lat11","f64",dimids1(3:4))
    call var%setData(lat11)
    call var%setAttribute("units","degerees_north")
    call var%setAttribute("long_name","latitude")
    call var%setAttribute("missing_value","-9999.")

    ! global attributes
    call date_and_time(date=date, time=time)
    write(datetime,"(a4,'-',a2,'-',a2,1x,a2,':',a2,':',a2)") date(1:4), &
         date(5:6), date(7:8), time(1:2), time(3:4), time(5:6)

    call ncn%setAttribute("title","WQMv"//trim(Version_WQM)//" simulation outputs")
    call ncn%setAttribute("creation_date",datetime)
    call ncn%setAttribute("institution",&
         "Helmholtz Centre for Environmental Research - UFZ, "// &
         "Department Aquatic Ecosystem Analysis and Management, " // &
         "Hydrological and river water quality modelling Group" )

  end function createOutputFile

  !------------------------------------------------------------------
  !     NAME
  !         writeVariableAttributes
  !
  !     PURPOSE
  !>        \brief Write output variable attributes
  !
  !     CALLING SEQUENCE
  !         call writeVariableAttributes(var, long_name, unit) 
  !
  !     INTENT(IN)
  !>        \param[in] "type(OutputVariablewqm) :: var" 
  !>        \param[in] "character(*)         :: long_name"    -> variable name
  !>        \param[in] "character(*)         :: unit"         -> physical unit
  !
  !     INTENT(INOUT)
  !         None
  !
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
  !        None
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author David Schaefer
  !>        \date June 2015
  subroutine writeVariableAttributes(var, long_name, unit)
    type(OutputVariablewqm), intent(in) :: var
    character(*)        , intent(in) :: long_name, unit

    call var%ncn%setAttribute("_FillValue",nodata_dp)
    call var%ncn%setAttribute("long_name",long_name)
    call var%ncn%setAttribute("unit",unit)
    call var%ncn%setAttribute("scale_factor",1.0_dp)
    call var%ncn%setAttribute("missing_value", "-9999.")
    call var%ncn%setAttribute("coordinates","lat lon")

  end subroutine writeVariableAttributes
  

  !------------------------------------------------------------------
  !     NAME
  !         mapCoordinates
  !
  !     PURPOSE
  !>        \brief Generate map coordinates
  !>        \details Generate map coordinate arrays for given basin and level
  !
  !     CALLING SEQUENCE
  !         call mapCoordinates(ibasin, level, y, x)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4)      :: iBasin" -> basin number
  !>        \param[in] "type(geoGridRef) :: level"  -> grid reference
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !>        \param[out] "real(:)  :: y(:)"          -> y-coordinates
  !>        \param[out] "real(dp) :: x(:)"          -> x-coorindates
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
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Apr 2013
  !         Modified:
  !             Stephan Thober, Nov 2013 - removed fproj dependency
  !             David Schaefer, Jun 2015 - refactored the former subroutine CoordSystem
  subroutine mapCoordinates(ibasin, level, y, x)

    implicit none

    integer(i4),      intent(in)               :: iBasin
    type(gridGeoRef), intent(in)               :: level
    real(dp),         intent(out), allocatable :: x(:), y(:)
    integer(i4)                                :: ii, ncols, nrows
    real(dp)                                   :: cellsize

    cellsize = level%cellsize(ibasin)
    nrows    = level%nrows(ibasin)
    ncols    = level%ncols(ibasin)

    allocate(x(nrows), y(ncols))

    x(1) =  level%xllcorner(ibasin) + 0.5_dp * cellsize
    do ii = 2, nrows
       x(ii)   =  x(ii-1) + cellsize
    end do

    ! inverse for Panoply, ncview display
    y(ncols) =  level%yllcorner(ibasin) + 0.5_dp * cellsize
    do ii = ncols-1,1,-1
       y(ii)   =  y(ii+1) + cellsize
    end do

  end subroutine mapCoordinates

  !------------------------------------------------------------------
  !     NAME
  !         geoCoordinates
  !
  !     PURPOSE
  !>        \brief Generate geographic coordinates
  !>        \details Generate geographic coordinate arrays for given basin and level
  !
  !     CALLING SEQUENCE
  !         call mapCoordinates(ibasin, level, y, x)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4)      :: iBasin"    -> basin number
  !>        \param[in] "type(gridGeoRef) :: level"     -> grid reference 
  !
  !     INTENT(INOUT)
  !         None
  !
  !     INTENT(OUT)
  !>        \param[out] "real(dp) :: lat(:,:)"         -> lat-coordinates
  !>        \param[out] "real(dp) :: lon(:,:)"         -> lon-coorindates
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
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author Matthias Zink
  !>        \date Apr 2013
  !         Modified:
  !             Stephan Thober, Nov 2013 - removed fproj dependency
  !             David Schaefer, Jun 2015 - refactored the former subroutine CoordSystem
  !             Stephan Thober, Sep 2015 - using mask to unpack coordinates
  !             Stephan Thober, Oct 2015 - writing full lat/lon again
  !             Stephan Thober, Oct 2015 - adapted to mRM
  !             Xiaoqiang Yang, Jul 2017 - Modified for WQM and be universal for both L1&L11 
  subroutine geoCoordinates(ibasin, level, lat, lon)

    use mo_mrm_global_variables, only :  level1  !level11,L11_rect_latitude, L11_rect_longitude,
    use mo_global_variables, only : L1_rect_latitude, L1_rect_longitude

    implicit none

    integer(i4),      intent(in)               :: iBasin
    type(gridGeoRef), intent(in)               :: level
    real(dp),         intent(out), allocatable :: lat(:,:), lon(:,:)
    integer(i4)                                :: ncols, nrows
    integer(i4)                                :: ii, pos 


    nrows = level%nrows(ibasin) 
    ncols = level%ncols(ibasin) 
	
    if (allocated(lat)) deallocate(lat)
    if (allocated(lon)) deallocate(lon)    
      
    pos = 1 
    if ( ibasin .gt. 1 ) then 
       do ii = 1, ibasin -1 
          pos = pos + level%ncols(ii) * level%nrows(ii) 
       end do
    end if

    if ( nrows .eq. level1%nrows(ibasin) ) then   
    lat = reshape(L1_rect_latitude(pos:pos+nrows*ncols-1),  (/nrows, ncols/)) 
    lon = reshape(L1_rect_longitude(pos:pos+nrows*ncols-1), (/nrows, ncols/)) 
    else
    !*****************************************************************************************************
    !constrained by mo_mrm_read_latlon.f90 which assign L11_rect_latitude&L11_rect_longitude
    !shape according to the input file: latlon.nc.(line:217,241)
    !However, when the routing resolution is set unequally with model resolution, a new latlon.nc needs to
    !be generated...Here should be further improved!!! ---X. Yang 20-07-2017
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !lat = reshape(L11_rect_latitude(pos:pos+nrows*ncols-1),  (/nrows, ncols/)) 
    !lon = reshape(L11_rect_longitude(pos:pos+nrows*ncols-1), (/nrows, ncols/))
    lat = reshape((/(ii, ii=1,nrows*ncols, 1)/), (/nrows, ncols/))
    lon = reshape((/(ii, ii=1,nrows*ncols, 1)/), (/nrows, ncols/))
    end if
    
  end subroutine geoCoordinates


  !------------------------------------------------------------------
  !     NAME
  !         fluxesUnit
  !
  !     PURPOSE
  !>        \brief Generate a unit string
  !>        \details Generate the unit string for the output variable
  !>                 netcdf attribute based on modeling timestep
  !
  !     CALLING SEQUENCE
  !         unit = fluxesUnit(iBasin)
  !
  !     INTENT(IN)
  !>        \param[in] "integer(i4) :: iBasin"    -> basin id
  !
  !     INTENT(INOUT)
  !         None
  !
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
  !         \return character(16)
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !         None
  !
  !     LITERATURE
  !         None
  !
  !     HISTORY
  !>        \author David Schaefer
  !>        \date June 2015
  function fluxesUnit(ibasin)

    use mo_mrm_global_variables,  only : simPer, NTSTEPDAY, timestep, timeStep_model_outputs_mrm

    integer(i4), intent(in) :: ibasin
    character(16)           :: fluxesUnit
    real(dp)                :: ntsteps

    if ( timestep*timestep_model_outputs_mrm .eq. 1 ) then
       fluxesUnit = 'mm h-1'
    else if (timestep_model_outputs_mrm > 1) then
       fluxesUnit = 'mm '//trim(adjustl(num2str(timestep)))//'h-1'
    else if (timestep_model_outputs_mrm .eq. 0) then
       ntsteps = ( simPer(iBasin)%julEnd - simPer(iBasin)%julStart + 1 ) * NTSTEPDAY
       fluxesUnit = 'mm '//trim(adjustl(num2str(nint(ntsteps))))//'h-1'
    else if (timestep_model_outputs_mrm .eq. -1) then
       fluxesUnit = 'mm d-1'
    else if (timestep_model_outputs_mrm .eq. -2) then
       fluxesUnit = 'mm month-1'
    else if (timestep_model_outputs_mrm .eq. -3) then
       fluxesUnit = 'mm a-1'
    else
       fluxesUnit = ''
    endif

  end function fluxesUnit

end module mo_wqm_write_fluxes_states
