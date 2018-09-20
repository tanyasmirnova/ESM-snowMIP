program simple_driver

  !
  ! Through MODULE_IO, we access most of our top-level input/output routines.  This
  ! includes routines OPEN_FORCING_FILE, READ_FORCING_TEXT, INITIALIZE_OUTPUT,
  ! OUTPUT_TIME, OUTPUT_LEVELS, OUTPUT_VAR, FINISH_OUTPUT_FOR_TIME, and OUTPUT_CLOSE.
  !
  ! Module MODULE_IO uses lower-level modules MODULE_ASCII_IO and MODULE_NETCDF_IO.
  !

  use module_io

  !
  ! Module MODULE_RUCLSM_UTILITY contains a few utility routines needed for setting
  ! up some data for RUC LSM.  

!This driver program makes use of module subroutines
  ! CALTMP and CALHUM.
  !

  use module_model_constants
  use module_Noahlsm_utility

  !
  ! MODULE_SF_RUCLSM contains the RUC LSM physics code.  This driver program
  ! makes use of module subroutines LSMRUC, SOILVEGIN, and SFCDIF_off.
  !

  use module_sf_ruclsm

  !
  ! MODULE_SFCDIF_WRF contains the MYJSFC version of SFCDIF.  This driver program
  ! makes use of subroutines MYJSFCINIT and SFCDIF_MYJ.
  !

  use module_sfcdif_wrf

  !
  ! KWM_DATE_UTILITIES contains handy subroutines for manipulating and
  ! computing date/time information.
  !

  use kwm_date_utilities

  implicit none

  integer, parameter  :: ime = 1
  integer, parameter  :: jme = 1

!  INTEGER, PARAMETER :: NSLTYPE=30
  !
  ! Command-line:
  !

  character(len=256) :: executable_name   ! The name of the executable, as found by Fortran library 
  !                                       ! function GETARG

  character(len=256) :: forcing_filename  ! The name of the initial/forcing conditions file, as found
  !                                       ! by Fortran library function GETARG

  character(len=256) :: FPEpara_filename  ! The name of the FPE parameters file, as found
  !                                       ! by Fortran library function GETARG
  
  ! Loop increment variables:
  !

  integer            :: ktime   ! A counter for the timesteps in the main loop TIMELOOP.
  character(len=12)  :: nowdate ! The date of each time step, ( YYYYMMDDHHmm ), updated in each step 
  !                             ! of the main loop TIMELOOP

  character(len=9) :: teststr
  character(len=256) :: cmd
  integer :: reloop_count = 0

  !
  ! Useful data attributes describing the data in the initial/forcing conditions file
  !

  character(len=4096) :: infotext         ! Character string returned by subroutine OPEN_FORCING_FILE, 
                                          ! giving some possibly useful information for the user.

  real                :: latitude         ! Latitude of the point ( Degrees North )
  real                :: longitude        ! Longitude of the point ( Degrees East )
  logical             :: loop_for_a_while ! Whether to loop the same year ad infinitum
  character(len=12)   :: startdate        ! Starting date of the data ( YYYYMMDDHHmm ) 
  character(len=12)   :: enddate          ! Ending date of the data ( YYYYMMDDHHmm )
  integer             :: forcing_timestep ! The time interval ( seconds ) of the data in the forcing file
  integer             :: ruclsm_timestep  ! The timestep ( seconds ) to use when integrating the RUC g LSM
  real, dimension(12) :: albedo_monthly   ! Monthly values of background (i.e., snow-free) albedo ( Fraction [0.0-1.0] )
  real, dimension(12) :: shdfac_monthly   ! Monthly values for green vegetation fraction ( Fraction [0.0-1.0] )
  real, dimension(12) :: z0brd_monthly    ! Monthly values for background (i.e., snow-free) roughness length ( m )
  real, dimension(12) :: lai_monthly      ! Monthly values for Leaf Area Index ( dimensionless )

  !
  ! Various arguments to subroutine LSMRUC:
  !

  logical :: frpcpn     ! Flag: .true. - mixed phase precipitation
  real, dimension (1:ime,1:jme) :: FFROZP     ! Fraction of precip which is frozen (0.0 - 1.0).
  integer, dimension (1:ime,1:jme) :: ICE        ! Flag for sea-ice (1) or land (0).
  real, dimension (1:ime,1:jme)    :: XICE       ! Flag for sea-ice (1) or land (0).
  logical :: fractional_seaice  ! 0 no fractional, 1 - yes
  real    :: xice_threshold  ! threshold below which no ice
  integer :: ISURBAN    ! Vegetation category for urban land class.
  real    :: DT         ! Time step (seconds).
  !real    :: ZLVL       ! Height at which atmospheric forcing variables are taken to be valid (m)
  real    :: zlvl_wind  ! Height at which the wind forcing variable is taken to be valid (m)
  real, pointer, dimension(:) :: SLDPTH ! Thicknesses of each soil level
  integer :: NSOIL      ! Number of soil levels.
  logical :: LOCAL      ! 
  logical :: myj        ! use MYJ scheme to compute exchange coefficients
  character(len=256) :: LLANDUSE  ! Land-use dataset.  Valid values are :
  !                               ! "USGS" (USGS 24/27 category dataset) and
  !                               ! "MODIFIED_IGBP_MODIS_NOAH" (MODIS 20-category dataset)
  character(len=256) :: LSOIL     ! Soil-category dateset.  Only "STAS" (STATSGO dataset) supported.
  real, dimension (1:ime,1:jme)    :: ZLVL       ! Downward longwave radiation flux at surface (W m-2) [Forcing]
  real, dimension (1:ime,1:jme)    :: LWDN       ! Downward longwave radiation flux at surface (W m-2) [Forcing]
  real, dimension (1:ime,1:jme)    :: SOLDN      ! Downward shortwave radiation flux at surface (W m-2) [Forcing]
  real, dimension (1:ime,1:jme)    :: SOLNET     ! Net downward shortwave radiation flux at the surface (W m-2)
  real, dimension (1:ime,1:jme)    :: SFCPRS     ! Surface atmospheric pressure (Pa) [Forcing]
  real, dimension (1:ime,1:jme)    :: PRCP       ! RAIN rate (kg m-2 s-1) [Forcing]
  real, dimension (1:ime,1:jme)    :: SNFALL     ! SNOW rate (kg m-2 s-1) [Forcing]
  real, dimension (1:ime,1:jme)    :: TABS       ! Air temperature (K) [Forcing]
  real, dimension (1:ime,1:jme)    :: Q2         ! Surface specific humidity (kg kg-1) [Forcing]
  real, dimension (1:ime,1:jme)    :: SFCSPD     ! Surface wind speed (m s-1) [Forcing]
  real    :: SFCU       ! West-to-east component of the surface wind (m s-1)
  real    :: SFCV       ! South-to-north component of the surface wind (m s-1)
  real    :: COSZ       ! Unused if we're not using urban canopy model.
  real, dimension (1:ime,1:jme)    :: PRCPRAIN   ! Unused.
  real, dimension (1:ime,1:jme)    :: SOLARDIRECT! Unused.
  real, dimension (1:ime,1:jme)    :: TH2        ! Potential temperature at level ZLVL (K)
  real, dimension (1:ime,1:jme)    :: t1V       ! Virtual skin temperature (K).  Used in SFCDIF_off for computing CM and CH, but not passed to LSMRUC
  real, dimension (1:ime,1:jme)    :: TH2V       ! Virtual potential temperature at level ZLVL (K).  Used in SFCDIF_off
  !                     ! for computing CM and CH, but not passed to LSMRUC
  real, dimension (1:ime,1:jme)    :: RHO        ! Air density (dummy value output from CALTMP, needs to be passed to LSMRUC).
  real, dimension (1:ime,1:jme)    :: Q2SAT      ! Saturated specific humidity (kg kg-1)
  real, dimension (1:ime,1:jme)    :: DQSDT2     ! Slope of the Saturated specific humidity curve W.R.T. Temperature.
  integer, dimension (1:ime,1:jme) :: VEGTYP     ! Vegetation category.
  integer, dimension (1:ime,1:jme) :: SOILTYP    ! Soil category.
  integer, dimension (1:ime,1:jme) :: SLOPETYP   ! Slope category.
  real, dimension (1:ime,1:jme)    :: SHDFAC     ! Shade factor (0.0-1.0).
  real, dimension (1:ime,1:jme)    :: SHDMIN     ! Minimum shade factor (0.0-1.0).
  real, dimension (1:ime,1:jme)    :: SHDMAX     ! Maximum shade factor (0.0-1.0).
  real, dimension (1:ime,1:jme)    :: ALBBCK     ! Background snow-free albedo (0.0-1.0).
  real, dimension (1:ime,1:jme)    :: SNOALB     ! Maximum snow albedo over deep snow (0.0-1.0)
  real, dimension (1:ime,1:jme)    :: TBOT       ! Deep-soil time-invariant temperature (K).  Representing sort of a mean annual air temperature.
  real, dimension (1:ime,1:jme)    :: Z0BRD      ! Background Z0 value (m).
  real, dimension (1:ime,1:jme)    :: Z0         ! Roughness length (m)
  real, dimension (1:ime,1:jme)    :: EMISSI     ! Surface emissivity (0.0 - 1.0).  This includes the snow-cover effect.
  real, dimension (1:ime,1:jme)    :: EMBRD      ! Background value (i.e., not including snow-cover effect) of surface emissivity (0.0 - 1.0)
  real, dimension (1:ime,1:jme)    :: CMC        ! Canopy moisture content (kg m-2)
  real, dimension (1:ime,1:jme)    :: t1         ! Skin temperature (K)
  real, pointer,     dimension(:) :: STC  ! Soil temperature (K)
  real, pointer,     dimension(:) :: SMC  ! Total soil moisture content (m3 m-3)
  real, pointer,     dimension(:) :: SH2O ! Liquid soil moisture content (m3 m-3)
  real, allocatable, dimension(:) :: ET   ! Plant transpiration from each soil level.
  real, allocatable, dimension(:) :: SMAV ! Soil Moisture Availability at each level, fraction between
                                          ! SMCWLT (SMAV=0.0) and SMCMAX (SMAV=1.0)

  real, allocatable, dimension(:,:,:) :: TSLB ! 3-D Soil temperature (K)
  real, allocatable, dimension(:,:,:) :: SOILMOIS ! 3-D Soil moisture
  real, allocatable, dimension(:,:,:) :: SLIQW ! 3-D Liquid soil moisture
  real, allocatable, dimension(:,:,:) :: SMFR3D  ! Soil temperature (K)
  real, allocatable, dimension(:,:,:) :: KEEPFR3DFLAG  ! Total soil moisture content (m3 m-3)
  real, allocatable, dimension(:,:,:) :: smcmass
  real, allocatable, dimension(:,:,:) :: smfrfrac
  real, allocatable, dimension(:,:,:) :: smliqfrac
  real, allocatable, dimension(:) :: zsmain
  real, allocatable, dimension(:) :: zshalf
  real :: smcfrmass
  real :: smcliqmass

  real, dimension (1:ime,1:jme)    :: RHOSNF
  real, dimension (1:ime,1:jme)    :: snowfallac
  real, dimension (1:ime,1:jme)    :: graupelncv
  real, dimension (1:ime,1:jme)    :: snowncv
  real, dimension (1:ime,1:jme)    :: rainncv
  real, dimension (1:ime,1:jme)    :: rainbl
  real, dimension (1:ime,1:jme)    :: precipfr
  real, dimension (1:ime,1:jme)    :: qcatm
  real, dimension (1:ime,1:28,1:jme) :: landusef 
  real, dimension (1:ime,1:19,1:jme) :: soilctop 
  integer :: nlcat, nscat  
  integer    :: mosaic_lu
  integer    :: mosaic_soil
  real, dimension (1:ime,1:jme)    :: hfrain
  real, dimension (1:ime,1:jme)    :: snowcan
  real, dimension (1:ime,1:jme)    :: lwup
  real, dimension (1:ime,1:jme)    :: swup
  real, dimension (1:ime,1:jme)    :: qsg
  real, dimension (1:ime,1:jme)    :: soilt1
  real, dimension (1:ime,1:jme)    :: tsnav
  real, dimension (1:ime,1:jme)    :: xland
  real, dimension (1:ime,1:jme)    :: sfcevp
  real, dimension (1:ime,1:jme)    :: acsnow
  real, dimension (1:ime,1:jme)    :: mavail
  real, dimension (1:ime,1:jme)    :: fltot
  real, dimension (1:ime,1:jme)    :: snflx
  integer :: iswater, isice
  real, dimension (1:ime,1:jme)    :: SNOWH      ! Physical snow depth.
  real, dimension (1:ime,1:jme)    :: SNOW       ! Snow water equivalent in kg/m^2/s
  real, dimension (1:ime,1:jme)    :: Hfus       ! energy of fusion ( W m{-2} )
  real, dimension (1:ime,1:jme)    :: SNEQV      ! Water equivalent of accumulated snow depth (m).
  real, dimension (1:ime,1:jme)    :: ALBEDO     ! Surface albedo including possible snow-cover effect.  This is set in LSMRUC,
  !                     ! overriding any value given; it should perhaps be INTENT(OUT) from LSMRUC.
  real, dimension (1:ime,1:jme)    :: CH         ! Exchange coefficient for head and moisture (m s-1).  An initial value is needed for SFCDIF_off.
  real, dimension (1:ime,1:jme)    :: CM         ! Exchange coefficient for momentum (m s-1).  An initial value is needed for SFCDIF_off.
  real, dimension (1:ime,1:jme)    :: ETA        ! Latent heat flux (evapotranspiration) ( W m{-2} )
  real, dimension (1:ime,1:jme)    :: SHEAT      ! Sensible heat flux ( W m{-2} )
  real, dimension (1:ime,1:jme)    :: ETAKIN     ! Latent heat flux (evapotranspiration) ( kg m{-2} s{-1} )
  real, dimension (1:ime,1:jme)    :: FDOWN      ! Radiation forcing at the surface ( W m{-2} )
  real, dimension (1:ime,1:jme)    :: EC         ! Latent heat flux component: canopy water evaporation ( W m{-2} )
  real, dimension (1:ime,1:jme)    :: EDIR       ! Latent heat flux component: direct soil evaporation ( W m{-2} )
  real, dimension (1:ime,1:jme)    :: ETT        ! Latent heat flux component: total plant transpiration ( W m{-2} )
  real, dimension (1:ime,1:jme)    :: ESNOW      ! Latent heat flux component: sublimation from (or deposition to) snowpack ( W m{-2} )
  real, dimension (1:ime,1:jme)    :: DRIP       ! Precipitation or dew falling through canopy, in excess of canopy holding capacity ( m )
  real, dimension (1:ime,1:jme)    :: DEW        ! Dewfall (or frostfall for T<273.15) ( m )
  real    :: BETA       ! Ratio of actual to potential evapotranspiration ( Fraction [0.0-1.0] )
  real, dimension (1:ime,1:jme)    :: ETP        ! Potential evapotranspiration ( W m{-2} )
  real, dimension (1:ime,1:jme)    :: SSOIL      ! Soil heat flux ( W m{-2} )
  real    :: FLX1       ! Latent heat flux from precipitation accumulating as snow ( W m{-2} )
  real    :: FLX2       ! Latent heat flux from freezing rain converting to ice ( W m{-2} )
  real    :: FLX3       ! Latent heat flux from melting snow ( W m{-2} )
  real, dimension (1:ime,1:jme)    :: SNOMLT     ! Snow melt water ( m )
  real, dimension (1:ime,1:jme)    :: SNCOVR     ! Fractional snow cover ( Fraction [0.0-1.0] )
  real, dimension (1:ime,1:jme)    :: RUNOFF1    ! Surface runoff, not infiltrating the soil ( m s{-1} )
  real, dimension (1:ime,1:jme)    :: RUNOFF2    ! Subsurface runoff, drainage out the bottom of the last soil layer ( m s{-1} )
  real    :: RUNOFF3    ! Internal soil layer runoff ( m s{-1} )
  real    :: RC         ! Canopy resistance ( s m{-1} )
  real, dimension (1:ime,1:jme)    :: PC         ! Plant coefficient, where PC * ETP = ETA ( Fraction [0.0-1.0] )
  real    :: RSMIN      ! Minimum canopy resistance ( s m{-1} )
  real, dimension (1:ime,1:jme)    :: XLAI       ! Leaf area index ( dimensionless )
  real    :: RCS        ! Incoming solar RC factor ( dimensionless )
  real    :: RCT        ! Air temperature RC factor ( dimensionless )
  real    :: RCQ        ! Atmospheric water vapor deficit RC factor ( dimensionless )
  real    :: RCSOIL     ! Soil moisture RC factor ( dimensionless )
  real, dimension (1:ime,1:jme)    :: SOILW      ! Available soil moisture in the root zone ( Fraction [SMCWLT-SMCMAX] )
  real, dimension (1:ime,1:jme)    :: SOILM      ! Total soil column moisture content, frozen and unfrozen ( m )
  real, dimension (1:ime,1:jme)    :: qsfc       ! Effective specific humidity at the surface ( kg kg{-1} )
  real, dimension (1:ime,1:jme)    :: QCG        ! Effective cloud water mixing ratio at the surface ( kg kg{-1} )
  real, dimension (1:ime,1:jme)    :: Q1         ! Effective water vapor mixing ratio at the surface ( kg kg{-1} )
  logical :: RDLAI2D    ! If RDLAI2D == .TRUE., then the XLAI value that we pass to LSMRUC will be used.
  !                     ! If RDLAI2d == .FALSE., then XLAI will be computed within LSMRUC, from table 
  !                     ! minimum and maximum values in VEGPARM.TBL, and the current Green Vegetation Fraction.
  logical :: USEMONALB  ! If USEMONALB == .TRUE., then the ALB value passed to LSMRUC will be used as the background
  !                     ! snow-free albedo term.  If USEMONALB == .FALSE., then ALB will be computed within LSMRUC
  !                     ! from minimum and maximum values in VEGPARM.TBL, and the current Green Vegetation Fraction.
  real    :: SNOTIME1   ! Age of the snow on the ground.
  real    :: RIBB       ! Bulk Richardson number used to limit the dew/frost.
!  real, dimension (1:NSLTYPE)    :: SMCWLT     ! oilting point ( m{3} m{-3} )
!  real, dimension (1:NSLTYPE)    :: SMCDRY     ! Dry soil moisture threshold where direct evaporation from the top layer ends ( m{3} m{-3} )
!  real, dimension (1:NSLTYPE)    :: SMCREF     ! Soil moisture threshold where transpiration begins to stress ( m{3} m{-3} )
!  real, dimension (1:NSLTYPE)    :: SMCMAX     ! Porosity, i.e., saturated value of soil moisture ( m{3} m{-3} )
!  integer :: NROOT      ! Number of root layers ( count )

  integer :: iz0tlnd    ! Option to turn on (IZ0TLND=1) or off (IZ0TLND=0) the vegetation-category-dependent
                        ! calculation of the Zilitinkivich coefficient CZIL in the SFCDIF subroutines.

  integer :: sfcdif_option ! Option to use previous (SFCDIF_OPTION=0) or updated (SFCDIF_OPTION=1) version of 
                           ! SFCDIF subroutine.

  !
  ! Some diagnostics computed from the output of subroutine LSMRUC
  !
  real, dimension (1:ime,1:jme) :: QFX       ! Evapotranspiration ( W m{-2} )  the sum of 1) direct evaporation
  !                 ! from soil; 2) evaporation from canopy; 3) total plant transpiration; 
  !                 ! and 4) evaporation from snowpack.  Mostly, this should be the 
  !                 ! same as ETA

  real, dimension (1:ime,1:jme) :: RES       ! Residual of the surface energy balance equation ( W m{-2} )
  real :: FUP       ! Upward longwave radiation flux from the surface ( W m{-2} )
  real :: F         ! Incoming shortwave and longwave radiation flux  ( W m{-2} )

  !
  ! Miscellaneous declarations
  !
  integer            :: ierr             ! Error flag returned from read routines.
  integer, parameter :: iunit = 10       ! Fortran unit number for reading initial/forcing conditions file.
  logical            :: use_urban_module ! Flag, set to TRUE in the initial/forcing conditions file, if the 
  !                                      ! user wants to use the urban canopy model.  Since this code does not
  !                                      ! include the urban canopy model, a TRUE value of this flag will simply
  !                                      ! stop the execution.
  real, external     :: month_d          ! External function (follows this main program):  given an array (dimension 12)
  !                                      ! representing monthly values for some parameter, return a value for 
  !                                      ! a specified date.
  real                :: CZIL            ! Zilitinkevich constant, read from GENPARM.TBL and used to compute surface
  !                                      ! exchange coefficients

  real, dimension (1:ime,1:jme) :: LONGWAVE  ! Longwave radiation as read from the forcing data, which is immediately 
  !                                          ! adjusted (by the emissivity factor) to set variable LWDN.

  integer :: i,j,k, kt
  !
  ! Get the command-line arguments
  !

  call getarg(0, executable_name)
  call getarg(1, forcing_filename)
  call getarg(2, FPEpara_filename)

  if (forcing_filename == " ") then
     write(*,'(/," ***** Problem:  Program expects a command-line argument *****")')
     write(*,'(" ***** Please specify the forcing filename on the command-line.")')
     write(*,'(" ***** E.g.:  ''",A,1x,A,"''",/)') trim(executable_name), "sodanskyla.dat"
     stop ":  ERROR EXIT"
  endif

  if (FPEpara_filename == " ") then
     write(*,'(/," ***** Problem:  Program expects a command-line argument *****")')
     write(*,'(" ***** Please specify the FPE parameters filename on the command-line.")')
     write(*,'(" ***** E.g.:  ''",A,1x,A,"''",/)') trim(executable_name), "sodanskyla.dat FPE.para"
     stop ":  ERROR EXIT"
  endif
  !
  ! Some defaults
  !

  snotime1  = 0.0
  RIBB      = 0.0

     DO I = 1,IME
       DO J = 1,JME

  sheat(i,j) = badval
  etakin(i,j) = badval
  fdown(i,j) = badval
  ec(i,j) = badval
  edir(i,j) = badval
  ett(i,j) = badval
  esnow(i,j) = badval
  drip(i,j) = badval
  dew(i,j) = badval
  beta=badval
  t1(i,j) = badval
  snowh(i,j) = badval
  sneqv(i,j) = badval
  etp(i,j) = badval
  ssoil(i,j) = badval
  flx1 = badval
  flx2 = badval
  flx3 = badval
  snomlt(i,j) = badval
  sncovr(i,j) = badval
  runoff1(i,j) = badval
  runoff2(i,j) = badval
  runoff3 = badval
  rc = badval
  pc(i,j) = badval
  rcs = badval
  rct = badval
  rcq = badval
  rcsoil = badval
  soilw(i,j) = badval
  soilm(i,j) = badval
  qsfc(i,j) = badval
  q1(i,j) = badval
!  wilt(i,j) = badval
!  qmin(i,j) = badval
!  ref(i,j) = badval
!  qmax(i,j) = badval
  rsmin = badval
!  nroot = -999999

soilt1 = 263.69089
tsnav = -9.4591
!tsnav = -9.309082
sfcevp = 0.0! - run total evaporation flux - kg/m^2
snowfallac =0.0 ! - run total snowfall accumulation [mm]
acsnow =0.0 ! - run total frozen precipitation

  !
  ! Read initial conditions
  !

  ! NSOIL            -- Number of soil levels
  ! STARTDATE        -- Starting date ("YYYYMMDDHHmm") of the data in the file
  ! ENDDATE          -- Ending date ("YYYYMMDDHHmm") of the data in the file
  ! LOOP_FOR_A_WHILE -- Whether to loop the same year over and over
  ! LATITUDE         -- Degrees N
  ! LONGITUDE        -- Degrees E
  ! FORCING_TIMESTEP -- Time interval (s) between data records (s) in the forcing file
  ! RUCLSM_TIMESTEP  -- Time step (s) for the RUC LSM integration
  ! ICE              -- Whether this is a sea-ice point (ICE==1) a glacial land ice point (ICE==-1) or a non-glacial land point (ICE==0)
  ! t1               -- Skin temperature (K)
  ! STC              -- Soil temperatures in the soil layers (K).  A pointer array allocated within subroutine open_forcing_file
  ! SMC              -- Soil moisture in the soil layers (m3 m{-3}). A pointer array allocated within subroutine open_forcing_file
  ! SH2O             -- Liquid soil moisture content (m3 m{-3}).  A pointer array allocated within subroutine open_forcing_file
  ! SLDPTH           -- The thicknesses of each soil layer.  A pointer array allocated within subroutine open_forcing_file
  ! CMC              -- Canopy moisture content (kg m-2)
  ! SNEQV            -- Water equivalent accumulated snow depth (m)
  ! TBOT             -- Deep soil temperature (K), a time invariant value
  ! VEGTYP           -- Vegetation category
  ! SOILTYP          -- Soil category
  ! SLOPETYP         -- Slope category
  ! SNOALB           -- Maximum snow albedo -- the albedo of the point when covered by deep snow
  ! ZLVL             -- The level (m AGL) at which the atmospheric thermodynamic forcing fields are considered to be valid
  ! ZLVL_WIND        -- The level (m AGL) at which the atmospheric momentum forcing fields are considered to be valid
  ! ALBEDO_MONTHLY   -- Appropriate background (i.e., snow-free)albedo values for each month of the year
  ! SHDFAC_MONTHLY   -- Appropriate green vegetation fraction values for each month of the year
  ! Z0BRD_MONTHLY    -- Appropriate background (i.e, snow-free) roughness-length values for each month of the year
  ! LAI_MONTHLY      -- Appropriate Leaf Area Index values for each month of the year
  ! USE_URBAN_MODULE -- Whether to call the Urban Canopy Model.  Must be .FALSE. for this code.
  ! ISURBAN          -- Vegetation index that refers to the urban category in the selected landuse dataset
  ! SHDMIN           -- Minimum green vegetation fraction through the year
  ! SHDMAX           -- Maximum green vegetation fraction through the year
  ! USEMONALB        -- Whether to use the provided monthly albedo values
  ! RDLAI2D          -- Whether to use the provided monthly LAI values
  ! LLANDUSE         -- Landuse dataset; either "USGS" or "MODIFIED_IGBP_MODIS_NOAH"


  !
  ! Open the forcing file, and read some metadata and the initial conditions.
  !)

  call open_forcing_file(iunit, forcing_filename, infotext, nsoil, startdate, enddate, loop_for_a_while,    &
       latitude, longitude, &
       forcing_timestep, ruclsm_timestep, ice(i,j), t1(i,j), stc, smc, &
       sh2o, sldpth, cmc(i,j), snowh(i,j), &
       sneqv(i,j), tbot(i,j),vegtyp(i,j), soiltyp(i,j), slopetyp(i,j), snoalb(i,j), zlvl(i,j), zlvl_wind,        &
       albedo_monthly, shdfac_monthly,z0brd_monthly, lai_monthly, use_urban_module, isurban, shdmin(i,j),   &
       shdmax(i,j), usemonalb, rdlai2d, llanduse,iz0tlnd, sfcdif_option)

  CALL open_FPEpara_file(FPEpara_filename)

  allocate( soilmois (1:ime, 1:nsoil ,1:jme) )
  allocate( sliqw (1:ime, 1:nsoil ,1:jme) )
  allocate( smfr3d (1:ime, 1:nsoil ,1:jme) )
  allocate( keepfr3dflag (1:ime, 1:nsoil ,1:jme) )
  allocate( tslb (1:ime, 1:nsoil ,1:jme) )
  allocate( smcmass (1:ime, 1:nsoil ,1:jme) )
  allocate( smliqfrac (1:ime, 1:nsoil ,1:jme) )
  allocate( smfrfrac(1:ime, 1:nsoil ,1:jme) )
  allocate( zsmain(1:nsoil))
  allocate( zshalf(1:nsoil))

  keepfr3dflag=0

       do k=1,nsoil
        soilmois(i,k,j)=smc(k)
        sliqw(i,k,j)=sh2o(k)
        smfr3d(i,k,j)=(smc(k)-sh2o(k))/0.9
        tslb(i,k,j)=stc(k)
       enddo

        sncovr(i,j)=MIN(1.,sneqv(i,j)/0.016)

!tgs - For RUC LSM set use_urban_module = .false. (?)

  dt = real(ruclsm_timestep)

  if (use_urban_module) STOP "This is not urban code."

  !
  ! Allocate additonal arrays (dimensioned by the number of soil levels) which we will need for LSMRUC.
  !

  allocate( et ( nsoil ) )
  et = -1.E36

  allocate( smav ( nsoil ) )
  smav = -1.E36

  !
  ! Set up some input variables for LSMRUC.
  !

  !
  ! LLANDUSE:  Currently only the USGS vegetation dataset as used in WRF is supported.
  !

  LLANDUSE = "USGS-RUC"

  !
  ! LSOIL:  Currently, only the STATSGO soil dataset as used in WRF is supported.
  !

  LSOIL = "STAS-RUC"

  !
  ! Read our lookup tables and parameter tables:  VEGPARM.TBL, SOILPARM.TBL, GENPARM.TBL
  !

  call ruclsm_soilvegprm( LLANDUSE, LSOIL )

  !
  ! COSZ is unused if we're not using the urban canopy model.  If we implement the
  ! urban canopy model for this simple point driver, we will need to compute a COSZ
  ! somewhere.
  !

  COSZ = badval

  !
  ! RCPRAIN is unused.
  !

  PRCPRAIN(i,j) = badval

  !
  ! SOLARDIRECT is unused.
  !

  SOLARDIRECT(i,j) = badval

  !
  ! Set EMISSI for our first time step.  Just a guess, but it's only for the 
  ! first step.  Later time steps get EMISSI from what was set in the prior
  ! time step by LSMRUC.
  !

!tgs! in RUC LSM emissivity is defined from the VEGPARM.TBL
  EMISSI(i,j) = 0.95

  !
  ! For the initial ALBEDO value used in computing SOLNET, just use our
  ! snow-free value.  Subsequent timesteps use the value computed in the
  ! previous call to LSMRUC:
  !

  ALBBCK(i,j) = month_d(albedo_monthly, startdate)
  ALBEDO(i,j) = month_d(albedo_monthly, startdate)

  !
  ! For the initial value of Z0 (used in SFCDIF_off to compute CH and CM),
  ! just use a snow-free background value.  Subsequent timesteps use this
  ! value as computed in the previous call to LSMRUC:
  !

!tgs! in RUC LSM z0 is defined from VEGPARM.TBL
  Z0(i,j) = month_d(z0brd_monthly, startdate)

  !
  ! Z0BRD is computed within LSMRUC.  But we need an initial value, so the call to 
  ! SFCDIF_MYJ can do its thing.  Subsequent timesteps will recycle the Z0BRD 
  ! value as returned from LSMRUC in the previous timestep.
  !

  if ( sfcdif_option == 1 ) then
     z0brd  = z0
  else
     z0brd = badval
  endif

  !
  ! CZIL is needed for the SFCDIF_OFF step.  This comes from CZIL_DATA, as read
  ! from the GENPARM.TBL file, which is how REDPRM ultimately gets it as well:
  !

  CZIL = CZIL_DATA

  !
  !  CM and CH, computed in subroutine SFCDIF_OFF, need initial values.  Values are
  !  subsequently updated for each time step.  So, just take a guess at reasonable
  !  initial values:
  !

  CH(i,j) = 1.E-4
  CM(i,j) = 1.E-4

  if ( sfcdif_option == 1 ) then
     call MYJSFCINIT()
  endif

  !
  ! Enter time loop:
  !

  nowdate = startdate
  ktime = 0

  TIMELOOP : do while ( nowdate <= enddate) 

     !
     ! Increment our counter KTIME and our time variable NOWDATE
     !

     call geth_newdate(nowdate, startdate, ktime*(ruclsm_timestep/60))
     ktime = ktime + 1
  
     print *,'nowdate=',nowdate,'enddate=',enddate
     print *,'loop_for_a_while=',loop_for_a_while

     !
     ! Check if we need to cycle back to our starting data
     !

     if ((loop_for_a_while) .and. (nowdate == enddate)) then
        print*, 'Nowdate: '//nowdate//"  Switching to startdate: "//startdate
        call output_close()
   
        print *,'nowdate=',nowdate
        write(teststr, '(I4.4)') reloop_count
        print *,'reloop_count=',reloop_count,'teststr=',teststr
        cmd =  "mv OUTPUT.nc OUTPUT."//trim(teststr)//".nc"
        reloop_count = reloop_count + 1
        call system(trim(cmd))
   
        ktime = 1
        call read_forcing_text(iunit, nowdate, forcing_timestep,      &
          sfcspd(i,j), sfcu, sfcv, tabs(i,j), q2(i,j),             &
          sfcprs(i,j), soldn(i,j), longwave(i,j), prcp(i,j), snfall(i,j),kt,ierr)
   
        if (ierr == 0) stop "Wrong input for looping a year."
        nowdate = startdate
        rewind(iunit)

     endif


     !
     ! Read the forcing fields, updated from external data every time step:
     ! SFCSPD, SFCU, SFCV, TABS, Q2, SFCPRS, SOLDN, LONGWAVE, PRCP.
     !

     print *,'nowdate=',nowdate, ' at reloop_count=',reloop_count
     call read_forcing_text(iunit, nowdate, forcing_timestep,      &
          sfcspd(i,j), sfcu, sfcv, TABS(i,j), q2(i,j), sfcprs(i,j),&
          soldn(i,j), longwave(i,j), prcp(i,j), snfall(i,j),kt,ierr)
     if (ierr /= 0) then
        exit TIMELOOP
        stop ":  FORCING DATA READ PROBLEM"
     endif

     !
     ! Update FFROZP for each time step, depending on the air temperature in the forcing data.
     ! FFROZP indicates the fraction of the total precipitation which is considered to be
     ! frozen.
     !
!tgs - any possibility to get mixed precipitation from the forcing?
     FFROZP(i,j) = 0.0
     if ( (PRCP(i,j)+snfall(i,j)) > 0.) then
        FFROZP(i,j)=snfall(i,j)/(snfall(i,j)+prcp(i,j))
     endif

!     if ( (PRCP(i,j) > 0) .and. (TABS(i,j) < 273.15) ) then
!        FFROZP(i,j) = 1.0
!     else
!        FFROZP(i,j) = 0.0
!     endif

     !
     ! At each time step, using the forcing fields (and t1, the skin temperature, which
     ! gets updated by LSMRUC), we need to compute a few additional thermodynamic variables.
     ! Ultimately, TH2, Q2SAT and DQSDT2 get passed to SFLX;
     !             t1V and TH2V get used in SFCDIF_off but are not used by LSMRUC.
     !             RHO is used in LSMRUC.
     !

     CALL CALTMP(t1(i,j), TABS(i,j), SFCPRS(i,j), ZLVL(i,j), Q2(i,j), TH2(i,j), t1V(i,j), TH2V(i,j), RHO(i,j)) ! Returns TH2, t1V, TH2V, RHO
     CALL CALHUM(TABS(i,j), SFCPRS(i,j), Q2SAT(i,j), DQSDT2(i,j)) ! Returns Q2SAT, DQSDT2

     !
     ! If the USEMONALB flag is .TRUE., we want to provide ALB from the user-specified
     ! trend through the year, rather than let LSMRUC calculate it for us.
     !
     if (USEMONALB) then
        albbck(i,j) = month_d(albedo_monthly, nowdate)
     else
        albbck(i,j) = 0.18
     endif

     Z0(i,j) = month_d(z0brd_monthly, nowdate)
     if(sncovr(i,j) > 0.99) Z0(i,j) = 0.011

     if ( sfcdif_option == 1 ) then
        z0brd  = z0
     else
        z0brd = badval
     endif

     !
     ! If the RDLAI2D flag is .TRUE., we want to provide XLAI from the user-specified
     ! trend through the year, rather than let LSMRUC calculate it for us.
     !

     if (RDLAI2D) then
        xlai(i,j) = month_d(lai_monthly, nowdate)
     else
        xlai(i,j) = badval
     endif

     !
     ! SHDFAC comes from the user-specified trend through the year.  No other option
     ! at the moment
     !

     shdfac(i,j) = month_d(shdfac_monthly, nowdate)
     shdmin(i,j) = minval(shdfac_monthly)
     shdmax(i,j) = maxval(shdfac_monthly)

     !
     ! qsfc is computed within LSMRUC.  But we need an initial value, so the call to 
     ! SFCDIF_MYJ can do its thing.  Subsequent timesteps will recycle the qsfc
     ! value from as returned from LSMRUC in the previous time step.
     !

!     if (qsfc(i,j) == badval) then
!        qsfc(i,j) = q2(i,j)
!     endif

     if (q1(i,j) == badval) then
        q1(i,j) = 1.6286864E-03
        qsfc(1,1) = q1(i,j)/(1.+q1(i,j))
     endif

     !
     ! SFCDIF_OFF computes mixing lengths for momentum and heat, CM and CH.
     ! Z0 is needed for SFCDIF_OFF.  We use the Z0 as computed in the previous
     ! timestep of LSMRUC, but on the first time step, we need a value of Z0.  This
     ! is set above from our background value.  The initial value may not be quite
     ! what we want, but for that one timestep, it should be OK.  Additionally,
     ! CH and CM need some values for the initial timestep.  These values are
     ! set above.
     !

     if ( SFCDIF_OPTION == 0 ) then

!        CALL SFCDIF_OFF ( ZLVL, ZLVL_WIND , Z0(i,j) , t1V(i,j) , TH2V(i,j) , SFCSPD(i,j) , CZIL , CM(i,j) , CH(i,j) , &
!             VEGTYP(i,j) , ISURBAN , IZ0TLND ) ! Out:  CM, CH

     else if ( SFCDIF_OPTION == 1 ) then

        CALL SFCDIF_MYJ ( ZLVL(i,j), ZLVL_WIND , Z0(i,j) , Z0BRD(i,j) , SFCPRS(i,j) , t1(i,j) , TABS(i,j) , qsfc(i,j) , &
             Q2(i,j) , SFCSPD(i,j) , CZIL , RIBB , CM(i,j) , CH(i,j) , VEGTYP(i,j) , ISURBAN , IZ0TLND )

     endif

     !
     ! SOLNET is an additional forcing field, created by applying the albedo to SOLDN.
     ! ALBEDO is returned each time step from LSMRUC.  The initial value is perhaps
     ! not quite what we want, but each subsequent timestep should be OK.
     !

     SOLNET(i,j) = SOLDN(i,j) * (1.0-ALBEDO(i,j))

     !
     ! Apply the emissivity factor to the given longwave radiation.
     ! This takes the EMISSI value from the previous time step, except
     ! for the first time through the loop, when EMISSI is set above.
     !

     LWDN(i,j) = LONGWAVE(i,j) * EMISSI(i,j)

     !
!tgs As default, let's use FRPCPN == .false.. Reset FRPCPN to .true. if there is data
! about amounts of different precipitation types (rain, snow, graupel, ice, etc.)

     frpcpn = .true.

     !  Call the RUC g LSM routine for a single time step.
     !

     !
     ! Input:
     !
     !    FFROZP      -- Fraction of total precipitation which is frozen ( Fraction [0.0-1.0] )
     !    FRPCPN      -- .true. if mixed phase precipitation available
     !    ICE         -- Land point (ICE==0) or sea-ice point (ICE==1) ( Integer flag 0 or 1 )
     !    ISURBAN     -- The vegetation category for Urban points
     !    DT          -- Time step ( seconds )
     !    ZLVL        -- Height of atmospheric forcing variables ( m AGL )
     !    NSOIL       -- Number of soil layers ( count )
     !    SLDPTH      -- Thickness of each soil layer ( m )
     !    LOCAL       -- Logical flag, .TRUE. to use table values for ALBBCK, SHDFAC, and Z0BRD
     !                   .FALSE. to use values for ALBBCK, SHDFAC, and Z0BRD as set in this driver routine
     !    LLANDUSE    -- Land-use dataset we're using.  "USGS" is the only dataset supported
     !    LSOIL       -- Soil dataset we're using.  "STAS" (for STATSGO) is the only dataset supported
     !    LWDN        -- Longwave downward radiation flux ( W m{-2} )
     !    SOLDN       -- Shortwave downward radiation flux ( W m{-2} )
     !    SOLNET      -- Shortwave net radiation flux ( W m{-2} )
     !    SFCPRS      -- Atmospheric pressure at height ZLVL m AGL ( Pa )
     !    PRCP        -- RAIN rate ( kg m{-2} s{-1} )
     !    PSNFALL     -- Snow rate ( kg m{-2} s{-1} )
     !    LSMRUC      -- Air temperature at height ZLVL m AGL ( K )
     !    Q2          -- Atmospheric mixing ratio at height ZLVL m AGL ( kg kg{-1} )
     !    SFCSPD      -- Wind speed at height ZLVL m AGL ( m s{-1} )
     !    COSZ        -- Cosine of the Solar Zenith Angle (unused in LSMRUC)
     !    PRCPRAIN    -- Liquid precipitation rate ( kg m{-2} s{-1} ) (unused)
     !    SOLARDIRECT -- Direct component of downward solar radiation ( W m{-2} ) (unused)
     !    TH2         -- Air potential temperature at height ZLVL m AGL ( K )
     !    Q2SAT       -- Saturation specific humidity at height ZLVL m AGL ( kg kg{-1} )
     !    DQSDT2      -- Slope of the Saturation specific humidity curve at temperature LSMRUC ( kg kg{-1} K{-1} )
     !    VEGTYP      -- Vegetation category ( index )
     !    SOILTYP     -- Soil category ( index )
     !    SLOPETYP    -- Slope category ( index )
     !    SHDFAC      -- Areal fractional coverage of green vegetation ( fraction [0.0-1.0] ).
     !                   SHDFAC will be set by REDPRM if (LOCAL == .TRUE.)
     !    SHDMIN      -- Minimum areal fractional coverage of green vegetation ( fraction [0.0-1.0] )
     !    SHDMAX      -- Maximum areal fractional coverage of green vegetation ( fraction [0.0-1.0] )
     !    ALBBCK      -- Surface background snow-free albedo (fraction [0.0-1.0]).  ALB will
     !                   be set by REDPRM if (LOCAL == .TRUE.).
     !    SNOALB      -- Maximum deep-snow albedo. ( fraction [0.0-1.0] )
     !    TBOT        -- Constant deep-soil temperature ( K )
     !    Z0BRD       -- Background (i.e., without snow-cover effects) roughness length ( M )
     !    Z0          -- Roughness length, including snow-cover effects ( M )
     !    EMBRD       -- Background emissivity (i.e., not including snow-cover effects) ( fraction [0.0-1.0] )
     !
     ! Updated:
     !
     !    EMISSI      -- Emissivity ( fraction )
     !    CMC         -- Canopy moisture content ( kg m{-2} )
     !    t1          -- Skin temperature ( K )
     !    STC         -- Soil temperature at NSOIL levels ( K )
     !    SMC         -- Volumetric soil moisture content at NSOIL levels ( m{3} m{-3} )
     !    SH2O        -- Liquid portion of the volumetric soil moisture content at NSOIL levels ( m{3} m{-3} )
     !    SNOWH       -- Snow depth ( m )
     !    SNEQV       -- Water equivalent snow depth ( m )
     !    ALBEDO      -- Surface albedo, including any snow-cover effects ( Fraction [0.0-1.0] )
     !    CH          -- Surface exchange coefficient for heat and moisture ( m s{-1} )
     !    CM          -- Surface exchange coefficient for momentum, unused in this code ( m s{-1} )
     !    ETA         -- Latent heat flux (evapotranspiration) ( W m{-2} )
     !    SHEAT       -- Sensible heat flux ( W m{-2} )
     !    ETAKIN      -- Latent heat flux (evapotranspiration) ( kg m{-2} s{-1} )
     !    FDOWN       -- Radiation forcing at the surface ( W m{-2} )
     !    EC          -- Latent heat flux component: canopy water evaporation ( W m{-2} )
     !    EDIR        -- Latent heat flux component: direct soil evaporation ( W m{-2} )
     !    ET          -- Latent heat flux component: plant transpiration from each of NSOIL levels ( W m{-2} )
     !    ETT         -- Latent heat flux component: total plant transpiration ( W m{-2} )
     !    ESNOW       -- Latent heat flux component: sublimation from (or deposition to) snowpack ( W m{-2} )
     !    DRIP        -- Precipitation or dew falling through canopy, in excess of canopy holding capacity ( m )
     !    DEW         -- Dewfall (or frostfall for T<273.15) ( m )
     !    BETA        -- Ratio of actual to potential evapotranspiration ( Fraction [0.0-1.0] )
     !    ETP         -- Potential evapotranspiration ( W m{-2} )
     !    SSOIL       -- Soil heat flux ( W m{-2} )
     !    FLX1        -- Latent heat flux from precipitation accumulating as snow ( W m{-2} )
     !    FLX2        -- Latent heat flux from freezing rain converting to ice ( W m{-2} )
     !    FLX3        -- Latent heat flux from melting snow ( W m{-2} )
     !    SNOMLT      -- Snow melt water ( m )
     !    SNCOVR      -- Fractional snow cover ( Fraction [0.0-1.0] )
     !    RUNOFF1     -- Surface runoff, not infiltrating the soil ( m s{-1} )
     !    RUNOFF2     -- Subsurface runoff, drainage out the bottom of the last soil layer ( m s{-1} )
     !    RUNOFF3     -- Internal soil layer runoff ( m s{-1} )
     !    RC          -- Canopy resistance ( s m{-1} )
     !    PC          -- Plant coefficient, where PC * ETP = ETA ( Fraction [0.0-1.0] )
     !    RSMIN       -- Minimum canopy resistance ( s m{-1} )
     !    XLAI        -- Leaf area index ( dimensionless )
     !    RCS         -- Incoming solar RC factor ( dimensionless )
     !    RCT         -- Air temperature RC factor ( dimensionless )
     !    RCQ         -- Atmospheric water vapor deficit RC factor ( dimensionless )
     !    RCSOIL      -- Soil moisture RC factor ( dimensionless )
     !    SOILW       -- Available soil moisture in the root zone ( Fraction [SMCWLT-SMCMAX] )
     !    SOILM       -- Total soil column moisture content, frozen and unfrozen ( m )
     !    q1          -- Effective mixing ratio at the surface ( kg kg{-1} )
     !    SMAV        -- Soil Moisture Availability at each level, fraction between SMCWLT (SMAV=0.0) and SMCMAX (SMAV=1.0)
     !    SMCWLT      -- Wilting point ( m{3} m{-3} )
     !    SMCDRY      -- Dry soil moisture threshold where direct evaporation from the top layer ends ( m{3} m{-3} )
     !    SMCREF      -- Soil moisture threshold where transpiration begins to stress ( m{3} m{-3} )
     !    SMCMAX      -- Porosity, i.e., saturated value of soil moisture ( m{3} m{-3} )
     !    NROOT       -- Number of root layers ( count )
     !    precipfr    -- Time-step frozen precipitation, output from RUC LSM
     !


#define _DEBUG_PRINT_ 0
#if _DEBUG_PRINT_
              print*, "Before LSMRUC time step=",ktime
              print*, 'b  FFROZP = ', FFROZP(i,j)
              print*, 'b  ICE = ', ICE(i,j)
!              print*, 'b  ISURBAN = ', ISURBAN
              print*, 'b  DT = ', DT
              print*, 'b  ZLVL = ', ZLVL
              print*, 'b  NSOIL = ', NSOIL
              print*, 'b  SLDPTH = ', SLDPTH
              print*, 'b  LLANDUSE = ', trim(LLANDUSE)
              print*, 'b  LSOIL = ', trim(LSOIL)
              print*, 'b  LWDN = ', LWDN(i,j)
              print*, 'b  SOLDN = ', SOLDN(i,j)
              print*, 'b  SOLNET = ', SOLNET(i,j)
              print*, 'b  SFCPRS = ', SFCPRS(i,j)
              print*, 'b  RHO = ', RHO(i,j)
              print*, 'b  PRCP = ', PRCP(i,j)
              print*, 'b  SNFALL = ', SNFALL(i,j)
              print*, 'b  TABS = ', TABS(i,j)
              print*, 'b  Q2 = ', Q2(i,j)
              print*, 'b  SFCSPD = ', SFCSPD(i,j)
!              print*, 'b  COSZ = ', COSZ
!              print*, 'b  PRCPRAIN = ', PRCPRAIN(i,j)
!              print*, 'b  SOLARDIRECT = ', SOLARDIRECT(i,j)
              print*, 'b  TH2 = ', TH2(i,j)
              print*, 'b  Q2SAT = ', Q2SAT(i,j)
              print*, 'b  DQSDT2 = ', DQSDT2(i,j)
              print*, 'b  VEGTYP = ', VEGTYP(i,j)
!              print*, 'b  XLAI = ', xlai(i,j)
              print*, 'b  SOILTYP = ', SOILTYP(i,j)
              print*, 'b  SLOPETYP = ', SLOPETYP(i,j)
              print*, 'b  SHDFAC = ', SHDFAC(i,j)
              print*, 'b  SHDMIN = ', SHDMIN(i,j)
              print*, 'b  SHDMAX = ', SHDMAX(i,j)
              print*, 'b  ALBBCK = ', ALBBCK(i,j)
              print*, 'b  SNOALB = ', SNOALB(i,j)
              print*, 'b  TBOT = ', TBOT(i,j)
!              print*, 'b  Z0BRD = ', Z0BRD
              print*, 'b  Z0 = ', Z0(i,j)
              print*, 'b  EMISSI = ', EMISSI(i,j)
!              print*, 'b  EMBRD = ', EMBRD
              print*, 'b  CMC = ', CMC(i,j)
              print*, 'b  t1 = ', t1(i,j)
              print*, 'b  STC = ', STC
              print*, 'b  SMC = ', SMC
              print*, 'b  SH2O = ', SH2O
              print*, 'b  SNOWH = ', SNOWH(i,j)
              print*, 'b  SNEQV = ', SNEQV(i,j)
              print*, 'b  ALBEDO = ', ALBEDO(i,j)
              print*, 'b  CH = ', CH(i,j)
              print*, 'b  CM = ', CM(i,j)
              print*, 'b  ETA = ', ETA(i,j)
              print*, 'b  SHEAT = ', SHEAT(i,j)
              print*, 'b  ETAKIN = ', ETAKIN(i,j)
!              print*, 'b  FDOWN = ', FDOWN(i,j)
              print*, 'b  EC = ', EC(i,j)
              print*, 'b  EDIR = ', EDIR(i,j)
!              print*, 'b  ET = ', ET
              print*, 'b  ETT = ', ETT(i,j)
              print*, 'b  ESNOW = ', ESNOW(i,j)
              print*, 'b  DRIP = ', DRIP(i,j)
              print*, 'b  DEW = ', DEW(i,j)
!              print*, 'b  BETA = ', BETA
!              print*, 'b  ETP = ', ETP(i,j)
              print*, 'b  SSOIL = ', SSOIL(i,j)
!              print*, 'b  FLX1 = ', FLX1
!              print*, 'b  FLX2 = ', FLX2
!              print*, 'b  FLX3 = ', FLX3
              print*, 'b  SNOMLT = ', SNOMLT(i,j)
              print*, 'b  SNCOVR = ', SNCOVR(i,j)
              print*, 'b  RUNOFF1 = ', RUNOFF1(i,j)
              print*, 'b  RUNOFF2 = ', RUNOFF2(i,j)
!              print*, 'b  RUNOFF3 = ', RUNOFF3
!              print*, 'b  RC = ', RC
!              print*, 'b  PC = ', PC(i,j)
!              print*, 'b  RSMIN = ', RSMIN
!              write(*, '("XLAI = ", F7.5)') XLAI
!              print*, 'b  RCS = ', RCS
!              print*, 'b  RCT = ', RCT
!              print*, 'b  RCQ = ', RCQ
!              print*, 'b  RCSOIL = ', RCSOIL
              print*, 'b  SOILW = ', SOILW
              print*, 'b  SOILM = ', SOILM
              print*, 'b  q1 = ', q1(i,j)
#if 0
              print*, 'b  RDLAI2D = ', RDLAI2D
              print*, 'b  USEMONALB = ', USEMONALB
              print*, 'b  SNOTIME1 = ', SNOTIME1
              print*, 'b  RIBB = ', RIBB
#endif
!              print*, 'b  SMCWLT = ', SMCWLT
!              print*, 'b  SMCDRY = ', SMCDRY
!              print*, 'b  SMCREF = ', SMCREF
!              print*, 'b  SMCMAX = ', SMCMAX
!              print*, 'b  NROOT = ', NROOT
#endif

! RUC LSM additional variables
! -----------------------------
! rhosnf - density of frozen precipitation, output from RUC LSM
! graupelncv - rate of graupel fallin on the ground
! snowncv - rate of snowfall
! rainncv - rate of rain
! precipfr - time step frozen precipitation in [mm]
! zlvl is multiplied by 2., because inside RUC LSM surface layer CONFLX=0.5*DZ8W
! qcatm - cloud water at the first atm. level
! chklowq - can be removed
! flqc and flhc could be removed
! landusef - fraction of each land use category
! nlcat (lucats) - number of landuse categories
! mosaic_lu - 0 or 1
! mosaic_soil - 0 or 1
! soilctop - fractions of soil categories in the top portion of soil domain
! nscat - number of soil categories
! qsfc - specific humidity at the surface
! qsg - surface water vapor mixing ratio at satration
! q1 - mixing ratio at the surface
! soilt1 - snow temperature at the interface between 2 snow layers in K (if exists)
! tsnav - average snow temperathre in C
! xland - land/water mask: 1 - land, 2 - water.
! iswater - number of water category in LLANDUSE classification
! isice - number of ice category in LLANDUSE classification
! ice - fraction of sea ice, could be binary 0/1 flag.
! t1 - skin temperature
! sfcevp - run total evaporation flux - kg/m^2
! snowfallac - run total snowfall accumulation [mm]
! acsnow - run total frozen precipitation
! smfr3d - ice in soil
! keepfr3dflag - flag for frozen soil physics: 0. or 1.
! add initialization here 
! temporary set fractional_seaice=0
!       fractional_seaice = 0
! if ( fractional_seaice == 0 ) then
!    xice_threshold = 0.5
! else if ( fractional_seaice == 1 ) then
!    xice_threshold = 0.02
! endif
  
  fractional_seaice = .FALSE.
  if ( fractional_seaice) then
     xice_threshold = 0.02
  else 
     xice_threshold = 0.5
  endif

    xice(i,j) = ice(i,j)

     myj = .false.
  if ( sfcdif_option == 1) then
     myj = .true.
  endif

    xland(i,j) = 1.
    nlcat=lucats
    nscat=slcats
    mosaic_lu=0
    mosaic_soil=0
    local = .true.
    iswater = 16  ! USGS
    isice = 24  ! USGS

! snow - [mm]
    snow(i,j) = sneqv(i,j)*1.e3
    graupelncv(i,j) = 0.
    snowncv(i,j) = snfall(i,j)*dt
    rainncv(i,j) = prcp(i,j)*dt

    rainbl(i,j)=snowncv(i,j)+rainncv(i,j)

!   if(tabs(i,j) > 273.15) then
!    snowncv(i,j) = 0.
!    rainncv(i,j) = prcp(i,j)
!   else
!    snowncv(i,j) = prcp(i,j)
!    rainncv(i,j) = 0.
!   endif
    mavail(i,j) = soilmois(i,1,j)/maxsmc(nscat)
    landusef(1,1:28,1) = 0.0
    soilctop(1,1:19,1) = 0.0 
    q2(i,j)=q2(i,j)/(1.-q2(i,j))
    qcatm(i,j)=0.
    hfus(i,j) = 0.
    hfrain(i,j) = 0.
    esnow(i,j) = 0.
    snomlt(i,j) = 0.
    CALL LSMRUC(nowdate,local,dt,ktime,nsoil,                            &
                sldpth,rainbl,snow,snowh,sncovr,ffrozp,frpcpn,   &
                rhosnf,graupelncv,snowncv,rainncv,precipfr,sfcspd,  &
                zlvl_wind, zlvl*2.,sfcprs,tabs,q2,qcatm,rho,     & !p_phy in [pa]
                longwave,solnet,emissi,ch,                       & !flqc,flhc,
                mavail,cmc,shdfac*100.,albedo,z0,                &
                z0brd,snoalb, albbck, xlai,                      &   !new
                llanduse, landusef, nlcat, mosaic_lu,            &
                mosaic_soil, soilctop, nscat,                    &   !new
                qsfc,qsg,q1,qcg,dew,soilt1,tsnav,                &
                tbot,vegtyp,soiltyp,xland,                       &
                iswater,isice,xice,xice_threshold,               &
                soilmois,sliqw,soilm,soilw,                      &
                tslb,t1,sheat,etakin,eta,                        &
                hfus,hfrain,snowcan, &
                edir,ec,ett,et,snflx,fltot,                      &
                runoff1,runoff2,drip,esnow,                      &
                sfcevp,ssoil,snowfallac,acsnow,snomlt,           &
                smfr3d,keepfr3dflag,                             &
                myj,shdmin*100.,shdmax*100.,rdlai2d,             &
                1,1, 1,1, 1,1, 1,1, 1,1, 1,1, 1,1, 1,1, 1,1 )

    sneqv(i,j) = snow(i,j)*1.e-3
    q2(i,j)=q2(i,j)/(1.+q2(i,j))
    lwup(i,j)= emissi(i,j) * STBOLT * (t1(i,j)**4)
    swup(i,j)=soldn(i,j)*albedo(i,j)

 if (snowh(i,j) > 0.)  print *,'ktime, rhosn=',ktime,snow(i,j)/snowh(i,j)
             zsmain(1)=0.
             zshalf(1)=0.
          do k=2,nsoil
             zsmain(k)= sldpth(k)
             zshalf(k)=0.5*(zsmain(k-1) + zsmain(k))
          enddo
          do k=1,nsoil-1
           smcmass(i,k,j)=(zshalf(k+1)-zshalf(k))*soilmois(i,k,j)*1.e3
           smcfrmass=(zshalf(k+1)-zshalf(k))*smfr3d(i,k,j)*1.e3
           smcliqmass=(zshalf(k+1)-zshalf(k))*sliqw(i,k,j)*1.e3
           smliqfrac(i,k,j)=min(1.,smcliqmass/smcmass(i,k,j))
           smfrfrac(i,k,j)=min(1.,smcfrmass/smcmass(i,k,j))
          enddo
           smcmass(i,nsoil,j)=(zsmain(nsoil)-zshalf(nsoil))*soilmois(i,nsoil,j)*1.e3
           smcfrmass=(zsmain(nsoil)-zshalf(nsoil))*smfr3d(i,nsoil,j)*1.e3
           smcliqmass=(zsmain(nsoil)-zshalf(nsoil))*sliqw(i,nsoil,j)*1.e3
           smliqfrac(i,nsoil,j)=min(1.,smcliqmass/smcmass(i,nsoil,j))
           smfrfrac(i,nsoil,j)=min(1.,smcfrmass/smcmass(i,nsoil,j))


#if _DEBUG_PRINT_
              print*, "After LSMRUC"
              print*, 'a  FFROZP = ', FFROZP(i,j)
              print*, 'a  ICE = ', ICE(i,j)
!              print*, 'a  ISURBAN = ', ISURBAN
              print*, 'a  DT = ', DT
              print*, 'a  ZLVL = ', ZLVL
              print*, 'a  NSOIL = ', NSOIL
              print*, 'a  SLDPTH = ', SLDPTH
              print*, 'a  LLANDUSE = ', trim(LLANDUSE)
              print*, 'a  LSOIL = ', trim(LSOIL)
              print*, 'a  LWDN = ', LWDN(i,j)
              print*, 'a  SOLDN = ', SOLDN(i,j)
              print*, 'a  SOLNET = ', SOLNET(i,j)
              print*, 'a  SFCPRS = ', SFCPRS(i,j)
              print*, 'a  RHO = ', RHO(i,j)
              print*, 'a  PRCP = ', PRCP(i,j)
              print*, 'a  SNFALL = ', SNFALL(i,j)
              print*, 'a  TABS = ', TABS(i,j)
              print*, 'a  Q2 = ', Q2(i,j)
              print*, 'a  SFCSPD = ', SFCSPD(i,j)
!              print*, 'a  COSZ = ', COSZ
!              print*, 'a  PRCPRAIN = ', PRCPRAIN(i,j)
!              print*, 'a  SOLARDIRECT = ', SOLARDIRECT
              print*, 'a  TH2 = ', TH2(i,j)
              print*, 'a  Q2SAT = ', Q2SAT(i,j)
              print*, 'a  DQSDT2 = ', DQSDT2(i,j)
              print*, 'a  VEGTYP = ', VEGTYP(i,j)
              print*, 'a  SOILTYP = ', SOILTYP(i,j)
              print*, 'a  SlOPETYP = ', SlOPETYP(i,j)
              print*, 'a  SHDFAC = ', SHDFAC(i,j)
              print*, 'a  SHDMIN = ', SHDMIN(i,j)
              print*, 'a  SHDMAX = ', SHDMAX(i,j)
              print*, 'a  ALBBCK = ', ALBBCK(i,j)
              print*, 'a  SNOALB = ', SNOALB(i,j)
              print*, 'a  TBOT = ', TBOT(i,j)
              print*, 'a  Z0BRD = ', Z0BRD(i,j)
              print*, 'a  Z0 = ', Z0(i,j)
              print*, 'a  EMISSI = ', EMISSI(i,j)
!              print*, 'a  EMBRD = ', EMBRD
              print*, 'a  CMC = ', CMC(i,j)
!              print*, 'a  t1 = ', t1
              print*, 'a  t1 = ', t1(i,j)
              print*, 'a  STC = ', STC
              print*, 'a  SMC = ', SMC
              print*, 'a  SH2O = ', SH2O
              print*, 'a  SNOWH = ', SNOWH(i,j)
              print*, 'a  SNEQV = ', SNEQV(i,j)
              print*, 'a  ALBEDO = ', ALBEDO(i,j)
              print*, 'a  CH = ', CH(i,j)
              print*, 'a  CM = ', CM(i,j)
              print*, 'a  ETA = ', eta(i,j)
              print*, 'a  ET = ', et
              print*, 'a  SHEAT = ', SHEAT(i,j)
              print*, 'a  ETAKIN = ', ETAKIN(i,j)
!              print*, 'a  FDOWN = ', FDOWN(i,j)
              print*, 'a  EC = ', EC(i,j)*2.5e6
              print*, 'a  EDIR = ', EDIR(i,j)*2.5e6
!              print*, 'a  ET = ', ET
              print*, 'a  ETT = ', ETT(i,j)*2.5e6
              print*, 'a  ESNOW = ', ESNOW(i,j)*2.835e6
              print*, 'a  FLTOT = ', fltot(i,j)
              print*, 'a  DRIP = ', DRIP(i,j)
              print*, 'a  DEW = ', DEW(i,j)
!              print*, 'a  BETA = ', BETA
!              print*, 'a  ETP = ', ETP(i,j)
              print*, 'a  SSOIL = ', SSOIL(i,j)
!              print*, 'a  FLX1 = ', FLX1
!              print*, 'a  FLX2 = ', FLX2
!              print*, 'a  FLX3 = ', FLX3
              print*, 'a  SNOMLT = ', SNOMLT(i,j)
              print*, 'a  SNCOVR = ', SNCOVR(i,j)
              print*, 'a  RUNOFF1 = ', RUNOFF1(i,j)
              print*, 'a  RUNOFF2 = ', RUNOFF2(i,j)
!              print*, 'a  RUNOFF3 = ', RUNOFF3
!              print*, 'a  RC = ', RC
!              print*, 'a  PC = ', PC(i,j)
!              print*, 'a  RSMIN = ', RSMIN
              print*, 'a  XLAI = ', XLAI(i,j)
!              print*, 'a  RCS = ', RCS
!              print*, 'a  RCT = ', RCT
!              print*, 'a  RCQ = ', RCQ
!              print*, 'a  RCSOIL = ', RCSOIL
              print*, 'a  SOILW = ', SOILW
              print*, 'a  SOILM = ', SOILM
              print*, 'a  q1 = ', q1(i,j)
              print*, 'a  qsfc = ', qsfc(i,j)
#if 0
              print*, 'a  RDLAI2D = ', RDLAI2D
              print*, 'a  USEMONALB = ', USEMONALB
!              print*, 'a  SNOTIME1 = ', SNOTIME1
!              print*, 'a  RIBB = ', RIBB
#endif
!              print*, 'a  SMCWLT = ', SMCWLT
!              print*, 'a  SMCDRY = ', SMCDRY
!              print*, 'a  SMCREF = ', SMCREF
!              print*, 'a  SMCMAX = ', SMCMAX
!              print*, 'a  NROOT = ', NROOT
#endif

       do k=1,nsoil
        smc(k)=soilmois(i,k,j)
        sh2o(k)=sliqw(i,k,j)
        stc(k)=tslb(i,k,j)
       enddo

     !
     ! Compute some diagnostics for output.
     !

!tgs QFX is already computed inside RUC LSM
!     qfx = edir + ec + ett + esnow

     !
     ! Residual of surface energy balance equation terms
     !

!     f = solnet(i,j) + lwdn(i,j)
!     fup = emissi(i,j) * STBOLT * (t1(i,j)**4)
     res(i,j)=fltot(i,j)
!tgs     res = f - sheat + ssoil - eta - fup - flx1 - flx2 - flx3

     !
     ! Write the output data for this timestep.
     !

     if (ktime == 1) call initialize_output(nsoil, 0, 0, 0, dt, iz0tlnd, sfcdif_option)

     ! Time variable
     call output_time(ktime, nowdate,   "Times",   "UTC time of data output",                           "YYYYMMDD HH:mm"  )

! ESM-SnowMIP output
! Table 2:
! hfds - downward heat flux at ground surface
! hfdsn  - downward heat flux into snowpack
! hfls - surface upward latent heat flux
! hfmlt - energy of fusion solid to liquid
! hfrs - heat transferred to snowpack by rain downward
! hfsbl - energy of sublimation solid to vapour
! hfss - surface upward sensible heat flux
! rlus - surface upwelling longwave radiation upward
! rsus - surface upwelling shortwave radiation upward
     call output_var(ktime, -ssoil(i,j),       "hfds",   "downward heat flux at ground surface",              "W m{-2}"       )
! print *,'ktime,-ssoil(i,j)',ktime,-ssoil(i,j),i,j
  if(snow(i,j) > 0.) then
     call output_var(ktime, snflx(i,j),        "hfdsn",  "downward heat flux into snowpack",                  "W m{-2}"       )
  else
     call output_var(ktime, 0.,                "hfdsn",  "downward heat flux into snowpack",                  "W m{-2}"       )
  endif
     call output_var(ktime, etakin(i,j)*XLV,   "hfls",   "surface upward latent heat flux",                   "W m{-2}"       )
  print*,'hfus(i,j)=',ktime,hfus(i,j)
     call output_var(ktime, hfus(i,j),         "hfmlt",  "energy of fusion",                                  "W m{-2}"       )
     call output_var(ktime, hfrain(i,j),       "hfrs",   "heat transferred to snowpack by rain",              "W m{-2}"       )
     call output_var(ktime, 2.835e6*esnow(i,j),"hfsbl",  "energy of sublimation",                             "W m{-2}"       )
     call output_var(ktime, sheat(i,j),        "hfss",   "surface upward sensible heat flux",                 "W m{-2}"       )
     call output_var(ktime, lwup(i,j),         "rlus",   "surface upwelling longwave radiation",              "W m{-2}"       )
     call output_var(ktime, swup(i,j),         "rsus",   "surface upwelling shortwave radiation",             "W m{-2}"       )
! Table 3 - kg m{-2} s{-1}
! esn - liquid water evaporation from snowpack
! evspsbl - total water vapour flux from the surface to the atmosphere
! evspsblsoi - evaporation and sublimation from soil
! evspsblveg - evaporation and sublimation from canopy
! mrrob - subsurface runoff
! mrros - surface runoff
! sbl - sublimation of snow
! snm - surface snow melt
! snmsl - water flowing out of snowpack
! tran - transpiration
     call output_var(ktime, 999.,              "esn",    "liquid water evaporation from snowpack",     "kg m{-2} s{-1}"      )
     call output_var(ktime, etakin(i,j),       "evspsbl","total water vapour flux from the surface to the atmosphere",  "kg m{-2} s{-1}"      )
     call output_var(ktime, edir(i,j),         "evspsblsoi","evaporation and sublimation from soil",   "kg m{-2} s{-1}"      )
     call output_var(ktime, ec(i,j),           "evspsblveg","evaporation and sublimation from canopy", "kg m{-2} s{-1}"      )
     call output_var(ktime, runoff2(i,j)*1.e-3,"mrrob",  "subsurface runoff",                          "kg m{-2} s{-1}"      )
     call output_var(ktime, runoff1(i,j)*1.e-3,"mrros",  "surface runoff",                             "kg m{-2} s{-1}"      )
 print *,'RUNOFF1=',ktime,runoff1(i,j)*1.e-3
     call output_var(ktime, esnow(i,j),        "sbl",    "sublimation of snow",                        "kg m{-2} s{-1}"      )
     call output_var(ktime, snomlt(i,j)/dt,    "snm",    "surface snow melt",                          "kg m{-2} s{-1}"      )
 print *,'snomlt(i,j)/dt',ktime,snomlt(i,j)/dt
     call output_var(ktime, 999.,              "snmsl",  "water flowing out of snowpack",              "kg m{-2} s{-1}"      )
     call output_var(ktime, ett(i,j),          "tran",   "transpiration",                              "kg m{-2} s{-1}"      )
! Table 3:
! albs surface albedo –
! albsn snow albedo –
! cw total canopy water storage kg m–2
! lqsn mass fraction of liquid water in snowpack –
! lwsnl liquid water content of snowpack kg m–2
! mrfsofr mass fractions of frozen water in soil layers –
! mrlqso mass fractions of unfrozen water in soil layers –
! mrlsl masses of frozen and unfrozen moisture in soil layers kg m–2
! snc snow area fraction –
! snd snowdepth m
! snw mass of snowpack kg m–2
! snwc mass of snow intercepted by vegetation kg m–2
! tcs vegetation canopy temperature K
! tgs temperature of bare soil K
! ts surface temperature K
! tsl temperatures of soil layers K
! tsn snow internal temperature K
! tsns snow surface temperature K
     call output_var(ktime, albbck(i,j),    "albs",  "surface albedo",                             "fraction"      )
     call output_var(ktime, albedo(i,j),    "albsn", "snow albedo",                                "fraction"      )
     call output_var(ktime, cmc(i,j),       "cw",    "total canopy water storage",                 "kg m{–2}"      )
     call output_var(ktime, 999.,           "lqsn",  "mass fraction of liquid water in snowpack",  "fraction"      )
     call output_var(ktime, snow(i,j),      "lwsnl", "liquid water content of snowpack",           "kg m{–2}"      )
     call output_levels(ktime, nsoil, "num_soil_layers", smfrfrac, "mrfsofr","mass fractions of frozen water in soil layers",         "fraction"  )
     call output_levels(ktime, nsoil, "num_soil_layers", smliqfrac,"mrlqso", "mass fractions of unfrozen water in soil layers",       "fraction"  )
     call output_levels(ktime, nsoil, "num_soil_layers", smcmass,"mrlsl",  "masses of frozen and unfrozen moisture in soil layers", "kg m{–2}"  )
     call output_var(ktime, sncovr(i,j),    "snc",   "snow area fraction",                         "fraction"      )
     call output_var(ktime, snowh(i,j),     "snd",   "snowdepth",                                  "m"             )
     call output_var(ktime, snow(i,j),      "snw",   "mass of snowpack",                           "kg m{–2}"      )
 print *,'SWE=',ktime,snow(i,j)
     call output_var(ktime, snowcan(i,j),   "snwc",  "mass of snow intercepted by vegetation",     "kg m{–2}"      )
     call output_var(ktime, 999.,           "tcs",   "vegetation canopy temperature",              "K"             )
     call output_var(ktime, 999.,           "tgs",   "temperature of bare soil",                   "K"             )
     call output_var(ktime, t1(i,j),        "ts",    "surface temperature",                        "K"             )
 print *,'TSK=',ktime,t1(i,j)
     call output_levels(ktime, nsoil, "num_soil_layers", stc,  "tsl",   "temperature of soil layers", "K"          )
     call output_var(ktime, min(273.,soilt1(i,j)),    "tsn",   "snow internal temperature",                  "K"             )
     call output_var(ktime, min(273.,t1(i,j)),        "tsns",  "snow surface temperature",                   "K"             )
! End of ESM-SnowMIP output`

if(1==2) then
     ! Multi-layer variables
     call output_levels(ktime, nsoil, "num_soil_layers", 1.e3*2.501e6*et, "ET",  "Plant transpiration from a particular root layer",  "W m{-2}")

     ! Single-layer variables

     call output_var(ktime, mavail(i,j),    "MAVAIL",  "Top layer Soil Moisture Availability",              "fraction"      )
     call output_var(ktime, soldn(i,j),     "SOLDN",   "Downward short-wave radiation flux at the surface", "W m{-2}"       )
     call output_var(ktime, emissi(i,j),    "EMISSI",  "Emissivity",                                        "fraction"      )
     call output_var(ktime, z0brd(i,j),     "Z0BRD",   "Background roughness length (not including snow-cover effect)", "m" )
     call output_var(ktime, z0(i,j),        "Z0",      "Roughness length (including snow-cover effect)",    "m"             )
     call output_var(ktime, sfcprs(i,j),    "SFCPRS",  "Atmospheric pressure at ZLVL m AGL",                "Pa"            )
     call output_var(ktime, prcp(i,j),      "PRCP",    "Liquid Precipitation rate",                         "kg m{-2} s{-1}")
     call output_var(ktime, snfall(i,j),    "SNFALL",  "Snow rate",                                "kg m{-2} s{-1}")
     call output_var(ktime, TABS(i,j),      "TABS",    "Air temperature at ZLVL m AGL",                     "K"             )
     call output_var(ktime, q2(i,j),        "Q2",      "Mixing ratio at ZLVL m AGL",                        "kg kg{-1}"     )
     call output_var(ktime, sfcspd(i,j),    "SFCSPD",  "Wind speed",                                        "m s{-1}"       )
     call output_var(ktime, snowh(i,j),     "SNOWH",   "Snow depth",                                        "m"             )
     call output_var(ktime, sneqv(i,j),     "SNEQV",   "Liquid equivalent of accumulated snow depth",       "m"             )
     call output_var(ktime, shdfac(i,j),    "SHDFAC",  "shdfac",                                            "fraction"      )
     call output_var(ktime, ch(i,j),        "CH",      "Surface exchange coefficient for heat and moisture","m s{-1}"       )
     call output_var(ktime, cm(i,j),        "CM",      "Surface exchange coefficient for momentum",         "m s{-1}"       )
     call output_var(ktime, res(i,j),       "RES",     "Residual of surface energy balance equation",       "W m{-2}"       )
     call output_var(ktime, etakin(i,j),    "ETAKIN",  "Actual latent heat flux",                           "kg m{-2} s{-1}")
!     call output_var(ktime, qfx(i,j),       "QFX",     "Latent heat flux",                                  "W m{-2}"       )
!     call output_var(ktime, fdown(i,j),     "FDOWN",   "Radiation forcing at the surface",                  "W m{-2}"       )
     call output_var(ktime, 1.e3*2.5e6*ec(i,j),"EC",      "Canopy water evaporation",                          "W m{-2}"       )
     call output_var(ktime, 1.e3*2.5e6*edir(i,j),"EDIR",  "Direct soil evaporation",                           "W m{-2}"       )
     call output_var(ktime, 1.e3*2.5e6*ett(i,j),"ETT",    "Total plant transpiration",                         "W m{-2}"       )
     ! Convert DRIP from m/timestep to kg m{-2} s{-1} (mm/s)
     drip(i,j) = 1.E3 * drip(i,j) / dt
     call output_var(ktime, drip(i,j),      "DRIP",    "Throughfall of precipitation from canopy",          "kg m{-2} s{-1}")
     ! Convert DEW from m s{-1} to kg m{-2} s{-1}
     call output_var(ktime, dew(i,j)*rho(i,j),"DEW",   "Dewfall (or frostfall for T < 273.15)",             "kg m{-2} s{-1}")
!     call output_var(ktime, beta,      "BETA",    "Ratio of actual to potential evaporation",          "dimensionless" )
!     call output_var(ktime, etp(i,j),       "ETP",     "Potential evaporation",                             "W m{-2}"       )
!     call output_var(ktime, flx1,      "FLX1",    "Heat flux from snow surface to accumulating precip","W m{-2}"       )
!     call output_var(ktime, flx2,      "FLX2",    "Freezing rain latent heat flux",                    "W m{-2}"       )
     call output_var(ktime, snomlt(i,j)*1.e-3,"SNOMLT",  "Water equivalent snow melt",                      "m"             )
     call output_var(ktime, sncovr(i,j),    "SNCOVR",    "Fractional snow cover",                           "fraction"      )
     call output_var(ktime, runoff2(i,j),   "UDRUNOFF",  "Underground runoff",                              "m s{-1}"       )
!     call output_var(ktime, rc,        "RC",      "Canopy resistance",                                 "s m{-2}"       )
!     call output_var(ktime, pc(i,j),        "PC",      "Plant coefficient (PC * ETP = TRANSP)",             "fraction"      )
!     call output_var(ktime, rsmin,     "RSMIN",   "Minimum canopy resistance",                         "s m{-1}"       )
     call output_var(ktime, xlai(i,j),      "LAI",     "Leaf Area Index",                                   "dimensionless" )
!     call output_var(ktime, rcs,       "RCS",     "Incoming solar RC factor",                          "dimensionless" )
!     call output_var(ktime, rct,       "RCT",     "Air temperature RC factor",                         "dimensionless" )
!     call output_var(ktime, rcq,       "RCQ",     "Atmospheric vapor deficit RC factor",               "dimensionless" )
!     call output_var(ktime, rcsoil,    "RCSOIL",  "Soil moisture RC factor",                           "dimensionless" )
     call output_var(ktime, soilw(i,j),      "SOILW",  "Available soil moisture in root zone",              "fraction"      )
     call output_var(ktime, soilm(i,j)*1.e-3,"SOILM",  "Total column moisture content",                     "m"             )
     call output_var(ktime, q1(i,j),         "Q1",     "Effective mixing ratio at surface",                 "kg kg{-1}"     )
!     call output_var(ktime, smcwlt,    "SMCWLT",  "Wilting-point soil moisture threshold",             "m{3} m{-3}"    )
!     call output_var(ktime, smcdry,    "SMCDRY",  "Dry-soil soil moisture threshold",                  "m{3} m{-3}"    )
!     call output_var(ktime, smcref,    "SMCREF",  "Soil moisture threshold for transpiration stress",  "m{3} m{-3}"    )
!     call output_var(ktime, smcmax,    "SMCMAX",  "Saturated value of soil moisture (Porosity)",       "m{3} m{-3}"    )
endif
   print *,'finish_output_for_time(ktime)',ktime
     call finish_output_for_time(ktime)

   print *,'nowdate(7:10) =',nowdate(7:10)
     if (nowdate(7:10) == "0101") then
     print '(I10, 3x, A4,"-", A2, "-", A2, " ", A2, ":", A2, 20(2x,F9.4))', ktime, &
          nowdate(1:4), nowdate(5:6), nowdate(7:8), nowdate(9:10), nowdate(11:12),  &
          stc(1:nsoil),smc(1:nsoil),sh2o(1:nsoil)
     endif
!   print*, nowdate
!   print*, 'XXXXXXXXXXX', teststr
!  if (trim(teststr) > '0010') print*, trim(teststr)
!      loop_for_a_while = .true.
!    if (trim(teststr) > '0010') then
!        print *,'Spin-up and Main runs are finished'
!        exit TIMELOOP
!        print *,'Spin-up and Main runs are finished'
!        stop
!    endif

  enddo TIMELOOP


    ENDDO
  ENDDO
  !
  ! Shut down the output streams
  !

  call output_close()

  write(*, FMT=trim(infotext))

  ! All done.

end program simple_driver

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

real function month_d(a12, nowdate) result (nowval)
  !
  ! Given a set of 12 values, taken to be valid on the fifteenth of each month (Jan through Dec)
  ! and a date in the form <YYYYMMDD[HHmmss]> ....
  ! 
  ! Return a value valid for the day given in <nowdate>, as an interpolation from the 12
  ! monthly values.
  !
  use kwm_date_utilities
  implicit none
  real, dimension(12), intent(in) :: a12 ! 12 monthly values, taken to be valid on the 15th of
  !                                      ! the month
  character(len=10), intent(in) :: nowdate ! Date, in the form <YYYYMMDD[HHmmss]>
  integer :: nowy, nowm, nowd
  integer :: prevm, postm
  real    :: factor
  integer, dimension(12) :: ndays = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  !
  ! Handle leap year by setting the number of days in February for the year in question.
  !
  read(nowdate(1:8),'(I4,I2,I2)') nowy, nowm, nowd
  ndays(2) = nfeb(nowy)

  !
  ! Do interpolation between the fifteenth of two successive months.
  !
  if (nowd == 15) then
     nowval = a12(nowm)
     return
  else if (nowd < 15) then
     postm = nowm
     prevm = nowm - 1
     if (prevm == 0) prevm = 12
     factor = real(ndays(prevm)-15+nowd)/real(ndays(prevm))
  else if (nowd > 15) then
     prevm = nowm
     postm = nowm + 1
     if (postm == 13) postm = 1
     factor = real(nowd-15)/real(ndays(prevm))
  endif

  nowval = a12(prevm)*(1.0-factor) + a12(postm)*factor

end function month_d

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
! Adapted from the WRF subroutine in module_sf_ruclsm.F:
!-----------------------------------------------------------------
SUBROUTINE RUCLSM_SOILVEGPRM( MMINLU, MMINSL)
!-----------------------------------------------------------------

  USE module_sf_ruclsm
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: MMINLU, MMINSL
  integer :: LUMATCH, IINDEX, LC, NUM_SLOPE
  integer :: ierr
  INTEGER , PARAMETER :: OPEN_OK = 0

  character*128 :: mess , message
  logical, external :: wrf_dm_on_monitor


!-----SPECIFY VEGETATION RELATED CHARACTERISTICS :
!             ALBBCK: SFC albedo (in percentage)
!       note: recommend to use monthly climatological albedo 
!                 Z0: Roughness length (m)
!            LEMITBL: Emmisivity
!             PCTBL : Plant coefficient
!            SHDTBL : Shading factor
!       note: recommend to use monthly climatological greeness fractioin 
!            IFORTBL: parameter to define rootiung depth
!              RSMIN: Mimimum stomatal resistance (s m-1) - not used
!                RGL: Parameters used in radiation stress function - not used
!                 HS: Parameter used in vapor pressure deficit functio
!               SNUP: Threshold snow depth (in water equivalent m) that - not used
!             LAITBL: Leaf area index (dimensionless)
!       note: recommend to use monthly climatological LAI
!             MAXALB: Upper bound on maximum albedo over deep snow
!       note: recommend to use monthly climatological data of max snow albedo
!               TOPT: Optimum transpiration air temperature. (K) - not used
!             CMCMAX: Maximum canopy water capacity
!             CMXTBL: MAX CNPY Capacity (m)
!             CFACTR: Parameter used in the canopy inteception calculati
!                     implies 100% snow cover - not used
!              RSMAX: Max. stomatal resistance (s m-1) - not used
!               BARE: category number for bare soil
!            NATURAL: category number for grassland
!               CROP: category number for crops
!
!-----READ IN VEGETAION PROPERTIES FROM VEGPARM.TBL
!

  IF ( wrf_dm_on_monitor() ) THEN

     OPEN(19, FILE='VEGPARM.TBL',FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
     IF(ierr .NE. OPEN_OK ) THEN
        WRITE(message,FMT='(A)') &
             'module_sf_ruclsm.F: soil_veg_gen_parm: failure opening VEGPARM.TBL'
        CALL wrf_error_fatal ( message )
     END IF


     LUMATCH=0

     FIND_LUTYPE : DO WHILE (LUMATCH == 0)
        READ (19,*,END=2002)
        READ (19,*,END=2002)LUTYPE
        READ (19,*)LUCATS,IINDEX

        IF(LUTYPE.EQ.MMINLU)THEN
           WRITE( mess , * ) 'LANDUSE TYPE = ' // TRIM ( LUTYPE ) // ' FOUND', LUCATS,' CATEGORIES'
           CALL wrf_message( mess )
           LUMATCH=1
        ELSE
           call wrf_message ( "Skipping over LUTYPE = " // TRIM ( LUTYPE ) )
           DO LC = 1, LUCATS+14
              read(19,*)
           ENDDO
        ENDIF
     ENDDO FIND_LUTYPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
     IF ( SIZE(SHDTBL)       < LUCATS .OR. &
          SIZE(RSTBL)        < LUCATS .OR. &
          SIZE(RGLTBL)       < LUCATS .OR. &
          SIZE(HSTBL)        < LUCATS .OR. &
          SIZE(SNUPTBL)      < LUCATS .OR. &
          SIZE(MAXALB)       < LUCATS .OR. &
          SIZE(LAITBL)       < LUCATS .OR. & 
          SIZE(Z0TBL)        < LUCATS .OR. &
          SIZE(ALBTBL)       < LUCATS .OR. &
          SIZE(LEMITBL )     < LUCATS ) THEN
        CALL wrf_error_fatal('Table sizes too small for value of LUCATS in module_sf_ruclsm.F')
     ENDIF

     IF(LUTYPE.EQ.MMINLU)THEN
        DO LC=1,LUCATS
           READ (19,*)IINDEX,ALBTBL(LC),Z0TBL(LC),         &
                 LEMITBL(LC),PCTBL(LC), SHDTBL(LC),        &
                 IFORTBL(LC),RSTBL(LC),RGLTBL(LC),         &
                 HSTBL(LC),SNUPTBL(LC),LAITBL(LC),MAXALB(LC)

        ENDDO
!
        READ (19,*)
        READ (19,*)TOPT_DATA
        READ (19,*)
        READ (19,*)CMCMAX_DATA
        READ (19,*)
        READ (19,*)CFACTR_DATA
        READ (19,*)
        READ (19,*)RSMAX_DATA
        READ (19,*)
        READ (19,*)BARE
        READ (19,*)
        READ (19,*)NATURAL
        READ (19,*) 
        READ (19,*)CROP
     ENDIF
!
2002 CONTINUE

     CLOSE (19)
     IF (LUMATCH == 0) then
        CALL wrf_error_fatal ("Land Use Dataset '"//MMINLU//"' not found in VEGPARM.TBL.")
     ENDIF
  ENDIF

#ifndef _SIMPLE_DRIVER_
      CALL wrf_dm_bcast_string  ( LUTYPE  , 8 )
      CALL wrf_dm_bcast_integer ( LUCATS  , 1 )
      CALL wrf_dm_bcast_integer ( IINDEX  , 1 )
      CALL wrf_dm_bcast_integer ( LUMATCH , 1 )
      CALL wrf_dm_bcast_real    ( ALBTBL  , NLUS )
      CALL wrf_dm_bcast_real    ( Z0TBL   , NLUS )
      CALL wrf_dm_bcast_real    ( LEMITBL , NLUS )
      CALL wrf_dm_bcast_real    ( PCTBL   , NLUS )
      CALL wrf_dm_bcast_real    ( SHDTBL  , NLUS )
      CALL wrf_dm_bcast_real    ( IFORTBL , NLUS )
      CALL wrf_dm_bcast_real    ( RSTBL   , NLUS )
      CALL wrf_dm_bcast_real    ( RGLTBL  , NLUS )
      CALL wrf_dm_bcast_real    ( HSTBL   , NLUS )
      CALL wrf_dm_bcast_real    ( SNUPTBL , NLUS )
      CALL wrf_dm_bcast_real    ( LAITBL  , NLUS )
      CALL wrf_dm_bcast_real    ( MAXALB  , NLUS )
      CALL wrf_dm_bcast_real    ( TOPT_DATA    , 1 )
      CALL wrf_dm_bcast_real    ( CMCMAX_DATA  , 1 )
      CALL wrf_dm_bcast_real    ( CFACTR_DATA  , 1 )
      CALL wrf_dm_bcast_real    ( RSMAX_DATA  , 1 )
      CALL wrf_dm_bcast_integer ( BARE    , 1 )
      CALL wrf_dm_bcast_integer ( NATURAL    , 1 )
      CALL wrf_dm_bcast_integer ( CROP    , 1 )
#endif

!
!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
!
  IF ( wrf_dm_on_monitor() ) THEN
     OPEN(19, FILE='SOILPARM.TBL',FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
     IF(ierr .NE. OPEN_OK ) THEN
        WRITE(message,FMT='(A)') &
             'module_sf_ruclsm.F: soil_veg_gen_parm: failure opening SOILPARM.TBL'
        CALL wrf_error_fatal ( message )
     END IF

     WRITE(mess,*) 'INPUT SOIL TEXTURE CLASSIFICAION = ', TRIM ( MMINSL )
     CALL wrf_message( mess )

     LUMATCH=0

     READ (19,*)
     READ (19,2000,END=2003)SLTYPE
2000 FORMAT (A8)
     READ (19,*)SLCATS,IINDEX
        IF(SLTYPE.NE.MMINSL)THEN
          DO LC=1,SLCATS
              READ (19,*) IINDEX,BB(LC),DRYSMC(LC),HC(LC),MAXSMC(LC),&
                        REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
                        WLTSMC(LC), QTZ(LC)
          ENDDO
        ENDIF

     READ (19,*)
     READ (19,2000,END=2003)SLTYPE
     READ (19,*)SLCATS,IINDEX

print *,'SLTYPE=',SLTYPE,'slcats=',slcats
 
     IF(SLTYPE.EQ.MMINSL)THEN
        WRITE( mess , * ) 'SOIL TEXTURE CLASSIFICATION = ', TRIM ( SLTYPE ) , ' FOUND', &
             SLCATS,' CATEGORIES'
        CALL wrf_message ( mess )
        LUMATCH=1
     ENDIF
!   ENDDO FIND_LUTYPE

      print *,'soil lumatch=',lumatch
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
     IF ( SIZE(BB    ) < SLCATS .OR. &
          SIZE(DRYSMC) < SLCATS .OR. &
          SIZE(HC   )  < SLCATS .OR. &
          SIZE(MAXSMC) < SLCATS .OR. &
          SIZE(REFSMC) < SLCATS .OR. &
          SIZE(SATPSI) < SLCATS .OR. &
          SIZE(SATDK ) < SLCATS .OR. &
          SIZE(SATDW ) < SLCATS .OR. &
          SIZE(WLTSMC) < SLCATS .OR. &
          SIZE(QTZ   ) < SLCATS  ) THEN
        CALL wrf_error_fatal('Table sizes too small for value of SLCATS in module_sf_rucdrv.F')
     ENDIF
     IF(SLTYPE.EQ.MMINSL)THEN
        DO LC=1,SLCATS
           READ (19,*) IINDEX,BB(LC),DRYSMC(LC),HC(LC),MAXSMC(LC),   &
                       REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
                       WLTSMC(LC), QTZ(LC)
        ENDDO
     ENDIF

2003 CONTINUE

     CLOSE (19)
  ENDIF

#ifndef _SIMPLE_DRIVER_
  CALL wrf_dm_bcast_integer ( LUMATCH , 1 )
  CALL wrf_dm_bcast_string  ( SLTYPE  , 8 )
  CALL wrf_dm_bcast_string  ( MMINSL  , 8 )  ! since this is reset above, see oct2 ^
  CALL wrf_dm_bcast_integer ( SLCATS  , 1 )
  CALL wrf_dm_bcast_integer ( IINDEX  , 1 )
  CALL wrf_dm_bcast_real    ( BB      , NSLTYPE )
  CALL wrf_dm_bcast_real    ( DRYSMC  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( HC      , NSLTYPE )
  CALL wrf_dm_bcast_real    ( MAXSMC  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( REFSMC  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( SATPSI  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( SATDK   , NSLTYPE )
  CALL wrf_dm_bcast_real    ( SATDW   , NSLTYPE )
  CALL wrf_dm_bcast_real    ( WLTSMC  , NSLTYPE )
  CALL wrf_dm_bcast_real    ( QTZ     , NSLTYPE )
#endif

  IF(LUMATCH.EQ.0)THEN
     CALL wrf_message( 'SOIl TEXTURE IN INPUT FILE DOES NOT ' )
     CALL wrf_message( 'MATCH SOILPARM TABLE'                 )
     CALL wrf_error_fatal ( 'INCONSISTENT OR MISSING SOILPARM FILE' )
  ENDIF

!
!-----READ IN GENERAL PARAMETERS FROM GENPARM.TBL
!
  IF ( wrf_dm_on_monitor() ) THEN
     OPEN(19, FILE='GENPARM.TBL',FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
     IF(ierr .NE. OPEN_OK ) THEN
        WRITE(message,FMT='(A)') &
             'module_sf_ruclsm.F: soil_veg_gen_parm: failure opening GENPARM.TBL'
        CALL wrf_error_fatal ( message )
     END IF

     READ (19,*)
     READ (19,*)
     READ (19,*) NUM_SLOPE

     SLPCATS=NUM_SLOPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
     IF ( SIZE(slope_data) < NUM_SLOPE ) THEN
        CALL wrf_error_fatal('NUM_SLOPE too large for slope_data array in module_sf_ruclsm')
     ENDIF

     DO LC=1,SLPCATS
        READ (19,*)SLOPE_DATA(LC)
     ENDDO

     READ (19,*)
     READ (19,*)SBETA_DATA
     READ (19,*)
     READ (19,*)FXEXP_DATA
     READ (19,*)
     READ (19,*)CSOIL_DATA
     READ (19,*)
     READ (19,*)SALP_DATA
     READ (19,*)
     READ (19,*)REFDK_DATA
     READ (19,*)
     READ (19,*)REFKDT_DATA
     READ (19,*)
     READ (19,*)FRZK_DATA
     READ (19,*)
     READ (19,*)ZBOT_DATA
     READ (19,*)
     READ (19,*)CZIL_DATA
     READ (19,*)
     READ (19,*)SMLOW_DATA
     READ (19,*)
     READ (19,*)SMHIGH_DATA
     READ (19,*)
     READ (19,*)LVCOEF_DATA
     CLOSE (19)
  ENDIF

  CALL wrf_dm_bcast_integer ( NUM_SLOPE    ,  1 )
  CALL wrf_dm_bcast_integer ( SLPCATS      ,  1 )
  CALL wrf_dm_bcast_real    ( SLOPE_DATA   ,  NSLOPE )
  CALL wrf_dm_bcast_real    ( SBETA_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( FXEXP_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( CSOIL_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( SALP_DATA    ,  1 )
  CALL wrf_dm_bcast_real    ( REFDK_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( REFKDT_DATA  ,  1 )
  CALL wrf_dm_bcast_real    ( FRZK_DATA    ,  1 )
  CALL wrf_dm_bcast_real    ( ZBOT_DATA    ,  1 )
  CALL wrf_dm_bcast_real    ( CZIL_DATA    ,  1 )
  CALL wrf_dm_bcast_real    ( SMLOW_DATA   ,  1 )
  CALL wrf_dm_bcast_real    ( SMHIGH_DATA  ,  1 )
  CALL wrf_dm_bcast_real    ( LVCOEF_DATA  ,  1 )


!-----------------------------------------------------------------
END SUBROUTINE RUCLSM_SOILVEGPRM
!-----------------------------------------------------------------
