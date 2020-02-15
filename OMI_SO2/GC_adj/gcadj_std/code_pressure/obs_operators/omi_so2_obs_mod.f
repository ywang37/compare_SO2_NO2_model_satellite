! $Id: omi_so2_obs_mod.f
      MODULE OMI_SO2_OBS_MOD
!
!*****************************************************************************
! MODULE OMI_SO2_OBS_MOD contians subroutines necessary to
! 1. Read OMI L2 swath SO2 observations
! 2. Compute OMI-GEOS-Chem difference, cost function, and adjoint
! forcing
!
!    (ywang, 07/16/04)
!
!  Module Variables:
!  ===========================================================================
!
!  Module Routines:
!  ===========================================================================
!  ( 1) CALC_OMI_SO2_FORCE
!  ( 2) CHECK                : Check status for calling netCDF
!  ( 3) INIT_OMI_SO2_OBS     : Initialize OMI SO2 observation operator
!  ( 4) WRITE_GC_OMI_SO2_OBS
!  ( 5) READ_GC_OMI_SO2_OBS
!  ( 6) MAKE_CURRENT_OMI_SO2 :
!  ( 7) MAKE_AVERAGE_OMI_SO2 :
!
!  ===========================================================================
!  NOTES:
!  ( 1) The original OMI L2 swath SO2 HDFEOS5 data are preprocessed by
!  IDL code OMI_SO2_L2_preprocess.pro. Quality control is done through
!  the IDL code.
!  ( 2) Local air mass factors are considered. (ywang, 08/04/15)
!
!*****************************************************************************
!

      IMPLICIT NONE

      ! Header files
#     include "define.h"
#     include "CMN_SIZE"       ! Size parameters

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep centain internal variables
      ! and routines from being seen outside "omi_so2_obs_mod.f"
      !=================================================================

      ! Make everything PRIVATE...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CALC_OMI_SO2_FORCE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Parameters
      INTEGER, PARAMETER :: TIME_WINDOW = 60 ! (units: min)

      ! VLIDORT parameters
      INTEGER, PARAMETER :: N_SZA = 12
      INTEGER, PARAMETER :: N_VZA = 8
      INTEGER, PARAMETER :: N_LEV = 48

      ! Variables

      ! VLIDORT varibales
      REAL*4,  ALLOCATABLE :: SZA(:),      VZA(:)
      REAL*4               :: O3_C0,       O3_C1,       O3_C2
      REAL*4               :: Factor
      REAL*4,  ALLOCATABLE :: I0(:,:),     I1(:,:),     I2(:,:)
      REAL*4,  ALLOCATABLE :: Ir(:,:)
      REAL*4               :: Sb
      REAL*4,  ALLOCATABLE :: dI0(:,:,:),  dI1(:,:,:),  dI2(:,:,:)
      REAL*4,  ALLOCATABLE :: dIr(:,:,:)
      REAL*4,  ALLOCATABLE :: AIR_P(:),    ALTITUDE(:), O3_P(:)
      REAL*4,  ALLOCATABLE :: TEMP_P(:)
      REAL*8,  ALLOCATABLE :: VLIDORT_PRESS(:) ! (unit: hPa)

      ! Recored to store OMI SO2 observations
      TYPE OMI_SO2_OBS
         REAL    :: LON           ! Latitude (units: degree)
         REAL    :: LAT           ! Longitude (units: dgeree)
         REAL*8  :: TIME          ! Time at start of scan (TAI93) (units: s)
         REAL*8  :: SZA           ! Solar zenith angle (units: degree)
         REAL*8  :: VZA           ! View zenith angle (units: dgeree)
         REAL*8  :: V_OMI_SO2     ! OMI vertical SO2 (units: D.U.)
         REAL*8  :: S_OMI_SO2     ! OMI slant column SO2 (units: D.U.)
         REAL*8  :: AMF           ! New air mass factor
         REAL*8  :: V_SO2         ! OMI vertical SO2 correted by new AMF (units: D.U.)
         REAL*8  :: V_GC_SO2      ! GEOS-Chem vertical column SO2 (units: D.U.)
         REAL*8  :: S_GC_SO2      ! GEOS-Chem slant column SO2 (units: D.U.)
         REAL*8  :: V_ERR         ! V_OMI_SO2 error (units: D.U.)
         REAL*8  :: S_ERR         ! S_OMI_SO2 error (units: D.U.)
      END TYPE OMI_SO2_OBS

      TYPE(OMI_SO2_OBS), ALLOCATABLE :: OMI_SO2(:)

      INTEGER                        :: N_SO2
      INTEGER                        :: N_CURR

      LOGICAL,           ALLOCATABLE :: FLAGS(:)

      REAL*8                         :: CURR_GC_V_SO2(IIPAR,JJPAR)
      REAL*8                         :: CURR_GC_S_SO2(IIPAR,JJPAR)
      REAL*8                         :: CURR_AMF(IIPAR,JJPAR)
      REAL*8                         :: CURR_OMI_V_SO2(IIPAR,JJPAR)
      REAL*8                         :: CURR_OMI_S_SO2(IIPAR,JJPAR)
      REAL*8                         :: CURR_V_SO2(IIPAR,JJPAR)
      INTEGER                        :: CURR_COUNT(IIPAR,JJPAR)
      REAL*8                         :: CURR_DIFF_SO2(IIPAR,JJPAR)
      REAL*8                         :: CURR_FORCING(IIPAR,JJPAR)
      REAL*8                         :: CURR_COST(IIPAR,JJPAR)

      REAL*8                         :: ALL_GC_V_SO2(IIPAR,JJPAR)  = 0D0
      REAL*8                         :: ALL_GC_S_SO2(IIPAR,JJPAR)  = 0D0
      REAL*8                         :: ALL_AMF(IIPAR,JJPAR)       = 0D0
      REAL*8                         :: ALL_OMI_V_SO2(IIPAR,JJPAR) = 0D0
      REAL*8                         :: ALL_OMI_S_SO2(IIPAR,JJPAR) = 0D0
      REAL*8                         :: ALL_V_SO2(IIPAR,JJPAR)     = 0D0
      INTEGER                        :: ALL_COUNT(IIPAR,JJPAR)     = 0
      REAL*8                         :: ALL_DIFF_SO2(IIPAR,JJPAR)  = 0D0
      REAL*8                         :: ALL_FORCING(IIPAR,JJPAR)   = 0D0
      REAL*8                         :: ALL_COST(IIPAR,JJPAR)      = 0D0

      REAL*8                         :: OMI_SZA, OMI_VZA
      INTEGER                        :: SZA_MIN, SZA_MAX
      INTEGER                        :: VZA_MIN, VZA_MAX

      ! for interpolation
      REAL*8                         :: INTER_I0,  INTER_I1
      REAL*8                         :: INTER_I2,  INTER_Ir
      REAL*8                         :: INTER_dI0(N_LEV)
      REAL*8                         :: INTER_dI1(N_LEV)
      REAL*8                         :: INTER_dI2(N_LEV)
      REAL*8                         :: INTER_dIr(N_LEV)


      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAIN" statement
      !=================================================================
      CONTAINS
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE INIT_OMI_SO2_OBS
!
!*****************************************************************************
!  Subroutine INIT_OMI_SO2_OBS initialize the OMI SO2 observation 
!   ( 1) Check if SO2 observation is specified by the adjoint input
!   files
!  (ywang, 07/16/14)
!*****************************************************************************
!
      ! Reference to f90 modules
      USE TRACER_MOD,         ONLY : TRACER_NAME
      USE TRACERID_MOD,       ONLY : IDTSO2
      USE ADJ_ARRAYS_MOD,     ONLY : OBS_THIS_TRACER
      USE ERROR_MOD,          ONLY : ERROR_STOP

      !=================================================================
      ! INIT_OMI_SO2_OBS begins here!
      !=================================================================

      IF ( OBS_THIS_TRACER(IDTSO2) ) THEN

         WRITE( 6, 100 ) IDTSO2, TRACER_NAME(IDTSO2)

      ELSE

         CALL ERROR_STOP( 'Error: Obs SO2 tracer is not specified in
     &OBSERVATION MENU in input.gcadj', 
     &                    'INIT_OMI_SO2_OBS (omi_so2_obs_mod.f)' ) 

      END IF

 100  FORMAT( 3X, 'Tracer ID: ', I4, ' Use OMI ', A6 ) 

      ! Return to the calling routine
      END SUBROUTINE INIT_OMI_SO2_OBS
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE READ_OMI_SO2_OBS( YYYYMMDD )
!
!*****************************************************************************
!  Subroutine READ_OMI_SO2_OBS read OMI SO2 data 
!  (ywang, 07/09/15)
!
!  Arguements as Input:
!  ===========================================================================
!  ( 1) YYYYMMDD    (INTEGER) : Current year-month-day
!
!  Module variable as Output:
!  ===========================================================================
!  ( 1) N_SO2       (INTEGER) : Number of OMI SO2 observations for
!  current day
!  ( 2) OMI_SO2 (OMI_SO2_OBS) : OMI SO2 observations
!
!*****************************************************************************
!
      ! Reference to f90 modules
      USE NETCDF
      USE TIME_MOD,         ONLY : EXPAND_DATE, GET_NYMD

      ! Arguements
      INTEGER, INTENT( IN)   :: YYYYMMDD

      ! Local variables
      INTEGER                :: FID,          N_ID
      INTEGER                :: LON_ID,       LAT_ID,     TIME_ID
      INTEGER                :: SZA_ID,       VZA_ID
      INTEGER                :: V_OMI_SO2_ID, ERR_ID
      CHARACTER(LEN=255)     :: READ_FILENAME
      CHARACTER(LEN=4)       :: TMP
      INTEGER, ALLOCATABLE   :: TMPINT(:)
      REAL*4,  ALLOCATABLE   :: TMP4(:)
      REAL*8,  ALLOCATABLE   :: TMP8(:)
      LOGICAL                :: LF

      !=================================================================
      ! READ_OMI_SO2_OBS begins here!
      !================================================================= 

      READ_FILENAME = 'OMI_SO2_YYYYMMDD.nc'

      ! Expand date tokens in filename
      CALL EXPAND_DATE( READ_FILENAME, YYYYMMDD, 9999 )

      ! Construct complete filename 
      READ_FILENAME = "./data/OMI_SO2/" // TRIM( READ_FILENAME )

      ! Does data file exist? If not, it means no data in the day.
      INQUIRE( FILE = TRIM( READ_FILENAME ), EXIST = LF )
      IF ( .NOT. LF ) THEN

         ! No data
         N_SO2 = 0

         WRITE(6, 120) GET_NYMD()

         RETURN

      END IF

 120  FORMAT(' - READ_OMI_SO2_OBS: No data file (warning) in', I10)

      ! Print to screen
      WRITE(6, 100) TRIM( READ_FILENAME )
 100  FORMAT(' - READ_OMI_SO2_OBS: reading file: ', A)

      ! Open file and assign file id (FID)
      CALL CHECK( NF90_OPEN( READ_FILENAME, NF90_NOWRITE, FID ), 0 )

      !-----------------------------------
      ! Get data record IDs
      !-----------------------------------
      CALL CHECK( NF90_INQ_DIMID( FID, "time",    N_ID         ), 100 )

      CALL CHECK( NF90_INQ_VARID( FID, "lon",     LON_ID       ), 101 )
      CALL CHECK( NF90_INQ_VARID( FID, "lat",     LAT_ID       ), 102 )
      CALL CHECK( NF90_INQ_VARID( FID, "time",    TIME_ID      ), 103 )
      CALL CHECK( NF90_INQ_VARID( FID, "SZA",     SZA_ID       ), 104 )
      CALL CHECK( NF90_INQ_VARID( FID, "VZA",     VZA_ID       ), 105 )
      CALL CHECK( NF90_INQ_VARID( FID, "SO2",     V_OMI_SO2_ID ), 106 )
      CALL CHECK( NF90_INQ_VARID( FID, "ERR",     ERR_ID       ), 107 )

      !------------------------------------
      ! Read dimensions
      !------------------------------------

      ! Read number of observations, N_SO2
      CALL CHECK( NF90_INQUIRE_DIMENSION( FID, N_ID, TMP, N_SO2 ), 200 )

      ! Print to screen
      WRITE(6, 110) N_SO2, GET_NYMD()
 110  FORMAT('      Number of OMI SO2 observations: ' I10, ' in ' I10)

      !-------------------------------------
      ! Read 1D data
      !-------------------------------------

      ! Allocate temporal arrays for 1D data
      ALLOCATE( TMPINT(N_SO2) )
      ALLOCATE( TMP4(N_SO2)   )
      ALLOCATE( TMP8(N_SO2)   )
      TMPINT = 0
      TMP4   = 0E0
      TMP8   = 0D0

      ! Allocate OMI SO2 observations array
      IF ( ALLOCATED( OMI_SO2 ) ) DEALLOCATE( OMI_SO2 )
      ALLOCATE( OMI_SO2(N_SO2) )

      IF ( ALLOCATED( FLAGS) ) DEALLOCATE( FLAGS )
      ALLOCATE( FLAGS(N_SO2) )

      ! Read longitude
      CALL CHECK( NF90_GET_VAR( FID, LON_ID, TMP4 ), 301 )
      OMI_SO2(1:N_SO2)%LON = TMP4(1:N_SO2)

      ! Read latitude
      CALL CHECK( NF90_GET_VAR( FID, LAT_ID, TMP4 ), 302 )
      OMI_SO2(1:N_SO2)%LAT = TMP4(1:N_SO2)

      ! Read time
      CALL CHECK( NF90_GET_VAR( FID, TIME_ID, TMP8 ), 303 )
      OMI_SO2(1:N_SO2)%TIME = TMP8(1:N_SO2)

      ! Read SZA
      CALL CHECK( NF90_GET_VAR( FID, SZA_ID, TMP4 ), 304 )
      OMI_SO2(1:N_SO2)%SZA = TMP4(1:N_SO2)

      ! Read VZA
      CALL CHECK( NF90_GET_VAR( FID, VZA_ID, TMP4 ), 305 )
      OMI_SO2(1:N_SO2)%VZA = TMP4(1:N_SO2)

      ! Read ColumnAmountSO2_PBL
      CALL CHECK( NF90_GET_VAR( FID, V_OMI_SO2_ID, TMP4 ), 306 )
      OMI_SO2(1:N_SO2)%V_OMI_SO2 = TMP4(1:N_SO2)

      ! Read ERR
      CALL CHECK( NF90_GET_VAR( FID, ERR_ID, TMP4       ), 306 )
      OMI_SO2(1:N_SO2)%V_ERR = TMP4(1:N_SO2)

      ! Close the file
      CALL CHECK( NF90_CLOSE( FID ), 9999 )

      DEALLOCATE( TMPINT )
      DEALLOCATE( TMP4   )
      DEALLOCATE( TMP8   )

      ! Return to the calling routines
      END SUBROUTINE READ_OMI_SO2_OBS
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE CALC_OMI_SO2_FORCE( COST_FUNC )
!
!*****************************************************************************
!  Subroutine CALC_OMI_SO2_FORCE calculate the adjoint forcing from OMI
!  L2 swath SO2 observation and updates the cost function.
!  (ywang, 07/16/14)
!
!  Arguments as Input/Output:
!  ===========================================================================
!  ( 1) COST_FUNC (REAL*8) : Cost funtion                     [unitless]
!
!  NOTES:
!*****************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,     ONLY : OBS_FREQ
      USE CHECKPT_MOD,        ONLY : CHK_STT
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE DAO_MOD,            ONLY : AIRVOL
      USE DAO_MOD,            ONLY : BXHEIGHT
      USE DAO_MOD,            ONLY : AD, T
      USE GRID_MOD,           ONLY : GET_IJ
      USE PRESSURE_MOD,       ONLY : GET_PCENTER
      USE TIME_MOD,           ONLY : GET_NYMD, GET_NHMS
      USE TIME_MOD,           ONLY : GET_TAUb, GET_TAU, GET_TAUe
      USE TRACER_MOD,         ONLY : XNUMOL, TCVV
      USE TRACERID_MOD,       ONLY : IDTSO2
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP

      ! Arguements
      REAL*8, INTENT(INOUT)       :: COST_FUNC

      ! Parameters
      REAL*8, PARAMETER           :: ALBEDO  = 0.05
      REAL*8, PARAMETER           :: OMI_AMF = 0.36
      REAL*8, PARAMETER           :: M2DU    = 2.69D20 ! 1 DU = 2.69D20 molec/m2

      ! Local variables
      INTEGER                     :: NT
      INTEGER                     :: NC
      INTEGER                     :: IIJJ(2), I, J, L
      INTEGER                     :: K
      INTEGER                     :: IL, L1, L2
      INTEGER, ALLOCATABLE        :: STORE_L1(:,:), STORE_L2(:,:)

      REAL*8                      :: OLD_COST
      REAL*8                      :: TMP_COST
      REAL*8, ALLOCATABLE         :: NEW_COST(:)

      REAL*8                      :: GC_SO2_CONC(LLPAR)  ! GEOS-Chem SO2 at each layer (units: kg/m3)
      REAL*8                      :: LAYER_GC_SO2(LLPAR) ! GEOS-Chem SO2 at each layer (units: DU)
      REAL*8, ALLOCATABLE         :: COLUMN_GC_SO2(:,:)  ! GEOS-Chem vertical column SO2 (units: DU)
      REAL*8                      :: DIFF
      REAL*8                      :: FORCING

      ! These variables are reversed (from top to bottom)
      REAL*8                      :: REV_GC_SO2(LLPAR)         ! v/v
      REAL*8                      :: REV_GC_TMPU(LLPAR)        ! K
      REAL*8                      :: REV_GC_PCEN(LLPAR)        ! hPa
      REAL*8, ALLOCATABLE         :: STORE_REV_GC_SO2(:,:,:)   ! v/v
      REAL*8, ALLOCATABLE         :: STORE_REV_GC_TMPU(:,:,:)  ! K
      REAL*8, ALLOCATABLE         :: STORE_REV_GC_PCEN(:,:,:)  ! hPa
      REAL*8                      :: REV_INTER_GC_SO2(N_LEV)   ! v/v
      REAL*8                      :: REV_INTER_GC_TMPU(N_LEV)  ! K
      REAL*8                      :: REV_INTER_GC_T(N_LEV)     ! Celsius degree
      REAL*8                      :: REV_INTER_GC_TT(N_LEV)    ! (Celsius degree) ^ 2
      REAL*8, ALLOCATABLE         :: STORE_REV_INTER_GC_SO2(:,:,:)   ! v/v
      REAL*8, ALLOCATABLE         :: STORE_REV_INTER_GC_TMPU(:,:,:)  ! K

      REAL*8                      :: RADIANCE          ! radiance
      REAL*8                      :: PROFILE_WF(N_LEV) ! profile weighting function
      REAL*8                      :: O3_XS(N_LEV)
      REAL*8                      :: ABS_XS_O3(N_LEV)  ! ozone absorption cross section
      REAL*8                      :: SCATW(N_LEV)      ! scattering weight

      REAL*8                      :: MIN_GC_PRESS, MAX_GC_PRESS


      LOGICAL                     :: FIRST_DEBUG ! for debug
      LOGICAL                     :: FIRST = .TRUE.
      LOGICAL, ALLOCATABLE        :: PROFILE_FLAG(:,:)
   
      !=================================================================
      ! CALC_OMI_SO2_FORCE begins here!
      !=================================================================

      PRINT*, ' - CALC_OMI_SO2_FORCE: OMI SO2 forcing '

      ! Initialize
      IF ( FIRST ) THEN

         CALL INIT_OMI_SO2_OBS

         CALL INIT_LOOKUP

         CALL GET_LOOKUP

         FIRST = .FALSE.

      END IF

      ! Save a value of the cost function first
      OLD_COST = COST_FUNC

      ! Check if it is the last time of a day
      IF ( OBS_FREQ > 60 ) THEN
         PRINT*, '236000 - OBS_FREQ * 100 is not valid'
         STOP
      END IF
      IF ( GET_NHMS() == 236000 - OBS_FREQ * 100 ) THEN
!      IF ( (GET_NHMS() == 236000 - OBS_FREQ * 100) .OR.
!     &     ( ABS(GET_TAU() -GET_TAUe()) < 1e-6  ) ) THEN

         CALL READ_OMI_SO2_OBS( GET_NYMD() )

      END IF

      ! No observations for current day
      IF ( N_SO2 == 0 ) THEN

         PRINT*, '    - CALC_OMI_SO2_FORCE: No OMI SO2 obsevations for
     &current day'

         IF ( ABS( GET_TAUb() - GET_TAU() ) < 1E-6) THEN

            CALL MAKE_AVERAGE_OMI_SO2

         END IF

         RETURN

      END IF

      ! GET observations in time window
      CALL GET_OBS

      ! No observations for time window
      IF ( N_CURR == 0 ) THEN

         PRINT*, '    - CALC_OMI_SO2_FORCE: No OMI SO2 obsevations for
     &current time window'

         IF ( ABS( GET_TAUb() - GET_TAU() ) < 1E-6) THEN

            CALL MAKE_AVERAGE_OMI_SO2

         END IF

         RETURN

      END IF

      ALLOCATE( STORE_L1(IIPAR,JJPAR) )
      ALLOCATE( STORE_L2(IIPAR,JJPAR) )

      ALLOCATE( STORE_REV_GC_SO2(IIPAR,JJPAR,LLPAR)  )
      ALLOCATE( STORE_REV_GC_TMPU(IIPAR,JJPAR,LLPAR) )
      ALLOCATE( STORE_REV_GC_PCEN(IIPAR,JJPAR,LLPAR) )

      ALLOCATE( STORE_REV_INTER_GC_SO2(IIPAR,JJPAR,N_LEV ) )
      ALLOCATE( STORE_REV_INTER_GC_TMPU(IIPAR,JJPAR,N_LEV) )

      ! Reset
      CURR_GC_V_SO2  = 0D0
      CURR_GC_S_SO2  = 0D0
      CURR_AMF       = 0D0
      CURR_OMI_V_SO2 = 0D0
      CURR_OMI_S_SO2 = 0D0
      CURR_V_SO2     = 0D0
      CURR_COUNT     = 0
      CURR_DIFF_SO2  = 0D0
      CURR_FORCING   = 0D0
      CURR_COST      = 0D0

      ! Reset
      IF ( ALLOCATED( NEW_COST ) ) DEALLOCATE( NEW_COST )
      ALLOCATE ( NEW_COST(N_CURR) )
      NEW_COST = 0D0

      ! Reset
      IF ( ALLOCATED( PROFILE_FLAG ) ) DEALLOCATE( PROFILE_FLAG )
      ALLOCATE( PROFILE_FLAG(IIPAR,JJPAR) )
      PROFILE_FLAG = .FALSE.

      IF ( ALLOCATED( COLUMN_GC_SO2 ) ) DEALLOCATE( COLUMN_GC_SO2 )
      ALLOCATE( COLUMN_GC_SO2(IIPAR,JJPAR) )
      COLUMN_GC_SO2 = 1D20

      FIRST_DEBUG = .TRUE.

      NC = 0
      ! Loop for all observations
      DO NT = 1, N_SO2, 1

         ! Observations in time window and simulation area
         IF ( FLAGS(NT) ) THEN

            ! Get grid box of current record
            IIJJ = GET_IJ( REAL(OMI_SO2(NT)%LON ,4),
     &                     REAL(OMI_SO2(NT)%LAT ,4) )

            I    = IIJJ(1)
            J    = IIJJ(2)

            ! There are some OMI pixels located in the same GC grid box,
            ! thus we only calculte GC variables once.
            IF ( .NOT. PROFILE_FLAG(I,J) ) THEN

               !----------------------------------------------------
               ! Calculate GC vertial column SO2 density (unit: DU)
               !----------------------------------------------------
               ! SO2 outside troposphere is set to 0
               GC_SO2_CONC  = 0D0
               LAYER_GC_SO2 = 0D0
               DO L = 1, LLPAR, 1

!                  IF ( ITS_IN_THE_TROP(I,J,L) ) THEN

                     ! Units conversion [ kg/gridbox => DU ]
                     ! Units of LAYER_GC_SO2: DU
                     !           GC_SO2_CONC: kg/m3
                     !              CHK_STT : kg/gridbox
                     !               AIRVOL : m3/gridbox
                     !             BXHEIGHT : m 
                     !               XNUMOL : molec/kg
                     ! M2DU = 2.69D20 , 1 DU = 2.69D20 molec/m2
                
                     GC_SO2_CONC(L)  = CHK_STT(I,J,L,IDTSO2) /
     &                                 AIRVOL(I,J,L)

                     LAYER_GC_SO2(L) = GC_SO2_CONC(L)        *
     &                                 BXHEIGHT(I,J,L)       *
     &                                 XNUMOL(IDTSO2)        /
     &                                 M2DU
           
!                  END IF

               END DO
               COLUMN_GC_SO2(I,J)    = SUM( LAYER_GC_SO2 )

!         PRINT*,'*******************************DEBUG'
!         PRINT*, 'M2DU = ', M2DU
!         PRINT*, 'XNUMOL = ', XNUMOL(IDTSO2)
!         DO L = 1, LLPAR, 1
!         WRITE(6, '(5E13.5)') CHK_STT(I,J,L,IDTSO2), AIRVOL(I,J,L), 
!     &     GC_SO2_CONC(L), BXHEIGHT(I,J,L), LAYER_GC_SO2(L)
!         END DO
!         PRINT*, 'COLUMN_GC_SO2(I,J) = ', COLUMN_GC_SO2(I,J)
!         PRINT*,'*******************************'

               !---------------------------------------------------
               ! reversed (from top to bottom)
               !---------------------------------------------------
               DO L = 1, LLPAR, 1

                  K = LLPAR - L + 1

               !  The conversion is as follows:
               !
               !   kg tracer(N)       1        Air mol wt     
               !   -----------  * -------- *  -------------   
               !        1          kg air     tracer mol wt   
               !
               !       moles tracer     volume tracer
               !   =   ------------  =  -------------
               !        moles air        volume air
               !
               ! Since the volume of a gas depends on the number of moles.
               ! Therefore, with:
               !
               !  TCMASS(N) = mol. wt. of tracer (AMU)
               !  TCVV(N)   = 28.97 / TCMASS(N)
               !            = mol. wt. of air (AMU) / mol. wt. of tracer (AMU)
               !  AD(I,J,L) = mass of air (kg) in grid box (I,J,L)
               !     
               ! the conversion is:
               ! 
               !  CHK_STT(I,J,L,N) [kg] * TCVV(N) / AD(I,J,L) = [v/v]
                  REV_GC_SO2(L)  = CHK_STT(I,J,K,IDTSO2) * 
     &                             TCVV(IDTSO2)          /
     &                             AD(I,J,K)


                  REV_GC_TMPU(L) = T(I,J,K)
                  REV_GC_PCEN(L) = GET_PCENTER(I,J,K)

                  IF ( REV_GC_SO2(L) < 0D0 ) THEN

                     CALL ERROR_STOP( 'REV_GC_SO2(L) < 0.0', 
     &                                'omi_so2_obs_mod.f' )

                  END IF

               END DO

               !--------------------------------------------
               ! Interpolate REV_GC_SO2 and REV_GC_TMPU at
               ! REV_GC_TMPU into VLIDORT_PRESS
               !--------------------------------------------
               CALL SPLINE( N_LEV, LLPAR, REV_GC_PCEN,   REV_GC_SO2,
     &                                VLIDORT_PRESS, REV_INTER_GC_SO2  )
               CALL SPLINE( N_LEV, LLPAR, REV_GC_PCEN,   REV_GC_TMPU,
     &                              VLIDORT_PRESS, REV_INTER_GC_TMPU )

               MIN_GC_PRESS = REV_GC_PCEN(1)
               MAX_GC_PRESS = REV_GC_PCEN(LLPAR)

               L1 = MINLOC( ABS( VLIDORT_PRESS - MIN_GC_PRESS ), 1 )
               L2 = MINLOC( ABS( VLIDORT_PRESS - MAX_GC_PRESS ), 1 )

               DO IL = 1, L1-1, 1
                  REV_INTER_GC_SO2(IL)  = REV_GC_SO2(1)
                  REV_INTER_GC_TMPU(IL) = REV_GC_TMPU(1)
               END DO

               IF ( VLIDORT_PRESS(L1) < MIN_GC_PRESS ) THEN
                  REV_INTER_GC_SO2(L1)  = REV_GC_SO2(1)
                  REV_INTER_GC_TMPU(L1) = REV_GC_TMPU(1)
               END IF

               DO IL = L2+1, N_LEV, 1
                  REV_INTER_GC_SO2(IL)  = REV_GC_SO2(LLPAR)
                  REV_INTER_GC_TMPU(IL) = REV_GC_TMPU(LLPAR)
               END DO

               IF ( VLIDORT_PRESS(L2) > MAX_GC_PRESS ) THEN
                  REV_INTER_GC_SO2(L2)  = REV_GC_SO2(LLPAR)
                  REV_INTER_GC_TMPU(L2) = REV_GC_TMPU(LLPAR)
               END IF

               DO IL = 1, N_LEV, 1

                  IF ( REV_INTER_GC_SO2(IL) < 0D0 ) THEN
                     REV_INTER_GC_SO2(IL) = MINVAL( REV_GC_SO2 )
                  END IF

                  IF ( REV_INTER_GC_TMPU(IL) < 0D0 ) THEN
                     CALL ERROR_STOP( 'REV_INTER_GC_TMPU(IL) < 0D0',
     &                                'omi_so2_obs_mod.f' )
                  END IF

               END DO

               !----------------------------------------------
               ! Store L1 and L2
               !----------------------------------------------
               STORE_L1(I,J) = L1
               STORE_L2(I,J) = L2

               !----------------------------------------------
               ! Store REV_GC_SO2, REV_GC_TMPU, and
               !       REV_GC_PCEN
               !----------------------------------------------
               DO L = 1, LLPAR, 1
                  STORE_REV_GC_SO2(I,J,L)  = REV_GC_SO2(L)
                  STORE_REV_GC_TMPU(I,J,L) = REV_GC_TMPU(L)
                  STORE_REV_GC_PCEN(I,J,L) = REV_GC_PCEN(L)
               END DO

               !----------------------------------------------
               ! Store REV_INTER_GC_SO2 and REV_INTER_GC_TMPU
               !----------------------------------------------
               DO IL = 1, N_LEV, 1
                  STORE_REV_INTER_GC_SO2(I,J,IL) = REV_INTER_GC_SO2(IL)
                  STORE_REV_INTER_GC_TMPU(I,J,IL) =
     &                                      REV_INTER_GC_TMPU(IL)
               END DO

               PROFILE_FLAG(I,J) = .TRUE.

            END IF ! .NOT. PROFILE_FLAG(I,J)

            !---------------------------------
            ! Get SZA an VZA for interpolation
            !---------------------------------
            OMI_SZA = REAL( OMI_SO2(NT)%SZA, 4 )
            OMI_VZA = REAL( OMI_SO2(NT)%VZA, 4 )
            CALL INTERPOLATION

            !---------------------------------
            ! Calculate scattering weight
            !---------------------------------

            ! radiance
            RADIANCE = INTER_I0 + INTER_I1 + INTER_I2 +
     &                 ( ALBEDO * INTER_Ir ) / ( 1.0 - ( ALBEDO * Sb ) )

            ! profile weighting function
            PROFILE_WF(:) = (
     &                       INTER_dI0(:) + INTER_dI1(:) + 
     &                       INTER_dI2(:) +
     &                       ( ALBEDO * INTER_dIr(:) ) /
     &                       (1.0 - ( ALBEDO * Sb )  )
     &                      ) / Factor

            ! ozone absorption cross section
            DO IL = 1, N_LEV, 1
               REV_INTER_GC_TMPU(IL) = STORE_REV_INTER_GC_TMPU(I,J,IL)
            END DO
            REV_INTER_GC_T(:)  = REV_INTER_GC_TMPU(:) - 273.15
            REV_INTER_GC_TT(:) = REV_INTER_GC_T(:) ** 2
            O3_XS(:) =   O3_C0                        +
     &                 ( O3_C1 * REV_INTER_GC_T(:)  ) +
     &                 ( O3_C2 * REV_INTER_GC_TT(:) )
            ABS_XS_O3(:) =  O3_XS(:) * O3_P(:)

            ! scattering weigth
            SCATW(:) = (-PROFILE_WF(:)) / ( RADIANCE * ABS_XS_O3(:) )

            !----------------------------------------------
            ! Retrieve L1 and L2
            !----------------------------------------------
            L1 = STORE_L1(I,J)
            L2 = STORE_L2(I,J)

            !----------------------------------------------
            ! Retrieve REV_GC_SO2, REV_GC_TMPU, and
            !          REV_GC_PCEN
            !----------------------------------------------
            DO L = 1, LLPAR, 1
               REV_GC_SO2(L)  = STORE_REV_GC_SO2(I,J,L)
               REV_GC_TMPU(L) = STORE_REV_GC_TMPU(I,J,L)
               REV_GC_PCEN(L) = STORE_REV_GC_PCEN(I,J,L)
            END DO

            !--------------------
            ! air mass factor
            !--------------------
            DO IL = 1, N_LEV, 1
               REV_INTER_GC_SO2(IL) = STORE_REV_INTER_GC_SO2(I,J,IL)
            END DO
            OMI_SO2(NT)%AMF =
     &                   SUM( REV_INTER_GC_SO2(L1:L2) * SCATW(L1:l2) ) /
     &                   SUM( REV_INTER_GC_SO2(L1:L2) )

            !---------------------
            ! OMI_SO2(NT)
            !---------------------
            OMI_SO2(NT)%S_OMI_SO2 = OMI_SO2(NT)%V_OMI_SO2 * OMI_AMF

            OMI_SO2(NT)%V_SO2     = OMI_SO2(NT)%S_OMI_SO2 / 
     &                              OMI_SO2(NT)%AMF

            OMI_SO2(NT)%V_GC_SO2  = COLUMN_GC_SO2(I,J)

            OMI_SO2(NT)%S_GC_SO2  = OMI_SO2(NT)%V_GC_SO2 * 
     &                              OMI_SO2(NT)%AMF
!            OMI_SO2(NT)%V_ERR     = 0.7
            OMI_SO2(NT)%S_ERR     = OMI_SO2(NT)%V_ERR * OMI_AMF

            !----------------------
            ! current statistics
            !----------------------
            CURR_GC_V_SO2(I,J)  = CURR_GC_V_SO2(I,J)  +
     &                            OMI_SO2(NT)%V_GC_SO2
            CURR_GC_S_SO2(I,J)  = CURR_GC_S_SO2(I,J)  +
     &                            OMI_SO2(NT)%S_GC_SO2
            CURR_AMF(I,J)       = CURR_AMF(I,J)       +
     &                            OMI_SO2(NT)%AMF
            CURR_OMI_V_SO2(I,J) = CURR_OMI_V_SO2(I,J) +
     &                            OMI_SO2(NT)%V_OMI_SO2
            CURR_OMI_S_SO2(I,J) = CURR_OMI_S_SO2(I,J) +
     &                            OMI_SO2(NT)%S_OMI_SO2
            CURR_V_SO2(I,J)     = CURR_V_SO2(I,J)     +
     &                            OMI_SO2(NT)%V_SO2
            CURR_COUNT(I,J)     = CURR_COUNT(I,J)     + 1

!            !----------------------
!            ! all statistics
!            !----------------------
!            ALL_GC_V_SO2(I,J)  = ALL_GC_V_SO2(I,J)  +
!     &                            OMI_SO2(NT)%V_GC_SO2
!            ALL_GC_S_SO2(I,J)  = ALL_GC_S_SO2(I,J)  +
!     &                            OMI_SO2(NT)%S_GC_SO2
!            ALL_AMF(I,J)       = ALL_AMF(I,J)       +
!     &                            OMI_SO2(NT)%AMF
!            ALL_OMI_V_SO2(I,J) = ALL_OMI_V_SO2(I,J) +
!     &                            OMI_SO2(NT)%V_OMI_SO2
!            ALL_OMI_S_SO2(I,J) = ALL_OMI_S_SO2(I,J) +
!     &                            OMI_SO2(NT)%S_OMI_SO2
!            ALL_V_SO2(I,J)     = ALL_V_SO2(I,J)     +
!     &                            OMI_SO2(NT)%V_SO2
!            ALL_COUNT(I,J)     = ALL_COUNT(I,J)     + 1

            !-------------------
            ! debug
            !-------------------
            IF ( FIRST_DEBUG ) THEN

               WRITE(6, 100) NT, OMI_SO2(NT)%TIME
               WRITE(6, 101) OMI_SO2(NT)%LON, OMI_SO2(NT)%LAT
               WRITE(6, 102) I, J
               WRITE(6, '(A16,F10.5,A9)')
     &                  ' --- Radiance = ', RADIANCE, ' W/(cm^2)'
               WRITE(6, 110) OMI_SZA, OMI_VZA
               WRITE(6, 120) SZA_MIN, SZA_MAX
               WRITE(6, 130) VZA_MIN, VZA_MAX
               WRITE(6, 140) INTER_I0, INTER_I1, INTER_I2, INTER_Ir
               WRITE(6, 150)
               DO IL = 1, N_LEV, 1
               WRITE(6, 160) IL, INTER_dI0(IL), INTER_dI1(IL),
     &                           INTER_dI2(IL), INTER_dIr(IL)
               END DO
               WRITE(6, 103)
               DO L = 1, LLPAR, 1
               WRITE(6, 104) L, REV_GC_PCEN(L),   REV_GC_TMPU(L),
     &                          REV_GC_SO2(L) * 1D9
               END DO
               WRITE(6, 161)
               DO IL = 1, N_LEV, 1
               IF ( IL == L1 ) THEN
                  WRITE(6, '(A)') REPEAT(' ', 5) //  REPEAT('-', 199)
               END IF
               WRITE(6, 162) IL, ALTITUDE(IL), VLIDORT_PRESS(IL),
     &                       REV_INTER_GC_TMPU(IL),
     &                       REV_INTER_GC_SO2(IL) * 1D9,
     &                       O3_XS(IL), O3_P(IL), ABS_XS_O3(IL),
     &                       PROFILE_WF(IL), SCATW(IL)
               IF ( IL == L2 ) THEN
                  WRITE(6, '(A)') REPEAT(' ', 5) //  REPEAT('-', 199)
               END IF
               END DO
               WRITE(6, '(A11,F10.5)') ' --- AMF = ', OMI_SO2(NT)%AMF
               WRITE(6, '(A38,F10.5)') ' --- OMI_SO2(NT)%V_OMI_SO2 = ',
     &                                       OMI_SO2(NT)%V_OMI_SO2
               WRITE(6, '(A38,F10.5)') ' --- OMI_SO2(NT)%S_OMI_SO2 = ',
     &                                       OMI_SO2(NT)%S_OMI_SO2
               WRITE(6, '(A38,F10.5)') ' ---     OMI_SO2(NT)%V_SO2 = ',
     &                                           OMI_SO2(NT)%V_SO2
               WRITE(6, '(A38,F10.5)') ' ---  OMI_SO2(NT)%V_GC_SO2 = ',
     &                                        OMI_SO2(NT)%V_GC_SO2
               WRITE(6, '(A38,F10.5)') ' ---  OMI_SO2(NT)%S_GC_SO2 = ',
     &                                        OMI_SO2(NT)%S_GC_SO2
               WRITE(6, '(A38,F10.5)') ' ---     OMI_SO2(NT)%V_ERR = ',
     &                                           OMI_SO2(NT)%V_ERR
               WRITE(6, '(A38,F10.5)') ' ---     OMI_SO2(NT)%S_ERR = ',
     &                                           OMI_SO2(NT)%S_ERR


            END IF
            FIRST_DEBUG = .FALSE.

            !-----------------------------
            ! Calculate adjoint forcing
            !-----------------------------

            ! The difference between GEOS-Chem SO2 and OMI SO2
            DIFF = OMI_SO2(NT)%S_GC_SO2 - OMI_SO2(NT)%S_OMI_SO2
            CURR_DIFF_SO2(I,J) = CURR_DIFF_SO2(I,J) + DIFF
!            ALL_DIFF_SO2(I,J)  = ALL_DIFF_SO2(I,J)  + DIFF

            ! S_{obs}^{-1} * DIFF
            FORCING           = DIFF / (OMI_SO2(NT)%S_ERR ** 2)
            CURR_FORCING(I,J) = CURR_FORCING(I,J) + FORCING
            ALL_FORCING(I,J)  = ALL_FORCING(I,J)  + FORCING

            ! Contribution to the cost function
            TMP_COST       = 0.5D0 * DIFF * FORCING
            NC             = NC + 1
            NEW_COST(NC)   = NEW_COST(NC)   + TMP_COST
            CURR_COST(I,J) = CURR_COST(I,J) + TMP_COST
            ALL_COST(I,J)  = ALL_COST(I,J)  + TMP_COST

            ! Now pass the adjoint back to the adjoint tracer array
            DO L = 1, LLPAR, 1

               IF ( ITS_IN_THE_TROP(I,J,L) ) THEN

                  STT_ADJ(I,J,L,IDTSO2) = STT_ADJ(I,J,L,IDTSO2)   +
     &                                    FORCING         * 
     &                                    OMI_SO2(NT)%AMF /
     &                                    AIRVOL(I,J,L)   *
     &                                    BXHEIGHT(I,J,L) *
     &                                    XNUMOL(IDTSO2)  / 
     &                                    M2DU

               END IF

            END DO

         END IF ! FLAGS(NT)

      END DO ! NT

 100  FORMAT( ' --- NT = ' I10 ' TAU93 = ' F18.2 )
 101  FORMAT( ' --- LON = ' F6.2 ' LAT = ' F6.2 )
 102  FORMAT( ' --- I = ' I6 ' J = ' I6 )
 110  FORMAT( ' --- OMI_SZA = ' F6.2 ', OMI_VZA = ' F6.2 )
 120  FORMAT( ' --- SZA_MIN = ' I6,  ', SZA_MAX = ' I6   )
 130  FORMAT( ' --- VZA_MIN = ' I6,  ', VZA_MAX = ' I6   )
 140  FORMAT( ' --- Interpolation I0, Ir, I2, and Ir: ', 4E13.5,
     &'W/(cm^2)' )
 150  FORMAT( ' --- Interpolation',2X,'dI0',7X,'dI1',9X,'dI2',6X,
     &'dIr (1E-21 W mol)')
 160  FORMAT( ' ---    ' I4, 4E12.3 )
 103  FORMAT( ' --- Reversed GC:  Pressure (hPa)' 3x 'Temperature (K)'
     &3x 'SO2 (ppbv)' )
 104  FORMAT( ' --- '  I6, 6X, F10.3, 6X, F10.2, 5x, F15.5 )
 161  FORMAT( ' --- Interpolation Altitude', 3X, 'Pressure (hPa)', 3X,
     &   'Temperature (K)', 3X, 'SO2 (ppbv)', 3X, 'O3_XS (cm^2/mol)',
     &   3X, 'Ozone profile (mol/cm^2)', 3X, 'ABS_XS_O3 (unitless)',
     &   3x, 'Weighting function (W mol)',
     &   3X, 'Scattering weight (unitless)')
 162  FORMAT( 7X, I5, 4X, F10.3, 5X, F10.3, 6X, F10.2, 4x, F12.5, 3x,
     &        E16.5, 3X, E20.5, 4x, E20.5, 6X, E20.5, 12X, F10.5 )

      DEALLOCATE( STORE_L1 )
      DEALLOCATE( STORE_L2 )

      DEALLOCATE( STORE_REV_GC_SO2  )
      DEALLOCATE( STORE_REV_GC_TMPU )
      DEALLOCATE( STORE_REV_GC_PCEN )

      DEALLOCATE( STORE_REV_INTER_GC_SO2  )
      DEALLOCATE( STORE_REV_INTER_GC_TMPU )

      ! Update cost function
      COST_FUNC = COST_FUNC + SUM( NEW_COST )
      PRINT*, ' Update value of COST_FUNC = ', COST_FUNC
      PRINT*, ' OMI SO2 contribution      = ', COST_FUNC - OLD_COST

      PRINT*, ' MIN/MAX STT_ADJ  = ', MINVAL(STT_ADJ), MAXVAL(STT_ADJ)
      PRINT*, ' MIN/MAX in       = ', MINLOC(STT_ADJ), MAXLOC(STT_ADJ)
      PRINT*, ' MIN/MAX NEW_COST = ', MINVAL(NEW_COST), MAXVAL(NEW_COST)
      PRINT*, ' MIN/MAX cost in  = ', MINLOC(NEW_COST),MAXLOC(NEW_COST)

      CALL MAKE_CURRENT_OMI_SO2

      IF ( ABS( GET_TAUb() - GET_TAU() ) < 1E-6) THEN

         CALL MAKE_AVERAGE_OMI_SO2

      END IF

      ! Return to the calling routines
      END SUBROUTINE CALC_OMI_SO2_FORCE
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE MAKE_CURRENT_OMI_SO2
!
!*****************************************************************************
!  Subroutine MAKE_CURRENT_OMI_SO2 output some dignostic data for
!  current assimilation time window
!  (ywang, 07/19/14)
!
!*****************************************************************************
!
      ! Reference to f90 module
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,    ONLY : EXPAND_NAME
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
      USE FILE_MOD,          ONLY : IU_OMSO2,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TIME_MOD,          ONLY : GET_NYMD, GET_NHMS
      USE TIME_MOD,          ONLY : GET_TAU

      REAL*4, PARAMETER    :: UNDEF = -999.0

      ! Local Variables
      INTEGER              :: I,    I0, J,  J0, L
      REAL*4               :: SO2(IIPAR,JJPAR,10)
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      REAL*4               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      CHARACTER(LEN=255)   :: OUTPUT_FILE
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE

      !=================================================================
      ! MAKE_CURRENT_OMI_SO2 begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_FILE = 'gctm.omi.so2.YYYYMMDD.hhmm.NN'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM diag File: ' //
     &           'OMI SO2'
      UNIT     = 'DU'
      CATEGORY = 'IJ-AVG-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_FILE )

      ! Replace YYMMDD and hhmmss token w/ actucl value
      CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )

      ! Replace NN token w/ actual value
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add DIAGADJ_DIR prefix to FILENAME
      FILENAME = TRIM( DIAGADJ_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_CURRENT_OMI_SO2:  Writing ', a )

      ! Open file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_OMSO2, FILENAME, TITLE )

      SO2 = 0.0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         IF ( CURR_COUNT(I,J) > 0 ) THEN
            SO2(I,J, 1) = REAL( CURR_GC_V_SO2(I,J),  4 ) /
     &                    REAL( CURR_COUNT(I,J),     4 )
            SO2(I,J, 2) = REAL( CURR_GC_S_SO2(I,J),  4 ) /
     &                    REAL( CURR_COUNT(I,J),     4 )
            SO2(I,J, 3) = REAL( CURR_AMF(I,J),       4 ) /
     &                    REAL( CURR_COUNT(I,J),     4 )
            SO2(I,J, 4) = REAL( CURR_OMI_V_SO2(I,J), 4 ) /
     &                    REAL( CURR_COUNT(I,J),     4 )
            SO2(I,J, 5) = REAL( CURR_OMI_S_SO2(I,J), 4 ) /
     &                    REAL( CURR_COUNT(I,J),     4 )
            SO2(I,J, 6) = REAL( CURR_V_SO2(I,J),     4 ) /
     &                    REAL( CURR_COUNT(I,J),     4 )
            SO2(I,J, 7) = REAL( CURR_DIFF_SO2(I,J),  4 ) /
     &                    REAL( CURR_COUNT(I,J),     4 )
            SO2(I,J, 8) = REAL( CURR_FORCING(I,J),   4 )
            SO2(I,J, 9) = REAL( CURR_COST(I,J),      4 )
            SO2(I,J,10) = REAL( CURR_COUNT(I,J),     4 )

            ALL_GC_V_SO2(I,J)  = ALL_GC_V_SO2(I,J)  + SO2(I,J, 1)
            ALL_GC_S_SO2(I,J)  = ALL_GC_S_SO2(I,J)  + SO2(I,J, 2)
            ALL_AMF(I,J)       = ALL_AMF(I,J)       + SO2(I,J, 3)
            ALL_OMI_V_SO2(I,J) = ALL_OMI_V_SO2(I,J) + SO2(I,J, 4)
            ALL_OMI_S_SO2(I,J) = ALL_OMI_S_SO2(I,J) + SO2(I,J, 5)
            ALL_V_SO2(I,J)     = ALL_V_SO2(I,J)     + SO2(I,J, 6)
            ALL_DIFF_SO2(I,J)  = ALL_DIFF_SO2(I,J)  + SO2(I,J, 7)
            ALL_COUNT(I,J)     = ALL_COUNT(I,J)     + 1
         END IF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL BPCH2( IU_OMSO2,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  26,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     10,        I0+1,
     &            J0+1,      1,         SO2 )

      ! Close file
      CLOSE( IU_OMSO2 )

      ! Return to the calling routines
      END SUBROUTINE MAKE_CURRENT_OMI_SO2
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE MAKE_AVERAGE_OMI_SO2
!
!*****************************************************************************
!  Subroutine MAKE_AVERAGE_OMI_SO2 output some dignostic data for
!  simulation time
!  (ywang, 07/19/14)
!
!*****************************************************************************
!
      ! Reference to f90 module
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,    ONLY : EXPAND_NAME
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
      USE FILE_MOD,          ONLY : IU_OMSO2,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TIME_MOD,          ONLY : GET_NYMD
      USE TIME_MOD,          ONLY : GET_TAUb, GET_TAUe

      REAL*4, PARAMETER    :: UNDEF = -999.0

      ! Local Variables
      INTEGER              :: I,    I0, J,  J0, L
      REAL*4               :: SO2(IIPAR,JJPAR,10)
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      REAL*4               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      CHARACTER(LEN=255)   :: OUTPUT_FILE
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE

      !=================================================================
      ! MAKE_AVERAGE_OMI_SO2 begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_FILE = 'gctm.omi.so2.ave.YYYYMMDD.NN'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM diag File: ' //
     &           'OMI SO2'
      UNIT     = 'DU'
      CATEGORY = 'IJ-AVG-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_FILE )

      ! Replace YYMMDD token w/ actucl value
      CALL EXPAND_DATE( FILENAME, GET_NYMD(), 9999 )

      ! Replace NN token w/ actual value
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add DIAGADJ_DIR prefix to FILENAME
      FILENAME = TRIM( DIAGADJ_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_AVERAGE_OMI_SO2:  Writing ', a )

      ! Open file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_OMSO2, FILENAME, TITLE )

      SO2 = 0.0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         IF ( ALL_COUNT(I,J) > 0 ) THEN
            SO2(I,J, 1) = REAL( ALL_GC_V_SO2(I,J),  4 ) /
     &                    REAL( ALL_COUNT(I,J),     4 )
            SO2(I,J, 2) = REAL( ALL_GC_S_SO2(I,J),  4 ) /
     &                    REAL( ALL_COUNT(I,J),     4 )
            SO2(I,J, 3) = REAL( ALL_AMF(I,J),       4 ) /
     &                    REAL( ALL_COUNT(I,J),     4 )
            SO2(I,J, 4) = REAL( ALL_OMI_V_SO2(I,J), 4 ) /
     &                    REAL( ALL_COUNT(I,J),     4 )
            SO2(I,J, 5) = REAL( ALL_OMI_S_SO2(I,J), 4 ) /
     &                    REAL( ALL_COUNT(I,J),     4 )
            SO2(I,J, 6) = REAL( ALL_V_SO2(I,J),     4 ) /
     &                    REAL( ALL_COUNT(I,J),     4 )
            SO2(I,J, 7) = REAL( ALL_DIFF_SO2(I,J),  4 ) /
     &                    REAL( ALL_COUNT(I,J),     4 )
            SO2(I,J, 8) = REAL( ALL_FORCING(I,J),   4 )
            SO2(I,J, 9) = REAL( ALL_COST(I,J),      4 )
            SO2(I,J,10) = REAL( ALL_COUNT(I,J),     4 )
         END IF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL BPCH2( IU_OMSO2,    MODELNAME,  LONRES,     LATRES,
     &            HALFPOLAR, CENTER180,  CATEGORY,   26,
     &            UNIT,      GET_TAUb(), GET_TAUe(), RESERVED,
     &            IIPAR,     JJPAR,      10,         I0+1,
     &            J0+1,      1,          SO2 )

      ! Close file
      CLOSE( IU_OMSO2 )

      ! Return to the calling routines
      END SUBROUTINE MAKE_AVERAGE_OMI_SO2
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE GET_OBS
!
!*****************************************************************************
!  Subroutine GET_OBS finds all obsevations in the current time window
!  and simulation area.
!  (ywang, 07/17/14)
!
!  Module variable  as Input:
!  ===========================================================================
!  ( 1) N_SO2     (INTEGER) : Number of observation in current day
!  Arguements as Output:
!
!  Module variable as Output:
!  ===========================================================================
!  ( 1) FLAGS     (LOGICAL) : Whether or not a speicific obsevation is
!  in current time window and simulation area.
!  ( 2) N_CURR    (INTEGER) : Number of obsevations in the current time
!  window and simulation area.
!
!*****************************************************************************
!
      ! Reference to f90 module
      USE GRID_MOD,      ONLY : GET_XEDGE, GET_YEDGE
      USE TIME_MOD,      ONLY : GET_JD,   GET_TAU
      USE TIME_MOD,      ONLY : GET_NYMD, GET_NHMS

      ! Local variables
      REAL*8                 :: HALF_TIME_WINDOW
      REAL*8                 :: WINDOW_BEGIN,    WINDOW_END
      REAL*8                 :: JD85,            JD93,      JD93_85
      REAL*8                 :: CURRENT_TAU
      INTEGER                :: NT

#if defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_SD )
      REAL*8, SAVE           :: XEDGE_MIN, XEDGE_MAX
      REAL*8, SAVE           :: YEDGE_MIN, YEDGE_MAX
#endif

      !=================================================================
      ! GET_OBS begins here!
      !=================================================================

      ! Get the difference between JD93 and JD85
      ! In GEOS-Chem, it is since 1/1/1985, while in OMI, it is since
      ! 1/1/1993
      JD85    = GET_JD( 19850000, 000000 )
      JD93    = GET_JD( 19930000, 000000 )
      JD93_85 = JD93 - JD85
      JD93_85 = JD93_85 * 24D0 ! days => hours


      ! Get current GEOS-Chem TAU
      CURRENT_TAU = GET_TAU()

      ! Change current GEOS-Chem TAU into OMI TAU
      CURRENT_TAU = CURRENT_TAU - JD93_85

      ! Change TAU units ( hours => second )
      CURRENT_TAU = CURRENT_TAU * 3600D0

      ! Get half time window
      HALF_TIME_WINDOW = TIME_WINDOW / 2D0
      HALF_TIME_WINDOW = HALF_TIME_WINDOW * 60D0 ! ( minute => second )

      ! Get current time window
      WINDOW_BEGIN = CURRENT_TAU - HALF_TIME_WINDOW
      WINDOW_END   = CURRENT_TAU + HALF_TIME_WINDOW 

#if defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_SD )
         XEDGE_MIN = GET_XEDGE( 1 )
         XEDGE_MAX = GET_XEDGE( IIPAR+1 )
         YEDGE_MIN = GET_YEDGE( 1 )
         YEDGE_MAX = GET_YEDGE( JJPAR+1 )
         PRINT*, 'Nested region edge limit'
         PRINT*, 'XEDGE_MIN: ', XEDGE_MIN
         PRINT*, 'XEDGE_MAX: ', XEDGE_MAX
         PRINT*, 'YEDGE_MIN: ', YEDGE_MIN
         PRINT*, 'YEDGE_MAX: ', YEDGE_MAX
#endif
      N_CURR = 0
      DO NT = 1, N_SO2, 1

         IF (      ( OMI_SO2(NT)%TIME >= WINDOW_BEGIN )
     &       .AND. ( OMI_SO2(NT)%TIME <  WINDOW_END   )
#if defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_SD )
     &       .AND. ( OMI_SO2(NT)%LON  >= XEDGE_MIN    )
     &       .AND. ( OMI_SO2(NT)%LON  <= XEDGE_MAX    )
     &       .AND. ( OMI_SO2(NT)%LAT  >= YEDGE_MIN    )
     &       .AND. ( OMI_SO2(NT)%LAT  <= YEDGE_MAX    )
#endif 
     &                                                  ) THEN


            FLAGS(NT) = .TRUE.

            N_CURR    = N_CURR + 1

         ELSE

            FLAGS(NT) = .FALSE.

         END IF

      END DO

      WRITE(6, 100) N_CURR, GET_NHMS()
 100  FORMAT('      Number of OMI SO2 observations: ' I10 ' at ' I10.6)

      ! Return to calling program
      END SUBROUTINE GET_OBS
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE CHECK( STATUS, LOCATION )
!
!*****************************************************************************
!  Subroutine CHECK checks the status of calls to netCDF libraries
!  routines
!  (dkh, 02/15/09) 
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) STATUS    (INTEGER) : Completion status of netCDF library call    
!  (2 ) LOCATION  (INTEGER) : Location at which netCDF library call was
!  made   
!     
!  NOTES:
!
!*****************************************************************************
!
      ! Reference to f90 modules 
      USE ERROR_MOD,    ONLY  : ERROR_STOP
      USE NETCDF

      ! Arguments
      INTEGER, INTENT(IN)    :: STATUS
      INTEGER, INTENT(IN)    :: LOCATION

      !=================================================================
      ! CHECK begins here!
      !=================================================================

      IF ( STATUS /= NF90_NOERR ) THEN
        WRITE(6,*) TRIM( NF90_STRERROR( STATUS ) )
        WRITE(6,*) 'At location = ', LOCATION
        CALL ERROR_STOP('netCDF error', 'omi_so2_mod')
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHECK
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE NCIO_1D (NCID, VAR1D, VARNAME, LONGNAME, VARUNIT,
     &                    DIMID, DIMV, INDEX)

      ! References to F90 modules
      USE NETCDF

      ! Arguments
      INTEGER, INTENT(IN)          :: NCID, DIMID, DIMV, INDEX
      REAL*4,  INTENT(IN)          :: VAR1D(DIMV)
      CHARACTER(LEN=*), INTENT(IN) :: VARNAME, LONGNAME, VARUNIT

      ! Local variable            
      INTEGER                     :: DIMS(1)
      INTEGER                     :: VAR_ID

      CALL CHECK( NF90_REDEF(NCID), INDEX )

      DIMS(1) = DIMID
      CALL CHECK( NF90_DEF_VAR(NCID, TRIM(VARNAME),
     &            NF90_FLOAT, DIMS, VAR_ID), INDEX)
      CALL CHECK( NF90_PUT_ATT(NCID,VAR_ID,"name",TRIM(LONGNAME)),INDEX)
      CALL CHECK( NF90_PUT_ATT(NCID,VAR_ID,"units",TRIM(VARUNIT)),INDEX)
      CALL CHECK( NF90_ENDDEF(NCID), INDEX )
      CALL CHECK( NF90_PUT_VAR(NCID, VAR_ID, VAR1D ), INDEX )

      END SUBROUTINE NCIO_1D
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE NCIO_1D_DBL (NCID, VAR1D, VARNAME, LONGNAME, VARUNIT,
     &                        DIMID, DIMV, INDEX)

      ! References to F90 modules
      USE NETCDF

      ! Arguments
      INTEGER, INTENT(IN)          :: NCID, DIMID, DIMV, INDEX
      REAL*8,  INTENT(IN)          :: VAR1D(DIMV)
      CHARACTER(LEN=*), INTENT(IN) :: VARNAME, LONGNAME, VARUNIT

      ! Local variable            
      INTEGER                     :: DIMS(1)
      INTEGER                     :: VAR_ID

      CALL CHECK( NF90_REDEF(NCID), INDEX )

      DIMS(1) = DIMID
      CALL CHECK( NF90_DEF_VAR(NCID, TRIM(VARNAME),
     &            NF90_DOUBLE, DIMS, VAR_ID), INDEX)
      CALL CHECK( NF90_PUT_ATT(NCID,VAR_ID,"name",TRIM(LONGNAME)),INDEX)
      CALL CHECK( NF90_PUT_ATT(NCID,VAR_ID,"units",TRIM(VARUNIT)),INDEX)
      CALL CHECK( NF90_ENDDEF(NCID), INDEX )
      CALL CHECK( NF90_PUT_VAR(NCID, VAR_ID, VAR1D ), INDEX )

      END SUBROUTINE NCIO_1D_DBL
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE NCIO_1D_INT (NCID, VAR1D, VARNAME, LONGNAME, VARUNIT,
     &                        DIMID, DIMV, INDEX)

      ! References to F90 modules
      USE NETCDF

      ! Arguments
      INTEGER, INTENT(IN)          :: NCID, DIMID, DIMV, INDEX
      INTEGER, INTENT(IN)          :: VAR1D(DIMV)
      CHARACTER(LEN=*), INTENT(IN) :: VARNAME, LONGNAME, VARUNIT

      ! Local variable            
      INTEGER                     :: DIMS(1)
      INTEGER                     :: VAR_ID

      CALL CHECK( NF90_REDEF(NCID), INDEX )

      DIMS(1) = DIMID
      CALL CHECK( NF90_DEF_VAR(NCID, TRIM(VARNAME),
     &            NF90_INT, DIMS, VAR_ID), INDEX)
      CALL CHECK( NF90_PUT_ATT(NCID,VAR_ID,"name",TRIM(LONGNAME)),INDEX)
      CALL CHECK( NF90_PUT_ATT(NCID,VAR_ID,"units",TRIM(VARUNIT)),INDEX)
      CALL CHECK( NF90_ENDDEF(NCID), INDEX )
      CALL CHECK( NF90_PUT_VAR(NCID, VAR_ID, VAR1D ), INDEX )

      END SUBROUTINE NCIO_1D_INT
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE GET_LOOKUP
!
!*****************************************************************************
! Subroutine GET_LOOKUP read data from "AMFs_lookup_310_320nm.he5"
! (ywang, 07/07/15)
!
!*****************************************************************************
!
      ! References to F90 modules
      USE HDF5

      ! Parameters
      INTEGER, PARAMETER :: TOMS_OFFSET = 10 ! TOMS ozone (m325)
      INTEGER, PARAMETER :: LLP_OFFSET  =  0 ! Lower level pressure (1.0)
      INTEGER, PARAMETER :: WL_OFFSET   =  3 ! Wavelength (313 nm)

      ! Local variables
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=255) :: DSETNAME
      INTEGER            :: DSETRANK
      INTEGER(HSIZE_T)   :: OFFSET_1(1), OFFSET_3(3), OFFSET_5(5)
      INTEGER(HSIZE_T)   :: OFFSET_6(6)
      INTEGER(HSIZE_T)   :: COUNT_1(1),  COUNT_3(3),  COUNT_5(5)
      INTEGER(HSIZE_T)   :: COUNT_6(6)
      INTEGER            :: MEMRANK
      INTEGER(HSIZE_T)   :: DIMSM_1(1),  DIMSM_2(2),  DIMSM_3(3)

      INTEGER            :: I, J, K
      REAL*4             :: TMP(1)

      FILENAME = "./data/VLIDORT/AMFs_lookup_310_320nm.he5"

      ! Ozone C0
      DSETNAME = "/Cross sections/Ozone C0"
      DSETRANK    = 1
      OFFSET_1(1) = WL_OFFSET
      COUNT_1(1)  = 1
      MEMRANK     = 1
      DIMSM_1(1)  = 1
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_1,
     &                     COUNT_1,  MEMRANK,  DIMSM_1,  TMP      )
      O3_C0 = TMP(1)

      ! Ozone C1
      DSETNAME = "/Cross sections/Ozone C1"
      DSETRANK    = 1
      OFFSET_1(1) = WL_OFFSET
      COUNT_1(1)  = 1
      MEMRANK     = 1
      DIMSM_1(1)  = 1
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_1,
     &                     COUNT_1,  MEMRANK,  DIMSM_1,  TMP      )
      O3_C1 = TMP(1)

      ! Ozone C2
      DSETNAME = "/Cross sections/Ozone C2"
      DSETRANK    = 1
      OFFSET_1(1) = WL_OFFSET
      COUNT_1(1)  = 1
      MEMRANK     = 1
      DIMSM_1(1)  = 1
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_1,
     &                     COUNT_1,  MEMRANK,  DIMSM_1,  TMP      )
      O3_C2 = TMP(1)

      ! SZA
      DSETNAME    = "/Grid/SZA"
      DSETRANK    = 1
      OFFSET_1(1) = 0
      COUNT_1(1)  = N_SZA
      MEMRANK     = 1
      DIMSM_1     = N_SZA
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_1,
     &                     COUNT_1,  MEMRANK,  DIMSM_1,  SZA      )

      ! VZA
      DSETNAME    = "/Grid/VZA"
      DSETRANK    = 1
      OFFSET_1(1) = 0
      COUNT_1(1)  = N_VZA
      MEMRANK     = 1
      DIMSM_1     = N_VZA
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_1,
     &                     COUNT_1,  MEMRANK,  DIMSM_1,  VZA      )

      !-------------------------
      ! I0, I1, I2, Ir
      !-------------------------
      DSETRANK    = 5
      OFFSET_5    = 0
      OFFSET_5(1) = TOMS_OFFSET
      OFFSET_5(2) = LLP_OFFSET
      OFFSET_5(5) = WL_OFFSET
      COUNT_5     = 1
      COUNT_5(3)  = N_SZA
      COUNT_5(4)  = N_VZA
      MEMRANK     = 2
      DIMSM_2(1)  = N_SZA
      DIMSM_2(2)  = N_VZA
      ! I0
      DSETNAME = "/Intensity/I0"
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_5,
     &                     COUNT_5,  MEMRANK,  DIMSM_2,  I0       )
      ! I1
      DSETNAME = "/Intensity/I1"
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_5,
     &                     COUNT_5,  MEMRANK,  DIMSM_2,  I1       )
      ! I2
      DSETNAME = "/Intensity/I2"
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_5,
     &                     COUNT_5,  MEMRANK,  DIMSM_2,  I2       )
      ! Ir
      DSETNAME = "/Intensity/Ir"
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_5,
     &                     COUNT_5,  MEMRANK,  DIMSM_2,  Ir       )

      ! Sb
      DSETNAME = "/Intensity/Sb"
      DSETRANK    = 3
      OFFSET_3(1) = TOMS_OFFSET
      OFFSET_3(2) = LLP_OFFSET
      OFFSET_3(3) = WL_OFFSET
      COUNT_3     = 1
      MEMRANK     = 1
      DIMSM_1(1)  = 1
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_3,
     &                     COUNT_3,  MEMRANK,  DIMSM_1,  TMP      )
      Sb = TMP(1)

      ! Factor
      DSETNAME = "/Jacobians/Factor"
      DSETRANK    = 1
      OFFSET_1(1) = 0
      COUNT_1(1)  = 1
      MEMRANK     = 1
      DIMSM_1(1)  = 1
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_1,
     &                     COUNT_1,  MEMRANK,  DIMSM_1,  TMP      )
      Factor = TMP(1)

      !-------------------------
      ! dI0, dI1, dI2, dIr
      !-------------------------
      DSETRANK    = 6
      OFFSET_6    = 0
      OFFSET_6(1) = TOMS_OFFSET
      OFFSET_6(2) = LLP_OFFSET
      OFFSET_6(5) = WL_OFFSET
      COUNT_6     = 1
      COUNT_6(3)  = N_SZA
      COUNT_6(4)  = N_VZA
      COUNT_6(6)  = N_LEV
      MEMRANK     = 3
      DIMSM_3(1)  = N_SZA
      DIMSM_3(2)  = N_VZA
      DIMSM_3(3)  = N_LEV
      ! dI0
      DSETNAME = "/Jacobians/dI0"
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_6,
     &                     COUNT_6,  MEMRANK,  DIMSM_3,  dI0      )
      ! dI1
      DSETNAME = "/Jacobians/dI1"
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_6,
     &                     COUNT_6,  MEMRANK,  DIMSM_3,  dI1      )
      ! dI2
      DSETNAME = "/Jacobians/dI2"
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_6,
     &                     COUNT_6,  MEMRANK,  DIMSM_3,  dI2      )
      ! dIr
      DSETNAME = "/Jacobians/dIr"
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_6,
     &                     COUNT_6,  MEMRANK,  DIMSM_3,  dIr      )

      !-------------------------
      ! Profiles
      !-------------------------
      DSETRANK    = 3
      OFFSET_3(1) = TOMS_OFFSET
      OFFSET_3(2) = LLP_OFFSET
      OFFSET_3(3) = 0
      COUNT_3     = 1
      COUNT_3(3)  = N_LEV
      MEMRANK     = 1
      DIMSM_1(1)  = N_LEV
      ! Air profile
      DSETNAME = "/Profiles/Air profile"
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_3,
     &                     COUNT_3,  MEMRANK,  DIMSM_1,  AIR_P    )
      ! Altitude
      DSETNAME = "/Profiles/Altitude"
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_3,
     &                     COUNT_3,  MEMRANK,  DIMSM_1,  ALTITUDE )
      ! Ozone profile
      DSETNAME = "/Profiles/Ozone profile"
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_3,
     &                     COUNT_3,  MEMRANK,  DIMSM_1,  O3_P     )
      ! Temperature profile
      DSETNAME = "/Profiles/Temperature profile"
      CALL READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET_3,
     &                     COUNT_3,  MEMRANK,  DIMSM_1,  TEMP_P   )

      VLIDORT_PRESS =
     &       (/   0.00772193,    0.0224936,    0.112247,    0.356718,
     &            0.972622,      2.44691,      5.87901,    13.3036,
     &           23.5377,       33.5917,      47.2958,     66.0496,
     &           84.3372,       99.0540,     116.339,     136.856,
     &          161.245,       189.980,      224.898,     265.986,
     &          313.881,       358.581,      396.921,     434.870,
     &          473.106,       510.561,      548.918,     586.483,
     &          624.478,       662.736,      693.982,     720.011,
     &          744.951,       770.585,      795.940,     818.928,
     &          837.292,       851.804,      867.578,     882.513,
     &          898.746,       914.113,      928.573,     944.351,
     &          960.346,       975.394,      990.632,    1006.06      /)


      IF ( .TRUE. ) THEN

         WRITE(6, '(A)') REPEAT( '-', 79 )
         WRITE(6, '(3X,A11,E15.5,X,A8)')
     &            'Ozone C0 = ', O3_C0, 'cm^2/mol'
         WRITE(6, '(3X,A11,E15.5,X,A12)')
     &             'Ozone C1 = ', O3_C1, 'cm^2/(mol K)'
         WRITE(6, '(3X,A11,E15.5,X,A14)')
     &             'Ozone C2 = ', O3_C2, 'cm^2/(mol K^2)'

         !--------------------
         ! I0, I1, I2, Ir
         !--------------------
         CALL PRINT_INTENSITY( I0, "I0" )
         CALL PRINT_INTENSITY( I1, "I1" )
         CALL PRINT_INTENSITY( I2, "I2" )
         CALL PRINT_INTENSITY( Ir, "Ir" )

         ! Sb
         WRITE(6, '(A)') REPEAT( '-', 79 )
         WRITE(6, '(20X,A5,F10.6)') "Sb = ", Sb

         ! Factor
         WRITE(6, '(A)') REPEAT( '-', 79 )
         WRITE(6, '(20X,A9,E10.3)') "Factor = ", Factor

         !---------------------
         ! dI0, dI1, dI2, dIr
         !---------------------
         CALL PRINT_JACOBIANS( dI0, "dI0" )
         CALL PRINT_JACOBIANS( dI1, "dI1" )
         CALL PRINT_JACOBIANS( dI2, "dI2" )
         CALL PRINT_JACOBIANS( dIr, "dIr" )

         !---------------------
         ! Profiles
         !---------------------
         WRITE(6, '(A)') REPEAT( '-', 79 )
         WRITE(6, 100)
         DO K = 1, N_LEV, 1
            WRITE(6, 110) K, AIR_P(K), ALTITUDE(K), VLIDORT_PRESS(K),
     &                       O3_P(K),  TEMP_P(K)
         END DO

 100  FORMAT(7X, 'Air profile', 3X, 'Altitude', 3X, 'Pressure',
     &3X, 'Ozone profile', 3X, 'Temperature profile')
 110  FORMAT(I5, E13.5, F10.3, F10.3, 3X, E13.5, 6X, F10.2)
      END IF

      ! Return to the calling routine
      END SUBROUTINE GET_LOOKUP
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE PRINT_INTENSITY( INTENSITY, LABEL )

      REAL*4,           INTENT(IN) :: INTENSITY(N_SZA,N_VZA)
      CHARACTER(LEN=2), INTENT(IN) :: LABEL
      INTEGER                      :: I, J

      WRITE(6, '(A)') REPEAT( '-', 79 )
      WRITE(6, '(35X,A2)') LABEL
      WRITE(6, 100) "VZA", (VZA(J), J = 1, N_VZA, 1)
      WRITE(6, '(2X,3A)') "SZA"
      DO I = 1, N_SZA, 1
         WRITE(6, 110) SZA(I), (INTENSITY(I,J), J = 1, N_VZA, 1)
      END DO

 100  FORMAT(2X, A, 8F13.1)
 110  FORMAT(F6.1, 8E13.5)

      END SUBROUTINE PRINT_INTENSITY
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE PRINT_JACOBIANS( JACOBIANS, LABEL )

      REAL*4,           INTENT(IN) :: JACOBIANS(N_SZA, N_VZA, N_LEV)
      CHARACTER(LEN=3), INTENT(IN) :: LABEL
      INTEGER                      :: I, J, K

      WRITE(6, '(A)') REPEAT( '-', 79 )
      DO K = 1, N_LEV, 1
!      DO K = 1, 1, 1
         WRITE(6, '(20X,A13,I3)') TRIM(LABEL)//" at level ", K
         WRITE(6, 100) "VZA", (VZA(J), J = 1, N_VZA, 1)
         WRITE(6, '(2X,3A)') "SZA"
         DO I = 1, N_SZA, 1
            WRITE(6, 110) SZA(I), (JACOBIANS(I,J,K), J = 1, N_VZA, 1)
         END DO
         WRITE(6, *)
      END DO

 100  FORMAT(2X, A, 8F13.1)
 110  FORMAT(F6.1, 8E13.5)

      END SUBROUTINE PRINT_JACOBIANS
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE READ_HYPERSLAB( FILENAME, DSETNAME, DSETRANK, OFFSET,
     &                           COUNT,    MEMRANK,  DIMSM,    OUT_DATA
     &                         )
!
!*****************************************************************************
! Subroutine READ_HYPERSLAB read hyperslab from a HDF5 file.
! (ywang, 07/06/15)
!
! Arguements as Input:
! ============================================================================
! ( 1) FILENAME    CHARACTER(LEN=255) : File name
! ( 2) DSETNAME    CHARACTER(LEN=255) : Dataset name
! ( 3) DSETRANK    INTEGER            : Dataset rank in the file
! ( 4) OFFSET      INTEGER(HSIZE_T)   : Hyperslab offset in the file
! ( 5) COUNT       INTEGER(HSIZE_T)   : Size of the hyperslab in the
!                                       file
! ( 6) MEMRANK     INTEGER            : Dataset rank in memory
! ( 7) DIMSM       INTEGER            : Dataset dimensions in memory
!
! Arguements as Output:
! ===========================================================================
! ( 1) OUT_DATA    REAL*4             : Output buffer
!
!*****************************************************************************
!
      ! References to F90 modules
      USE HDF5

      ! Parameters
      CHARACTER(LEN=255), INTENT(IN ) :: FILENAME
      CHARACTER(LEN=255), INTENT(IN ) :: DSETNAME
      INTEGER,            INTENT(IN ) :: DSETRANK
      INTEGER(HSIZE_T),   INTENT(IN ) :: OFFSET(DSETRANK)
      INTEGER(HSIZE_T),   INTENT(IN ) :: COUNT(DSETRANK)
      INTEGER,            INTENT(IN ) :: MEMRANK
      INTEGER(HSIZE_T),   INTENT(IN ) :: DIMSM(MEMRANK)

      REAL*4,             INTENT(OUT) :: OUT_DATA(*)

      ! Local variables
      INTEGER            :: IERR
      INTEGER(HID_T)     :: FILE_ID
      INTEGER(HID_T)     :: DSET_ID
      INTEGER(HID_T)     :: DATASPACE
      INTEGER(HID_T)     :: MEMSPACE

      !=================================================================
      ! READ_HYPERSLAB begins here!
      !=================================================================

      PRINT*, " - READ_HYPERSLAB: reading ", TRIM(DSETNAME)

      ! Initialize HDF5 interface
      CALL H5OPEN_F( IERR )

      ! Open HDF5 file
      CALL H5FOPEN_F( FILENAME, H5F_ACC_RDONLY_F, FILE_ID, IERR )

      ! Open dataset
      CALL H5DOPEN_F( FILE_ID, DSETNAME, DSET_ID, IERR )

      ! Get dataspace identifier
      CALL H5DGET_SPACE_F( DSET_ID, DATASPACE, IERR )

      ! Select hyperslab in the dataset
      CALL H5SSELECT_HYPERSLAB_F( DATASPACE, H5S_SELECT_SET_F, OFFSET,
     &                           COUNT, IERR )

      ! Create memory dataspace
      CALL H5SCREATE_SIMPLE_F( MEMRANK, DIMSM, MEMSPACE, IERR )

      ! Read data form hyperslab in the file
      CALL H5DREAD_F( DSET_ID,  H5T_NATIVE_REAL, OUT_DATA, DIMSM, IERR,
     &                MEMSPACE, DATASPACE )

      ! Close the dataspace for the dataset
      CALL H5SCLOSE_F( DATASPACE, IERR )

      ! Close the memoryspace
      CALL H5SCLOSE_F( MEMSPACE, IERR )

      ! Close the dataset
      CALL H5DCLOSE_F( DSET_ID, IERR )

      ! Close the file
      CALL H5FCLOSE_F( FILE_ID, IERR )

      ! Close HD5 interface
      CALL H5CLOSE_F( IERR )

      ! Return to the calling routine
      END SUBROUTINE READ_HYPERSLAB
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE INIT_LOOKUP( RESET_IS_INIT )
!
!*****************************************************************************
! Subroutine INIT_LOOKUP allocate arrays in this module.
! (ywang, 07/06/15)
!*****************************************************************************
!
      ! Reference to F90 modules
      USE ERROR_MOD,     ONLY : ALLOC_ERR

      ! Arguements
      LOGICAL, OPTIONAL    :: RESET_IS_INIT

      ! Local variables
      INTEGER              :: AS
      LOGICAL, SAVE        :: IS_INIT = .FALSE.

      !=================================================================
      ! INIT_LOOKUP begins here
      !=================================================================

      ! Reset IS_INIT as false when clean up
      IF ( PRESENT( RESET_IS_INIT ) ) THEN
         IS_INIT = .FALSE.
         WRITE(6, '(A)') ' - INIT_LOOKUP: IS_INIT is reset to FALSE'
         RETURN
      END IF

      ! Return if we have already initialized
      IF ( IS_INIT ) RETURN

      ! SZA
      ALLOCATE ( SZA(N_SZA), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SZA' )

      ! VZA
      ALLOCATE ( VZA(N_VZA), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VZA' )

      ! I0
      ALLOCATE ( I0(N_SZA,N_VZA), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'I0' )

      ! I1
      ALLOCATE ( I1(N_SZA,N_VZA), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'I1' )

      ! I2
      ALLOCATE ( I2(N_SZA,N_VZA), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'I2' )

      ! Ir
      ALLOCATE ( Ir(N_SZA,N_VZA), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Ir' )

      ! dI0
      ALLOCATE ( dI0(N_SZA,N_VZA,N_LEV), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'dI0' )

      ! dI1
      ALLOCATE ( dI1(N_SZA,N_VZA,N_LEV), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'dI1' )

      ! dI2
      ALLOCATE ( dI2(N_SZA,N_VZA,N_LEV), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'dI2' )

      ! dIr
      ALLOCATE ( dIr(N_SZA,N_VZA,N_LEV), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'dIr' )

      ! AIR_P
      ALLOCATE ( AIR_P(N_LEV), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AIR_P' )

      ! ALTITUDE
      ALLOCATE ( ALTITUDE(N_LEV), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALTITUDE' )

      ! O3_P
      ALLOCATE ( O3_P(N_LEV), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'O3_P' )

      ! TEMP_P
      ALLOCATE ( TEMP_P(N_LEV), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TEMP_P' )

      ! VLIDORT_PRESS
      ALLOCATE ( VLIDORT_PRESS(N_LEV), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VLIDORT_PRESS' )

      ! Reset IS_INIT so we do not allocate arrays again
      IS_INIT = .TRUE.

      WRITE(6, '(A)')
     &     ' - INIT_LOOKUP: Initialized the lookup table arrays'

      ! Return to the calling routine
      END SUBROUTINE INIT_LOOKUP
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE CLEANUP_LOOKUP

      ! Local varibales
      LOGICAL         :: RESET_IS_INIT = .TRUE.

      ! Clean up the output arrays of this module
      IF ( ALLOCATED( SZA      ) ) DEALLOCATE( SZA      )
      IF ( ALLOCATED( VZA      ) ) DEALLOCATE( VZA      )
      IF ( ALLOCATED( I0       ) ) DEALLOCATE( I0       )
      IF ( ALLOCATED( I1       ) ) DEALLOCATE( I1       )
      IF ( ALLOCATED( I2       ) ) DEALLOCATE( I2       )
      IF ( ALLOCATED( Ir       ) ) DEALLOCATE( Ir       )
      IF ( ALLOCATED( dI0      ) ) DEALLOCATE( dI0      )
      IF ( ALLOCATED( dI1      ) ) DEALLOCATE( dI1      )
      IF ( ALLOCATED( dI2      ) ) DEALLOCATE( dI2      )
      IF ( ALLOCATED( dIr      ) ) DEALLOCATE( dIr      )
      IF ( ALLOCATED( AIR_P    ) ) DEALLOCATE( AIR_P    )
      IF ( ALLOCATED( ALTITUDE ) ) DEALLOCATE( ALTITUDE )
      IF ( ALLOCATED( O3_P     ) ) DEALLOCATE( O3_P     )
      IF ( ALLOCATED( TEMP_P   ) ) DEALLOCATE( TEMP_P   )
      IF ( ALLOCATED( VLIDORT_PRESS )) DEALLOCATE( VLIDORT_PRESS )

      ! Reset IS_INIT in routine INIT_ND49
      CALL INIT_LOOKUP( RESET_IS_INIT )

      END SUBROUTINE CLEANUP_LOOKUP
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE SPLINE(M,N,X,Y,T,SS)
      !***************************************************************** 
      !*       THIS SUBROUTINE IS VERTICAL INTERPOLATIONAL PROGRAM.    * 
      !*****************************************************************
      !Program name: SPLINE 
      !Input variables: 
      !M: vertical levels to be interpolated 
      !N: vertical levels of observed data. 
      !X(N): 1-dimensional array storing height or pressure or Ln(p)
      !etc. 
      !    of observed data 
      !Y(N): 1-dimensional array storing observed data (say temperature) 
      !T(M): 1-dimensional array storing height or pressure or Ln(p)
      !etc. 
      !    of levels to be interpolated 
      !Output variables: 
      !SS(M): 1-dimensional array storing interpolated data (say
      !temperature) 
      IMPLICIT NONE

      INTEGER, INTENT(IN ) :: M, N
      REAL*8,  INTENT(IN ) :: X(N),Y(N)
      REAL*8,  INTENT(IN ) :: T(M)
      REAL*8,  INTENT(OUT) :: SS(M)
      REAL*8               :: S2(N),H(N),H1(N),B(N),C(N)
      REAL*8               :: DELY(N),S3(N),DELSQY(N)
      INTEGER              :: N1, I, J, JJ
      REAL*8               :: HT1, HT2, PROD, DELSQS

      N1=N-1
      DO I=1,N1
      H(I)=X(I+1)-X(I)
      DELY(I)=(Y(I+1)-Y(I))/H(I)
      enddo
      DO I=2,N1
      H1(I)=H(I-1)+H(I)
      B(I)=0.5*H(I-1)/H1(I)
      DELSQY(I)=(DELY(I)-DELY(I-1))/H1(I)
      S2(I)=2.0*DELSQY(I)
      C(I)=3.0*DELSQY(I)
      enddo
      S2(1)=0.0
      S2(N)=0.0
      DO 30 JJ=1,26
      DO 30 I=2,N1
      S2(I)=(C(I)-B(I)*S2(I-1)-(0.5-B(I))*S2(I+1)-S2(I))*1.0717968+S2(I) 
30      CONTINUE
      DO 40 I=1,N1
40      S3(I)=(S2(I+1)-S2(I))/H(I)
      DO 50 J=1,M
      I=1
      IF((T(J)-X(I)).LE.0.0) GOTO 17
      IF((T(J)-X(N)).LT.0.0) GOTO 57
      GOTO 59
56      IF((T(J)-X(I)).LT.0.0) GOTO 60
      IF((T(J)-X(I)).EQ.0.0) GOTO 17
57      I=I+1
      GOTO 56
59      I=N
60      I=I-1
17      HT1=T(J)-X(I)
        HT2=T(J)-X(I+1)
        PROD=HT1*HT2
        DELSQS=(2.0*S2(I)+S2(I+1)+HT1*S3(I))/6.0
        SS(J)=Y(I)+HT1*DELY(I)+PROD*DELSQS
50      CONTINUE

      RETURN
      END SUBROUTINE SPLINE
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE INTERPOLATION
!
!*****************************************************************************
! Suboutine INTERPOLATION do bilinear interpolation for I0, I1, I2, Ir,
! dI0, dI1, dIr, dIr accoring to SZA and VZA
! (ywang, 07/15/18)
!*****************************************************************************
!
      ! Reference to F90 module
      USE ERROR_MOD,         ONLY  : ERROR_STOP

      REAL*4          :: DENOM
      REAL*4          :: C11, C12, C21, C22

      CALL DETERMINE_SZA_VZA

      DENOM = ( SZA(SZA_MAX) - SZA(SZA_MIN) ) *
     &        ( VZA(VZA_MAX) - VZA(VZA_MIN) )

      C11 = ( SZA(SZA_MAX) - OMI_SZA      ) *
     &      ( VZA(VZA_MAX) - OMI_VZA      ) /
     &        DENOM
      C12 = ( SZA(SZA_MAX) - OMI_SZA      ) *
     &      ( OMI_VZA      - VZA(VZA_MIN) ) /
     &        DENOM
      C21 = ( OMI_SZA      - SZA(SZA_MIN) ) *
     &      ( VZA(VZA_MAX) - OMI_VZA      ) /
     &        DENOM
      C22 = ( OMI_SZA      - SZA(SZA_MIN) ) *
     &      ( OMI_VZA      - VZA(VZA_MIN) ) /
     &        DENOM

      ! for debug
      IF ( (C11 < 0.0) .OR. (C12 < 0.0) .OR.
     &     (C21 < 0.0) .OR. (C22 < 0.0)     ) THEN
         CALL ERROR_STOP( 'Bilinear interpolation is wrong',
     &                    'INTERPOLATION' )
      END IF

      ! I0
      INTER_I0 = I0(SZA_MIN,VZA_MIN) * C11 + I0(SZA_MIN,VZA_MAX) * C12 +
     &           I0(SZA_MAX,VZA_MIN) * C21 + I0(SZA_MAX,VZA_MAX) * C22
      ! I1
      INTER_I1 = I1(SZA_MIN,VZA_MIN) * C11 + I1(SZA_MIN,VZA_MAX) * C12 +
     &           I1(SZA_MAX,VZA_MIN) * C21 + I1(SZA_MAX,VZA_MAX) * C22
      ! I2
      INTER_I2 = I2(SZA_MIN,VZA_MIN) * C11 + I2(SZA_MIN,VZA_MAX) * C12 +
     &           I2(SZA_MAX,VZA_MIN) * C21 + I2(SZA_MAX,VZA_MAX) * C22
      ! Ir
      INTER_Ir = Ir(SZA_MIN,VZA_MIN) * C11 + Ir(SZA_MIN,VZA_MAX) * C12 +
     &           Ir(SZA_MAX,VZA_MIN) * C21 + Ir(SZA_MAX,VZA_MAX) * C22

      ! dI0
      INTER_dI0(:) = dI0(SZA_MIN,VZA_MIN,:) * C11 +
     &               dI0(SZA_MIN,VZA_MAX,:) * C12 +
     &               dI0(SZA_MAX,VZA_MIN,:) * C21 +
     &               dI0(SZA_MAX,VZA_MAX,:) * C22
      ! dI1
      INTER_dI1(:) = dI1(SZA_MIN,VZA_MIN,:) * C11 +
     &               dI1(SZA_MIN,VZA_MAX,:) * C12 +
     &               dI1(SZA_MAX,VZA_MIN,:) * C21 +
     &               dI1(SZA_MAX,VZA_MAX,:) * C22
      ! dI2
      INTER_dI2(:) = dI2(SZA_MIN,VZA_MIN,:) * C11 +
     &               dI2(SZA_MIN,VZA_MAX,:) * C12 +
     &               dI2(SZA_MAX,VZA_MIN,:) * C21 +
     &               dI2(SZA_MAX,VZA_MAX,:) * C22
      ! dIr
      INTER_dIr(:) = dIr(SZA_MIN,VZA_MIN,:) * C11 +
     &               dIr(SZA_MIN,VZA_MAX,:) * C12 +
     &               dIr(SZA_MAX,VZA_MIN,:) * C21 +
     &               dIr(SZA_MAX,VZA_MAX,:) * C22

      ! Return to the calling routine
      END SUBROUTINE INTERPOLATION
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE DETERMINE_SZA_VZA
!
!*****************************************************************************
! Find index for interpolation
!*****************************************************************************
!
      ! Local variables
      INTEGER               :: I

      ! SZA
      DO I = 2, N_SZA, 1
      IF ( OMI_SZA >= SZA(I-1) .AND.
     &     OMI_SZA <= SZA(I)         ) THEN
         SZA_MIN = I - 1
         SZA_MAX = I
      END IF
      END DO

      ! VZA
      DO I = 2, N_VZA, 1
      IF ( OMI_VZA >= VZA(I-1) .AND.
     &     OMI_VZA <= VZA(I)         ) THEN
         VZA_MIN = I - 1
         VZA_MAX = I
      END IF
      END DO

      END SUBROUTINE DETERMINE_SZA_VZA
!
!-----------------------------------------------------------------------------
!
      END MODULE OMI_SO2_OBS_MOD
