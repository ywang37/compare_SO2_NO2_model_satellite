! $Id: omps_so2_obs_mod.f
      MODULE OMPS_SO2_OBS_MOD
!
!*****************************************************************************
! MODULE OMPS_SO2_OBS_MOD contians subroutines necessary to
! 1. Read OMPS L2 SO2 observations
! 2. Compute OMPS-GEOS-Chem difference, cost function, and adjoint
! forcing
!
!    (ywang (yi-wang-4@uiowa.edu), 08/07/2017)
!
!  Module Variables:
!  ===========================================================================
!
!  Module Routines:
!  ===========================================================================
!  ( 1) CALC_OMPS_SO2_FORCE
!  ( 2) CHECK                : Check status for calling netCDF
!  ( 1) INIT_OMPS_SO2_OBS     : Initialize OMPS SO2 observation operator
!  ( 3) WRITE_GC_OMPS_SO2_OBS
!  ( 4) READ_GC_OMPS_SO2_OBS
!  ( 5) MAKE_CURRENT_OMPS_SO2 :
!  (  ) MAKE_AVERAGE_OMPS_SO2 :
!
!  ===========================================================================
!  NOTES:
!  ( 1) The original OMPS L2 SO2 HDFEOS5 data are preprocessed by
!  python code.
!
!*****************************************************************************
!

      IMPLICIT NONE

      ! Header files
#     include "define.h"
#     include "CMN_SIZE"       ! Size parameters

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep centain internal variables
      ! and routines from being seen outside "omps_so2_obs_mod.f"
      !=================================================================

      ! Make everything PRIVATE...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CALC_OMPS_SO2_FORCE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Parameters
      INTEGER, PARAMETER :: TIME_WINDOW = 60 ! (units: min)

      ! Variables

      ! Recored to store OMPS SO2 observations
      TYPE OMPS_SO2_OBS
         REAL    :: LON           ! Latitude (units: degree)
         REAL    :: LAT           ! Longitude (units: dgeree)
         REAL*8  :: TIME          ! Time at start of scan (TAI93) (units: s)
         REAL*8  :: SO2           ! OMPS SO2 whose CMA is closest to GEOS-Chem CMA (units: DU)
         REAL*8  :: ERR           ! SO2 error (units: DU)
      ENDTYPE OMPS_SO2_OBS

      TYPE(OMPS_SO2_OBS), ALLOCATABLE :: OMPS_SO2(:)

      INTEGER                        :: N_SO2
      INTEGER                        :: N_CURR

      LOGICAL,           ALLOCATABLE :: FLAGS(:)

      REAL*8                      :: CURR_GC_SO2(IIPAR, JJPAR)
      REAL*8                      :: CURR_SO2(IIPAR, JJPAR)
      REAL*8                      :: CURR_DIFF_SO2(IIPAR, JJPAR)
      REAL*8                      :: CURR_FORCING(IIPAR, JJPAR)
      REAL*8                      :: CURR_COST(IIPAR, JJPAR)
      INTEGER                     :: CURR_COUNT(IIPAR, JJPAR) ! number of observations in current time window

      REAL*8                      :: ALL_GC_SO2(IIPAR, JJPAR)   = 0D0
      REAL*8                      :: ALL_SO2(IIPAR, JJPAR)      = 0D0
      REAL*8                      :: ALL_DIFF_SO2(IIPAR, JJPAR) = 0D0
      REAL*8                      :: ALL_FORCING(IIPAR, JJPAR)  = 0D0
      REAL*8                      :: ALL_COST(IIPAR, JJPAR)     = 0D0
      INTEGER                     :: ALL_COUNT(IIPAR, JJPAR)    = 0 ! number of observations in simulation time

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAIN" statement
      !=================================================================
      CONTAINS
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE INIT_OMPS_SO2_OBS
!
!*****************************************************************************
!  Subroutine INIT_OMPS_SO2_OBS initialize the OMPS SO2 observation
!   ( 1) Check if SO2 observation is specified by the adjoint input
!   files
!  (ywang, 08/07/17)
!*****************************************************************************
!
      ! Reference to f90 modules
      USE TRACER_MOD,         ONLY : TRACER_NAME
      USE TRACERID_MOD,       ONLY : IDTSO2
      USE ADJ_ARRAYS_MOD,     ONLY : OBS_THIS_TRACER
      USE ERROR_MOD,          ONLY : ERROR_STOP

      !=================================================================
      ! INIT_OMPS_SO2_OBS begins here!
      !=================================================================

      IF ( OBS_THIS_TRACER(IDTSO2) ) THEN

         WRITE( 6, 100 ) IDTSO2, TRACER_NAME(IDTSO2)

      ELSE

         CALL ERROR_STOP( 'Error: Obs SO2 tracer is not specified in
     &OBSERVATION MENU in input.gcadj',
     &                    'INIT_OMPS_SO2_OBS (omps_so2_obs_mod.f)' )

      END IF

 100  FORMAT( 3X, 'Tracer ID: ', I4, 'Use OMPS L2 ', A6 )

      ! Return to the calling routine
      END SUBROUTINE INIT_OMPS_SO2_OBS
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE READ_OMPS_SO2_OBS( YYYYMMDD )
!
!*****************************************************************************
!  Subroutine READ_OMPS_SO2_OBS read OMPS SO2 data preprocessed by
!  OMPS_SO2_L2_preprocess.py if N_CALC > 1 (ref: note 1)
!  (ywang, 08/07/17)
!
!  Arguements as Input:
!  ===========================================================================
!  ( 1) YYYYMMDD    (INTEGER) : Current year-month-day
!  Arguements as Output:
!
!  Module variable as Output:
!  ===========================================================================
!  ( 1) N_SO2       (INTEGER) : Number of OMPS SO2 observations for
!  current day
!  ( 2) OMPS_SO2 (OMPS_SO2_OBS) : OMPS SO2 observations
!
!*****************************************************************************
!
      ! Reference to f90 modules
      USE NETCDF
      USE TIME_MOD,         ONLY : EXPAND_DATE, GET_NYMD

      ! Arguements
      INTEGER, INTENT( IN)   :: YYYYMMDD

      ! Local variables
      INTEGER                :: FID,        N_ID
      INTEGER                :: LON_ID,     LAT_ID,     TIME_ID
      INTEGER                :: SO2_ID
      INTEGER                :: ERR_ID
      INTEGER                :: OPT_ID  ! (ywang, 09/12/14)
      INTEGER                :: FLAG_ID ! (ywang, 09/12/14)
      CHARACTER(LEN=255)     :: DIR
      CHARACTER(LEN=255)     :: READ_FILENAME
      CHARACTER(LEN=4)       :: TMP
      INTEGER, ALLOCATABLE   :: TMPINT(:)
      REAL*4,  ALLOCATABLE   :: TMP4(:)
      REAL*8,  ALLOCATABLE   :: TMP8(:)
      LOGICAL                :: LF

      !=================================================================
      ! READ_OMPS_SO2_OBS begins here!
      !=================================================================

      ! Filename root
      READ_FILENAME = 'OMPS_GC_adj_SO2_YYYYMMDD.nc'

      ! Expand date tokens in filename
      CALL EXPAND_DATE( READ_FILENAME, YYYYMMDD, 9999 )

      ! Construct complete filename
      DIR           = './data/OMPS_SO2/'
      READ_FILENAME = TRIM( DIR ) // TRIM( READ_FILENAME )

      ! Does data file exist? If not, it means no data in the day.
      INQUIRE( FILE = TRIM( READ_FILENAME ), EXIST = LF )
      IF ( .NOT. LF ) THEN

         ! No data
         N_SO2 = 0

         WRITE(6, 120) GET_NYMD()

         RETURN

      END IF

 120  FORMAT(' - READ_OMPS_SO2_OBS: No data file (warning) in', I10)

      ! Print to screen
      WRITE(6, 100) TRIM( READ_FILENAME )
 100  FORMAT(' - READ_OMPS_SO2_OBS: reading file: ', A)


      ! Open file and assign file id (FID)
      CALL CHECK( NF90_OPEN( READ_FILENAME, NF90_NOWRITE, FID ), 0 )

      !-----------------------------------
      ! Get data record IDs
      !-----------------------------------
      CALL CHECK( NF90_INQ_DIMID( FID, "n_obs",   N_ID       ), 100 )

      CALL CHECK( NF90_INQ_VARID( FID, "lon",     LON_ID     ), 101 )
      CALL CHECK( NF90_INQ_VARID( FID, "lat",     LAT_ID     ), 102 )
      CALL CHECK( NF90_INQ_VARID( FID, "TAI93",   TIME_ID    ), 103 )
      CALL CHECK( NF90_INQ_VARID( FID, "SO2",     SO2_ID     ), 104 )
      CALL CHECK( NF90_INQ_VARID( FID, "ERR",     ERR_ID     ), 105 )

      !------------------------------------
      ! Read dimensions
      !------------------------------------

      ! Read number of observations, N_SO2
      CALL CHECK( NF90_INQUIRE_DIMENSION( FID, N_ID, TMP, N_SO2 ), 200 )

      ! Print to screen
      WRITE(6, 110) N_SO2, GET_NYMD()
 110  FORMAT('      Number of OMPS SO2 observations: ' I10, ' in ' I10)
      IF (N_SO2 == 0) THEN
         ! Close the file
         CALL CHECK( NF90_CLOSE( FID ), 9999 )
         RETURN
      END IF

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

      ! Allocate OMPS SO2 observations array
      IF ( ALLOCATED( OMPS_SO2 ) ) DEALLOCATE( OMPS_SO2 )
      ALLOCATE( OMPS_SO2(N_SO2) )

      IF ( ALLOCATED( FLAGS) ) DEALLOCATE( FLAGS )
      ALLOCATE( FLAGS(N_SO2) )

      ! Read longitude
      CALL CHECK( NF90_GET_VAR( FID, LON_ID, TMP4 ), 301 )
      OMPS_SO2(1:N_SO2)%LON = TMP4(1:N_SO2)

      ! Read latitude
      CALL CHECK( NF90_GET_VAR( FID, LAT_ID, TMP4 ), 302 )
      OMPS_SO2(1:N_SO2)%LAT = TMP4(1:N_SO2)

      ! Read time
      CALL CHECK( NF90_GET_VAR( FID, TIME_ID, TMP8 ), 303 )
      OMPS_SO2(1:N_SO2)%TIME = TMP8(1:N_SO2)

      ! Read column SO2
      CALL CHECK( NF90_GET_VAR( FID, SO2_ID,  TMP4 ), 304 )
      OMPS_SO2(1:N_SO2)%SO2  = TMP4(1:N_SO2)

      ! Read SO2 error
      CALL CHECK( NF90_GET_VAR( FID, ERR_ID,  TMP4 ), 305 )
      OMPS_SO2(1:N_SO2)%ERR  = TMP4(1:N_SO2)

      ! Close the file
      CALL CHECK( NF90_CLOSE( FID ), 9999 )

      DEALLOCATE( TMPINT )
      DEALLOCATE( TMP4   )
      DEALLOCATE( TMP8   )

      ! Return to the calling routines
      END SUBROUTINE READ_OMPS_SO2_OBS
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE CALC_OMPS_SO2_FORCE( COST_FUNC )
!
!*****************************************************************************
!  Subroutine CALC_OMPS_SO2_FORCE calculate the adjoint forcing from
!  OMPS L2 SO2 observation and updates the cost function.
!  (ywang, 08/07/14)
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
      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
      USE GRID_MOD,           ONLY : GET_IJ
      USE TIME_MOD,           ONLY : GET_NYMD, GET_NHMS
      USE TIME_MOD,           ONLY : GET_TAUb, GET_TAU, GET_TAUe
      USE TRACER_MOD,         ONLY : XNUMOL
      USE TRACERID_MOD,       ONLY : IDTSO2
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP

      ! Arguements
      REAL*8, INTENT(INOUT)       :: COST_FUNC

      ! Parameters
      REAL*8, PARAMETER           :: M2DU = 2.69D20 ! 1 DU = 2.69D20 molec/m2

      ! Local variables
      INTEGER                     :: NT
      INTEGER                     :: NC
      INTEGER                     :: IIJJ(2), I, J, L

      REAL*8                      :: OLD_COST
      REAL*8                      :: TMP_COST
      REAL*8, ALLOCATABLE         :: NEW_COST(:)

      REAL*8                      :: GC_SO2_CONC(LLPAR)  ! GEOS-Chem SO2 at each layer (units: kg/m3)
      REAL*8                      :: LAYER_GC_SO2(LLPAR) ! GEOS-Chem SO2 at each layer (units: DU)
      REAL*8                      :: COLUMN_GC_SO2       ! GEOS-Chem column SO2 (units: DU)
      REAL*8                      :: DIFF
      REAL*8                      :: FORCING

      LOGICAL                     :: FIRST = .TRUE.

      !=================================================================
      ! CALC_OMPS_SO2_FORCE begins here!
      !=================================================================

      PRINT*, ' - CALC_OMPS_SO2_FORCE: OMPS SO2 forcing '

      ! Initialize
      IF ( FIRST ) THEN

         CALL INIT_OMPS_SO2_OBS

         FIRST = .FALSE.

      END IF

      ! Save a value of the cost function first
      OLD_COST = COST_FUNC
!      PRINT*, 'OLD_COST', OLD_COST, COST_FUNC

      ! Check if it is the last time of a day
      IF ( OBS_FREQ > 60 ) THEN
         PRINT*, '236000 - OBS_FREQ * 100 is not valid'
         STOP
      END IF
      IF ( GET_NHMS() == 236000 - OBS_FREQ * 100 ) THEN
!      IF ( (GET_NHMS() == 236000 - OBS_FREQ * 100) .OR.
!     &     ( ABS(GET_TAU() -GET_TAUe()) < 1e-6  ) ) THEN

         CALL READ_OMPS_SO2_OBS( GET_NYMD() )

      END IF

      ! No observations for current day
      IF ( N_SO2 == 0 ) THEN

         PRINT*, '    - CALC_OMPS_SO2_FORCE: No OMPS SO2 obsevations for
     &current day'

         IF ( ABS( GET_TAUb() - GET_TAU() ) < 1E-6) THEN

            CALL MAKE_AVERAGE_OMPS_SO2

         END IF

         RETURN

      END IF

      ! GET observations in time window
      CALL GET_OBS

      ! No observations for time window
      IF ( N_CURR == 0 ) THEN

         PRINT*, '    - CALC_OMPS_SO2_FORCE: No OMPS SO2 obsevations for
     &current time window'

         IF ( ABS( GET_TAUb() - GET_TAU() ) < 1E-6) THEN

            CALL MAKE_AVERAGE_OMPS_SO2

         END IF

         RETURN

      END IF

      ! Reset
      CURR_GC_SO2   = 0D0
      CURR_SO2      = 0D0
      CURR_DIFF_SO2 = 0D0
      CURR_FORCING  = 0D0
      CURR_COST     = 0D0
      CURR_COUNT    = 0

      ! Reset
      IF ( ALLOCATED( NEW_COST ) ) DEALLOCATE( NEW_COST )
      ALLOCATE ( NEW_COST(N_CURR) )
      NEW_COST = 0D0

      NC = 0
      ! Loop for all observations
      DO NT = 1, N_SO2, 1

         ! Observations in time window and simulation area
         IF ( FLAGS(NT) ) THEN

            ! Get grid box of current record
            IIJJ = GET_IJ( REAL(OMPS_SO2(NT)%LON ,4),
     &                     REAL(OMPS_SO2(NT)%LAT ,4) )

            I    = IIJJ(1)
            J    = IIJJ(2)


            ! SO2 outside troposphere is set to 0
            GC_SO2_CONC  = 0D0
            LAYER_GC_SO2 = 0D0
            DO L = 1, LLPAR, 1

               IF ( ITS_IN_THE_TROP(I,J,L) ) THEN

                  ! Units conversion [ kg/gridbox => DU ]
                  ! Units of LAYER_GC_SO2: DU
                  !           GC_SO2_CONC: kg/m3
                  !              CHK_STT : kg/gridbox
                  !               AIRVOL : m3/gridbox
                  !             BXHEIGHT : m
                  !               XNUMOL : molec/kg
                  ! M2DU = 2.69D20 , 1 DU = 2.69D20 molec/m2

                  GC_SO2_CONC(L)  = CHK_STT(I,J,L,IDTSO2) /
     &                              AIRVOL(I,J,L)

                  LAYER_GC_SO2(L) = GC_SO2_CONC(L)        *
     &                              BXHEIGHT(I,J,L)       *
     &                              XNUMOL(IDTSO2)        /
     &                              M2DU

               END IF

            END DO
            COLUMN_GC_SO2    = SUM( LAYER_GC_SO2 )

            ! Count observations in gridbox
            CURR_COUNT(I,J) = CURR_COUNT(I,J) + 1
            !ALL_COUNT(I,J)  = ALL_COUNT(I,J)  + 1

            CURR_GC_SO2(I,J) = CURR_GC_SO2(I,J) + COLUMN_GC_SO2
            !ALL_GC_SO2(I,J)  = ALL_GC_SO2(I,J)  + COLUMN_GC_SO2

            ! Calculate sum here. In the end, average will be calculated
            CURR_SO2(I,J) = CURR_SO2(I,J) + OMPS_SO2(NT)%SO2
            !ALL_SO2(I,J)  = ALL_SO2(I,J)  + OMPS_SO2(NT)%SO2

            !-----------------------------
            ! Calculate adjoint forcing
            !-----------------------------

            ! The difference between GEOS-Chem SO2 and OMPS SO2
            DIFF               = COLUMN_GC_SO2      - OMPS_SO2(NT)%SO2
            CURR_DIFF_SO2(I,J) = CURR_DIFF_SO2(I,J) + DIFF
            !ALL_DIFF_SO2(I,J)  = ALL_DIFF_SO2(I,J)  + DIFF

            ! S_{obs}^{-1} * DIFF
            FORCING           = DIFF / (OMPS_SO2(NT)%ERR ** 2)
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
     &                                    FORCING / AIRVOL(I,J,L) *
     &                                    BXHEIGHT(I,J,L)         *
     &                                    XNUMOL(IDTSO2) / M2DU

               END IF

            END DO

         END IF ! FLAGS(NT)

      END DO ! NT

      ! Update cost function
      COST_FUNC = COST_FUNC + SUM( NEW_COST )
      PRINT*, ' Update value of COST_FUNC = ', COST_FUNC
      PRINT*, ' OMPS SO2 contribution      = ', COST_FUNC - OLD_COST

      PRINT*, ' MIN/MAX STT_ADJ  = ', MINVAL(STT_ADJ), MAXVAL(STT_ADJ)
      PRINT*, ' MIN/MAX in       = ', MINLOC(STT_ADJ), MAXLOC(STT_ADJ)
      PRINT*, ' MIN/MAX NEW_COST = ', MINVAL(NEW_COST), MAXVAL(NEW_COST)
      PRINT*, ' MIN/MAX cost in  = ', MINLOC(NEW_COST),MAXLOC(NEW_COST)


      CALL MAKE_CURRENT_OMPS_SO2


      IF ( ABS( GET_TAUb() - GET_TAU() ) < 1E-6) THEN

         CALL MAKE_AVERAGE_OMPS_SO2

      END IF

      ! Return to the calling routines
      END SUBROUTINE CALC_OMPS_SO2_FORCE
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE MAKE_CURRENT_OMPS_SO2
!
!*****************************************************************************
!  Subroutine MAKE_CURRENT_OMPS_SO2 output some dignostic data for
!  current assimilation time window
!  (ywang, 08/07/17)
!
!*****************************************************************************
!
      ! Reference to f90 module
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,    ONLY : EXPAND_NAME
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
      USE FILE_MOD,          ONLY : IU_OMPSSO2,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TIME_MOD,          ONLY : GET_NYMD, GET_NHMS
      USE TIME_MOD,          ONLY : GET_TAU

      REAL*4, PARAMETER    :: UNDEF = -999.0

      ! Local Variables
      INTEGER              :: I,    I0, J,  J0, L
      REAL*4               :: SO2(IIPAR,JJPAR,6)
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
      ! MAKE_CURRENT_OMPS_SO2 begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_FILE = 'gctm.omps.so2.YYYYMMDD.hhmm.NN'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM diag File: ' //
     &           'OMPS SO2'
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
 100  FORMAT( '     - MAKE_CURRENT_OMPS_SO2:  Writing ', a )

      ! Open file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_OMPSSO2, FILENAME, TITLE )

      SO2 = 0.0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         IF ( CURR_COUNT(I,J) > 0 ) THEN
            SO2(I,J,1) = REAL(CURR_GC_SO2(I,J))   / CURR_COUNT(I,J)
            SO2(I,J,2) = REAL(CURR_SO2(I,J))      / CURR_COUNT(I,J)
            SO2(I,J,3) = REAL(CURR_DIFF_SO2(I,J)) / CURR_COUNT(I,J)
            SO2(I,J,4) = REAL(CURR_FORCING(I,J))
            SO2(I,J,5) = REAL(CURR_COST(I,J))
            SO2(I,J,6) = REAL(CURR_COUNT(I,J))

            ALL_GC_SO2(I,J)   = ALL_GC_SO2(I,J)   + SO2(I,J,1)
            ALL_SO2(I,J)      = ALL_SO2(I,J)      + SO2(I,J,2)
            ALL_DIFF_SO2(I,J) = ALL_DIFF_SO2(I,J) + SO2(I,J,3)
            ALL_COUNT(I,J)    = ALL_COUNT(I,J)    + 1
         END IF

      ENDDO
      ENDDO

!$OMP END PARALLEL DO

      CALL BPCH2( IU_OMPSSO2,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  26,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     6,         I0+1,
     &            J0+1,      1,         SO2 )

      ! Close file
      CLOSE( IU_OMPSSO2 )

      ! Return to the calling routines
      END SUBROUTINE MAKE_CURRENT_OMPS_SO2
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE MAKE_AVERAGE_OMPS_SO2
!
!*****************************************************************************
!  Subroutine MAKE_AVERAGE_OMPS_SO2 output some dignostic data for
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
      USE FILE_MOD,          ONLY : IU_OMPSSO2,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TIME_MOD,          ONLY : GET_NYMD
      USE TIME_MOD,          ONLY : GET_TAUb, GET_TAUe

      REAL*4, PARAMETER    :: UNDEF = -999.0

      ! Local Variables
      INTEGER              :: I,    I0, J,  J0, L
      REAL*4               :: SO2(IIPAR,JJPAR,6)
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
      ! MAKE_AVERAGE_OMPS_SO2 begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_FILE = 'gctm.omps.so2.ave.YYYYMMDD.NN'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM diag File: ' //
     &           'OMPS SO2'
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
 100  FORMAT( '     - MAKE_AVERAGE_OMPS_SO2:  Writing ', a )

      ! Open file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_OMPSSO2, FILENAME, TITLE )

      SO2 = 0.0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         IF ( ALL_COUNT(I,J) > 0 ) THEN
            SO2(I,J,1) = REAL(ALL_GC_SO2(I,J))   / ALL_COUNT(I,J)
            SO2(I,J,2) = REAL(ALL_SO2(I,J))      / ALL_COUNT(I,J)
            SO2(I,J,3) = REAL(ALL_DIFF_SO2(I,J)) / ALL_COUNT(I,J)
            SO2(I,J,4) = REAL(ALL_FORCING(I,J))
            SO2(I,J,5) = REAL(ALL_COST(I,J))
            SO2(I,J,6) = REAL(ALL_COUNT(I,J))
         END IF

      ENDDO
      ENDDO

!$OMP END PARALLEL DO

      CALL BPCH2( IU_OMPSSO2,    MODELNAME,  LONRES,     LATRES,
     &            HALFPOLAR, CENTER180,  CATEGORY,   26,
     &            UNIT,      GET_TAUb(), GET_TAUe(), RESERVED,
     &            IIPAR,     JJPAR,      6,         I0+1,
     &            J0+1,      1,          SO2 )

      ! Close file
      CLOSE( IU_OMPSSO2 )

      ! Return to the calling routines
      END SUBROUTINE MAKE_AVERAGE_OMPS_SO2
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE GET_OBS
!
!*****************************************************************************
!  Subroutine GET_OBS finds all obsevations in the current time window
!  and simulation area.
!  (ywang, 08/07/17)
!
!  Module variable as Input:
!  ===========================================================================
!  ( 1) N_SO2     (INTEGER) : Number of observation in current day
!  Arguements as Output:
!
!  Module variable as Output:
!  ===========================================================================
!  ( 1) FLAGS     (LOGICAL) : Whether or not a speicific obsevation is
!  in current time window and simulation area.
!  ( 2) N_CURR    (INTEGER) : Number of obsevations in the current time
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
      ! In GEOS-Chem, it is since 1/1/1985, while in OMPS, it is since
      ! 1/1/1993
      JD85    = GET_JD( 19850000, 000000 )
      JD93    = GET_JD( 19930000, 000000 )
      JD93_85 = JD93 - JD85
      JD93_85 = JD93_85 * 24D0 ! days => hours

!      write(*, '(A15, f30.15)' ) 'JD93_85', JD93_85

      ! Get current GEOS-Chem TAU
      CURRENT_TAU = GET_TAU()
!      write(*,'(A15, f30.15)') 'CURRENT_TAU',CURRENT_TAU

      ! Change current GEOS-Chem TAU into OMPS TAU
      CURRENT_TAU = CURRENT_TAU - JD93_85
!      write(*,'(A15, f30.15)') 'CURRENT_TAU',CURRENT_TAU

      ! Change TAU units ( hours => second )
      CURRENT_TAU = CURRENT_TAU * 3600D0
!      write(*,'(A15, f30.15)') 'CURRENT_TAU',CURRENT_TAU

      ! Get half time window
      HALF_TIME_WINDOW = TIME_WINDOW / 2D0
      HALF_TIME_WINDOW = HALF_TIME_WINDOW * 60D0 ! ( minute => second )
!      write(*, '(A15, f30.15)') 'HALF_TIME_WINDOW',HALF_TIME_WINDOW

      ! Get current time window
      WINDOW_BEGIN = CURRENT_TAU - HALF_TIME_WINDOW
      WINDOW_END   = CURRENT_TAU + HALF_TIME_WINDOW

!      write(*, '(A15,f30.15)') 'WINDOW_BEGIN',WINDOW_BEGIN
!      write(*, '(A15,f30.15)') 'WINDOW_END',WINDOW_END

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
!      Write(*, '(f30.15)') OMPS_SO2(1:10)%TIME
      N_CURR = 0
      ! Find observations in current time window and simulation area
      DO NT = 1, N_SO2, 1

         IF (      ( OMPS_SO2(NT)%TIME >= WINDOW_BEGIN )
     &       .AND. ( OMPS_SO2(NT)%TIME <  WINDOW_END   )
#if defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_SD )
     &       .AND. ( OMPS_SO2(NT)%LON  >= XEDGE_MIN    )
     &       .AND. ( OMPS_SO2(NT)%LON  <= XEDGE_MAX    )
     &       .AND. ( OMPS_SO2(NT)%LAT  >= YEDGE_MIN    )
     &       .AND. ( OMPS_SO2(NT)%LAT  <= YEDGE_MAX    )
#endif
     &                                                  ) THEN


            FLAGS(NT) = .TRUE.

            N_CURR    = N_CURR + 1

         ELSE

            FLAGS(NT) = .FALSE.

         END IF

      END DO

      WRITE(6, 100) N_CURR, GET_NHMS()
 100  FORMAT('      Number of OMPS SO2 observations: ' I10 ' at ' I10.6)

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
        CALL ERROR_STOP('netCDF error', 'omps_so2_mod')
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
      END MODULE OMPS_SO2_OBS_MOD
