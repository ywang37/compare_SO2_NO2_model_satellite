MODULE omi_no2_obs_mod
!
!******************************************************************************
! Module OMI_NO2_OBS_MOD contains all subroutine to 
! (1) Read observations
! (2) Compute cost function and adjoint forcing
!
! (ywang, 11/07/2017)
! The code is modified from Martin's code
!
! Module Varibales:
! =============================================================================
!
! Module Routines:
! =============================================================================
! (1) 
!
! =============================================================================
! NOTES:
!
!******************************************************************************
!

 IMPLICIT NONE

 ! Header files
#include "define.h"
#include "CMN_SIZE"       ! size parameters

 ! Make everything private ...
 PRIVATE

 ! ... except these routines
 PUBLIC :: calc_omi_no2_force

 !======================================================================
 ! Module variables
 !======================================================================

 ! Parameters
 INTEGER, PARAMETER :: time_window = 60 ! (units: min)
 INTEGER, PARAMETER :: NPres = 35      ! # of OMI pressure levels
 LOGICAL, PARAMETER :: super_obs = .TRUE.

 ! Record to store OMI NO2 observations
 TYPE omi_no2_obs
    REAL*8 :: time
    REAL*8 :: lat
    REAL*8 :: lon
    REAL*8 :: sza
    REAL*8 :: vza
    REAL*8 :: raa
    REAL*8 :: sw(NPres)    ! scattering weight
    REAL*8 :: swp(NPres)   ! scattering weight pressure
    REAL*8 :: omi_no2_trop !
    REAL*8 :: omi_amf_trop !
    REAL*8 :: err_s        ! slant column error [molec/cm2]
 END TYPE omi_no2_obs

 TYPE(omi_no2_obs), ALLOCATABLE :: omi_no2(:)

 INTEGER                        :: n_no2
 INTEGER                        :: n_curr

 LOGICAL,           ALLOCATABLE :: flags(:)

 REAL*8                         :: curr_gc_v_no2(iipar, jjpar)
 REAL*8                         :: curr_gc_s_no2(iipar, jjpar)
 REAL*8                         :: curr_gc_amf(iipar, jjpar)
 REAL*8                         :: curr_omi_v_no2(iipar, jjpar)
 REAL*8                         :: curr_omi_s_no2(iipar, jjpar)
 REAL*8                         :: curr_v_no2(iipar, jjpar)
 REAL*8                         :: curr_omi_amf(iipar, jjpar)
 INTEGER                        :: curr_count(iipar, jjpar)
 REAL*8                         :: curr_diff_no2(iipar, jjpar)
 REAL*8                         :: curr_forcing(iipar, jjpar)
 REAL*8                         :: curr_cost(iipar, jjpar)

 REAL*8                         :: all_gc_v_no2(iipar, jjpar)  = 0D0
 REAL*8                         :: all_gc_s_no2(iipar, jjpar)  = 0D0
 REAL*8                         :: all_gc_amf(iipar, jjpar)    = 0D0
 REAL*8                         :: all_omi_v_no2(iipar, jjpar) = 0D0
 REAL*8                         :: all_omi_s_no2(iipar, jjpar) = 0D0
 REAL*8                         :: all_v_no2(iipar, jjpar)     = 0D0
 REAL*8                         :: all_omi_amf(iipar, jjpar)   = 0D0
 INTEGER                        :: all_count(iipar, jjpar)     = 0
 REAL*8                         :: all_diff_no2(iipar, jjpar)  = 0D0
 REAL*8                         :: all_forcing(iipar, jjpar)   = 0D0
 REAL*8                         :: all_cost(iipar, jjpar)      = 0D0

 !======================================================================
 ! Module routines -- follow below the "CONTAIN" statement
 !======================================================================
 CONTAINS

!
!-------------------------------------------------------------------------------
!
 SUBROUTINE init_omi_no2_obs
!
!*******************************************************************************
! Subroutine init_omi_no2_obs initialze the OMI NO2 observation
! (1) Check if NO2 observation is specified by the adjoint input files
! (ywang, 11/07/2017)
!*******************************************************************************
!
 ! Reference to F90 modules
 USE tracer_mod,          ONLY : tracer_name
 USE tracerid_mod,        ONLY : idtno2
 USE adj_arrays_mod,      ONLY : obs_this_tracer
 USE error_mod,           ONLY : error_stop

 !======================================================================
 ! init_omi_no2_obs begins here
 !======================================================================

 IF ( obs_this_tracer(idtno2) ) THEN

    WRITE(6, 100) idtno2, tracer_name(idtno2)

 ELSE

    CALL error_stop('Error: Obs NO2 tracer is not specified in OBSERVATION MENU in input.gcadj', &
                    'init_omi_no2_obs (omi_no2_obs_mod.f)')

 ENDIF

 100 FORMAT(3X, 'Tracer ID: ', I4, 'Use OMI ', A6)

 END SUBROUTINE init_omi_no2_obs
!
!-------------------------------------------------------------------------------
!
 SUBROUTINE read_omi_no2_obs( yyyymmdd )
!
!*******************************************************************************
!
!*******************************************************************************
!
 USE netcdf
 USE time_mod,            ONLY : expand_date, get_nymd

 ! arguement
 INTEGER, INTENT(IN)   :: YYYYMMDD

 ! local variables
 INTEGER               :: fid,             n_id,            npres_id
 INTEGER               :: lon_id,          lat_id,          time_id
 INTEGER               :: sza_id,          vza_id,          raa_id
 INTEGER               :: sw_id,           swp_id
 INTEGER               :: omi_no2_trop_id, omi_amf_trop_id, err_s_id
 CHARACTER(LEN=255)    :: read_filename
 CHARACTER(LEN=4)      :: tmp
 INTEGER, ALLOCATABLE  :: tmpint(:)
 REAL*4,  ALLOCATABLE  :: tmp4(:), tmp4_(:)
 REAL*4,  ALLOCATABLE  :: tmp8(:)
 REAL*4,  ALLOCATABLE  :: tmp4_2d(:,:)
 LOGICAL               :: lf
 INTEGER               :: i

 !----------------------------------------------------------------------
 ! read_omi_no2_obs begins here
 !----------------------------------------------------------------------
 
 read_filename = 'OMI_GC_adj_NO2_YYYYMMDD.nc'

 ! expand date tokens in filename
 CALL expand_date( read_filename, YYYYMMDD, 9999 )
 
 ! construct complete filename
 read_filename = "./data/OMI_NO2/" // TRIM( read_filename )

 ! Does data file exist? If not, it means no data in the day.
 INQUIRE( FILE = TRIM( read_filename), EXIST = lf )
 IF ( .NOT. lf ) THEN

    ! no data
    n_no2 = 0

    WRITE(6, 120) get_nymd()

    RETURN

 ENDIF

 120 FORMAT(' - read_omi_no2_obs: No data file (warning) in ', I10)

 ! print to screen
 WRITE(6, 100) TRIM( read_filename )
 100 FORMAT(' - read_omi_so2_obs: reading file: ', A)

 ! open file and assign file id
 CALL check( nf90_open( read_filename, nf90_nowrite, fid ), 0 )

 !------------------------------------
 ! get data record IDs
 !------------------------------------
 CALL check( nf90_inq_dimid( fid, "time",  n_id    ), 100 )
 CALL check( nf90_inq_dimid( fid, "npres", npres_id), 101 )

 CALL check( nf90_inq_varid( fid, "lon",   lon_id          ), 103 )
 CALL check( nf90_inq_varid( fid, "lat",   lat_id          ), 104 )
 CALL check( nf90_inq_varid( fid, "TAI93", time_id         ), 105 )
! CALL check( nf90_inq_varid( fid, "sza",          sza_id          ), 106 )
! CALL check( nf90_inq_varid( fid, "vza",          vza_id          ), 107 )
! CALL check( nf90_inq_varid( fid, "raa",          raa_id          ), 108 )
 CALL check( nf90_inq_varid( fid, "SW",    sw_id           ), 109 )
 CALL check( nf90_inq_varid( fid, "SWP",   swp_id          ), 110 )
 CALL check( nf90_inq_varid( fid, "NO2",   omi_no2_trop_id ), 111 )
 CALL check( nf90_inq_varid( fid, "amf",   omi_amf_trop_id ), 112 )
 CALL check( nf90_inq_varid( fid, "ERR",   err_s_id        ), 113 )

 !-------------------------------------
 ! read dimensions
 !-------------------------------------

 ! read number of observation, n_no2
 CALL check( nf90_inquire_dimension( fid, n_id, tmp, n_no2 ), 200 ) 

 ! print to screen
 WRITE(6, 110) n_no2, get_nymd()
 110 FORMAT( '      Number of OMI NO2 observation: ', I10, ' in ', I10 )

 !----------------------------
 ! read 1D data
 !----------------------------
 ALLOCATE( tmpint(n_no2) )
 ALLOCATE( tmp4(n_no2)   )
 ALLOCATE( tmp4_(n_no2)  )
 ALLOCATE( tmp8(n_no2)   )
 ALLOCATE( tmp4_2d(NPres, n_no2) )
 tmpint = 0
 tmp4   = 0E0
 tmp4_  = 0E0
 tmp8   = 0D0
 tmp4_2d = 0E0

 ! allocate OMI NO2 observation array
 IF ( ALLOCATED( omi_no2 ) ) DEALLOCATE( omi_no2 )
 ALLOCATE( omi_no2(n_no2) )

 IF ( ALLOCATED( flags ) ) DEALLOCATE( flags )
 ALLOCATE( flags(n_no2) )

 ! read longitude
 CALL check( nf90_get_var( fid, lon_id, tmp4 ), 301 )
 omi_no2(1:n_no2)%lon = tmp4(1:n_no2)
 
 ! read latitude
 CALL check( nf90_get_var( fid, lat_id, tmp4 ), 302 )
 omi_no2(1:n_no2)%lat = tmp4(1:n_no2)

 ! read time
 CALL check( nf90_get_var( fid, time_id, tmp8 ), 303 )
 omi_no2(1:n_no2)%time = tmp8(1:n_no2)

! ! read sza
! CALL check( nf90_get_var( fid, sza_id, tmp4 ), 304 )
! omi_no2(1:n_no2)%sza = tmp4(1:n_no2)

! ! read vza
! CALL check( nf90_get_var( fid, vza_id, tmp4 ), 305 )
! omi_no2(1:n_no2)%vza = tmp4(1:n_no2)

! ! read raa
! CALL check( nf90_get_var( fid, raa_id, tmp4 ), 306 )
! omi_no2(1:n_no2)%raa = tmp4(1:n_no2)

 ! read sw
 CALL check( nf90_get_var( fid, sw_id, tmp4_2d ), 307 )
 DO i = 1, n_no2, 1
    omi_no2(i)%sw(:) = tmp4_2d(:,i)
 ENDDO

 ! read swp
 CALL check( nf90_get_var( fid, swp_id, tmp4_), 308 )
 DO i = 1, n_no2, 1
    omi_no2(i)%swp(:) = tmp4_(:)
 ENDDO

 ! read omi_no2_trop
 CALL check( nf90_get_var( fid, omi_no2_trop_id, tmp4 ), 309 )
 omi_no2(1:n_no2)%omi_no2_trop = tmp4(1:n_no2)

 ! read omi_amf_trop
 CALL check( nf90_get_var( fid, omi_amf_trop_id, tmp4 ), 310 )
 omi_no2(1:n_no2)%omi_amf_trop = tmp4(1:n_no2)

 ! read err_s
 CALL check( nf90_get_var( fid, err_s_id, tmp4 ), 311 )
 omi_no2(1:n_no2)%err_s = tmp4(1:n_no2)

 ! close the file
 CALL check( nf90_close( fid ), 9999 )

 DEALLOCATE( tmpint  )
 DEALLOCATE( tmp4    )
 DEALLOCATE( tmp4_   )
 DEALLOCATE( tmp8    )
 DEALLOCATE( tmp4_2d )

 END SUBROUTINE read_omi_no2_obs
!
!-------------------------------------------------------------------------------
!
 SUBROUTINE calc_omi_no2_force( cost_func )
!
!*******************************************************************************
!
!*******************************************************************************
!
 USE adj_arrays_mod,      ONLY : obs_freq
 USE comode_mod,          ONLY : jlop
 USE comode_mod,          ONLY : cspec_after_chem
 USE comode_mod,          ONLY : cspec_after_chem_adj
 USE grid_mod,            ONLY : get_ij
 USE grid_mod,            ONLY : get_xmid, get_ymid
 USE dao_mod,             ONLY : bxheight
 USE dao_mod,             ONLY : t, airden
 USE tropopause_mod,      ONLY : its_in_the_trop
 USE adj_arrays_mod,      ONLY : id2c
 USE tracerid_mod,        ONLY : idno2
 USE pressure_mod,        ONLY : get_pcenter, get_pedge
 USE tracer_mod,          ONLY : xnumolair
 USE time_mod,            ONLY : get_nymd,    get_nhms
 USE time_mod,            ONLY : get_taub,    get_tau,   get_taue

 ! arguements
 REAL*8, INTENT(INOUT)    :: cost_func

 ! local variables
 INTEGER                  :: nt
 INTEGER                  :: nc
 INTEGER                  :: iijj(2), i, j, l
 INTEGER                  :: k
 INTEGER                  :: jloop

 REAL*8                   :: amf_gc
 REAL*8                   :: gc_no2(llpar)
 REAL*8                   :: gc_no2_col
 REAL*8                   :: sw(NPres)    ! scattering weight
 REAL*8                   :: swp(NPres)   ! scattering weight pressure
 REAL*8                   :: sw_gc(llpar), dp(llpar)
 REAL*8                   :: diff
 REAL*8                   :: forcing
 REAL*8                   :: one_costa
 REAL*8                   :: obs_err
 REAL*8                   :: omi_s_no2
 REAL*8                   :: v_no2

 REAL*8                   :: one_cost
 REAL*8                   :: old_cost


 REAL*8                   :: sobs_adj_force(iipar,jjpar,llpar)

 LOGICAL                  :: first = .TRUE.

 !======================================================================
 ! calc_omi_no2_force begins here.
 !====================================================================== 

 ! Initialize
 IF ( first ) THEN

    CALL init_omi_no2_obs

 ENDIF

 ! save a vaule of the cost function first
 old_cost = cost_func

 ! Check if it is the last time of a day
 IF ( obs_freq > 60 ) THEN
     PRINT*, '23600 - obs_freq * 100 is not valid'
     STOP
 ENDIF
 IF ( get_nhms() == 23600 - obs_freq * 100 ) THEN

    ! CALL read_omi_no2_obs( get_nymd() )

 ENDIF

 ! No observation for current day
 IF ( n_no2 == 0 ) THEN

    PRINT*, '   - calc_omi_no2_force: no omi no2 for current day'

    IF ( ABS( get_taub() - get_taub() ) < 1E6 ) THEN

       CALL make_average_omi_no2

    ENDIF

    RETURN

 ENDIF

 ! get observations in time window
 call get_obs

 ! no observations for time window
 IF ( n_curr == 0 ) THEN

    PRINT*, '   - calc_omi_no2_force: no omi no2 for current time window'

    IF ( ABS( get_taub() - get_tau() ) < 1E-6 ) THEN

       CALL make_average_omi_no2

    ENDIF

    RETURN

 ENDIF


  curr_gc_v_no2  = 0D0
  curr_gc_s_no2  = 0D0
  curr_gc_amf    = 0D0
  curr_omi_v_no2 = 0D0
  curr_omi_s_no2 = 0D0
  curr_v_no2     = 0D0
  curr_omi_amf   = 0D0
  curr_count     = 0
  curr_diff_no2  = 0D0
  curr_forcing   = 0D0
  curr_cost      = 0D0

  sobs_adj_force        = 0D0

  nc = 0
  ! loop for all observations
  DO nt = 1, n_no2, 1

     ! get model grid coordinate indices that correspond to the observation

     iijj = get_ij( REAL(omi_no2(nt)%lon, 4), REAL(omi_no2(nt)%lat, 4) )

     i    = iijj(1)
     j    = iijj(2)

     ! initialize variables & arrays
     gc_no2     = 0D0
     gc_no2_col = 0D0
     sw_gc      = 0D0
     dp         = 0D0

     ! Get GEOS-Chem NO2 values [#/cm3] 

     DO l = 1, llpar

        IF ( its_in_the_trop(i,j,l) ) THEN

           jloop = jlop(i,j,l)
           gc_no2(l) = cspec_after_chem(jloop, id2c(idno2))

        ENDIF

     ENDDO

     ! compute tropospheric no2 vertical column [#/cm2]
     gc_no2_col = SUM( gc_no2(:) * bxheight(i,j,:) * 100D0 )
     curr_gc_v_no2(i,j) = curr_gc_v_no2(i,j) + gc_no2_col

     ! interpolate scattering weighs to GEOS-Chem grid to compute GEOS-Chem air
     ! mass factors
     sw(:)  = omi_no2(nt)%sw(:)
     swp(:) = omi_no2(nt)%swp(:)
     DO l = 1, llpar, 1

        DO k = 2, NPres, 1

           IF ( (get_pcenter(i,j,l) <= swp(k-1)) .AND. &
                (get_pcenter(i,j,l) >  swp(k)  )  ) THEN

              sw_gc(l) = sw(k) + (sw(k-1) - sw(k)) * (get_pcenter(i,j,l)-swp(k)) / (swp(k-1)-swp(k))

              ! save pressure difference of edge pressures
              dp(l) = get_pedge(i,j,l) - get_pedge(i,j,l+1)

              ! apply temperature correction, as in Bucsela 2013, eq. (4)
              sw_gc(l) = sw_gc(l) * ( 1 - 0.003 * (t(i,j,l) - 220) )

              ! convert no2 concentrations from number density to mixing ratio,
              ! as required for the calculation of air mass factor from
              ! scattering weights
              gc_no2(l) = gc_no2(l) * 1D6 / (airden(l,i,j) * xnumolair)

              EXIT

           ENDIF

        ENDDO

     ENDDO

     ! use GEOS-Chem tropospheric air mass factor to convert vertical column to
     ! slant column
     amf_gc = SUM(gc_no2 * dp * sw_gc) / SUM (gc_no2 * dp)
     gc_no2_col = amf_gc * gc_no2_col
     omi_s_no2  = omi_no2(nt)%omi_no2_trop * omi_no2(nt)%omi_amf_trop
     v_no2      = omi_s_no2 / amf_gc
     curr_gc_s_no2(i,j)  = curr_gc_s_no2(i,j)  + gc_no2_col
     curr_gc_amf(i,j)    = curr_gc_amf(i,j)    + amf_gc
     curr_omi_amf(i,j)   = curr_omi_amf(i,j)   + omi_no2(nt)%omi_amf_trop
     curr_omi_v_no2(i,j) = curr_omi_v_no2(i,j) + omi_no2(nt)%omi_no2_trop
     curr_omi_s_no2(i,j) = curr_omi_s_no2(i,j) + omi_s_no2
     curr_v_no2(i,j)     = curr_v_no2(i,j)     + v_no2
     curr_count(i,j)     = curr_count(i,j)     + 1
     
     

     ! compute slant column difference
     diff = gc_no2_col - omi_s_no2
     curr_diff_no2(i,j) = curr_diff_no2(i,j) + diff

     obs_err = omi_no2(nt)%err_s

     forcing = diff / (obs_err**2)
     curr_forcing(i,j) = curr_forcing(i,j) + forcing
     all_forcing(i,j)  = all_forcing(i,j)  + forcing

     ! update adjoint no2 concentration
     DO l = 1, llpar, 1

         IF (its_in_the_trop(i,j,l)) THEN

            IF (super_obs) THEN

               sobs_adj_force(i,j,l) = sobs_adj_force(i,j,l) + &
                         forcing * bxheight(i,j,l) * 100D0 * amf_gc
                                       
            ELSE

               jloop = jlop(i,j,l)
               cspec_after_chem_adj(jloop,id2c(idno2)) = cspec_after_chem_adj(jloop,id2c(idno2)) + &
                         forcing * bxheight(i,j,l) * 100D0 * amf_gc

            ENDIF

         ENDIF

     ENDDO 

     ! cost function
     one_cost = 0.5 * ((diff / obs_err) ** 2)
     curr_cost(i,j) = curr_cost(i,j) + one_cost
     all_cost(i,j)  = all_cost(i,j)  + one_cost

     IF (.NOT. super_obs) THEN

        cost_func = cost_func + one_cost

     ENDIF

  ENDDO

  ! update adjoint forcing and cost function for super obs
  IF (super_obs) THEN

     DO j = 1, jjpar, 1
     DO i = 1, iipar, 1

        IF ( curr_count(i,j) > 0 ) THEN

           DO l = 1, llpar, 1

               IF ( its_in_the_trop(i,j,l) ) THEN

                  jloop = jlop(i,j,l)
                  cspec_after_chem_adj(jloop,id2c(idno2)) = cspec_after_chem_adj(jloop,id2c(idno2)) + &
                             sobs_adj_force(i,j,l) / curr_count(i,j)

               ENDIF

           ENDDO

           cost_func = cost_func + curr_cost(i,j) / curr_count(i,j) 

        ENDIF

     ENDDO
     ENDDO

  ENDIF

  PRINT*, ' Update valude of cost_func = ', cost_func
  PRINT*, ' OMI NO2 contribution = ', cost_func - old_cost

 
  CALL make_current_omi_no2

  IF ( ABS(get_taub() - get_tau()) < 1E-6 ) THEN

     CALL make_average_omi_no2

  ENDIF
 
 END SUBROUTINE calc_omi_no2_force
!
!-------------------------------------------------------------------------------
!
 SUBROUTINE make_current_omi_no2
!
!*******************************************************************************
!
! (ywang, 11/08/2017)
!*******************************************************************************
!
 USE adj_arrays_mod,      ONLY : n_calc
 USE adj_arrays_mod,      ONLY : expand_name
 USE directory_adj_mod,   ONLY : diagadj_dir
 USE netcdf
 USE time_mod,            ONLY : expand_date
 USE time_mod,            ONLY : get_nymd,    get_nhms

 ! Local variables
 INTEGER                     :: fid
 INTEGER                     :: lon_dim_id, lat_dim_id
 CHARACTER(LEN=255)          :: output_file
 CHARACTER(LEN=255)          :: filename

 REAL*4                      :: tmp(iipar, jjpar, 11)

 INTEGER                     :: i, j

 !======================================================================
 ! make_current_omi_no2 begins here
 !======================================================================
  
 output_file = 'gctm.omi.no2.YYYYMMDD.hhmm.NN'

 filename = TRIM( output_file )

 ! replace YYYYMMDD and hhmm token w/ actual value
 CALL expand_date( filename, get_nymd(), get_nhms() )

 ! replace NN token w/ actual value
 CALL expand_name( filename, n_calc )

 filename = TRIM( diagadj_dir ) // trim( filename ) // '.nc'

 WRITE(6, 100) TRIM( filename )
 100 FORMAT( '     - make_current_omi_no2:  Writing ', a )

 ! Open file
 CALL check( nf90_create( filename, nf90_clobber, fid), 0 )

 ! define dimensions
 CALL check( nf90_def_dim( fid, "lon", iipar, lon_dim_id ), 100 )
 CALL check( nf90_def_dim( fid, "lat", jjpar, lat_dim_id ), 101 )

 ! end define mode
 CALL check( nf90_enddef( fid ), 8888 )

 !-------------------------------
 ! put variables
 !-------------------------------

 tmp = 0.0
! !$OMP PARALLEL DO
! !$OMP+DEFAULT( SHARED )
! !$OMP+PRIVATE( I, J )
 DO j = 1, jjpar
 DO i = 1, iipar

     IF ( curr_count(i,j) > 0 ) THEN

        tmp(i,j, 1) = REAL(curr_gc_v_no2(i,j),  4) / REAL(curr_count(i,j), 4)
        tmp(i,j, 2) = REAL(curr_gc_s_no2(i,j),  4) / REAL(curr_count(i,j), 4)
        tmp(i,j, 3) = REAL(curr_gc_amf(i,j),    4) / REAL(curr_count(i,j), 4)
        tmp(i,j, 4) = REAL(curr_omi_v_no2(i,j), 4) / REAL(curr_count(i,j), 4)
        tmp(i,j, 5) = REAL(curr_omi_s_no2(i,j), 4) / REAL(curr_count(i,j), 4)
        tmp(i,j, 6) = REAL(curr_v_no2(i,j),     4) / REAL(curr_count(i,j), 4)
        tmp(i,j, 7) = REAL(curr_omi_amf(i,j),   4) / REAL(curr_count(i,j), 4)
        tmp(i,j, 8) = REAL(curr_count(i,j),     4)
        tmp(i,j, 9) = REAL(curr_diff_no2(i,j),  4) / REAL(curr_count(i,j), 4)
        tmp(i,j,10) = REAL(curr_forcing(i,j),   4)
        tmp(i,j,11) = REAL(curr_cost(i,j),      4)

        all_gc_v_no2(i,j)  = all_gc_v_no2(i,j)  + tmp(i,j, 1)
        all_gc_s_no2(i,j)  = all_gc_s_no2(i,j)  + tmp(i,j, 2)
        all_gc_amf(i,j)    = all_gc_amf(i,j)    + tmp(i,j, 3)
        all_omi_v_no2(i,j) = all_omi_v_no2(i,j) + tmp(i,j, 4)
        all_omi_s_no2(i,j) = all_omi_s_no2(i,j) + tmp(i,j, 5)
        all_v_no2(i,j)     = all_v_no2(i,j)     + tmp(i,j, 6)
        all_omi_amf(i,j)   = all_omi_amf(i,j)   + tmp(i,j, 7)
        all_count(i,j)     = all_count(i,j)     + 1
        all_diff_no2(i,j)  = all_diff_no2(i,j)  + tmp(i,j, 9)

     ENDIF

 ENDDO
 ENDDO
! !$OMP END PARALLEL DO


 ! curr_gc_v_no2
 CALL write_nc_2d_float( fid, "gc_v_no2", "GEOS-Chem vertical column NO2", &
         "molec/cm2", tmp(:,:,1), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! curr_gc_s_no2
 CALL write_nc_2d_float( fid, "gc_s_no2", "GEOS-Chem slant column NO2", &
         "molec/cm2", tmp(:,:,2), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! curr_gc_amf
 CALL write_nc_2d_float( fid, "gc_amf", "GEOS-Chem air mass factor", &
         "unitless", tmp(:,:,3), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! curr_omi_v_no2
 CALL write_nc_2d_float( fid, "omi_v_no2", "OMI vertical column NO2", &
         "molec/cm2", tmp(:,:,4), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! curr_omi_s_no2
 CALL write_nc_2d_float( fid, "omi_s_no2", "OMI slant column NO2", &
         "molec/cm2", tmp(:,:,5), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! curr_v_no2
 CALL write_nc_2d_float( fid, "v_no2", "OMI vertical column NO2 (corrected gc_amf)", &
         "molec/cm2", tmp(:,:,6), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! curr_omi_amf
 CALL write_nc_2d_float( fid, "omi_amf", "OMI air mass factor", &
         "molec/cm2", tmp(:,:,7), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! curr_count
 CALL write_nc_2d_float( fid, "count", "# of observations in a grid box", &
         "unitless", tmp(:,:,8), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! curr_diff
 CALL write_nc_2d_float( fid, "diff", "GEOS-Chem - OMI slant column NO2", &
         "molec/cm2", tmp(:,:,9), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! curr_forcing
 CALL write_nc_2d_float( fid, "forcing", "diff / (err ** 2)", &
         "1/(molec/cm2)", tmp(:,:,10), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! curr_cost
 CALL write_nc_2d_float( fid, "cost", "cost function", &
         "unitless", tmp(:,:,11), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! close file
 CALL CHECK( nf90_close( fid ), 9999 )

 END SUBROUTINE make_current_omi_no2
!
!-------------------------------------------------------------------------------
!
 SUBROUTINE make_average_omi_no2
!
!*******************************************************************************
!
! (ywang, 11/08/2017)
!*******************************************************************************
!
 USE adj_arrays_mod,      ONLY : n_calc
 USE adj_arrays_mod,      ONLY : expand_name
 USE directory_adj_mod,   ONLY : diagadj_dir
 USE netcdf
 USE time_mod,            ONLY : expand_date
 USE time_mod,            ONLY : get_nymd

 ! Local variables
 INTEGER                     :: fid
 INTEGER                     :: lon_dim_id, lat_dim_id
 CHARACTER(LEN=255)          :: output_file
 CHARACTER(LEN=255)          :: filename

 REAL*4                      :: tmp(iipar, jjpar, 11)

 INTEGER                     :: i, j

 !======================================================================
 ! make_average_omi_no2 begins here
 !======================================================================

 output_file = 'gctm.omi.no2.YYYYMMDD.NN'

 filename = TRIM( output_file )

 ! replace YYYYMMDD token w/ actual value
 CALL expand_date( filename, get_nymd(), 9999 )

 ! replace NN token w/ actual value
 CALL expand_name( filename, n_calc )

 filename = TRIM( diagadj_dir ) // trim( filename ) // '.nc'

 WRITE(6, 100) TRIM( filename )
 100 FORMAT( '     - make_average_omi_no2:  Writing ', a )

 ! Open file
 CALL check( nf90_create( filename, nf90_clobber, fid), 0 )

 ! define dimensions
 CALL check( nf90_def_dim( fid, "lon", iipar, lon_dim_id ), 100 )
 CALL check( nf90_def_dim( fid, "lat", jjpar, lat_dim_id ), 101 )

 ! end define mode
 CALL check( nf90_enddef( fid ), 8888 )

 !-------------------------------
 ! put variables
 !-------------------------------

 tmp = 0.0
! !$OMP PARALLEL DO
! !$OMP+DEFAULT( SHARED )
! !$OMP+PRIVATE( I, J )
 DO j = 1, jjpar
 DO i = 1, iipar

     IF ( all_count(i,j) > 0 ) THEN

        tmp(i,j, 1) = REAL(all_gc_v_no2(i,j),  4) / REAL(all_count(i,j), 4)
        tmp(i,j, 2) = REAL(all_gc_s_no2(i,j),  4) / REAL(all_count(i,j), 4)
        tmp(i,j, 3) = REAL(all_gc_amf(i,j),    4) / REAL(all_count(i,j), 4)
        tmp(i,j, 4) = REAL(all_omi_v_no2(i,j), 4) / REAL(all_count(i,j), 4)
        tmp(i,j, 5) = REAL(all_omi_s_no2(i,j), 4) / REAL(all_count(i,j), 4)
        tmp(i,j, 6) = REAL(all_v_no2(i,j),     4) / REAL(all_count(i,j), 4)
        tmp(i,j, 7) = REAL(all_omi_amf(i,j),   4) / REAL(all_count(i,j), 4)
        tmp(i,j, 8) = REAL(all_count(i,j),     4)
        tmp(i,j, 9) = REAL(all_diff_no2(i,j),  4) / REAL(all_count(i,j), 4)
        tmp(i,j,10) = REAL(all_forcing(i,j),   4)
        tmp(i,j,11) = REAL(all_cost(i,j),      4)

     ENDIF

 ENDDO
 ENDDO
! !$OMP END PARALLEL DO

 ! all_gc_v_no2
 CALL write_nc_2d_float( fid, "gc_v_no2", "GEOS-Chem vertical column NO2", &
         "molec/cm2", tmp(:,:,1), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! all_gc_s_no2
 CALL write_nc_2d_float( fid, "gc_s_no2", "GEOS-Chem slant column NO2", &
         "molec/cm2", tmp(:,:,2), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! all_gc_amf
 CALL write_nc_2d_float( fid, "gc_amf", "GEOS-Chem air mass factor", &
         "unitless", tmp(:,:,3), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! all_omi_v_no2
 CALL write_nc_2d_float( fid, "omi_v_no2", "OMI vertical column NO2", &
         "molec/cm2", tmp(:,:,4), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! all_omi_s_no2
 CALL write_nc_2d_float( fid, "omi_s_no2", "OMI slant column NO2", &
         "molec/cm2", tmp(:,:,5), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! all_v_no2
 CALL write_nc_2d_float( fid, "v_no2", "OMI vertical column NO2 (corrected gc_amf)", &
         "molec/cm2", tmp(:,:,6), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! all_omi_amf
 CALL write_nc_2d_float( fid, "omi_amf", "OMI air mass factor", &
         "molec/cm2", tmp(:,:,7), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! all_count
 CALL write_nc_2d_float( fid, "count", "# of observations in a grid box", &
         "unitless", tmp(:,:,8), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! all_diff
 CALL write_nc_2d_float( fid, "diff", "GEOS-Chem - OMI slant column NO2", &
         "molec/cm2", tmp(:,:,9), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! all_forcing
 CALL write_nc_2d_float( fid, "forcing", "diff / (err ** 2)", &
         "1/(molec/cm2)", tmp(:,:,10), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! all_cost
 CALL write_nc_2d_float( fid, "cost", "cost function", &
         "unitless", tmp(:,:,11), iipar, jjpar, lon_dim_id, lat_dim_id )

 ! close file
 CALL CHECK( nf90_close( fid ), 9999 )


 END SUBROUTINE make_average_omi_no2
!
!-------------------------------------------------------------------------------
!
 SUBROUTINE write_nc_2d_float( fid, names, longnames, units, var, iin, jjn, dim1, dim2 )

 USE NETCDF

 INTEGER,              INTENT(IN) :: iin, jjn
 INTEGER,              INTENT(IN) :: fid
 INTEGER,              INTENT(IN) :: dim1, dim2
 CHARACTER(*),         INTENT(IN) :: names
 CHARACTER(*),         INTENT(IN) :: longnames
 CHARACTER(*),         INTENT(IN) :: units
 REAL*4,               INTENT(IN) :: var(iin, jjn)

 INTEGER                          :: vid

 ! open def
 CALL check( nf90_redef(fid), 9000 )

 ! define variable
 CALL check( nf90_def_var(fid, TRIM(names), nf90_float, (/dim1,dim2/), vid), 9001 )

 ! define attributes
 CALL check( nf90_put_att(fid, vid, "longname", TRIM(longnames)), 9002 )
 CALL check( nf90_put_att(fid, vid, "unit",     TRIM(units)    ), 9003 )

 ! close def
 CALL check( nf90_enddef(fid), 9004 )

 ! put var
 CALL check( nf90_put_var(fid, vid, var), 9005)

 END SUBROUTINE
!
!-------------------------------------------------------------------------------
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
!  ( 1) N_NO2     (INTEGER) : Number of observation in current day
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
      DO NT = 1, N_NO2, 1

         IF (      ( OMI_NO2(NT)%TIME >= WINDOW_BEGIN ) &
             .AND. ( OMI_NO2(NT)%TIME <  WINDOW_END   ) &
#if defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_SD )
             .AND. ( OMI_NO2(NT)%LON  >= XEDGE_MIN    ) &
             .AND. ( OMI_NO2(NT)%LON  <= XEDGE_MAX    ) &
             .AND. ( OMI_NO2(NT)%LAT  >= YEDGE_MIN    ) &
             .AND. ( OMI_NO2(NT)%LAT  <= YEDGE_MAX    ) &
#endif 
                                                        ) THEN


            FLAGS(NT) = .TRUE.

            N_CURR    = N_CURR + 1

         ELSE

            FLAGS(NT) = .FALSE.

         END IF

      END DO
      WRITE(6, 100) N_CURR, GET_NHMS()
 100  FORMAT('      Number of OMI NO2 observations: ' I10 ' at ' I10.6)

      ! Return to calling program
      END SUBROUTINE GET_OBS
!
!-------------------------------------------------------------------------------
!
      SUBROUTINE CHECK( STATUS, LOCATION )
!
!******************************************************************************
!  Subroutine CHECK checks the status of calls to netCDF libraries routines
!  (dkh, 02/15/09) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) STATUS    (INTEGER) : Completion status of netCDF library call    
!  (2 ) LOCATION  (INTEGER) : Location at which netCDF library call was made   
!     
!  NOTES:
!
!******************************************************************************
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
        CALL ERROR_STOP('netCDF error', 'modis_aod_mod')
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHECK
!
!------------------------------------------------------------------------------
!
END MODULE omi_no2_obs_mod
