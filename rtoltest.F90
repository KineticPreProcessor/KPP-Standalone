program main


  USE GCKPP_GLOBAL
  USE GCKPP_JACOBIANSP
  USE GCKPP_PARAMETERS
  USE GCKPP_MONITOR
  USE INITIALIZE
  USE gckpp_StoichiomSP
  USE GCKPP_MODEL
  USE SETQUANTS

  IMPLICIT NONE

  INTEGER                :: Experiment, NFILES
  INTEGER                :: ICNTRL(20), IERR, I, II, III, N, filecount
  INTEGER                :: ISTATUS(20)
  INTEGER                :: fileTotSteps
  INTEGER                :: Exp2TotSteps
  INTEGER                :: StepsSavedExp2
  INTEGER                :: level
  REAL(dp)               :: RCNTRL(20)
  REAL(dp)               :: Hstart
  REAL(dp)               :: cosSZA
  REAL(dp)               :: RSTATE(20)
  REAL(dp)               :: T, TIN, TOUT, start, start2, end
  REAL :: full_sumtime, comp_sumtime, compact_avg, full_avg, setup_time

  ! COMPACTION
  INTEGER :: REMOVE(NVAR)    ! Vector of species indexes to be removed from full mechanism
  
  INTEGER :: idx, S
  INTEGER :: NSTEPSt ! Number of iterations in the timing averaging loop

  REAL(dp)             :: Vloc(NVAR), Cinit(NSPEC), Ctrue(NSPEC), Credux(NSPEC)

  ! Vars for reading batches of files
  character(len=256) :: filename
  integer :: iostat
  integer :: unit


  write(*,*) ' '
  write(*,*) 'The KPP boxmodel - testcase'
  write(*,*) ' '

   R     = 0._dp
   Cinit = 0._dp

  !
  ! Experiment 1: Run for one set of initial conditions
  !               and print end concentrations to file 

   filename = 'samples/BorneoTwilight_L23_20181001_1015.txt'
   ! Read the input
   read(unit, '(A)', iostat=iostat) filename
   call read_input(filename, R, Cinit, Hstart, cosSZA, level, fileTotSteps)
   ! Run the full mechanism
   call fullmech(.false.,1,RTOL_VALUE=1e-2_dp)

 CONTAINS

  subroutine fullmech( init , RELAX_RTOL_CASE, RTOL_VALUE)
    USE GCKPP_INTEGRATOR
    USE GCKPP_RATES
    USE GCKPP_INITIALIZE
    USE GCKPP_GLOBAL

    IMPLICIT NONE

    LOGICAL :: init
    INTEGER :: RELAX_RTOL_CASE
    REAL(dp)    :: RTOL_VALUE

    ! Set OPTIONS
    IERR      = 0                 ! Success or failure flag
    ISTATUS   = 0                 ! Rosenbrock output 
    RCNTRL    = 0.0_dp            ! Rosenbrock input
    RCNTRL(3) = Hstart
    RSTATE    = 0.0_dp            ! Rosenbrock output
    ICNTRL    = 0
    ICNTRL(1) = 1
    ICNTRL(2) = 0	
    ICNTRL(3) = 4
    ICNTRL(7) = 1
    ICNTRL(15) = -1
    
    ! Tolerances
    ATOL      = 1e-2_dp
    RTOL      = RTOL_VALUE ! default in GEOS-CF 2.0 is 0.5e-2_dp

    ! Set ENV
    T    = 0d0
    TIN  = T
    TOUT = T + 900._dp
    TEMP = 298.
    
    ! if (.not. init) write(*,*) 'Running the full mechanism'
    
    full_avg     = 0.
    full_sumtime = 0.
    start        = 0.
    end          = 0.

    ! C(1:NSPEC) = Cinit(1:NSPEC)
    ! --- INTEGRATION & TIMING LOOP
    ! Generate solution at tight tolerances
    C(1:NSPEC) = Cinit(1:NSPEC)
    call Update_RCONST()
    CALL Integrate( TIN,    TOUT,    ICNTRL,      &
         RCNTRL, ISTATUS, RSTATE, IERR )
    Ctrue = C
    NSTEPSt = ISTATUS(3)
    write(*,'(a,i5)') " Number of internal timesteps: ", ISTATUS(3)

    ! Run the integration with predicted RTOL
    call cpu_time(start)

    call Initialize()
    C(1:NSPEC) = Cinit(1:NSPEC)

    VAR(1:NVAR) => C(1:NVAR)
    FIX(1:NFIX) => C(NVAR+1:NSPEC)

    ! Set RTOL from predictor
    call setRTOLandC()

    ! Set RCONST
    call Update_RCONST()

    CALL Fun( C, FIX, RCONST, Vloc )

    ! Get a random RTOL
    CALL RANDOM_NUMBER(RTOL)
    RTOL = 10**(-2.*RTOL)

    ! Integrate
    CALL Integrate( TIN,    TOUT,    ICNTRL,      &
         RCNTRL, ISTATUS, RSTATE, IERR )

    write(*,*) 'NSTEPS: ', NSTEPSt, ISTATUS(3)
    write(*,*) C(ind_O3)
    write(*,*) Ctrue(ind_O3)
    write(*,*) 100*(C(ind_O3)-Ctrue(ind_O3))/Ctrue(ind_O3),"%"
  
    call cpu_time(end)
    full_sumtime = full_sumtime+end-start

    return
 end subroutine fullmech

end program main
