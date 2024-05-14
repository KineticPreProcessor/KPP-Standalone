program main


  USE GCKPP_GLOBAL
  USE GCKPP_JACOBIANSP
  USE GCKPP_PARAMETERS
  USE GCKPP_MONITOR
  USE INITIALIZE
  USE gckpp_StoichiomSP
  USE GCKPP_MODEL
!  USE SETQUANTS
  USE PREDICTIONS

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
  REAL(dp)               :: OperatorTimestep
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

  open(newunit=unit, file='testfiles.txt', status='old', action='read')
  read(unit, '(I10)', iostat=iostat) NFILES
  write(*,*) 'NFILES for testing: ', NFILES
  call set_rtol_preds()
  DO II = 1,NFILES
     ! Read the input
     read(unit, '(A)', iostat=iostat) filename
     write(*,*) 'Test file: ', trim(filename)!, ' ', trim(testfiles(II))
     call read_input(filename, R, Cinit, SPC_NAMES, Hstart, cosSZA, level, fileTotSteps, OperatorTimestep)
     call fullmech(II)
  ENDDO

  close(unit)

 CONTAINS

  subroutine fullmech( index )
    USE GCKPP_INTEGRATOR
    USE GCKPP_RATES
    USE GCKPP_INITIALIZE
    USE GCKPP_GLOBAL

    IMPLICIT NONE

    INTEGER :: index

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
    RTOL      = 1e-2_dp

    ! Set ENV
    T    = 0d0
    TIN  = T
    TOUT = T + OperatorTimestep
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
    RCONST = R
    CALL Integrate( TIN,    TOUT,    ICNTRL,      &
         RCNTRL, ISTATUS, RSTATE, IERR )
    Ctrue = C
    NSTEPSt = ISTATUS(3)
!    write(*,'(a,i5)') " Number of internal timesteps: ", ISTATUS(3)

    ! Run the integration with predicted RTOL
    call cpu_time(start)

    call Initialize()
    C(1:NSPEC) = Cinit(1:NSPEC)

    VAR(1:NVAR) => C(1:NVAR)
    FIX(1:NFIX) => C(NVAR+1:NSPEC)

    ! Set RCONST
!    call Update_RCONST()
    RCONST = R
    
    ! Set the RTOL vals
    RTOL = RTOL_p(II,:)

    ! Protect against negative tolerances. Not good!
    ! How to prevent this?
    where(RTOL <= 0.) RTOL = 0.5e-2
    where(RTOL >= 1.) RTOL = 0.999

    RCNTRL    = 0.0_dp            ! Rosenbrock input
    RCNTRL(3) = Hstart
    RSTATE    = 0.0_dp            ! Rosenbrock output

    ! Integrate
    CALL Integrate( TIN,    TOUT,    ICNTRL,      &
         RCNTRL, ISTATUS, RSTATE, IERR )

    write(*,*) 'NSTEPS: full', NSTEPSt, ' pred',ISTATUS(3), ' ', trim(filename)
    write(*,*) 100*(C(ind_O3)-Ctrue(ind_O3))/Ctrue(ind_O3),"%",&
         100*(C(ind_NO)-Ctrue(ind_NO))/Ctrue(ind_NO),"%",&
         100*(C(ind_NO2)-Ctrue(ind_NO2))/Ctrue(ind_NO2),"%"
    write(*,*) C(ind_O3), C(ind_NO), C(ind_NO2)
    write(*,*) Ctrue(ind_O3), Ctrue(ind_NO), Ctrue(ind_NO2)
  
    call cpu_time(end)
    full_sumtime = full_sumtime+end-start

    return
 end subroutine fullmech

end program main
