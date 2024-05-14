program main


  USE GCKPP_GLOBAL
  USE GCKPP_JACOBIANSP
  USE GCKPP_PARAMETERS
  USE GCKPP_MONITOR
  USE INITIALIZE
  USE gckpp_StoichiomSP
  USE GCKPP_MODEL


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
  REAL(dp)               :: cosSZA, OperatorTimestep
  REAL(dp)               :: RSTATE(20)
  REAL(dp)               :: T, TIN, TOUT, start, start2, end
  REAL :: full_sumtime, comp_sumtime, compact_avg, full_avg, setup_time

  ! COMPACTION
  INTEGER :: REMOVE(NVAR)    ! Vector of species indexes to be removed from full mechanism
  INTEGER :: NRMV                      ! Number of species removed (size(REMOVE))

  INTEGER :: idx, S
  INTEGER :: NRTOL, NSTEPSt ! Number of iterations in the timing averaging loop

  REAL(dp)             :: Vloc(NVAR), Cinit(NSPEC), Ctrue(NSPEC), Credux(NSPEC)

  ! Vars for reading batches of files
  character(len=256) :: filename
  integer :: iostat
  integer :: unit, testcsvunit

  write(*,*) ' '
  write(*,*) 'The KPP boxmodel'
  write(*,*) ' '

  ! Set up and populate the test data file
  open(testcsvunit,FILE='testdata.csv')
  write(testcsvunit,'(a)',advance='NO') 'filename'
  write(testcsvunit,'(a)',advance='NO') ',NSTEPS'
  do i=1,NSPEC
     write(testcsvunit,'(a)',advance='NO') ','//trim(spc_names(i))//'_i'
  end do
  do i=1,NVAR
     write(testcsvunit,'(a)',advance='NO') ',d'//trim(spc_names(i))//'_dt'
  end do
  write(testcsvunit,'(a)') ''
  open(newunit=unit, file='testfiles.txt', status='old', action='read')
  read(unit, '(I10)', iostat=iostat) NFILES
  write(*,*) 'NFILES for testing: ', NFILES
  DO II = 1,NFILES
     ! Read the input
     read(unit, '(A)', iostat=iostat) filename
     call read_input(filename, R, Cinit, SPC_NAMES, Hstart, cosSZA, level, fileTotSteps, OperatorTimestep)
     write(*,*) 'Reading '//trim(filename)//' to testdata.csv'
     write(testcsvunit,'(a,i3,a)',advance='NO') trim(filename)//',',fileTotSteps,','
     call dumptestdata()
  ENDDO

  close(unit)
  close(testcsvunit)

CONTAINS

  subroutine dumptestdata( )
    USE GCKPP_RATES
    USE GCKPP_INITIALIZE
    USE GCKPP_GLOBAL

    IMPLICIT NONE

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

    call Initialize()
    C(1:NSPEC) = Cinit(1:NSPEC)
    VAR(1:NVAR) => C(1:NVAR)
    FIX(1:NFIX) => C(NVAR+1:NSPEC)
    !call Update_RCONST()
    RCONST = R
    CALL Fun( VAR, FIX, RCONST, Vloc )
    write(testcsvunit,'(*(G0.18,:,","))') Cinit, Vloc
    VAR => null()
    FIX => null()

    return
  end subroutine dumptestdata

end program main
