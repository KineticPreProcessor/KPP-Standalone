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
  REAL(dp)               :: cosSZA
  REAL(dp)               :: RSTATE(20)
  REAL(dp)               :: T, TIN, TOUT, start, start2, end
  REAL :: full_sumtime, comp_sumtime, compact_avg, full_avg, setup_time

  ! COMPACTION
  INTEGER :: REMOVE(NVAR)    ! Vector of species indexes to be removed from full mechanism
  INTEGER :: NRMV                      ! Number of species removed (size(REMOVE))

  INTEGER :: idx, S
  INTEGER :: NRTOL, NSTEPSt ! Number of iterations in the timing averaging loop

  REAL(dp)             :: Vloc(NVAR), Cinit(NSPEC), Ctrue(NSPEC), Credux(NSPEC)

  LOGICAL              :: OUTPUT
  LOGICAL              :: ReInit
  ! Logical for getting the stoichiometry
  LOGICAL              :: STOICMAT

  ! Variable ATOL?
  INTEGER              :: RELAX_RTOL_CASE

  INTEGER              :: SCENARIO

  ! Formatting vars
  character(len=20) :: lunz, nv, clunz, cnv

  ! Vars for reading batches of files
  character(len=256) :: filename
  integer :: iostat
  integer :: unit, testcsvunit

  OUTPUT       = .false.
  REINIT       = .true.  ! Reset C every NITR,NRTOL iteration
  !  REINIT       = .false. ! Let C evolve over the NRTOL loop
  NRTOL         = 25
  RELAX_RTOL_CASE = 0 ! if 0, no P/L species ! 1 tight tolerances for P/L species ! 2 20 species ATOL to concentration ! 3 higher RTOL across the board
  STOICMAT     = .false. ! Print biadjacency matrix of species reaction graph, spc and rxn list

  write(*,*) ' '
  write(*,*) 'The KPP boxmodel'
  write(*,*) ' '
  write(*,*) '-- this program generates training data for various'
  write(*,*) '   integrator parameters and operations -- MSL, 2024'
  write(*,*) ' '

  ! Write biadjacency matrix of species reaction graph, spc and rxn list
  IF (STOICMAT) THEN
     open(990,FILE='BiadjacencyMatrix.csv')
     DO i = 1,size(STOICM)
        WRITE(990,*) IROW_STOICM(i), ",", ICOL_STOICM(i), ",", STOICM(i)
     END DO
     close(990)
     open(991,FILE='SPC_NAMES.txt')
     DO i = 1,size(SPC_NAMES)
        WRITE(991,*) SPC_NAMES(i)
     END DO
     close(991)
     open(992,FILE='EQN_NAMES.txt')
     DO i = 1,size(EQN_NAMES)
        WRITE(992,*) EQN_NAMES(i)
     END DO

     close(992)
  END IF

  R     = 0._dp
  Cinit = 0._dp

  ! Set up the training data file
  open(998,FILE='trainingdata.csv')
  write(998,'(a)',advance='NO') 'IERR'
  do i=1,NSPEC
     write(998,'(a)',advance='NO') ','//trim(spc_names(i))//'_i'
  end do
  do i=1,NSPEC
     write(998,'(a)',advance='NO') ','//trim(spc_names(i))//'_f'
  end do
  do i=1,NSPEC
     write(998,'(a)',advance='NO') ','//trim(spc_names(i))//'_t'
  end do
  write(998,'(a)',advance='NO') ',e_O3,e_NO,e_NO2'
  !   do i=1,NSPEC
  !      write(998,'(a)',advance='NO') ',e_'//trim(spc_names(i))
  !   end do
  do i=1,NVAR
     write(998,'(a)',advance='NO') ',d'//trim(spc_names(i))//'_dt'
  end do
  do i=1,NVAR
     write(998,'(a)',advance='NO') ',RTOL(ind_'//trim(spc_names(i))//')'
  end do
  write(998,'(a)',advance='NO') ',NSTEPSt'
  write(998,'(a)',advance='NO') ',NSTEPS'
  write(998,'(a)',advance='NO') ',Hend'
  write(998,'(a)',advance='NO') ',filename'
  write(998,'(a)') ''

  ! -------------------------------------------------------------------------- !
  ! Experiment 1: Run for one set of initial conditions
  !               and print end concentrations to file 
  ! Read the input file
  open(newunit=unit, file='filelist_twilight.txt', status='old', action='read')
  read(unit, '(I10)', iostat=iostat) NFILES
  write(*,*) 'NFILES: ', NFILES
  DO II = 1,NFILES
     ! Read the input
     read(unit, '(A)', iostat=iostat) filename
     call read_input(filename, R, Cinit, Hstart, cosSZA, level, fileTotSteps)

     ! Run the full mechanism
     call fullmech(.false.,1,RTOL_VALUE=1e-2_dp)
  ENDDO

  close(998)
  close(unit)

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

    ! Run the RTOL variation loop
    DO I=1,NRTOL
       call cpu_time(start)

       call Initialize()
       C(1:NSPEC) = Cinit(1:NSPEC)

       VAR(1:NVAR) => C(1:NVAR)
       FIX(1:NFIX) => C(NVAR+1:NSPEC)
       ! Set RCONST
       call Update_RCONST()

       CALL Fun( C, FIX, RCONST, Vloc )

       ! Get a random RTOL
       CALL RANDOM_NUMBER(RTOL)
       RTOL = 10**(-2.*RTOL)

       ! Integrate
       CALL Integrate( TIN,    TOUT,    ICNTRL,      &
            RCNTRL, ISTATUS, RSTATE, IERR )

       ! Diff against the "true" result

       ! Dump result to file
       write(998,'(*(G0.18,:,","))',advance='NO') IERR, Cinit, C, Ctrue, &
            abs((C(ind_O3)-Ctrue(ind_O3))/Ctrue(ind_O3)),abs((C(ind_NO)-Ctrue(ind_NO))/Ctrue(ind_NO)), &
            abs((C(ind_NO2)-Ctrue(ind_NO2))/Ctrue(ind_NO2)), Vloc, RTOL, NSTEPSt, ISTATUS(3), RSTATE(3)
       write(998,'(a)') ','//trim(filename)
       write(*,'(a,i5)') " Number of internal timesteps: ", ISTATUS(3)

       call cpu_time(end)
       full_sumtime = full_sumtime+end-start
    ENDDO


    return
  end subroutine fullmech

end program main
