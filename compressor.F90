program main


  USE GCKPP_GLOBAL
  USE GCKPP_JACOBIANSP
  USE GCKPP_PARAMETERS
  USE GCKPP_MONITOR
  USE INITIALIZE
  USE gckpp_StoichiomSP



  IMPLICIT NONE

  INTEGER                :: Experiment
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
  INTEGER :: NAVG ! Number of iterations in the timing averaging loop

  INTEGER, ALLOCATABLE :: rLU_IROW(:), rLU_ICOL(:) ! temporary, for display purposes
  LOGICAL, ALLOCATABLE :: tDO_FUN(:), tDO_SLV(:), tDO_JVS(:)

  REAL(dp)             :: dcdt(NVAR), RxR(NREACT), cinit(NSPEC), AR_threshold, Cfull(NSPEC),Credux(NSPEC), RRMS
  REAL(dp)             :: A(NREACT), Prod(NVAR), Loss(NVAR)

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
  integer :: unit

  OUTPUT       = .false.
  REINIT       = .true.  ! Reset C every NITR,NAVG iteration
!  REINIT       = .false. ! Let C evolve over the NAVG loop
  NAVG         = 1
  AR_threshold = 1e4 ! Threshold value for AR
  SCENARIO     = 1
  RELAX_RTOL_CASE = 0 ! if 0, no P/L species ! 1 tight tolerances for P/L species ! 2 20 species ATOL to concentration ! 3 higher RTOL across the board
  STOICMAT     = .false. ! Print biadjacency matrix of species reaction graph, spc and rxn list
  Experiment = 3 ! 1: Run for one initial condition 2: Get species-resolved errors for expensive twilight cases ! 3: Compare two runs with different tolerances

  write(*,*) ' '
  write(*,*) 'The KPP boxmodel'
  write(*,*) ' '
  write(*,*) ' '
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
   cNONZERO = 0 ! Initialize number of nonzero elements in reduced mechanism

  ! -------------------------------------------------------------------------- !
  ! Experiment 1: Run for one set of initial conditions
  !               and print end concentrations to file 
  IF (Experiment .eq. 1) THEN
    ! Read the input file
    open(newunit=unit, file='filelist_exp1.txt', status='old', action='read')
    read(unit, '(A)', iostat=iostat) filename
    call read_input(filename, R, Cinit, Hstart, cosSZA, level, fileTotSteps)
    ! Run the full mechanism
    call fullmech(.false.,1,RTOL_VALUE=0.5e-2_dp)
    ! Write the concentrations to file
    open(996,FILE ='final_concentrations_exp1.csv')
    DO i=1,NSPEC
        write(996,'(a,f30.6)')  trim(spc_names(i)//", "),C(i)
    ENDDO
    close(996)
  ENDIF

  ! -------------------------------------------------------------------------- !
  ! Experiment 2: Get species-resolved errors for expensive twilight cases
  IF (Experiment .eq. 2) THEN
    ! Open the list of files to be read
    open(newunit=unit, file='filelist_exp2.txt', status='old', action='read')
    ! Open the file to write the species-resolved errors
    open(998,FILE='TwilightSpeciesErrors.csv')
    write(998,'(a)',advance='NO') 'FileName,cosSZA,Level,NStp,Err,'//trim(spc_names(1))
    do i=2,NVAR
      write(998,'(a)',advance='NO') ','//trim(spc_names(i))
    end do
    write(998,'(a)') ''
    filecount = 0
    do ! Loop over the files in the list
      read(unit, '(A)', iostat=iostat) filename
      if (iostat /= 0) exit
      call read_input(trim(filename), R, Cinit, Hstart, cosSZA, level, fileTotSteps)
      ! Only focus on files with more than 10 internal timesteps
      if (fileTotSteps < 10) then
          cycle
      endif 
      filecount = filecount + 1
      write(*,*) 'Reading file: ', trim(filename)
      write(*,*) 'file ', filecount
      ! Run the full mechanism
      write(998,'(a)',advance='NO') trim(filename)//','
      write(998,'(e12.4,a,i2,a)',advance='NO') cosSZA,',',level,','
      call fullmech(.false.,RELAX_RTOL_CASE,RTOL_VALUE=0.5e-2_dp)
    enddo
    close(998)
  ENDIF

  ! -------------------------------------------------------------------------- !
  RRMS = 1.0_dp
  ! Experiment 3: Compare two runs with different tolerances
  IF (Experiment .eq. 3) THEN
    StepsSavedExp2 = 0
    Exp2TotSteps = 0
    open(newunit=unit, file='filelist_exp3.txt', status='old', action='read')
    filecount = 0
    open(998,FILE='Experiment2_Err.csv')
    open(888,FILE='Experiment2_Conc.csv')
    write(998,'(a)',advance='NO') 'FileName,ADJUST_RTOL,NStp,Err,'//trim(spc_names(1))
    do i=2,NVAR
      write(998,'(a)',advance='NO') ','//trim(spc_names(i))
    end do
    write(998,'(a)') ''
    write(888,'(a)',advance='NO') 'FileName,ADJUST_RTOL,NStp,'//trim(spc_names(1))
    do i=2,NVAR
      write(888,'(a)',advance='NO') ','//trim(spc_names(i))
    end do
    write(888,'(a)') ''
    do
      read(unit, '(A)', iostat=iostat) filename
      if (iostat /= 0) exit
      call read_input(trim(filename), R, Cinit, Hstart, cosSZA, level, fileTotSteps)
      if (fileTotSteps < 1.0) then
          cycle
      else
        filecount = filecount + 1
        write(*,*) "EXPERIMENT ", filecount
        write(*,*) "File name: ", trim(filename)
        ! Run the full mechanism
        write(*,*) 'Running the full mechanism'
        write(998,'(a)',advance='NO') trim(filename)//',0,'
        call fullmech(.false.,0,RTOL_VALUE=0.5e-2_dp) 
        Cfull = C
        ! Write the concentrations to 888 Experiment2_Conc.csv
        write(888,'(a,i0)',advance='NO') trim(filename)//',0,',ISTATUS(3)
        do i=1,NVAR
          write(888,'(a,e12.4)',advance='NO') ',',C(i)
        end do
        write(888,'(a)') ''
        StepsSavedExp2 = StepsSavedExp2 + ISTATUS(3)
        Exp2TotSteps = Exp2TotSteps + ISTATUS(3)
        if (fileTotSteps .ne. ISTATUS(3)) then
          write(*,*) "WARNING: fileTotSteps is not equal to ISTATUS(3)"
          write(*,*) "fileTotSteps: ", fileTotSteps
           write(*,*)"ISTATUS(3):   ", ISTATUS(3)
        end if
        ! Run the mechanism with relaxed tolerances
        C = Cinit
        write(*,*) 'Running with customized tolerances'
        write(998,'(a)',advance='NO') trim(filename)//',1,'
        ! call fullmech(.false.,12,RTOL_VALUE=0.5e-2_dp) 
        call fullmech(.false.,0,RTOL_VALUE=1e-1_dp) 

        Credux = C
        write(888,'(a,i0)',advance='NO') trim(filename)//',1,',ISTATUS(3)
        do i=1,NVAR
          write(888,'(a,e12.4)',advance='NO') ',',C(i)
        end do
        write(888,'(a)') ''
        StepsSavedExp2 = StepsSavedExp2 - ISTATUS(3)
        ! Calculate the error norm per Santillana et al. (2010) and Shen et al. (2020)
        RRMS = sqrt(sum(((Credux-Cfull)/Cfull)**2, MASK=Cfull.ne.0..and.Cfull.gt.1e6_dp)/dble(NSPEC))
        write(*,'(a,f6.2,a)') '          RRMS: ', 100.*RRMS,"%"
        write(*,'(a,f6.2,a)') '  NOx relative error: ', 100.*abs(Credux(ind_NO)+Credux(ind_NO2)- &
                                                                 Cfull(ind_NO)-Cfull(ind_NO2))/ &
                                                                 (Cfull(ind_NO)+Cfull(ind_NO2)), "%"
      endif
    enddo
    write(*,*) "Total steps saved in Experiment 2: ", StepsSavedExp2, " out of ", Exp2TotSteps, &
     " total steps, [", 100.*StepsSavedExp2/Exp2TotSteps, "%]"
  ENDIF


  ! call read_input("samples/IndianOcean_L1_20190702_0015.txt", R, Cinit, Hstart,fileTotSteps)

  ! potentially pass file string as an input argument

  ! IF (SCENARIO .eq. 1) &
  !      call initialize_sunrise(Cinit, R) ! Sunrise, surface, E. Indian ocean, July 1
  ! IF (SCENARIO .eq. 2) &
  !      call initialize_sunrisemidtrop(Cinit, R) ! Sunrise, surface, E. Indian ocean, July 1
  ! IF (SCENARIO .eq. 3) &
  !      call initialize_namericaday(Cinit, R) ! N. America, surface, day, July 1
  ! IF (SCENARIO .eq. 4) &
  !      call initialize_satlnight(Cinit, R) ! S. Mid Atlantic, midtrop, night, July 1

!  where (Cinit .eq. 0.d0) Cinit = 1e-20 ! Set min concentration, if needed


  ! -------------------------------------------------------------------------- !
  ! 1. Run the full mechanism
  ! initialization run. Priming memory. Otherwise, first pass is slow
  ! call fullmech(.true.,RELAX_RTOL_CASE, RTOL_VALUE=0.5e-2_dp) 
  ! call fullmech(.false.,RELAX_RTOL_CASE,  RTOL_VALUE=0.5e-2_dp) 



!  call fullmech(.false.)
  Cfull = C
  
!  call massbalance(Cfull, Cinit)

  ! -------------------------------------------------------------------------- !
  ! 2. Run the compacted mechanism

!  call compactedmech()
  ! Credux = C
  
  ! call massbalance(Credux, Cfull)

  ! -------------------------------------------------------------------------- !
  ! 3. Calculate the error norm per Santillana et al. (2010) and Shen et al. (2020)

  RRMS = sqrt(sum(((Credux(SPC_MAP(1:rNVAR))-Cfull(SPC_MAP(1:rNVAR)))/Cfull(SPC_MAP(1:rNVAR)))**2,&
       MASK=Cfull(SPC_MAP(1:rNVAR)).ne.0..and.Cfull(SPC_MAP(1:rNVAR)).gt.1e6_dp)/dble(rNVAR))

  ! -------------------------------------------------------------------------- !
  ! 5. Report timing comparison

  write(*,*) ' '
  write(*,*) ' '
  ! write(*,'(a,e9.1)')   '     threshold: ', AR_threshold
  ! write(*,'(a,f6.2,a)') '          RRMS: ', 100.*RRMS,"%"
  ! ! write(*,'(a,f6.2,a)') '  AR/full time: ', 100.*compact_avg/full_avg, "%" 
  ! write(*,'(a,f6.2,a)') '  problem size: ', 100.*(rNVAR)/(NVAR), "%"
  ! write(*,'(a,f6.2,a)') '  non-zero elm: ', 100.*(cNONZERO)/(LU_NONZERO), "%"
  
  ! -------------------------------------------------------------------------- !
  ! 6. Write concentrations
  ! open(997,FILE ='C.csv')
  ! DO i=1,NSPEC
  !    write(997,'(a,f30.6)')  trim(spc_names(i)//", "),C(i)
  ! ENDDO
  ! write(*,*) trim(spc_names(ind_NO2)),"  C(",ind_NO2,")=", C(ind_NO2)
  ! write(*,*) trim(spc_names(ind_O3)),"  C(",ind_NO2,")=", C(ind_O3)

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

    ! P/L Tolerances very high
    ATOL(ind_LBRO2H)    = 1e25_dp
    ATOL(ind_LBRO2N)    = 1e25_dp
    ATOL(ind_LCH4)      = 1e25_dp
    ATOL(ind_LCO)       = 1e25_dp
    ATOL(ind_LISOPNO3)  = 1e25_dp
    ATOL(ind_LISOPOH)   = 1e25_dp
    ATOL(ind_LNRO2H)    = 1e25_dp
    ATOL(ind_LNRO2N)    = 1e25_dp
    ATOL(ind_LOx)       = 1e25_dp
    ATOL(ind_LTRO2H)    = 1e25_dp
    ATOL(ind_LTRO2N)    = 1e25_dp
    ATOL(ind_LXRO2H)    = 1e25_dp
    ATOL(ind_LXRO2N)    = 1e25_dp
    ATOL(ind_PCO)       = 1e25_dp
    ATOL(ind_PH2O2)     = 1e25_dp
    ATOL(ind_POx)       = 1e25_dp 
    ATOL(ind_PSO4)      = 1e25_dp

    IF (RELAX_RTOL_CASE .eq. 1 ) THEN
      write(*,*) "Tight tolerance for all species also prod/loss species"
      ATOL      = 1e-2_dp

    ELSEIF (RELAX_RTOL_CASE .eq. 2) THEN
      ! Old set of stiff species, from a single grid cell in the Indian Ocean
      RTOL(ind_I2O2) = 8e-1_dp
      RTOL(ind_I2O4) = 8e-1_dp
      RTOL(ind_IO) = 8e-1_dp
      RTOL(ind_O) = 8e-1_dp
      RTOL(ind_HNO4) = 8e-1_dp
      RTOL(ind_OIO) = 8e-1_dp
      RTOL(ind_Cl) = 8e-1_dp
      RTOL(ind_ClOO) = 8e-1_dp
      RTOL(ind_I2O3) = 8e-1_dp
      RTOL(ind_HOBr) = 8e-1_dp
      RTOL(ind_HO2) = 8e-1_dp
      RTOL(ind_OH) = 8e-1_dp
      RTOL(ind_I) = 8e-1_dp
      RTOL(ind_Br) = 8e-1_dp
      RTOL(ind_CO2) = 8e-1_dp
      RTOL(ind_ClO) = 8e-1_dp
      RTOL(ind_IONO2) = 8e-1_dp
      RTOL(ind_N2O5) = 8e-1_dp
      RTOL(ind_IONO) = 8e-1_dp
      RTOL(ind_OClO) = 8e-1_dp


    ELSEIF (RELAX_RTOL_CASE .eq. 3) THEN
      RTOL(ind_BrO) = 0.8_dp
      RTOL(ind_RCHO) = 0.8_dp
      RTOL(ind_B3O2) = 0.8_dp
      RTOL(ind_ATO2) = 0.8_dp
      RTOL(ind_HOI) = 0.8_dp
      RTOL(ind_ALD2) = 0.8_dp
      RTOL(ind_PRPE) = 0.8_dp
      RTOL(ind_ClO) = 0.8_dp
      RTOL(ind_R4O2) = 0.8_dp
      RTOL(ind_HOCl) = 0.8_dp
      ! RTOL(ind_O3) = 0.8_dp
      RTOL(ind_MO2) = 0.8_dp
      RTOL(ind_H2O) = 0.8_dp
      RTOL(ind_CH2O) = 0.8_dp
      RTOL(ind_MCO3) = 0.8_dp
      RTOL(ind_MACR) = 0.8_dp
      RTOL(ind_MACRNO2) = 0.8_dp
      RTOL(ind_SALACL) = 0.8_dp
      RTOL(ind_IONO2) = 0.8_dp
      RTOL(ind_KO2) = 0.8_dp
      ! RTOL(ind_HNO3) = 0.8_dp
      RTOL(ind_HAC) = 0.8_dp

    ELSE IF (RELAX_RTOL_CASE .eq. 4) THEN
    RTOL(ind_BrO) = 0.8_dp
    RTOL(ind_RCHO) = 0.8_dp
    RTOL(ind_B3O2) = 0.8_dp
    RTOL(ind_ATO2) = 0.8_dp
    RTOL(ind_HOI) = 0.8_dp
    RTOL(ind_ALD2) = 0.8_dp
    RTOL(ind_PRPE) = 0.8_dp
    RTOL(ind_ClO) = 0.8_dp
    RTOL(ind_R4O2) = 0.8_dp
    RTOL(ind_HOCl) = 0.8_dp
    ! RTOL(ind_O3) = 0.8_dp
    RTOL(ind_MO2) = 0.8_dp
    RTOL(ind_H2O) = 0.8_dp
    RTOL(ind_CH2O) = 0.8_dp
    RTOL(ind_MCO3) = 0.8_dp
    RTOL(ind_MACR) = 0.8_dp
    RTOL(ind_MACRNO2) = 0.8_dp
    RTOL(ind_SALACL) = 0.8_dp
    RTOL(ind_IONO2) = 0.8_dp
    RTOL(ind_KO2) = 0.8_dp
    ! RTOL(ind_HNO3) = 0.8_dp
    RTOL(ind_HAC) = 0.8_dp
    RTOL(ind_BrSALC) = 0.8_dp
    RTOL(ind_NO2) = 0.8_dp
    RTOL(ind_INO2D) = 0.8_dp
    RTOL(ind_IHOO1) = 0.8_dp
    RTOL(ind_HO2) = 0.8_dp
    RTOL(ind_BrSALA) = 0.8_dp
    RTOL(ind_I) = 0.8_dp
    RTOL(ind_ClNO3) = 0.8_dp
    RTOL(ind_O1D) = 0.8_dp

    ELSEIF (RELAX_RTOL_CASE .eq. 5) THEN
      RTOL(ind_BrO) = 0.8_dp
    RTOL(ind_RCHO) = 0.8_dp
    RTOL(ind_B3O2) = 0.8_dp
    RTOL(ind_ATO2) = 0.8_dp
    RTOL(ind_HOI) = 0.8_dp
    RTOL(ind_ALD2) = 0.8_dp
    RTOL(ind_PRPE) = 0.8_dp
    RTOL(ind_ClO) = 0.8_dp
    RTOL(ind_R4O2) = 0.8_dp
    RTOL(ind_HOCl) = 0.8_dp
    ! RTOL(ind_O3) = 0.8_dp
    RTOL(ind_MO2) = 0.8_dp
    RTOL(ind_H2O) = 0.8_dp
    RTOL(ind_CH2O) = 0.8_dp
    RTOL(ind_MCO3) = 0.8_dp
    RTOL(ind_MACR) = 0.8_dp
    RTOL(ind_MACRNO2) = 0.8_dp
    RTOL(ind_SALACL) = 0.8_dp
    RTOL(ind_IONO2) = 0.8_dp
    RTOL(ind_KO2) = 0.8_dp
    ! RTOL(ind_HNO3) = 0.8_dp
    RTOL(ind_HAC) = 0.8_dp
    RTOL(ind_BrSALC) = 0.8_dp
    RTOL(ind_NO2) = 0.8_dp
    RTOL(ind_INO2D) = 0.8_dp
    RTOL(ind_IHOO1) = 0.8_dp
    RTOL(ind_HO2) = 0.8_dp
    RTOL(ind_BrSALA) = 0.8_dp
    RTOL(ind_I) = 0.8_dp
    RTOL(ind_ClNO3) = 0.8_dp
    RTOL(ind_O1D) = 0.8_dp
    RTOL(ind_R4N1) = 0.8_dp
    RTOL(ind_RCO3) = 0.8_dp
    ! RTOL(ind_NO) = 0.8_dp
    RTOL(ind_IHOO4) = 0.8_dp
    RTOL(ind_HCl) = 0.8_dp
    RTOL(ind_CO) = 0.8_dp
    RTOL(ind_HOBr) = 0.8_dp
    RTOL(ind_INO2B) = 0.8_dp
    RTOL(ind_ITHN) = 0.8_dp
    RTOL(ind_IDN) = 0.8_dp
    RTOL(ind_OH) = 0.8_dp

    ELSEIF (RELAX_RTOL_CASE .eq. 6) THEN
      ! top 20 species from all summer cases with 15 or more timesteps
      RTOL(ind_SALACL) = 0.8_dp
      RTOL(ind_OH) = 0.8_dp
      RTOL(ind_ClO) = 0.8_dp
      RTOL(ind_O) = 0.8_dp
      RTOL(ind_HO2) = 0.8_dp
      RTOL(ind_BrSALA) = 0.8_dp
      RTOL(ind_PRPE) = 0.8_dp
      RTOL(ind_NO) = 0.8_dp
      RTOL(ind_BrO) = 0.8_dp
      RTOL(ind_IDN) = 0.8_dp
      RTOL(ind_SALCCL) = 0.8_dp
      RTOL(ind_BrSALC) = 0.8_dp
      RTOL(ind_IHOO1) = 0.8_dp
      RTOL(ind_RCHO) = 0.8_dp
      RTOL(ind_NO3) = 0.8_dp
      RTOL(ind_ALD2) = 0.8_dp
      RTOL(ind_ITHN) = 0.8_dp
      RTOL(ind_INO2D) = 0.8_dp
      RTOL(ind_R4O2) = 0.8_dp
      RTOL(ind_CH2O) = 0.8_dp
      ELSEIF (RELAX_RTOL_CASE .eq. 7) THEN
        ! top 20 species from SumScaledErrorNorm of all summer cases with 10 or more timesteps
        RTOL(ind_INO) = 0.8_dp
        RTOL(ind_OIO) = 0.8_dp
        RTOL(ind_I2) = 0.8_dp
        RTOL(ind_I2O2) = 0.8_dp
        RTOL(ind_Br) = 0.8_dp
        RTOL(ind_I2O3) = 0.8_dp
        RTOL(ind_I2O4) = 0.8_dp
        RTOL(ind_OClO) = 0.8_dp
        RTOL(ind_CO2) = 0.8_dp
        RTOL(ind_I) = 0.8_dp
        RTOL(ind_OH) = 0.8_dp
        RTOL(ind_IONO) = 0.8_dp
        RTOL(ind_Cl) = 0.8_dp
        RTOL(ind_H) = 0.8_dp
        RTOL(ind_BrNO2) = 0.8_dp
        RTOL(ind_NO3) = 0.8_dp
        RTOL(ind_NO) = 0.8_dp
        RTOL(ind_IO) = 0.8_dp
        RTOL(ind_HOBr) = 0.8_dp
        RTOL(ind_HO2) = 0.8_dp
      ELSEIF (RELAX_RTOL_CASE .eq. 8) THEN
        ! top 20 species from SumScaledErrorNorm of all summer cases with 8 or more timesteps
        RTOL(ind_INO) = 0.8_dp
        RTOL(ind_I2) = 0.8_dp
        RTOL(ind_OIO) = 0.8_dp
        RTOL(ind_I2O2) = 0.8_dp
        RTOL(ind_H) = 0.8_dp
        RTOL(ind_Br) = 0.8_dp
        RTOL(ind_I2O3) = 0.8_dp
        RTOL(ind_HOBr) = 0.8_dp
        RTOL(ind_I2O4) = 0.8_dp
        RTOL(ind_IONO) = 0.8_dp
        RTOL(ind_BrNO2) = 0.8_dp
        RTOL(ind_CO2) = 0.8_dp
        RTOL(ind_OH) = 0.8_dp
        RTOL(ind_OClO) = 0.8_dp
        RTOL(ind_I) = 0.8_dp
        RTOL(ind_NO3) = 0.8_dp
        RTOL(ind_Cl) = 0.8_dp
        RTOL(ind_NO) = 0.8_dp
        RTOL(ind_Br2) = 0.8_dp
        RTOL(ind_HO2) = 0.8_dp

      ELSEIF (RELAX_RTOL_CASE .eq. 9) THEN
        ! top species from SumScaledErrorNorm of all twilight cases with 10 or more timesteps
        RTOL(ind_INO) = 0.8_dp
        RTOL(ind_I2) = 0.8_dp
        RTOL(ind_OIO) = 0.8_dp
        RTOL(ind_I2O2) = 0.8_dp
        RTOL(ind_I2O4) = 0.8_dp
        RTOL(ind_Br) = 0.8_dp
        RTOL(ind_I2O3) = 0.8_dp
        RTOL(ind_IONO) = 0.8_dp
        RTOL(ind_I) = 0.8_dp
        ! RTOL(ind_OH) = 0.8_dp
        RTOL(ind_CO2) = 0.8_dp
        RTOL(ind_OClO) = 0.8_dp
        RTOL(ind_Cl) = 0.8_dp
        RTOL(ind_BrNO2) = 0.8_dp
        RTOL(ind_H) = 0.8_dp
        RTOL(ind_NO3) = 0.8_dp
        ! RTOL(ind_NO) = 0.8_dp
        RTOL(ind_HOBr) = 0.8_dp
        RTOL(ind_IO) = 0.8_dp
        ! RTOL(ind_HO2) = 0.8_dp
        ! Another 10
        RTOL(ind_Br2) = 0.8
        RTOL(ind_BrSALA) = 0.8
        RTOL(ind_ICl) = 0.8
        RTOL(ind_MO2) = 0.8
        RTOL(ind_BrO) = 0.8
        RTOL(ind_ClOO) = 0.8
        RTOL(ind_HI) = 0.8
        RTOL(ind_Cl2O2) = 0.8
        RTOL(ind_R4N1) = 0.8
        RTOL(ind_IONO2) = 0.8
      ELSEIF (RELAX_RTOL_CASE .eq. 10) THEN
        ! top species from SumScaledErrorRanked of all twilight cases with 10 or more timesteps
        RTOL(ind_INO) = 0.8
        RTOL(ind_Br) = 0.8
        RTOL(ind_IONO) = 0.8
        RTOL(ind_I) = 0.8
        RTOL(ind_OIO) = 0.8
        ! !!RTOL(ind_OH) = 0.8
        RTOL(ind_CO2) = 0.8
        RTOL(ind_IO) = 0.8
        RTOL(ind_I2) = 0.8
        RTOL(ind_Cl) = 0.8
        RTOL(ind_NO3) = 0.8
        RTOL(ind_I2O2) = 0.8
        ! !!RTOL(ind_NO) = 0.8
        RTOL(ind_BrO) = 0.8
        RTOL(ind_OClO) = 0.8
        ! !!RTOL(ind_HO2) = 0.8
        RTOL(ind_ClO) = 0.8
        RTOL(ind_BrNO2) = 0.8
        RTOL(ind_ClOO) = 0.8
        RTOL(ind_I2O3) = 0.8
        ! Another 10
        RTOL(ind_HNO4) = 0.8
        RTOL(ind_N2O5) = 0.8
        RTOL(ind_IONO2) = 0.8
        RTOL(ind_R4N1) = 0.8
        RTOL(ind_MCO3) = 0.8
        RTOL(ind_MO2) = 0.8
        RTOL(ind_RCO3) = 0.8
        RTOL(ind_AROMRO2) = 0.8
        RTOL(ind_BENZO) = 0.8
        RTOL(ind_Br2) = 0.8
        
        ! Missing species:
        ! "I2O4"
        ! "H"
        ! "HOBr"
        ! "BrSALA"
        ! "ICl"
        ! "HI"
        ! "Cl2O2"
        ! RTOL(ind_I2O4) = 0.8
        ! RTOL(ind_H) = 0.8
        ! RTOL(ind_HOBr) = 0.8
        ! RTOL(ind_BrSALA) = 0.8
        ! RTOL(ind_ICl) = 0.8
        ! RTOL(ind_HI) = 0.8
        ! RTOL(ind_Cl2O2) = 0.8
      ELSE IF (RELAX_RTOL_CASE .eq. 11) THEN
          ! top species from SumScaledErrorNorm of all twilight cases with 10 or more timesteps
         write(*,*) "0.2 RTOL"
          RTOL(ind_INO) = 0.2_dp
          RTOL(ind_I2) = 0.2_dp
          RTOL(ind_OIO) = 0.2_dp
          RTOL(ind_I2O2) = 0.2_dp
          RTOL(ind_I2O4) = 0.2_dp
          RTOL(ind_Br) = 0.2_dp
          RTOL(ind_I2O3) = 0.2_dp
          RTOL(ind_IONO) = 0.2_dp
          RTOL(ind_I) = 0.2_dp
          ! RTOL(ind_OH) = 0.8_dp
          RTOL(ind_CO2) = 0.2_dp
          RTOL(ind_OClO) = 0.2_dp
          RTOL(ind_Cl) = 0.2_dp
          RTOL(ind_BrNO2) = 0.2_dp
          RTOL(ind_H) = 0.2_dp
          RTOL(ind_NO3) = 0.2_dp
          ! RTOL(ind_NO) = 0.8_dp
          RTOL(ind_HOBr) = 0.2_dp
          RTOL(ind_IO) = 0.2_dp
          ! RTOL(ind_HO2) = 0.8_dp
          ! Another 10
          RTOL(ind_Br2) = 0.2
          RTOL(ind_BrSALA) = 0.2
          RTOL(ind_ICl) = 0.2
          RTOL(ind_MO2) = 0.2
          RTOL(ind_BrO) = 0.2
          RTOL(ind_ClOO) = 0.2
          RTOL(ind_HI) = 0.2
          RTOL(ind_Cl2O2) = 0.2
          RTOL(ind_R4N1) = 0.2
          RTOL(ind_IONO2) = 0.2
        ELSE IF (RELAX_RTOL_CASE .eq. 12) THEN
          ! top species from SumScaledErrorNorm of all twilight cases with 10 or more timesteps
         write(*,*) "running with 0.02 RTOL"
         RTOL(ind_INO) = 0.02_dp
         RTOL(ind_I2) = 0.02_dp
         RTOL(ind_OIO) = 0.02_dp
         RTOL(ind_I2O2) = 0.02_dp
         RTOL(ind_I2O4) = 0.02_dp
         RTOL(ind_Br) = 0.02_dp
         RTOL(ind_I2O3) = 0.02_dp
         RTOL(ind_IONO) = 0.02_dp
         RTOL(ind_I) = 0.02_dp
         ! RTOL(ind_OH) = 0.02_dp
         RTOL(ind_CO2) = 0.02_dp
         RTOL(ind_OClO) = 0.02_dp
         RTOL(ind_Cl) = 0.02_dp
         RTOL(ind_BrNO2) = 0.02_dp
         RTOL(ind_H) = 0.02_dp
         RTOL(ind_NO3) = 0.02_dp
         ! RTOL(ind_NO) = 0.02_dp
         RTOL(ind_HOBr) = 0.02_dp
         RTOL(ind_IO) = 0.02_dp
         ! RTOL(ind_HO2) = 0.02_dp
         ! Another 10
         RTOL(ind_Br2) = 0.02_dp
         RTOL(ind_BrSALA) = 0.02_dp
         RTOL(ind_ICl) = 0.02_dp
         RTOL(ind_MO2) = 0.02_dp
         RTOL(ind_BrO) = 0.02_dp
         RTOL(ind_ClOO) = 0.02_dp
         RTOL(ind_HI) = 0.02_dp
         RTOL(ind_Cl2O2) = 0.02_dp
         RTOL(ind_R4N1) = 0.02_dp
         RTOL(ind_IONO2) = 0.02_dp

      END IF
    ! write(*,*) "ATOL(ind_OH) = ", ATOL(ind_OH)
      

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

    if (.not. reinit) C(1:NSPEC) = Cinit(1:NSPEC)
    ! C(1:NSPEC) = Cinit(1:NSPEC)
    ! --- INTEGRATION & TIMING LOOP
    ! open(998,FILE='data.csv')
    ! write(998,'(a)',advance='NO') 'NStp, Err, H, Fac, '//trim(spc_names(1))
    ! do i=2,NVAR
    !    write(998,'(a)',advance='NO') ','//trim(spc_names(i))
    ! end do
    ! write(998,'(a)') ''

    DO I=1,NAVG
       call cpu_time(start)
       
       call Initialize()
       if (reinit) C(1:NSPEC) = Cinit(1:NSPEC)

       VAR(1:NVAR) => C(1:NVAR)
       FIX(1:NFIX) => C(NVAR+1:NSPEC)
       ! Set RCONST
       call Update_RCONST()
       ! Integrate
       CALL Integrate( TIN,    TOUT,    ICNTRL,      &
            RCNTRL, ISTATUS, RSTATE, IERR )
      !  write(*,*) 'Line 280'
      !  C(1:NVAR)  = VAR(:)
      !  write(*,*) 'Line 282'
       call cpu_time(end)
       full_sumtime = full_sumtime+end-start
    ENDDO
!    write(*,*) 'full ISTATUS: ', ISTATUS(1:5)
    full_avg = full_sumtime/real(NAVG)
    ! if (.not. init) write(*,*) "Average integration time: ", full_avg
    ! if fileTotSteps is not equal to ISTATUS(3), print a warning

    write(*,'(a,i5)') " Number of internal timesteps: ", ISTATUS(3)

    ! close(998)

    return
 end subroutine fullmech

  subroutine compactedmech()
    USE gckpp_Integrator
    USE GCKPP_RATES
    USE GCKPP_INITIALIZE

    IMPLICIT NONE

    ! Set OPTIONS
    IERR      = 0                 ! Success or failure flag
    ISTATUS   = 0                 ! Rosenbrock output 
    RCNTRL    = 0.0_dp            ! Rosenbrock input
    RSTATE    = 0.0_dp            ! Rosenbrock output
    ICNTRL    = 0
    ICNTRL(1) = 1
    ICNTRL(2) = 0	
    ICNTRL(3) = 4                 ! Integrator: Rodas3
    ICNTRL(7) = 1
    
    ICNTRL(8) = 1
    RCNTRL(8) = AR_threshold !default is 1.d2

    ! Tolerances
    ATOL      = 1e-2_dp    
    RTOL      = 1e-2_dp
    
    ! Set ENV
    T    = 0d0
    TIN  = T
    TOUT = T + 600._dp
    TEMP = 298.
    
    write(*,*) ' '
    write(*,*) 'Running the reduced mechanism'
    
    compact_avg  = 0.
    comp_sumtime = 0.
    start        = 0.
    end          = 0.

    keepActive                    = .true.
    keepSpcActive(ind_ClNO2)      = .true.
    keepSpcActive(ind_ClOO)       = .true.
    keepSpcActive(ind_BrCl)       = .true.
    keepSpcActive(ind_Br2 )       = .true.
    keepSpcActive(ind_BrNO3)      = .true.
    keepSpcActive(ind_HOBr)       = .true.
    keepSpcActive(ind_HOCl)       = .true.
    keepSpcActive(ind_ClNO3)      = .true.
    keepSpcActive(ind_Cl  )       = .true.
    keepSpcActive(ind_HBr )       = .true.
    keepSpcActive(ind_ClO )       = .true.
    keepSpcActive(ind_HCl )       = .true.
    keepSpcActive(ind_I2O2)       = .true.
    keepSpcActive(ind_BrNO2)      = .true.
    keepSpcActive(ind_Cl2O2)      = .true.
    keepSpcActive(ind_IONO)       = .true.
    keepSpcActive(ind_OClO)       = .true.
    keepSpcActive(ind_HOI)        = .true.
    keepSpcActive(ind_IONO2)      = .true.
    keepSpcActive(ind_Cl2)        = .true.
    keepSpcActive(ind_I)          = .true.
    keepSpcActive(ind_IO)         = .true.
    keepSpcActive(ind_BrO)        = .true.
    keepSpcActive(ind_Br)         = .true.

    if (.not.reinit) C(1:NSPEC) = Cinit(1:NSPEC)
    ! --- INTEGRATION & TIMING LOOP
    DO I=1,NAVG ! Iterate to generate average comp time
       call cpu_time(start)
       if (reinit) C(1:NSPEC) = Cinit(1:NSPEC)
       ! Initialize
       call Initialize()
       VAR = C(1:NVAR)
       FIX = C(NVAR+1:NSPEC)
       ! Set RCONST
       call Update_RCONST()
       ! Integrate
       CALL Integrate( TIN,    TOUT,    ICNTRL,      &
            RCNTRL, ISTATUS, RSTATE, IERR )
       C(1:NVAR)       = VAR(:)
       call cpu_time(end)
       comp_sumtime = comp_sumtime+end-start
    ENDDO
    compact_avg = comp_sumtime/real(NAVG)
!    write(*,*) 'AR ISTATUS: ', ISTATUS(1:5)
    write(*,*) "Average integration time: ", compact_avg
    write(*,'(a,i5)') " Number of iterations: ", NAVG
    
    return
  end subroutine compactedmech


  subroutine showoutput()
    IMPLICIT none
  ! Initialize formatting values
  write(lunz,*) LU_NONZERO
  write(nv,*) NVAR+1
  write(clunz,*) cNONZERO
  write(cnv,*) rNVAR+1

  write(*,*) ' '
  write(*,*) 'Species Info:'
  write(*,*) '---------------'
!  do i=1,nvar
!     write(*,'(a,i3)') " Species "//trim(spc_names(i))//" has index: ", i
!  enddo
  write(*,*) '---------------'
  write(*,*) ' '
  write(*,*) 'Full-mech sparse data: '
  write(*,'(a,'//lunz//'i4)') ' LU_IROW:  ', LU_IROW
  write(*,'(a,'//lunz//'i4)') ' LU_ICOL:  ', LU_ICOL
  write(*,'(a,'//nv//'i4)') ' LU_CROW:  ', LU_CROW
  write(*,'(a,'//nv//'i4)') ' LU_DIAG:  ', LU_DIAG
  write(*,*) '---------------'
  write(*,*) ' '
!  write(*,*) 'Reduced-mech, uncompacted sparse data: '
!  write(*,*) '-- removes species ' // SPC_NAMES(REMOVE(:))! // ' with index ', REMOVE(:)
!  write(*,'(a,'//lunz//'i4)') ' rLU_IROW: ', rLU_IROW
!  write(*,'(a,'//lunz//'i4)') ' rLU_ICOL: ', rLU_ICOL
  write(*,*) '---------------'
  write(*,*) ' '
  write(*,*) 'Compacted sparse data: '
  write(*,'(a,'//clunz//'i4)') ' cLU_IROW: ', cLU_IROW
  write(*,'(a,'//clunz//'i4)') ' cLU_ICOL: ', cLU_ICOL
  write(*,'(a,'//cnv//'i4)') ' cLU_CROW: ', cLU_CROW
  write(*,'(a,'//cnv//'i4)') ' cLU_DIAG: ', cLU_DIAG
  write(*,*) '---------------'
  write(*,*) ' '
  write(*,*) 'JVS_MAP ensures that the right JVS values are indexed in the integration'
  write(*,'(a,'//clunz//'i4)') ' JVS_MAP:  ', JVS_MAP
  write(*,*) ' '
  write(*,*) 'SPC_MAP ensures that the right species values are indexed in the integration'
  write(*,*) '-- ', rNVAR, SPC_MAP(rNVAR)
  write(*,'(a,'//nv//'i4)') ' SPC_MAP:  ', SPC_MAP
  write(*,*) ' '
  write(*,*) 'DO_SLV controls the terms that will be computed in KppSolve(): 1=compute, 0=skip'
  write(*,'(a,'//nv//'l4)') ' DO_SLV:   ', DO_SLV
  write(*,*) ' '
  write(*,*) 'DO_FUN controls the terms that will be computed in Fun(): 1=compute, 0=skip'
  write(*,'(a,'//nv//'l4)') ' DO_FUN:   ', DO_FUN
  write(*,*) ' '
  write(*,*) 'DO_JVS controls the terms that will be computed in Jac_SP(): 1=compute, 0=skip'
  write(*,'(a,'//lunz//'l4)') ' DO_JVS:   ', DO_JVS

  end subroutine showoutput

  subroutine massbalance(C_f, C_o)

    USE GCKPP_Monitor
    USE GCKPP_Parameters

    IMPLICIT NONE

    real(dp) :: C_o(NSPEC), C_f(NSPEC)
    real(dp) :: Cl_o_total, Cl_f_total
    real(dp) :: Br_o_total, Br_f_total
    real(dp) :: Ox_o_total, Ox_f_total
    real(dp) :: N_o_total, N_f_total
    real(dp) :: NOx_o_total, NOx_f_total

    write(*,*) '---  Mass Balances --- '
    ! Check Cl mass balance
    Cl_o_total = C_o(ind_CH2ICl)+3.*C_o(ind_CH3CCl3)+4.*C_o(ind_CCl4)+     &
         2.*C_o(ind_CH2Cl2)+2.*C_o(ind_Cl2O2)+C_o(ind_ICl)+C_o(ind_CH3Cl)+ &
         C_o(ind_ClOO)+C_o(ind_OClO)+C_o(ind_BrCl)+2.*C_o(ind_Cl2)+        &
         C_o(ind_ClNO2)+C_o(ind_ClNO3)+C_o(ind_HOCl)+C_o(ind_HCl)+C_o(ind_Cl)+C_o(ind_ClO)
    
    Cl_f_total = C_f(ind_CH2ICl)+3.*C_f(ind_CH3CCl3)+4.*C_f(ind_CCl4)+     &
         2.*C_f(ind_CH2Cl2)+2.*C_f(ind_Cl2O2)+C_f(ind_ICl)+C_f(ind_CH3Cl)+ &
         C_f(ind_ClOO)+C_f(ind_OClO)+C_f(ind_BrCl)+2.*C_f(ind_Cl2)+        &
         C_f(ind_ClNO2)+C_f(ind_ClNO3)+C_f(ind_HOCl)+C_f(ind_HCl)+C_f(ind_Cl)+C_f(ind_ClO)

    write(*,*) 'Cl mass balance: ', Cl_f_total - Cl_o_total, 100.*(Cl_f_total - Cl_o_total)/Cl_o_total,'%'
    ! Check Br mass balance
    Br_o_total = C_o(ind_CH2IBr)+C_o(ind_BrNO2)+2.*C_o(ind_CH2Br2)+C_o(ind_IBr)+  &
         3.*C_o(ind_CHBr3)+C_o(ind_CH3Br)+2.*C_o(ind_Br2)+C_o(ind_BrCl)+          &
         C_o(ind_BrNO3)+C_o(ind_HOBr)+C_o(ind_BrO)+C_o(ind_BrSALA)+ &
         C_o(ind_BrSALC)+C_o(ind_Br)+C_o(ind_HBr)

    Br_f_total = C_f(ind_CH2IBr)+C_f(ind_BrNO2)+2.*C_f(ind_CH2Br2)+C_f(ind_IBr)+  &
         3.*C_f(ind_CHBr3)+C_f(ind_CH3Br)+2.*C_f(ind_Br2)+C_f(ind_BrCl)+          &
         C_f(ind_BrNO3)+C_f(ind_HOBr)+C_f(ind_BrO)+C_f(ind_BrSALA)+ &
         C_f(ind_BrSALC)+C_f(ind_Br)+C_f(ind_HBr)

    write(*,*) 'Br mass balance: ', Br_f_total - Br_o_total, 100.*(Br_f_total - Br_o_total)/Br_o_total,'%'
    ! Check N  mass balance

    ! Check Ox mass balance
    Ox_o_total = C_o(ind_O3) + C_o(ind_NO2) + 2.*C_o(ind_NO3) + C_o(ind_PAN) +   &
         C_o(ind_PPN) + C_o(ind_MPAN) + C_o(ind_HNO4) + 3.*C_o(ind_N2O5) +        &
         C_o(ind_HNO3) + C_o(ind_BrO) + C_o(ind_HOBr) + C_o(ind_BrNO2) +          &
         2.*C_o(ind_BrNO3) + C_o(ind_MPN) + C_o(ind_ETHLN) + C_o(ind_MVKN) +      &
         C_o(ind_MCRHN) + C_o(ind_MCRHNB) + C_o(ind_PROPNN) + C_o(ind_R4N2) +     &
         C_o(ind_PRN1) + C_o(ind_PRPN) + C_o(ind_R4N1) + C_o(ind_HONIT) +         &
         C_o(ind_MONITS) + C_o(ind_MONITU) + C_o(ind_OLND) + C_o(ind_OLNN) +      &
         C_o(ind_IHN1) + C_o(ind_IHN2) + C_o(ind_IHN3) + C_o(ind_IHN4) +          &
         C_o(ind_INPB) + C_o(ind_INPD) + C_o(ind_ICN) + 2.*C_o(ind_IDN) +         &
         C_o(ind_ITCN) + C_o(ind_ITHN) + C_o(ind_ISOPNOO1) + C_o(ind_ISOPNOO2) +  &
         C_o(ind_INO2B) + C_o(ind_INO2D) + C_o(ind_INA) + C_o(ind_IDHNBOO) +      &
         C_o(ind_IDHNDOO1) + C_o(ind_IDHNDOO2) + C_o(ind_IHPNBOO) +               &
         C_o(ind_IHPNDOO) + C_o(ind_ICNOO) + 2.*C_o(ind_IDNOO) + C_o(ind_MACRNO2) + &
         C_o(ind_ClO) + C_o(ind_HOCl) + C_o(ind_ClNO2) + 2.*C_o(ind_ClNO3) +      &
         2.*C_o(ind_Cl2O2) + 2.*C_o(ind_OClO) + C_o(ind_O) + C_o(ind_O1D) +       &
         C_o(ind_IO) + C_o(ind_HOI)+ C_o(ind_IONO) + 2.*C_o(ind_IONO2) +          &
         2.*C_o(ind_OIO) + 2.*C_o(ind_I2O2) + 3.*C_o(ind_I2O3) + 4.*C_o(ind_I2O4)

    Ox_f_total = C_f(ind_O3) + C_f(ind_NO2) + 2.*C_f(ind_NO3) + C_f(ind_PAN) +   &
         C_f(ind_PPN) + C_f(ind_MPAN) + C_f(ind_HNO4) + 3.*C_f(ind_N2O5) +        &
         C_f(ind_HNO3) + C_f(ind_BrO) + C_f(ind_HOBr) + C_f(ind_BrNO2) +          &
         2.*C_f(ind_BrNO3) + C_f(ind_MPN) + C_f(ind_ETHLN) + C_f(ind_MVKN) +      &
         C_f(ind_MCRHN) + C_f(ind_MCRHNB) + C_f(ind_PROPNN) + C_f(ind_R4N2) +     &
         C_f(ind_PRN1) + C_f(ind_PRPN) + C_f(ind_R4N1) + C_f(ind_HONIT) +         &
         C_f(ind_MONITS) + C_f(ind_MONITU) + C_f(ind_OLND) + C_f(ind_OLNN) +      &
         C_f(ind_IHN1) + C_f(ind_IHN2) + C_f(ind_IHN3) + C_f(ind_IHN4) +          &
         C_f(ind_INPB) + C_f(ind_INPD) + C_f(ind_ICN) + 2.*C_f(ind_IDN) +         &
         C_f(ind_ITCN) + C_f(ind_ITHN) + C_f(ind_ISOPNOO1) + C_f(ind_ISOPNOO2) +  &
         C_f(ind_INO2B) + C_f(ind_INO2D) + C_f(ind_INA) + C_f(ind_IDHNBOO) +      &
         C_f(ind_IDHNDOO1) + C_f(ind_IDHNDOO2) + C_f(ind_IHPNBOO) +               &
         C_f(ind_IHPNDOO) + C_f(ind_ICNOO) + 2.*C_f(ind_IDNOO) + C_f(ind_MACRNO2) + &
         C_f(ind_ClO) + C_f(ind_HOCl) + C_f(ind_ClNO2) + 2.*C_f(ind_ClNO3) +      &
         2.*C_f(ind_Cl2O2) + 2.*C_f(ind_OClO) + C_f(ind_O) + C_f(ind_O1D) +       &
         C_f(ind_IO) + C_f(ind_HOI)+ C_f(ind_IONO) + 2.*C_f(ind_IONO2) +          &
         2.*C_f(ind_OIO) + 2.*C_f(ind_I2O2) + 3.*C_f(ind_I2O3) + 4.*C_f(ind_I2O4)

    write(*,*) 'Ox mass balance: ', Ox_f_total - Ox_o_total, 100.*(Ox_f_total - Ox_o_total)/Ox_o_total,'%'

    N_o_total = C_o(ind_HNO2) + C_o(ind_HNO3) + C_o(ind_HNO4) + C_o(ind_MPAN) +    &
         C_o(ind_NIT) + C_o(ind_NO) + C_o(ind_NO2) + C_o(ind_NO3) +                &
         2.*C_o(ind_N2O5) + C_o(ind_MPN) + C_o(ind_PAN) + C_o(ind_PPN) +           &
         2.*C_o(ind_N2O) + C_o(ind_MENO3) + C_o(ind_ETNO3) + C_o(ind_IPRNO3) +     &
         C_o(ind_NPRNO3) + C_o(ind_AONITA) + C_o(ind_ETHN) + C_o(ind_BZPAN) + C_o(ind_NPHEN)

    N_f_total = C_f(ind_HNO2) + C_f(ind_HNO3) + C_f(ind_HNO4) + C_f(ind_MPAN) +    &
         C_f(ind_NIT) + C_f(ind_NO) + C_f(ind_NO2) + C_f(ind_NO3) +                &
         2.*C_f(ind_N2O5) + C_f(ind_MPN) + C_f(ind_PAN) + C_f(ind_PPN) +           &
         2.*C_f(ind_N2O) + C_f(ind_MENO3) + C_f(ind_ETNO3) + C_f(ind_IPRNO3) +     &
         C_f(ind_NPRNO3) + C_f(ind_AONITA) + C_f(ind_ETHN) + C_f(ind_BZPAN) + C_f(ind_NPHEN)

    write(*,*) 'N mass balance: ', N_f_total - N_o_total, 100.*(N_f_total - N_o_total)/N_o_total,'%'

    NOx_o_total = C_o(ind_NO) + C_o(ind_NO2) + C_o(ind_NO3)
    NOx_f_total = C_f(ind_NO) + C_f(ind_NO2) + C_f(ind_NO3)

    write(*,*) 'NOx mass balance: ', NOx_f_total - NOx_o_total, 100.*(NOx_f_total - NOx_o_total)/NOx_o_total,'%'

  end subroutine massbalance

end program main
