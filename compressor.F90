program main

  USE GCKPP_GLOBAL
  USE GCKPP_JACOBIANSP
  USE GCKPP_PARAMETERS
  USE GCKPP_MONITOR

  IMPLICIT NONE

  INTEGER                :: ICNTRL(20), IERR, I, II, III, N
  INTEGER                :: ISTATUS(20)
  REAL(dp)               :: RCNTRL(20)
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

  ! Formatting vars
  character(len=20) :: lunz, nv, clunz, cnv

  OUTPUT       = .false.
  REINIT       = .true.  ! Reset C every NITR,NAVG iteration
!  REINIT       = .false. ! Let C evolve over the NAVG loop
  NAVG         = 100
  AR_threshold = 1e-2 ! Threshold value for AR

  Cinit = 0._dp
  cNONZERO = 0 ! Initialize number of nonzero elements in reduced mechanism

  ! Initialize
  Cinit(1:NVAR)   = 1.e8_dp
  Cinit(ind_O2)   = 5.3e18_dp
  Cinit(ind_CH4)  = 4.2e13_dp
  Cinit(ind_CO)   = 1.0e12_dp
  Cinit(ind_NO)   = 1.25e9_dp
  Cinit(ind_NO2)  = 1.25e9_dp

  write(*,*) ' '
  write(*,*) 'The KPP Auto-reduction test boxmodel'
  write(*,*) ' '
  write(*,*) ' '
  write(*,*) 'Running the full mechanism for memory initialization'
  write(*,*) '... '
  write(*,*) ' '
  ! -------------------------------------------------------------------------- !
  ! 1. Run the full mechanism
  ! initialization run. Priming memory. Otherwise, first pass is slow
  call fullmech(.true.) 

  call fullmech(.false.)
  Cfull = C
  
  ! -------------------------------------------------------------------------- !
  ! 2. Run the compacted mechanism

  call compactedmech()
  Credux = C
  
  ! -------------------------------------------------------------------------- !
  ! 3. Calculate the error norm per Santillana et al. (2010) and Shen et al. (2020)

  RRMS = sqrt(sum(((Credux(SPC_MAP(1:rNVAR))-Cfull(SPC_MAP(1:rNVAR)))/Cfull(SPC_MAP(1:rNVAR)))**2,&
       MASK=Cfull(SPC_MAP(1:rNVAR)).ne.0..and.Cfull(SPC_MAP(1:rNVAR)).gt.1e6_dp)/dble(rNVAR))

  ! -------------------------------------------------------------------------- !
  ! 5. Report timing comparison

  write(*,*) ' '
  write(*,*) ' '
  write(*,'(a,e9.1)')   '     threshold: ', AR_threshold
  write(*,'(a,f6.2,a)') '          RRMS: ', 100.*RRMS,"%"
  write(*,'(a,f6.2,a)') '  AR/full time: ', 100.*compact_avg/full_avg, "%" 
  write(*,'(a,f6.2,a)') '  problem size: ', 100.*(rNVAR)/(NVAR), "%"
  write(*,'(a,f6.2,a)') '  non-zero elm: ', 100.*(cNONZERO)/(LU_NONZERO), "%"

!  DO i=1,NSPEC
!     write(*,*) SPC_NAMES(i)
!  ENDDO

CONTAINS

  subroutine fullmech( init )
    USE GCKPP_INTEGRATOR
    USE GCKPP_RATES
    USE GCKPP_INITIALIZE

    IMPLICIT NONE

    LOGICAL :: init

    ! Set OPTIONS
    IERR      = 0                 ! Success or failure flag
    ISTATUS   = 0                 ! Rosenbrock output 
    RCNTRL    = 0.0_dp            ! Rosenbrock input
    RSTATE    = 0.0_dp            ! Rosenbrock output
    ICNTRL    = 0
    ICNTRL(1) = 1
    ICNTRL(2) = 0	
    ICNTRL(3) = 4
    ICNTRL(7) = 1
    
    ! Tolerances
    ATOL      = 1e-2_dp    
    RTOL      = 1e-2_dp

    ! Set ENV
    T    = 0d0
    TIN  = T
    TOUT = T + 600._dp
    TEMP = 298.
    
    if (.not. init) write(*,*) 'Running the full mechanism'
    
    full_avg     = 0.
    full_sumtime = 0.
    start        = 0.
    end          = 0.

    if (.not. reinit) C(1:NSPEC) = Cinit(1:NSPEC)
    ! --- INTEGRATION & TIMING LOOP
    DO I=1,NAVG
       call cpu_time(start)
       if (reinit) C(1:NSPEC) = Cinit(1:NSPEC)
       ! Initialize
       call Initialize()
       VAR(1:NVAR) = C(1:NVAR)
       ! Set RCONST
       call Update_RCONST()
       ! Integrate
       CALL Integrate( TIN,    TOUT,    ICNTRL,      &
            RCNTRL, ISTATUS, RSTATE, IERR )
       C(1:NVAR)   = VAR(:)
       call cpu_time(end)
       full_sumtime = full_sumtime+end-start
    ENDDO
    full_avg = full_sumtime/real(NAVG)
    if (.not. init) write(*,*) "Average integration time: ", full_avg
    if (.not. init) write(*,'(a,i5)') " Number of iterations: ", NAVG

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
    ICNTRL(3) = 4
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

    if (.not.reinit) C(1:NSPEC) = Cinit(1:NSPEC)
    ! --- INTEGRATION & TIMING LOOP
    DO I=1,NAVG ! Iterate to generate average comp time
       call cpu_time(start)
       if (reinit) C(1:NSPEC) = Cinit(1:NSPEC)
       ! Initialize
       call Initialize()
       VAR = C(1:NVAR)
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

end program main
