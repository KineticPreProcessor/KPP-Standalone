program main

  ! GIT Branch fullchem_online
  ! 1) make the vectors non-allocatable by setting the bounding index
  !    -- this will allow not having to repeatedly allocate & deallocate
  !       vectors. We already know none of them will be longer than
  !       NVAR or LU_NONZERO

  USE GCKPP_GLOBAL
  USE GCKPP_JACOBIANSP
  USE GCKPP_PARAMETERS
  USE GCKPP_MONITOR
  USE SETQUANTS

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
  INTEGER :: NITR ! Number of integration iterations

  INTEGER, ALLOCATABLE :: rLU_IROW(:), rLU_ICOL(:) ! temporary, for display purposes
  LOGICAL, ALLOCATABLE :: tDO_FUN(:), tDO_SLV(:), tDO_JVS(:)

  REAL(dp)             :: dcdt(NVAR), RxR(NREACT), cinit(NSPEC), lim, Cfull(NSPEC),Credux(NSPEC), RRMS
  REAL(dp)             :: A(NREACT), Prod(NVAR), Loss(NVAR)

  LOGICAL              :: OUTPUT
  LOGICAL              :: FORCE_FULL ! Force cIntegrate() to use the full mechanim. Inelegantly done
  LOGICAL              :: ReInit
  
  ! Formatting vars
  character(len=20) :: lunz, nv, clunz, cnv

  OUTPUT     = .false.
  FORCE_FULL = .false.
  REINIT     = .true.  ! Reset C every NITR,NAVG iteration
!  REINIT     = .false. ! Let C evolve over the NITR, NAVG loop

  NITR = 1
  NAVG = 100

  lim = 5e2

  R     = 0._dp
  Cinit = 0._dp
  cNONZERO = 0 ! Initialize number of nonzero elements in reduced mechanism

!  call set_quantssfc(dcdt,Cinit,R)
!  call set_quants_uppertrop(dcdt,Cinit,R)
!  call set_quants_terminator(Cinit, R)
  call set_quants_terminator7(Cinit, R)
  where (Cinit .eq. 0.d0) Cinit = 1e-20
  ! -------------------------------------------------------------------------- !
  ! 1. Reconstruct the sparse data for a reduced mechanism
  ! e.g. compact the Jacobian
  ! -- remove row & column
  ! -- DO_SLV, DO_FUN, and DO_JVS will not change size (remain NVAR & NONZERO)
  !    But the appropriate elements are set to zero, so the appropriate terms
  !    in KppSolve(), Fun() and Jac_SP() are not computed
  ! -- -- Loop through vectors
  ! -- -- -- Allocate the new Jacobian elements based on new non-zero elements
  ! -- -- Set up SPC_MAP() to map full-mech Fcn to compressed Fcn
  ! -- -- -- recompute LU_IROW
  ! -- -- -- recompute LU_ICOL
  ! -- -- -- Count the number of nonzero elements in the reduced LU diag


  ! -------------------------------------------------------------------------- !
  ! 2. Run the full mechanism

  ! Make sure everything is calculated. This should be automatic if not running
  ! a compacted mech. (This is currently also done above, but repeated here for
  ! safety. -- MSL
  write(*,*) ' '
  write(*,*) 'Running the full mechanism for initialization'
  call fullmech() ! initialization run

  call fullmech()
  Cfull = C
  write(*,*) SPC_NAMES(ind_O3), C(ind_O3), Cinit(ind_O3)
  write(*,*) SPC_NAMES(ind_OH), C(ind_OH), Cinit(ind_OH)
  write(*,*) SPC_NAMES(ind_SO4), C(ind_SO4), Cinit(ind_SO4)
  write(*,*) SPC_NAMES(ind_NO2), C(ind_NO2), Cinit(ind_NO2)
  write(*,*) SPC_NAMES(ind_Br), C(ind_Br), Cinit(ind_Br)
  write(*,*) SPC_NAMES(ind_BrO), C(ind_BrO), Cinit(ind_BrO)
  write(*,*) SPC_NAMES(ind_Cl), C(ind_Cl), Cinit(ind_Cl)
  write(*,*) SPC_NAMES(ind_ClO), C(ind_ClO), Cinit(ind_ClO)
  
  ! -------------------------------------------------------------------------- !
  ! 3. Run the compacted mechanism
  ! - Need to somehow pass the compacted vectors to KPP
  ! - Should be pretty straight-forward
  ! - Will need to respond to new 'parameters' cNVAR, cNONZERO
  ! - the compacted mechanism will still process the full species vector, 
  !   C(NVAR), not c(rNVAR), only dC/dt of inactive species is zero

  call compactedmech()
  Credux = C
  write(*,*) SPC_NAMES(ind_O3), C(ind_O3), Cinit(ind_O3)
  write(*,*) SPC_NAMES(ind_OH), C(ind_OH), Cinit(ind_OH)
  write(*,*) SPC_NAMES(ind_SO4), C(ind_SO4), Cinit(ind_SO4)
  write(*,*) SPC_NAMES(ind_NO2), C(ind_NO2), Cinit(ind_NO2)
  write(*,*) SPC_NAMES(ind_Br), C(ind_Br), Cinit(ind_Br)
  write(*,*) SPC_NAMES(ind_BrO), C(ind_BrO), Cinit(ind_BrO)
  write(*,*) SPC_NAMES(ind_Cl), C(ind_Cl), Cinit(ind_Cl)
  write(*,*) SPC_NAMES(ind_ClO), C(ind_ClO), Cinit(ind_ClO)
  
  write(*,*) ''
!>>  write(*,'(a,f6.2,a)') '  '//trim(SPC_NAMES(ind_O3))//' redux-full/full ', &
!>>       100._dp*(Credux(ind_O3)-Cfull(ind_O3))/Cfull(ind_O3), '%'
!>>  write(*,'(a,f6.2,a)') '  '//trim(SPC_NAMES(ind_OH))//' redux-full/full ', &
!>>       100._dp*(Credux(ind_OH)-Cfull(ind_OH))/Cfull(ind_OH), '%'
!>>  write(*,'(a,f6.2,a)') '  '//trim(SPC_NAMES(ind_SO4))//' redux-full/full ', &
!>>       100._dp*(Credux(ind_SO4)-Cfull(ind_SO4))/Cfull(ind_SO4), '%'
!>>  write(*,'(a,f6.2,a)') '  '//trim(SPC_NAMES(ind_NO2))//' redux-full/full ', &
!>>       100._dp*(Credux(ind_NO2)-Cfull(ind_NO2))/Cfull(ind_NO2), '%'
!>>  write(*,'(a,f6.2,a)') '  '//trim(SPC_NAMES(ind_Br))//' redux-full/full ', &
!>>       100._dp*(Credux(ind_Br)-Cfull(ind_Br))/Cfull(ind_Br), '%'
!>>  write(*,'(a,f6.2,a)') '  '//trim(SPC_NAMES(ind_BrO))//' redux-full/full ', &
!>>       100._dp*(Credux(ind_BrO)-Cfull(ind_BrO))/Cfull(ind_BrO), '%'
!>>  write(*,'(a,f6.2,a)') '  '//trim(SPC_NAMES(ind_Cl))//' redux-full/full ', &
!>>       100._dp*(Credux(ind_Cl)-Cfull(ind_Cl))/Cfull(ind_Cl), '%'
!>>  write(*,'(a,f6.2,a)') '  '//trim(SPC_NAMES(ind_ClO))//' redux-full/full ', &
!>>       100._dp*(Credux(ind_ClO)-Cfull(ind_ClO))/Cfull(ind_ClO), '%'
!>>  write(*,'(a,f6.2,a)') '  '//trim(SPC_NAMES(ind_HCl))//' redux-full/full ', &
!>>       100._dp*(Credux(ind_HCl)-Cfull(ind_HCl))/Cfull(ind_HCl), '%'
!>>  write(*,'(a,f6.2,a)') '  '//trim(SPC_NAMES(ind_BrNO3))//' redux-full/full ', &
!>>       100._dp*(Credux(ind_BrNO3)-Cfull(ind_BrNO3))/Cfull(ind_BrNO3), '%'
  do i=1,NVAR
  write(*,'(a,f6.2,a)') '  '//trim(SPC_NAMES(i))//' redux-full/full ', &
       100._dp*(Credux(i)-Cfull(i))/Cfull(i), '%'
  enddo
  ! -------------------------------------------------------------------------- !
  ! 4. (optional) Calculate the error norm per Santillana et al. (2010) and
  !    Shen et al. (2020)
  !    --- this is valuable if the species removed are selected per
  !        a reduction criterion. Requires a properly posed mechanism.

  ! -------------------------------------------------------------------------- !

  ! 5. Report timing comparison

  write(*,'(a,f5.1,a)') '  compact/full: ', 100.*compact_avg/full_avg, "%" 
  write(*,'(a,f5.1,a)') '  problem size: ', 100.*(rNVAR**2)/(NVAR**2), "%"
  write(*,'(a,f5.1,a)') '  non-zero elm: ', 100.*(cNONZERO)/(LU_NONZERO), "%"

  IF (OUTPUT) THEN
  DO i=1,rNVAR
     ii = SPC_MAP(i)
!     write(*,*) SPC_NAMES(ii) // ': ', Cfull(ii), Credux(ii), (Credux(ii)-Cfull(ii))/Cfull(ii)
  ENDDO
  ENDIF

  RRMS = sqrt(sum(((Credux(SPC_MAP(1:rNVAR))-Cfull(SPC_MAP(1:rNVAR)))/Cfull(SPC_MAP(1:rNVAR)))**2,&
       MASK=Cfull(SPC_MAP(1:rNVAR)).ne.0..and.Cfull(SPC_MAP(1:rNVAR)).gt.1e6_dp)/dble(rNVAR))

  write(*,'(a,f8.1)') 'RRMS: ', 100.*RRMS

CONTAINS

  subroutine fullmech( )
    USE GCKPP_INTEGRATOR
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
    
    ! Tolerances
    ATOL      = 1e-2_dp    
    RTOL      = 1e-2_dp

    ! Set ENV
    T    = 0d0
    TIN  = T
    TOUT = T + 600._dp
    TEMP = 298.
    
    write(*,*) ' '
    write(*,*) 'Running the full mechanism'
    
    full_avg     = 0.
    full_sumtime = 0.
    start        = 0.
    end          = 0.

    if (.not. reinit) C(1:NSPEC) = Cinit(1:NSPEC)
    ! --- INTEGRATION & TIMING LOOP
    DO I=1,NAVG
       call cpu_time(start)
       DO N=1,NITR
          if (reinit) C(1:NSPEC) = Cinit(1:NSPEC)
          ! Initialize
          call Initialize()
          VAR(1:NVAR) = C(1:NVAR)
          FIX(1:NFIX) = C(NVAR+1:NSPEC)
          ! Set RCONST
          call Update_RCONST()
          ! Integrate
          CALL Integrate( TIN,    TOUT,    ICNTRL,      &
               RCNTRL, ISTATUS, RSTATE, IERR )
          C(1:NVAR)   = VAR(:)
!          write(*,*) ISTATUS
       ENDDO
       call cpu_time(end)
!       write(*,*) 'time: ', end-start
       full_sumtime = full_sumtime+end-start
    ENDDO
    full_avg = full_sumtime/real(NAVG)
    write(*,*) "Average integration time: ", full_avg
    write(*,*) '---------------'
!    write(*,*) 'fullmech ISTATUS: ', ISTATUS(:)
 
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
    RCNTRL(8) = lim !default is 1.d2

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

    keepActive = .true.
!    keepSpcActive(ind_HNO3) = .true.
!    keepSpcActive(ind_Cl2) = .true.
!    keepSpcActive(ind_BrCl) = .true.

    if (.not.reinit) C(1:NSPEC) = Cinit(1:NSPEC)
    ! --- INTEGRATION & TIMING LOOP
    DO I=1,NAVG ! Iterate to generate average comp time
       call cpu_time(start)
       DO N=1,NITR
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
!          write(*,*) ISTATUS
       ENDDO
       call cpu_time(end)
!       write(*,*) 'time: ', end-start
       comp_sumtime = comp_sumtime+end-start
    ENDDO
    compact_avg = comp_sumtime/real(NAVG)
    write(*,*) "Average integration time: ", compact_avg
    write(*,*) '---------------'
!    write(*,*) 'compmech ISTATUS: ', ISTATUS(:)
    
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
