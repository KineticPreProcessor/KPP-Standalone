program main

  use GCKPP_MODEL

  USE INITIALIZE

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
  integer :: sampleunit, outunit, wtf

  write(*,*) ' '
  write(*,*) 'THIS PROGRAM IS NOT REALLY VALUABLE <<MSL>>'
  write(*,*) ' '

  ! Set up and populate the test data file
  ! -- open the file to write filenames of samples to be used
  open(9992, file='twilightfiles.txt', action='write', status='replace')
  open(9993, file='strattwilightfiles.txt', action='write', status='replace')

  ! -- open the file will all the sample filenames in it
  open(998, file='sampleslist.txt', status='old', action='read')
  read(998, '(I10)', iostat=iostat) NFILES
  write(*,*) 'NFILES for testing: ', NFILES

!  write(*,*) '<<>> ', wtf, sampleunit

  DO II = 1,NFILES
!     write(*,*) II
     ! Read the input
     read(998, '(A)', iostat=iostat) filename
!     write(*,*) filename
     call read_input(filename, R, Cinit, SPC_NAMES, Hstart, cosSZA, level, fileTotSteps, OperatorTimestep)
     if (abs(cosSZA) .le. 0.139173101) write(9992, '(A)') trim(filename) ! All twilight files
     if (abs(cosSZA) .le. 0.139173101 .and. level .ge. 37 ) write(9993, '(A)') trim(filename) ! Strat-only twilight files
  ENDDO

  close(998)
  close(unit=9992)
  close(unit=9993)

end program main
