module initialize_final
  implicit none
  public
contains

  ! This will be passed by the $(CPP) to get a condition from the archived folder
  subroutine initialize_sunrise(C, R)
    USE gckpp_Parameters
    IMPLICIT NONE
    real(dp), intent(out) :: C(NSPEC)
    real(dp), intent(out) :: R(NREACT)
#include "archive_13.4_conditions/IOSR_109_45_1.txt"
  end subroutine initialize_sunrise

  ! This will be passed by the $(CPP) to get a condition from the archived folder
  subroutine initialize_sunrise_500(C, R)
    USE gckpp_Parameters
    IMPLICIT NONE
    real(dp), intent(out) :: C(NSPEC)
    real(dp), intent(out) :: R(NREACT)
#include "archive_13.4_conditions/IOSR_109_45_23.txt"
  end subroutine initialize_sunrise_500

  ! This will be passed by the $(CPP) to get a condition from the archived folder
  subroutine initialize_na_daytime(C, R)
    USE gckpp_Parameters
    IMPLICIT NONE
    real(dp), intent(out) :: C(NSPEC)
    real(dp), intent(out) :: R(NREACT)
#include "archive_13.4_conditions/NAsfc_40_66_1.txt"
  end subroutine initialize_na_daytime

  ! This will be passed by the $(CPP) to get a condition from the archived folder
  subroutine initialize_satl_nighttime(C, R)
    USE gckpp_Parameters
    IMPLICIT NONE
    real(dp), intent(out) :: C(NSPEC)
    real(dp), intent(out) :: R(NREACT)
#include "archive_13.4_conditions/SA500_66_34_23.txt"
  end subroutine initialize_satl_nighttime

end module initialize_final
! hplin 4/11/22
! IF (Scenario .eq. 1) call initialize_sunrise(Cinit, R)
! IF (Scenario .eq. 2) call initialize_sunrise_500(Cinit, R)
! IF (Scenario .eq. 3) call initialize_na_daytime(Cinit, R)
! IF (Scenario .eq. 4) call initialize_satl_nighttime(Cinit, R)