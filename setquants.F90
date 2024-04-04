module setquants
  implicit none
  public
contains

  subroutine setRTOLandC
    use gckpp_Model

RTOL(ind_I2O2) = 0.23842764
RTOL(ind_I2O4) = 0.35788715
RTOL(ind_IO) = 0.25153738
RTOL(ind_O) = 0.2532347
RTOL(ind_HNO4) = 0.26788014
RTOL(ind_OIO) = 0.2536648
RTOL(ind_Cl) = 0.2353785
RTOL(ind_ClOO) = 0.15982367
RTOL(ind_I2O3) = 0.22536644
RTOL(ind_HOBr) = 0.13655436
RTOL(ind_HO2) = 0.13075581
RTOL(ind_OH) = 0.39678
RTOL(ind_I) = 0.21552496
RTOL(ind_Br) = 0.17561801
RTOL(ind_CO2) = 0.18560833
RTOL(ind_ClO) = 0.12499675
RTOL(ind_IONO2) = 0.27637592
RTOL(ind_N2O5) = 0.24338068
RTOL(ind_IONO) = 0.39919004
RTOL(ind_OClO) = 0.20721224

  end subroutine setRTOLandC

end module setquants
