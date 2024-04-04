module setquants
  implicit none
  public
contains

  subroutine setRTOLandC
    use gckpp_Model

    RTOL(ind_I2O2) = -0.14912368
    RTOL(ind_I2O4) = 0.42927492
    RTOL(ind_IO) = 0.18233833
    RTOL(ind_O) = 0.43380463
    RTOL(ind_HNO4) = 0.004950523
    RTOL(ind_OIO) = 0.083717205
    RTOL(ind_Cl) = 0.96560603
    RTOL(ind_ClOO) = 0.737162
    RTOL(ind_I2O3) = 0.5185118
    RTOL(ind_HOBr) = 0.4525018
    RTOL(ind_HO2) = 0.2986274
    RTOL(ind_OH) = 0.6277464
    RTOL(ind_I) = 0.115892604
    RTOL(ind_Br) = 0.41638651
    RTOL(ind_CO2) = 0.2994907
    RTOL(ind_ClO) = -0.09867178
    RTOL(ind_IONO2) = 0.06640734
    RTOL(ind_N2O5) = 0.35829496
    RTOL(ind_IONO) = 0.5930198
    RTOL(ind_OClO) = 0.9927597

  end subroutine setRTOLandC

end module setquants
