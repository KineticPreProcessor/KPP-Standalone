#!/bin/bash

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: build_executable.sh
#
# !DESCRIPTION: Convenience script to build the KPP standalone executable.

# !AUTHOR:
#  Bob Yantosca (yantosca@seas.harvard.edu)
#
# !REMARKS:
#  (1) Requires KPP version 3.0.0 or later.
#
# !REVISION HISTORY:
#  See the subsequent Git history with the gitk browser!
#EOP
#------------------------------------------------------------------------------
#BOC


#============================================================================
# Remove prior files that have been built with KPP
# while leaving those files containing the chemistry mechanism specification
#
# NOTE: KPP-generated source code files for GEOS-Chem will have the
# prefix "gckpp".  This is necessary because a consistent naming scheme
# needs to used so that modules in GeosCore etc. can find KPP code.
#============================================================================

# Find the mechanism name (which is at the top of the .eqn file)
mechName=$(grep ".eqn" *.eqn)
mechName="${mechName/\{ /}"
mechName="${mechName/\.eqn/}"

# Exit if mechanism name isn't found
if [[ "x${mechName}" == "x" ]]; then
    echo "Could not find the mechanism name in the ${mechanismDir}/*.eqn file"
    exit 1
fi

# Remove these files, which will be will be regnerated by KPP
filesToRemove=(             \
    gckpp.log               \
    gckpp_Function.F90      \
    gckpp_Global.F90        \
    gckpp_Initialize.F90    \
    gckpp_Integrator.F90    \
    gckpp_Jacobian.F90      \
    gckpp_JacobianSP.F90    \
    gckpp_LinearAlgebra.F90 \
    gckpp_Model.F90         \
    gckpp_Monitor.F90       \
    gckpp_Parameters.F90    \
    gckpp_Precision.F90     \
    gckpp_Rates.F90         \
    gckpp_Util.F90          \
)
for f in ${filesToRemove[@]}; do 
    rm -f $f
done

# Also remove any generated files
rm -f *.o *.mod *.MOD *.a

#============================================================================
# Build the mechanism!
#============================================================================
if [[ -f gckpp.kpp ]]; then
    kpp gckpp.kpp
    # on discover, use:
    #/discover/nobackup/mslong1/KPP/KPP/bin/kpp gckpp.kpp
else
    echo "Could not find the 'gckpp.kpp' file... Aborting!"
    exit 1
fi

# If the gckpp_Rates.F90 file is not found, there was an error
if [[ ! -f gckpp_Rates.F90 ]]; then
  echo "KPP failed to build 'gckpp_Rates.F90'! Aborting."
  exit 1
fi

#============================================================================
# Build the executable
#============================================================================
make
if [[ ! -f kpp_standalone.exe ]]; then
    echo "Failed to build the 'kpp_executable.exe' file!  Aborting..."
    exit 1
fi

#============================================================================
# Exit and return status to shell
#============================================================================
exit $?
#EOC
