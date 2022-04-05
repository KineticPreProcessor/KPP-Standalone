KPP Auto-reduction test bed
Derived from KPP-Compressor box model
Citation: Lin, H., et al. ... 

M.S. Long - Feb 10, 2022 -- initial release version
H.P. Lin  - Apr 4, 2022 -- Add some documentation and update to v13.4

# Instructions to build
Build the mechanism using `kpp gckpp.kpp`.

Compile it using `make`, then `./gckpp.exe`


## Building based on a new GEOS-Chem mechanism
Copy over `fullchem.eqn` to target version, e.g., `v13.4.eqn`.

* To allow forcing of archived rate constants, do regex magic within `gckpp_Rates.F90`, replacing `RCONST\((\d+)\) =.*$` to `RCONST(\1) = R(\1)`

## Example code for archiving rates and concentrations from fullchem simulation
Add to `fullchem_mod.F90` after `Update_RConst()` call:

```fortran
IF (L .eq. 1 .and. I .eq. 134 .and. J .eq. 55) THEN
 OPEN(998,FILE="boxdebug_134_55_1.txt")
 DO N=1,NREACT
   write(998,*) 'R(',N,') = ', RCONST(N)
 ENDDO
 DO N=1,NSPEC
   write(998,*) 'C(',N,') = ', C(N)
 ENDDO
 CLOSE(998)
 write(*,*) 'Written...'
ENDIF
```

Paste this into `initialize.F90`