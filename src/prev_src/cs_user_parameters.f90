!-------------------------------------------------------------------------------

!                      Code_Saturne version 6.0-beta
!                      --------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!===============================================================================

!> \file cs_user_parameters.f90
!>
!> \brief User subroutines for input of calculation parameters (Fortran modules).
!>        These subroutines are called in all cases.
!>
!>  See \subpage f_parameters for examples.
!>
!>   If the Code_Saturne GUI is used, this file is not required (but may be
!>   used to override parameters entered through the GUI, and to set
!>   parameters not accessible through the GUI).
!>
!>   Several routines are present in the file, each destined to defined
!>   specific parameters.
!>
!>   To modify the default value of parameters which do not appear in the
!>   examples provided, code should be placed as follows:
!>   - usipsu   for numerical and physical options
!>   - usipes   for input-output related options
!>
!>   As a convention, "specific physics" defers to the following modules only:
!>   pulverized coal, gas combustion, electric arcs.
!>
!>   In addition, specific routines are provided for the definition of some
!>   "specific physics" options.
!>   These routines are described at the end of this file and will be activated
!>   when the corresponding option is selected in the usppmo routine.
!-------------------------------------------------------------------------------

!===============================================================================

!> \brief User subroutine for selection of specific physics module

!> Define the use of a specific physics amongst the following:
!>   - combustion with gas / coal / heavy fuel oil
!>   - compressible flows
!>   - electric arcs
!>   - atmospheric modelling
!>   - radiative transfer
!>   - cooling towers modelling
!>
!>    Only one specific physics module can be activated at once.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixmlpu        indicates if the XML file from the GUI is used
!>                              (1 : yes, 0 : no)
!______________________________________________________________________________!

subroutine usppmo &
 ( ixmlpu )

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use cstphy
use ppppar
use ppthch
use ppincl
use ppcpfu
use coincl
use radiat
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer ixmlpu

!===============================================================================
! 1.  Choice for a specific physics
!===============================================================================

! --- Soot model
! =================

!        if = -1   module not activated
!        if =  0   constant fraction of fuel Xsoot
!        if =  1   2 equations model of Moss et al.

isoot = 0

xsoot  = 0.1d0 ! ( if isoot = 0 )
rosoot = 2000.d0 ! kg/m3

!----
! End
!----

return
end subroutine usppmo

!===============================================================================

!> \brief User subroutine for the input of additional user parameters.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp         number of active specific physics models
!______________________________________________________________________________!

subroutine usipsu &
 ( nmodpp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use albase
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use field
use cavitation
use post
use rotation
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

integer       iscal, ifcvsl, isc
integer       f_id, idim1, itycat, ityloc
logical       ilved, inoprv

type(var_cal_opt) :: vcopt

!===============================================================================

!>  This subroutine allows setting parameters
!>  which do not already appear in the other subroutines of this file.
!>
!>  It is possible to add or remove parameters.
!>  The number of physical properties and variables is known here.

!===============================================================================


! Calculation options (optcal)
! ============================

! --- Algorithm to take into account the density variation in time
!
!     idilat = 0 : Boussinesq algorithm with constant density (not available)
!              1 : dilatable steady algorithm (default)
!              2 : dilatable unsteady algorithm
!              3 : low-Mach algorithm

idilat = 2

! --- Algorithm to take into account the thermodynamical pressure variation in time
!     (not used by default except if idilat = 3)

!     by default:
!     ----------
!      - the thermodynamic pressure (pther) is initialized with p0 = p_atmos
!      - the maximum thermodynamic pressure (pthermax) is initialized with -1
!        (no maximum by default, this term is used to model a venting effect when
!         a positive value is given by the user)
!      - a global leakage can be set through a leakage surface sleak with a head
!      loss kleak of 2.9 (Idelcick)

ipthrm = 1

sleak = 0.0025d0
kleak = 2.9d0

! --- Diffusive scheme

imvisf = 1

! Physical constants (cstphy)
! ===========================

! --- Reference fluid properties

!       t0         : reference temperature in Kelvin
!       Pther      : thermodynamic pressure, variable in time

t0    = 26.d0 + tkelvi
Pther = P0 - 6.23d0

! --- Variable diffusivity field id (ifcsl>=0) or constant
!     diffusivity (ifcvsl=1) for the thermal scalar and USER scalars.

do isc = 1, nscapp
  iscal = iscapp(isc)
  if (iscavr(iscal).le.0) then
    ifcvsl = 0
    call field_set_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
  endif
enddo

! Postprocessing-related fields
! =============================

itycat = FIELD_INTENSIVE + FIELD_PROPERTY
ityloc = 1 ! cells, 3 boundary faces
idim1 = 1 ! dimension
ilved = .true. ! interleaved
inoprv = .false. ! no previous time step values needed

call field_create('Xm_O2', itycat, ityloc, idim1, inoprv, f_id)
call field_set_key_int(f_id, keyvis, POST_ON_LOCATION + POST_MONITOR)
call field_set_key_int(f_id, keylog, 1)

call field_create('Xm_CO2', itycat, ityloc, idim1, inoprv, f_id)
call field_set_key_int(f_id, keyvis, POST_ON_LOCATION + POST_MONITOR)
call field_set_key_int(f_id, keylog, 1)

call field_create('Cm_Soot', itycat, ityloc, idim1, inoprv, f_id)
call field_set_key_int(f_id, keyvis, POST_ON_LOCATION + POST_MONITOR)
call field_set_key_int(f_id, keylog, 1)

!----
! Formats
!----

return
end subroutine usipsu


!===============================================================================

!> \brief User subroutine for the input of additional user parameters for
!>        input/output.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp       number of active specific physics models
!______________________________________________________________________________!

subroutine usipes &
 ( nmodpp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use field
use parall
use period
use post
use ppppar
use ppthch
use ppincl
use cs_c_bindings
use radiat
use user_module

!===============================================================================

implicit none

! Arguments

integer nmodpp

!===============================================================================

integer ii

!===============================================================================

!>     This subroutine allows setting parameters
!>     which do not already appear in the other subroutines of this file.
!>
!>     It is possible to add or remove parameters.
!>     The number of physical properties and variables is known here.

!===============================================================================

!===============================================================================
! Fine control of variables output
!===============================================================================

! Per variable output control.
! More examples are provided in cs_user_parameters-output.f90

! User scalar variables.

call field_set_key_int(itemp, keyvis, POST_ON_LOCATION + POST_MONITOR)
call field_set_key_int(itemp, keylog, 1)

do ii = 1, 3
  call field_set_key_int(iym(ii), keyvis, POST_ON_LOCATION + POST_MONITOR)
  call field_set_key_int(iym(ii), keylog, 1)
enddo

return
end subroutine usipes

!===============================================================================
! Purpose:
! -------
!
!> 1. Additional Calculation Options
!>    a. Density Relaxation
!>
!> 2. Physical Constants
!>    a.Dynamic Diffusion Coefficient
!>    b.Constants of the chosen model (EBU, Libby-Williams, ...)
!
!> This routine is called:
!>
!>
!>  - Eddy Break Up pre-mixed flame
!>  - Diffusion flame in the framework of ``3 points'' rapid complete chemistry
!>  - Libby-Williams pre-mixed flame
!>  - Lagrangian module coupled with pulverized coal:
!>    Eulerian combustion of pulverized coal and
!>    Lagrangian transport of coal particles
!>  - Pulverised coal combustion
!>  - Fuel (oil) combustion
!
!===============================================================================

subroutine cs_user_combustion

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use cs_c_bindings
use radiat

!===============================================================================

implicit none

integer ivar, iscal, isc

type(var_cal_opt) :: vcopt
double precision     turb_schmidt

!===============================================================================

!===============================================================================
! 1. Additional Calculation Options
!===============================================================================

! -->  Density Relaxation
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.d0

! --> Convective scheme

do iscal = 1,nscal
  ivar = isca(iscal)
  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
  vcopt%blencv = 0.d0
  vcopt%isstpc = 0
  vcopt%ischcv = 0
  call field_set_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
enddo

!===============================================================================
! 2. Physical Constants
!===============================================================================

! diftl0: Dynamic Diffusion Coefficient (kg/(m s))
diftl0 = 4.25d-5

! Reference temperature for fuel and oxydant (K)
tinfue = t0
tinoxy = t0 !- 5.d0

! Set a variable laminar viscosity defined in cs_user_physical properties
ivivar = 1

! Diffusion coefficient
do isc = 1, nscapp
  iscal = iscapp(isc)
  if (iscavr(iscal).le.0) then
    ! Constant diffusion coefficient
    call field_get_key_double(ivarfl(isca(iscal)), ksigmas, turb_schmidt)
    call field_set_key_double(ivarfl(isca(iscal)), kvisl0, 0.1d0)
  endif
enddo

!----
! End
!----

return
end subroutine cs_user_combustion
