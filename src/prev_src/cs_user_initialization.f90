!-------------------------------------------------------------------------------

!                      Code_Saturne version
!                      --------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

subroutine cs_user_f_initialization &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use ctincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field
use turbomachinery

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

integer          iel
double precision xkini, xeini, nut

double precision, dimension(:,:), pointer :: cvar_vel
double precision, dimension(:), pointer :: cvar_k, cvar_eps
double precision, dimension(:), pointer :: cvar_fm, cvar_fp2m, cvar_hm
double precision, dimension(:), pointer :: cvar_npm, cvar_fsm

!===============================================================================


!===============================================================================
! Initialization
!===============================================================================

if (isuite.eq.0) then

!  call field_get_val_s(ivarfl(isca(ifm))  , cvar_fm)
!  call field_get_val_s(ivarfl(isca(ifp2m)), cvar_fp2m)
  if (ippmod(icod3p).eq.1) then
    call field_get_val_s(ivarfl(isca(iscalt)), cvar_hm)
  endif
!  if (isoot.eq.1) then
!    call field_get_val_s(ivarfl(isca(inpm)), cvar_npm)
!    call field_get_val_s(ivarfl(isca(ifsm)), cvar_fsm)
!  endif

  ! Mean Mixture Fraction and variance
!  do iel = 1, ncel
!    cvar_fm  (iel) = 0.d0
!    cvar_fp2m(iel) = 0.d0
!  enddo

  ! Enthalpy
  if (ippmod(icod3p).eq.1) then
    do iel = 1, ncel
      cvar_hm(iel) = hinoxy
    enddo
  endif

  ! Soot
!  if (isoot.eq.1) then
!    do iel = 1, ncel
!      cvar_npm(iel) = 0.d0
!      cvar_fsm(iel) = 0.d0
!    enddo
!  endif

endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine cs_user_f_initialization
