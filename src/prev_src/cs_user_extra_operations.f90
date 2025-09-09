!-------------------------------------------------------------------------------

!                      Code_Saturne version 4.0.0-patch
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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_extra_operations.f90
!>
!> \brief This function is called at the end of each time step, and has a very
!>  general purpose
!>  (i.e. anything that does not have another dedicated user subroutine)
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine cs_f_user_extra_operations &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use lagran
use parall
use period
use ppppar
use ppthch
use ppincl
use coincl
use mesh
use field
use field_operator
use turbomachinery
use cs_c_bindings
use radiat
use user_module

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

character*40     chain, varnam
integer          impout, ivar  , ifac
integer          iel   , iel1  , ii    , jj
integer          irangv, irang1, npoint
integer          imom  , idtcm
integer          ipcmo1, ipcmo2, ipcmo3, ipcmo4, ipcmo5, ipcmo6
integer          ipcmo7, ipcmo8, ipcmo9

integer          n_wt

double precision xyz(3), xabs, xu, xv, xw
double precision xx(3), yy(3), zz(3), dtm
double precision xt, xz, xzp2, xyf, xyo, xyp

integer          nlelt, ilelt
double precision flumab, flumas, fluadm, fluext, flumb, flutot, qw
double precision hb, hf, hi, hj
double precision Qfuel, Qadm, Qext, Qwall, Qtot

integer          izone
double precision mfuel, madm, mext, rho, Sfuel

double precision wt, wall_temp, wf, wall_flux

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp

integer          iscal
integer          ipcvst
integer          nswrgp, imligp, iwarnp
integer          iccocg, inc
integer          iconvp, idiffp, ircflp
integer          ischcp, isstpp, ippvar
integer          ipcvsl, iflmas, iflmab
integer          imucpp, idftnp, ipcyf, ipcyfa

double precision hint, hext, pimp
double precision epsrgp, climgp, extrap
double precision blencp, relaxp, thetex

double precision en_masfl, en_o2_masfl
double precision so_masfl, so_o2_masfl

double precision Yfm, Yfa, Yfs, fm
double precision HRR, PCI, HRR2, HRR3

double precision rvoid(1)

double precision, allocatable, dimension(:) :: vistot, viscf, viscb, dYfdt
double precision, allocatable, dimension(:) :: coefay, cofafy

double precision, dimension(:), pointer :: crom, brom, bmasfl, cvar_wt, cpro_wt
double precision, dimension(:), pointer :: bpro_var, cvar_qinc, cvar_qconv

double precision, dimension(:), pointer :: cpro_Yf, cpro_Yo, cpro_Yp
double precision, dimension(:), pointer :: cpro_XO2, cpro_XCO2, cpro_CSoot

double precision, allocatable, dimension(:) :: yg, ye, xe

double precision, dimension(:), pointer :: bym2

double precision, dimension(:), pointer :: cvar_fm, bfconv, bhconv, bfqinc



!===============================================================================
! User operations
!===============================================================================

! BS affichage du temps
if ((ntcabs.eq.ntpabs+1.or.ntcabs.eq.ntmabs.or. &
    mod(ntcabs,ntlist).eq.0).and.irangp.le.0 ) then
    !print'(a10,i8,a10,f14.5)',"Iteration",ntcabs,"Time",ttcabs
endif

! Density
call field_get_val_s(ibrom, brom)
call field_get_val_s(icrom, crom)

! Mass fluxes
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmab, bmasfl)

!
call field_get_val_s(itempb, cpro_wt)

!===============================================================================
! Concentration postprocessing
!===============================================================================

call field_get_val_s_by_name("Xm_O2", cpro_XO2)
call field_get_val_s_by_name("Xm_CO2", cpro_XCO2)
call field_get_val_s_by_name("Cm_Soot", cpro_CSoot)

call field_get_val_s_by_name("rad_convective_flux", bfconv)
call field_get_val_s_by_name("rad_exchange_coefficient", bhconv)
call field_get_val_s_by_name("rad_incident_flux", bfqinc)

call field_get_val_s_by_name("wall_incident_flux", cvar_qinc)
call field_get_val_s_by_name("wall_convective_flux", cvar_qconv)
call field_get_val_s_by_name("wall_temp", cvar_wt)

call field_get_val_s(iym(1), cpro_Yf)
call field_get_val_s(iym(2), cpro_Yo)
call field_get_val_s(iym(3), cpro_Yp)

allocate (yg(ngazg), ye(ngaze), xe(ngaze))

do iel = 1, ncel

  yg(1) = cpro_Yf(iel)
  yg(2) = cpro_Yo(iel)
  yg(3) = cpro_Yp(iel)

  call yg2xye (yg, ye, xe)

  cpro_XO2(iel) = xe(iio2)
  cpro_XCO2(iel) = xe(iico2)
  cpro_Csoot(iel) = xe(iic) * crom(iel)

enddo

deallocate (yg, ye, xe)

call field_get_val_s(ibym(2), bym2)

call field_get_val_s(ivarfl(isca(ifm)), cvar_fm)

!===============================================================================
! Fuel mass flow rate, mass flow rates through the ventilation and room pressure
!===============================================================================

allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

if (ntcabs.eq.ntpabs+1.or.ntcabs.eq.ntmabs.or.mod(ntcabs,ntlist).eq.0) then

  mfuel  = 0.d0
  madm   = 0.d0
  mext   = 0.d0
  wall_temp = 0.d0
  wall_flux = 0.d0
  en_masfl = 0.d0
  en_o2_masfl = 0.d0
  so_masfl = 0.d0
  so_o2_masfl = 0.d0

  
  call getfbr('POOL_B', nlelt, lstelt)

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    flumab = bmasfl(ifac)
   
    mfuel = mfuel - flumab

  enddo

  call getfbr('INTAKE', nlelt, lstelt)

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    flumab = bmasfl(ifac)
   
    rho  = brom(ifac)

    madm = madm - 3600.d0 * flumab/rho

    en_masfl = en_masfl + bmasfl(ifac)

    en_o2_masfl = en_o2_masfl + bmasfl(ifac)*bym2(ifac)

    write(nfecra,*) "INTAKE - bmasfl, brom", bmasfl(ifac), brom(ifac)

  enddo

  call getfbr('EXHAUST', nlelt, lstelt)

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    iel  = ifabor(ifac)

    flumab = bmasfl(ifac)
   
    rho  = brom(ifac)

    mext = mext + 3600.d0 * flumab/rho

    so_masfl = so_masfl + bmasfl(ifac)

    so_o2_masfl = so_o2_masfl + bmasfl(ifac)*cvar_fm(iel)!bym2(ifac)

    write(nfecra,*) "EXHAUST - bmasfl, brom", bmasfl(ifac), brom(ifac)

  enddo

  n_wt = 0

  call getfbr('WALL OR CEILING', nlelt, lstelt)

  n_wt = nlelt

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    wt = cpro_wt(ifac)

    cvar_wt(ifac) = cpro_wt(ifac)

    cvar_qinc(ifac) = bfqinc(ifac)

    cvar_qconv(ifac) = bfconv(ifac)

    wall_temp = wall_temp + wt

    wall_flux = wall_flux + wf

  enddo
   
  if (irangp.ge.0) call parsom(mfuel)
  if (irangp.ge.0) call parsom(madm)
  if (irangp.ge.0) call parsom(mext)
  if (irangp.ge.0) call parsom(wall_temp)
  if (irangp.ge.0) call parsom(wall_flux)
  if (irangp.ge.0) call parcpt(n_wt)
  if (irangp.ge.0) call parsom(en_masfl)
  if (irangp.ge.0) call parsom(en_o2_masfl)
  if (irangp.ge.0) call parsom(so_masfl)
  if (irangp.ge.0) call parsom(so_o2_masfl)

  write(nfecra,*) "INTAKE - en_masfl en_o2_masfl", en_masfl, en_o2_masfl
  write(nfecra,*) "EXHAUST - so_masfl so_fm_masfl", so_masfl, so_o2_masfl

  if (irangp.le.0) then
    open(impusr(1),file="results.dat",position="append")
    if (ntcabs.eq.ntpabs+1) write(impusr(1),10) "TIME", "mfuel", "dP", "DADM_MOY", "DEXT_MOY", "WALL_TEMP", "WALL_FLUX"
    if (ntcabs.eq.ntpabs+1) write(impusr(1),10)  "(s)", "(kg/s)", "(Pa)", "(m3/h)"  , "(m3/h)", "deg.C", "W/m2"
    write(impusr(1),11) ttcabs, mfuel, pther-p0, madm, mext, wall_temp/n_wt, wall_flux/n_wt
    close(1)
  endif

10 format("#",a13,99(a14))
11 format(100(e14.5))

  ! !=============================================================================
  ! ! Power balance
  ! !=============================================================================

  ! call cs_user_power_balance (dt)

  ! !=============================================================================
  ! ! Heat release rate
  ! !=============================================================================

  ! call cs_user_heat_release_rate (dt)


endif

deallocate(lstelt)

!--------
! Formats
!--------
write(nfecra, *) 'Pther:', Pther
!===============================================================================
! free user array
!===============================================================================

if (ntcabs.eq.ntmabs) call finalize_user_module

return
end subroutine cs_f_user_extra_operations
