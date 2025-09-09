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
! Function:
! ---------

!> \file cs_user_boundary_conditions.f90
!>
!> \brief User subroutine which fills boundary conditions arrays
!> (\c icodcl, \c rcodcl) for solved variables.

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[out]    icodcl        boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!>                               - 13 Dirichlet for the advection operator and
!>                                    Neumann for the diffusion operator
!> \param[in]     itrifb        indirection for boundary faces ordering
!> \param[in,out] itypfb        boundary face types
!> \param[out]    izfppp        boundary face zone number
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradt \, \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!_______________________________________________________________________________

subroutine cs_f_user_boundary_conditions &
 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     ,                                                       &
   rcodcl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
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
use ppcpfu
use atchem
use atincl
use atsoil
use ctincl
!use cs_fuel_incl
use mesh
use field
!use turbomachinery
use iso_c_binding
use cs_c_bindings

!use user_module

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ifac, iel
integer          izone
integer          ilelt, nlelt
integer          iflmab

integer          iclke
double precision uref2, ustar2, d2s3
double precision xdh, xitur, viscl
double precision xkent, xeent, nut, rho

integer          mode, icg
double precision xx, yy, zz, rinlet, minlet, Qinlet, tinlet
double precision coefg(ngazgm)

double precision Kadm, madm, Padm, Qadm, dadm, dQadmdt, radm
double precision Kext, mext, Pext, Qext, dext, dQextdt, rext
double precision Ploc, S0, dp, tstart

double precision madd, s1


character*(200)  chain, fname
integer          input, ios, ipyro, ivent, ipth
double precision tt, t1, t2
double precision padm_temp, pext_temp, m1, m2
double precision padm1, pext1, padm2, pext2, dmdt
double precision dvent, Q1, Q2, P1, P2, padm_act, pext_act
double precision vadm, vext, dpyro

double precision, dimension(:,:), pointer :: cvar_vel, cvara_vel
double precision, dimension(:), pointer :: cpro_viscl
double precision, dimension(:), pointer :: crom, brom
double precision, dimension(:), pointer :: bmasfl

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

tstart = 100.d0


! Ventilation pressure
! **************

rinlet = 0.36d0

fname="pres_vent.csv"

if (ntcabs.eq.ntpabs+1) then

  ! Initialisation for allocation
  npyro = 0
  nvent = 0
  npth  = 0

  ! open file
  input = 1
  open(input,file=trim(fname),status="old",iostat=ios)
  if (ios.ne.0) then
    print'(a,1x,a,/)',"Problem opening file", trim(fname)
    call csexit(1)
  endif

  ! skip headfile
  !read(input,'(a)',iostat=ios) chain
  !read(input,'(a)',iostat=ios) chain
  !read(input,'(a)',iostat=ios) chain

  if (ios.ne.0) then
    print'(a,1x,a,/)', "Problem reading head of file", trim(fname)
    call csexit(1)
  endif

  ! compute number of lines
  ipyro = 0
  do while (ios.eq.0)
    read(input,*,iostat=ios) tt, padm_temp, pext_temp
    ipyro = ipyro + 1
  enddo
  nvent = ipyro-1
  if (nvent.eq.0) then
    print*,"Npyro = ", nvent," > STOP"
    call csexit(1)
  endif

  ! allocate memory, desallocate in cs_user_extra_operations
  call init_user_module (ncel, ncelet)

  ! fill data array
  rewind(input)

  !read(input,'(a)',iostat=ios) chain
  !read(input,'(a)',iostat=ios) chain
  !read(input,'(a)',iostat=ios) chain

  ipyro = 0
  do while (ipyro.lt.nvent)
    read(input,*,iostat=ios) tt, padm_temp, pext_temp
    ipyro = ipyro + 1
    if (ios.ne.0) then
      print'(a,1x,i4,a,1x,a/)',"Problem reading line",ipyro,"of file",trim(fname)
    else
      Qvent(ipyro,1) = tt
      Qvent(ipyro,2) = padm_temp ! donnee en kg/s
      Qvent(ipyro,3) = pext_temp ! donnee en kg/s
    endif
  enddo

  close(input)

endif

ipyro = 0
tt = ttcabs - tstart

do while (ipyro.lt.npyro-1)

  ipyro = ipyro + 1

  t1 = Qvent(ipyro  ,1)
  t2 = Qvent(ipyro+1,1)

  padm1 = Qvent(ipyro  ,2)
  padm2 = Qvent(ipyro+1,2)

  pext1 = Qvent(ipyro  ,3)
  pext2 = Qvent(ipyro+1,3)

  if (t1.le.tt.and.tt.le.t2) then
    padm_act = padm1 + (padm2-padm1)/(t2-t1)*(tt-t1)
    pext_act = pext1 + (pext2-pext1)/(t2-t1)*(tt-t1)
  endif

enddo



! Mass loss rate
! **************

rinlet = 0.36d0

fname="PRS_SI_D6_MLR.txt"

if (ntcabs.eq.ntpabs+1) then

  ! Initialisation for allocation
  npyro = 0
  nvent = 0
  npth  = 0

  ! open file
  input = 1
  open(input,file=trim(fname),status="old",iostat=ios)
  if (ios.ne.0) then
    print'(a,1x,a,/)',"Problem opening file", trim(fname)
    call csexit(1)
  endif

  ! skip headfile
  read(input,'(a)',iostat=ios) chain
  read(input,'(a)',iostat=ios) chain
  read(input,'(a)',iostat=ios) chain

  if (ios.ne.0) then
    print'(a,1x,a,/)', "Problem reading head of file", trim(fname)
    call csexit(1)
  endif

  ! compute number of lines
  ipyro = 0
  do while (ios.eq.0)
    read(input,*,iostat=ios) tt, dpyro
    ipyro = ipyro + 1
  enddo
  npyro = ipyro-1
  if (npyro.eq.0) then
    print*,"Npyro = ", npyro," > STOP"
    call csexit(1)
  endif

  ! allocate memory, desallocate in cs_user_extra_operations
  call init_user_module (ncel, ncelet)

  ! fill data array
  rewind(input)

  read(input,'(a)',iostat=ios) chain
  read(input,'(a)',iostat=ios) chain
  read(input,'(a)',iostat=ios) chain

  ipyro = 0
  do while (ipyro.lt.npyro)
    read(input,*,iostat=ios) tt, dpyro
    ipyro = ipyro + 1
    if (ios.ne.0) then
      print'(a,1x,i4,a,1x,a/)',"Problem reading line",ipyro,"of file",trim(fname)
    else
      mpyro(ipyro,1) = tt
      mpyro(ipyro,2) = dpyro ! donnee en kg/s
    endif
  enddo

  close(input)

endif
!
! compute current mass loss rate from ttpabs, last instant of the initial
! calculation used to establish the ventilation flow

ipyro = 0
tt = ttcabs - tstart

do while (ipyro.lt.npyro-1)

  ipyro = ipyro + 1

  t1 = mpyro(ipyro  ,1)
  t2 = mpyro(ipyro+1,1)

  m1 = mpyro(ipyro  ,2)
  m2 = mpyro(ipyro+1,2)

  if (t1.le.tt.and.tt.le.t2) then
    dmdt = m1 + (m2-m1)/(t2-t1)*(tt-t1)
  endif

  if (ttcabs.ge.2604.d0) then
    dmdt = 0
  endif

enddo

!if (tt.le.22) then
!  dmdt = 0
!else
!  dmdt = 0.007
!endif

! Pertes de charge dans le reseau de ventilation
! **********************************************

! Les conditions initiales de SI-D1 sont utilisees pour calculer le coefficient de pertes
! de charges de la conduite d'extraction

Ploc= -6.23d0 ! Pa
S0 = 0.18 ! m2
S1 = acos(-1.d0) * (75.d-3)**2 ! m2

radm = (P0+ploc)/(cs_physical_constants_r*tinoxy/wmolg(2))
Padm = 509.73d0 ! Pa
Qadm = 560.d0 ! m3/h

madm = Qadm*radm/3600.d0

Kadm = 2.d0*radm*(S0/madm)**2 * abs(Ploc-Padm) !- 50.d0

rext = (P0+ploc)/(cs_physical_constants_r*t0/wmolg(2))
Pext = -708.7d0 ! Pa
Qext = 560.d0 ! m3/h

mext = Qext*rext/3600.d0

Kext = 2.d0*rext*(S0/mext)**2 * abs(Ploc-Pext) !+ 50.d0

if (ntcabs.eq.1.and.irangp.le.0) then
  write(*,'(/,a,/)')  " Coefficients de pertes de charges : "
  write(*,'(a,E14.5)')"    Admission  K = ", Kadm
  write(*,'(a,E14.5,/)')"    Extraction K = ", Kext

  write(NFECRA,'(a,/)')  " RESEAU DE VENTILATION :"
  write(NFECRA,'(a,/)')  " Coefficients de pertes de charges : "
  write(NFECRA,'(a,E14.5)')"    Admission  K = ", Kadm
  write(NFECRA,'(a,E14.5)')"    Extraction K = ", Kext
endif

!===============================================================================
! Field mapping
!===============================================================================

! Velocity at previous time step
call field_get_val_v(ivarfl(iu), cvar_vel)
call field_get_val_prev_v(ivarfl(iu), cvara_vel)

! Cell viscosity
call field_get_val_s(iprpfl(iviscl), cpro_viscl)

! Density
call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

! Mass fluxes
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmab, bmasfl)

allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

d2s3 = 2.d0/3.d0

!===============================================================================
! Assign boundary conditions to boundary faces here

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - loop on faces from a subset
!   - set the boundary condition for each face
!===============================================================================

! BS needed by the combustion model when there is no oxydant inlet and for the wall enthalpy
! tinoxy = t0 - 5.d0
do icg = 1, ngazgm
  coefg(icg) = 0.d0
enddo
coefg(2) = 1.d0
mode    = -1
call cothht                                                   &
!==========
  ( mode   , ngazg , ngazgm  , coefg  ,                     &
    npo    , npot   , th     , ehgazg ,                     &
    hinoxy , tinoxy  )


! --- Wall by default

call getfbr('WALL OR CEILING OR FLOOR', nlelt, lstelt)

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ! Wall: zero flow (zero flux for pressure)
  !       friction for velocities (+ turbulent variables)
  !       zero flux for scalars

  itypfb(ifac)   = iparoi

  ! Zone number (arbitrary number between 1 and n)
  !izone = 5

  ! Allocation of the actual face to the zone
  !izfppp(ifac) = izone

enddo

! --- Fuel inlet

nlelt = 0
call getfbr('POOL_B', nlelt, lstelt)

write(NFECRA,*)  " MLR AT TIME", ttcabs, "S :", dmdt

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)
  iel = ifabor(ifac)

  xx = cdgfbo(1,ifac)
  yy = cdgfbo(2,ifac)
  zz = cdgfbo(3,ifac)

!  itypfb(ifac)   = iparoi


!  if ((xx**2.d0+yy**2.d0).lt.rinlet**2.d0) then

    itypfb(ifac) = i_convective_inlet

    ! Zone number (arbitrary number between 1 and n)
    !izone = 1

    ! Allocation of the actual face to the zone
    !izfppp(ifac) = izone
    izone = izfppp(ifac)

    ! Indicating the inlet as a fuel flow inlet
    ientfu(izone) = 1

    ! Inlet Temperature in K
    tinfue = t0

    ! The incoming fuel flow refers to:
    ! a) a massflow rate   -> iqimp()  = 1
    iqimp(izone)  = 1
    qimp(izone)   = dmdt

    ! b) an inlet velocity -> iqimp()  = 0

    rcodcl(ifac,iu,1) = 0.d0
    rcodcl(ifac,iv,1) = 0.d0
    rcodcl(ifac,iw,1) = dmdt/(brom(ifac)*0.4)

    !   Boundary conditions of turbulence
    icalke(izone) = 1

    ! If ICALKE = 0 the boundary conditions of turbulence
    ! are computed by the user
    if (icalke(izone).eq.0) then

      xkent = 1.d-10
      nut   = 1.d-10
      xeent = cmu*xkent**2/nut

      if (itytur.eq.2) then
        rcodcl(ifac,ik,1)  = xkent
        rcodcl(ifac,iep,1) = xeent
      endif

    ! If ICALKE = 1 the boundary conditions of turbulence at
    ! the inlet refer to both, a hydraulic diameter and a
    ! reference velocity given in usini1.f90.
    elseif (icalke(izone).eq.1) then

       dh(izone)     = 2.d0 * rinlet

    ! If ICALKE = 2 the boundary conditions of turbulence at
    ! the inlet refer to a turbulence intensity.
    elseif (icalke(izone).eq.2) then

      xintur(izone) = 0.05d0

    endif

!  endif

enddo

! ADMISSION

call getfbr('INTAKE', nlelt, lstelt)

madm = 0.d0
if (tt.le.22) then
  dp = (Pther - P0) - Padm
else
  dp   = (Pther - P0) - Padm!padm_act
endif

write(nfecra,*) "DP INTAKE: ", padm_act, dp

if (dp.gt.0) then
!  kadm = kadm*4
endif

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)
  iel  = ifabor(ifac)

  if (dp.gt.0.d0) then
    rho = crom(iel)
  else
    rho = radm ! brom(ifac)
  endif

  ! Flow rate : Attention pour une entree, le flux est oriente
  ! selon moins la normale.

  ! FIXME le debit est constant par zone, or je le calcule avec
  ! un rho dependant des faces. C'est le debit calcule sur la derniere
  ! face qui est conserve.

  madm = - sign(1.d0,dp) * s0 * sqrt(2.d0*rho/Kadm*abs(dp))
  vadm = - sqrt((2.d0*abs(dp))/(Kadm*rho))

  itypfb(ifac) = ientre !i_convective_inlet

!   Zone number (arbitrary number between 1 and n)
!  izone = 2

!   Allocation of the actual face to the zone
  izone = izfppp(ifac)

  ! Indicating the inlet as a fuel flow inlet
  ientox(izone) = 1

  ! Inlet Temperature in K
  ! tinoxy = t0 - 5.d0

  ! The incoming fuel flow refers to:
  ! a) a massflow rate   -> iqimp()  = 1
  iqimp(izone)  = 1
  qimp(izone)   = madm
  !write(nfecra,*) "madm: ", madm, izone, qimp(izone)

  ! b) an inlet velocity -> iqimp()  = 0

  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = -madm/(rho*s0)
  rcodcl(ifac,iw,1) = 0.d0

  uref2 = rcodcl(ifac,iu,1)**2  &
         +rcodcl(ifac,iv,1)**2  &
         +rcodcl(ifac,iw,1)**2
  uref2 = max(uref2,1.d-12)

  ! Boundary conditions of turbulence
  icalke(izone) = 1
  !
  !  - If ICALKE = 1 the boundary conditions of turbulence at
  !    the inlet refer to both, a hydraulic diameter and a
  !    reference velocity given in usini1.f90.
  !
  dh(izone)     = 0.15d0

  !  - If ICALKE = 2 the boundary conditions of turbulence at
  !    the inlet refer to a turbulence intensity.
  !
  xintur(izone) = 0.05d0

enddo




! EXTRACTION

call getfbr('EXHAUST', nlelt, lstelt)

mext = 0.d0
dp   = (Pther - P0) - Pext

write(nfecra,*) "DP EXHAUST: ", pext_act, dp

do ilelt = 1,nlelt

  ifac = lstelt(ilelt)
  iel = ifabor(ifac)

  if (dp.gt.0.d0) then
    rho = crom(iel)
  else
    rho = rext ! brom(ifac)
  endif

  ! Flow rate : Attention pour une entree, le flux est oriente
  ! selon moins la normale.

  ! FIXME le debit est constant par zone, or je le calcule avec
  ! un rho dependant des faces. C'est le debit calcule sur la derniere
  ! face qui est conserve.
  mext = - sign(1.d0,dp) * s0 * sqrt(2.d0*rho/Kext*abs(dp)) !Kext = 2.d0*rext*(S0/mext)**2 * abs(Ploc-Pext) !+ 50.d0
  vext = - sign(1.d0,dp) * sqrt((2.d0*abs(dp))/(Kext*rho))

  itypfb(ifac) = ientre !i_convective_inlet

!   Zone number (arbitrary number between 1 and n)
!  izone = 2

!   Allocation of the actual face to the zone
  izone = izfppp(ifac)

  ! Indicating the inlet as a fuel flow inlet
  ientox(izone) = 1

  ! Inlet Temperature in K
  ! tinoxy = t0 - 5.d0

  ! The incoming fuel flow refers to:
  ! a) a massflow rate   -> iqimp()  = 1
  iqimp(izone)  = 1
  qimp(izone)   = mext
  !write(nfecra,*) "mext: ", mext, izone, qimp(izone)

  ! b) an inlet velocity -> iqimp()  = 0

  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = mext/(rho*s0)
  rcodcl(ifac,iw,1) = 0.d0

  uref2 = rcodcl(ifac,iu,1)**2  &
         +rcodcl(ifac,iv,1)**2  &
         +rcodcl(ifac,iw,1)**2
  uref2 = max(uref2,1.d-12)

  ! Boundary conditions of turbulence
  icalke(izone) = 1
  !
  !  - If ICALKE = 1 the boundary conditions of turbulence at
  !    the inlet refer to both, a hydraulic diameter and a
  !    reference velocity given in usini1.f90.
  !
  dh(izone)     = 0.15d0

  !  - If ICALKE = 2 the boundary conditions of turbulence at
  !    the inlet refer to a turbulence intensity.
  !
  xintur(izone) = 0.05d0

enddo

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine cs_f_user_boundary_conditions
