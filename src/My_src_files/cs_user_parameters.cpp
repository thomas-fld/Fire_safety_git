/*============================================================================
 * cs_user_parameters.cpp - User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/
/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

#include <stdio.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_parameters-base.cpp
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select physical model options, including user fields.
 *
 * This function is called at the earliest stages of the data setup,
 * so field ids are not available yet.
 */
/*----------------------------------------------------------------------------*/


/*============================================================================
 * rief User subroutine for selection of specific physics module (usppmo(ixmlpu)
 *============================================================================*/

void
cs_user_model(void)
{
  /* Activate gas combustion model */
  cs_combustion_gas_model_t *gcm
    = cs_combustion_gas_set_model(CS_COMBUSTION_3PT_PERMEATIC);

  /* Soot model for gas combustion
   *  if = -1   module not activated
   *  if =  0   constant soot yield
   *  if =  1   2 equations model of Moss et al. */
  gcm->isoot = 0;
  gcm->xsoot = 0.1;    /* (only if isoot = 0 and soot yield is not
                          defined in the thermochemistry data file) */
  gcm->rosoot = 2000.; /* kg/m3 */

  /*----------------------------------cs_user_parameters-volume_mass_injection--------------------------------------*/
  /*! [user_property_addition] */

  /* Example: add a user property defined on boundary faces.
   *
   * parameters for cs_parameters_add_property():
   *   name        <-- name of property and associated field
   *   dim         <-- property dimension
   *   location_id <-- id of associated mesh location, which must be one of:
   *                     CS_MESH_LOCATION_CELLS
   *                     CS_MESH_LOCATION_INTERIOR_FACES
   *                     CS_MESH_LOCATION_BOUNDARY_FACES
   *                     CS_MESH_LOCATION_VERTICES
   */

  cs_parameters_add_property("wall_incident_flux",
                             1,
                             CS_MESH_LOCATION_BOUNDARY_FACES);
  cs_parameters_add_property("wall_total_flux",
                             1,
                             CS_MESH_LOCATION_BOUNDARY_FACES);
  cs_parameters_add_property("wall_convective_flux",
                             1,
                             CS_MESH_LOCATION_BOUNDARY_FACES);
  cs_parameters_add_property("wall_temp",
                             1,
                             CS_MESH_LOCATION_BOUNDARY_FACES);
}

void
cs_user_parameters(cs_domain_t *domain)
{
  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();

/* irovar, ivivar, icp: constant or variable density,
                              viscosity/diffusivity, and specific heat

     When a specific physics module is active
       (coal, combustion, electric arcs, compressible: see cs_user_model)
       we MUST NOT set variables 'irovar', 'ivivar', and 'icp' here, as
       they are defined automatically.
     Nonetheless, for the compressible case, ivivar may be modified
       in the uscfx2 user subroutine.

     When no specific physics module is active, we may specify if the
       density, specific heat, and the molecular viscosity
       are constant (irovar=0, ivivar=0, icp=-1), which is the default
       or variable (irovar=1, ivivar=1, icp=0)

     For those properties we choose as variable, the corresponding law
       must be defined in cs_user_physical_properties.
       if they are constant, they take values ro0, viscl0, and cp0.
    */

    fp->irovar = 1;
    fp->ivivar = 1;
    fp->icp = -1;


/* Account the thermodynamical pressure variation in time
       (only if cs_glob_velocity_pressure_param->idilat = 3)

      By default:
      ----------
      - The thermodynamic pressure (pther) is initialized with p0 = p_atmos.
      - The maximum thermodynamic pressure (pthermax) is initialized with -1
        (no maximum by default, this term is used to model a venting effect when
        a positive value is given by the user).
      - A global leak can be set through a leakage surface sleak with a head
       loss kleak of 2.9 (Idelcick) */

    fp->ipthrm = 1;

    fp->pthermax= -1.;

    fp->sleak = 0.0025;
    fp->kleak = 2.9;


/* Diffusive scheme */
    cs_equation_t *vbn = cs_equation_by_name("velocity");
    if (vbn != nullptr) {
      vbn->param->imvisf =1;
    }

/*
      ro0        : density in kg/m3
      viscl0     : dynamic viscosity in kg/(m s)
      cp0        : specific heat in J/(Kelvin kg)
      t0         : reference temperature in Kelvin
      p0         : total reference pressure in Pascal
                   the calculation is based on a
                   reduced pressure P*=Ptot-ro0*g.(x-xref)
                   (except in compressible case)
      xyzp0(3)   : coordinates of the reference point for
                   the total pressure (where it is equal to p0)

      In general, it is not necessary to furnish a reference point xyz0.
      If there are outlets, the code will take the center of the
      reference outlet face.
      On the other hand, if we plan to explicitly fix Dirichlet conditions
      for pressure, it is better to indicate to which reference the
      values relate (for a better resolution of reduced pressure).

      Other properties are given by default in all cases.

      Nonetheless, we may note that:

/*With gas combustion:
      --------------------
      ro0 is not useful (it is automatically recalculated by the
          law of ideal gases from t0 and p0).
      viscl0 is indispensable: it is the molecular dynamic viscosity,
          assumed constant for the fluid.
      cp0 is indispensable: it is the heat capacity, assumed constant,
          (modelization of source terms involving a local Nusselt in
          the Lagrangian module, reference value allowing the
          calculation of a radiative
          (temperature, exchange coefficient) couple).
      t0  is indispensible and must be in Kelvin (> 0).
      p0  is indispensable and must be in Pascal (> 0).*/

    fp->viscl0 = 1.83337e-5;
    fp->cp0    = 1017.24;
    fp->t0 = 20. + 273.15;
    fp->p0 = 1.01325e5;

    fp->pther = fp->p0 - 6.23;

    // ---- Laminar viscositysued to be defined i cs_user_finalize_setup ----
    fp->ivivar = 1;

    {
    const int kivisl = cs_field_key_id("diffusivity_id");
    const int k_scal = cs_field_key_id("scalar_id");
    const int kscavr = cs_field_key_id("first_moment_id");

    const int n_fields = cs_field_n_fields();

    for (int f_id = 0; f_id < n_fields; f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);

      if (!(f->type & (CS_FIELD_VARIABLE | CS_FIELD_USER)))
        continue;

      int s_num = cs_field_get_key_int(f, k_scal);
      int iscavr = cs_field_get_key_int(f, kscavr);

      if (s_num > 0 && iscavr <= 0)
        cs_field_set_key_int(f, kivisl, 0);
      }
    }

  /* Postprocessing-related fields
     ============================= */

  /* Example: enforce existence of 'yplus', 'tplus' and 'tstar' fields, so that
              yplus may be saved, or a local Nusselt number may be computed using
              the post_boundary_nusselt subroutine.
              When postprocessing of these quantities is activated, those fields
              are present, but if we need to compute them in the
              cs_user_extra_operations user subroutine without postprocessing
              them, forcing the definition of these fields to save the
              values computed for the boundary layer is necessary. */

  /*! [param_force_yplus] */
  {
    const int ityloc = CS_MESH_LOCATION_CELLS; // or CS_MESH_LOCATION_BOUNDARY_FACES if needed
    const int idim1 = 1;

    const int keyvis = cs_field_key_id("post_vis");
    const int keylog = cs_field_key_id("log");

    const int vis_mode = CS_POST_ON_LOCATION | CS_POST_MONITOR;

    // Liste des champs à créer
    const char *field_names[] = { "Xm_O2", "Xm_CO2", "Cm_Soot" };

    for (int i = 0; i < 3; ++i) {
      cs_field_t *f = cs_field_by_name_try(field_names[i]);

      if (f == nullptr) {
        cs_field_t *new_field = cs_field_create(field_names[i],
                                            CS_FIELD_INTENSIVE | CS_FIELD_PROPERTY,
                                            ityloc,
                                            idim1,
                                            false); // inoprv = false
        cs_field_set_key_int(new_field, keyvis, vis_mode);
        cs_field_set_key_int(new_field, keylog, 1);
  }
}

  }

/*----------------------------------cs_user_parameters- volume_mass_injection--------------------------------------*/

(void)domain;
cs_field_t *f_o2 = cs_field_by_name_try("ym_oxyd");

cs_field_set_key_double(f_o2,
                        cs_field_key_id("min_scalar_clipping"),
                        0.8);
cs_field_set_key_double(f_o2,
                        cs_field_key_id("max_scalar_clipping"),
                        1.);
}

/*----------------------------------------------------------------------------*/

void
cs_parameters_output_complete(void)
{
  const int key_log = cs_field_key_id("log");
  const int key_vis = cs_field_key_id("post_vis");

  const int post_flags = CS_POST_ON_LOCATION | CS_POST_MONITOR;

  // Température
  cs_field_t *f_temp = cs_field_by_name_try("temperature");

  if (f_temp != nullptr) {
    cs_field_set_key_int(f_temp, key_vis, post_flags);
    cs_field_set_key_int(f_temp, key_log, 1);
  }

  // Champs ym_* (fractions massiques)
  const char *ym_fields[] = { "ym_fuel", "ym_oxyd", "ym_prod" };

  for (int ii = 0; ii < 3; ++ii) {
    cs_field_t *f_ym = cs_field_by_name_try(ym_fields[ii]);
    if (f_ym != nullptr) {
      cs_field_set_key_int(f_ym, key_vis, post_flags);
      cs_field_set_key_int(f_ym, key_log, 1);
    }
  }
}

/*----------------------------------------------------------------------------*/

void
cs_user_finalize_setup(cs_domain_t     *domain)
{

}

END_C_DECLS
