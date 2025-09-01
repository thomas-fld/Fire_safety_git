/*============================================================================
 * User definition of boundary conditions - C version for Code_Saturne 8.3
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_boundary_conditions.cpp
 *
 * \brief User functions for input of calculation parameters.
 *
 * See \ref parameters for examples.
 */
/*----------------------------------------------------------------------------*/

#define MAX_LINES 10000 /* number of maximal lines read in the files */

static cs_lnum_t n_qvent = 0;                /* number of points from vent data */
static cs_lnum_t n_mpyro = 0;                /* number of points from mass loss data */
static cs_real_t Qvent[MAX_LINES][3];  /* time series for pressure and vent data */
static cs_real_t mpyro[MAX_LINES][2];  /* time series for mass loss rate (pyrolysis) */

/*----------------------------------------------------------------------------*/
static void
_load_time_series_data(void)
{
  FILE *f = fopen("PRS_SI_D6_MLR.txt", "r");
  if (f == NULL) {
    bft_printf("Erreur d'ouverture du fichier PRS_SI_D6_MLR.txt\n");
    cs_exit(1);
  }

  char line[256];
  /* read all lines, flexible format: time value  (two columns) */
  n_mpyro = 0;
  while (n_mpyro < MAX_LINES && fgets(line, sizeof(line), f) != NULL) {
    cs_real_t t, v;
    if (sscanf(line, "%lf %lf", &t, &v) == 2) {
      mpyro[n_mpyro][0] = (cs_real_t)t;
      mpyro[n_mpyro][1] = (cs_real_t)v;
      n_mpyro++;
    }
    /* otherwise ignore line (header or malformed) */
  }

  fclose(f);
}

/*----------------------------------------------------------------------------*/
static cs_real_t
_interpolate(const cs_real_t data[][2], cs_lnum_t n, cs_real_t time)
{
  /* simple linear interpolation (if time <= first entry returns first,
     if time >= last entry returns last). Protect n>=1. */
  if (n <= 0) return (cs_real_t)0.0;
  if (n == 1) return data[0][1];

  if (time <= data[0][0]) return data[0][1];
  if (time >= data[n-1][0]) return data[n-1][1];

  for (cs_lnum_t i = 0; i < n - 1; i++) {
    cs_real_t t1 = data[i][0];
    cs_real_t t2 = data[i+1][0];
    if (t1 <= time && time <= t2) {
      cs_real_t v1 = data[i][1];
      cs_real_t v2 = data[i+1][1];
      return v1 + (v2 - v1) * (time - t1) / (t2 - t1);
    }
  }
  /* fallback */
  return data[n-1][1];
}


/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* Static globals used by setup and runtime (initialized in setup) */
static cs_real_t Ploc = 0.0;
static cs_real_t S0   = 0.0;
static cs_real_t P0   = 0.0;
static cs_real_t t0_glob = 0.0;
static cs_real_t temp_oxyd = 0.0;
static cs_real_t M_glob = 0.0;
static cs_real_t radm = 0.0;
static cs_real_t rext = 0.0;
static cs_real_t Padm = 0.0;
static cs_real_t Pext = 0.0;
static cs_real_t Qadm = 0.0;
static cs_real_t Qext = 0.0;
static cs_real_t madm_setup = 0.0;
static cs_real_t mext_setup = 0.0;
static cs_real_t Kadm = 0.0;
static cs_real_t Kext = 0.0;

/*----------------------------------------------------------------------------*/
void
cs_user_boundary_conditions_setup(cs_domain_t *domain)
{
  /* load optional data files if present (pyro, vent pressures) */
  _load_time_series_data();

  const cs_fluid_properties_t *fp = cs_glob_fluid_properties;
  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  /* --- constants / reference values --- */
  Ploc = -6.23; /* Pa (local reference) */
  S0   = 0.18;  /* m2 (opening surface) */
  P0   = fp->p0; /* reference pressure (Pa) */
  t0_glob = fp->t0; /* reference temperature (K) */
  temp_oxyd = cm->tinoxy;
  M_glob = cm->wmolg[2]; /* molar mass (kg/mol) */

  /* compute densities used as fallbacks */
  radm = (P0 + Ploc) / (cs_physical_constants_r * (temp_oxyd / M_glob));
  rext = (P0 + Ploc) / (cs_physical_constants_r * (t0_glob / M_glob));

  /* set pressures / flows used as reference */
  Padm = 509.73;  /* Pa */
  Qadm  = 560.0;  /* m3/h */
  madm_setup = Qadm * radm / 3600.0; /* kg/s (approx) */

  Pext = -708.7;  /* Pa */
  Qext = 560.0;   /* m3/h */
  mext_setup = Qext * rext / 3600.0;

  /* coefficients (dimensionless in your formula) */
  Kadm = 2.0 * radm * (S0 / madm_setup) * (S0 / madm_setup) * fabs(Ploc - Padm);
  Kext = 2.0 * rext * (S0 / mext_setup) * (S0 / mext_setup) * fabs(Ploc - Pext);

  /* debug prints */
  bft_printf("BC setup: P0=%g Ploc=%g radm=%g rext=%g\n",
             (cs_real_t)P0, (cs_real_t)Ploc, (cs_real_t)radm, (cs_real_t)rext);
  bft_printf("BC setup: Padm=%g Pext=%g Qadm=%g Qext=%g\n",
             (cs_real_t)Padm, (cs_real_t)Pext, (cs_real_t)Qadm, (cs_real_t)Qext);
  bft_printf("BC setup: madm_setup=%g mext_setup=%g Kadm=%g Kext=%g\n",
             (cs_real_t)madm_setup, (cs_real_t)mext_setup, (cs_real_t)Kadm, (cs_real_t)Kext);

  /* === kPyrolysis (Fire_pan) initialization === */
  const cs_zone_t *z_fire = cs_boundary_zone_by_name("Fire_pan");
  cs_glob_combustion_gas_model->tinfue = cs_glob_fluid_properties->t0;

  if (z_fire != NULL) {
    cs_boundary_conditions_open_set_mass_flow_rate_by_value(z_fire, 0.0);
    cs_boundary_conditions_inlet_set_turbulence_hyd_diam(z_fire, sqrt(4.0 * 0.4 / M_PI));
  }

  /* === Inlet_All: initialize massflow to madm_setup (will be updated each step) === */
  {
    const cs_zone_t *z_inlet_All = cs_boundary_zone_by_name("Inlet_All");
    cs_real_t d_vent = sqrt(S0 / M_PI);
    if (z_inlet_All != NULL) {
      cs_boundary_conditions_open_set_mass_flow_rate_by_value(z_inlet_All, madm_setup);
      cs_boundary_conditions_inlet_set_turbulence_hyd_diam(z_inlet_All, d_vent);
    }
  }

  /* === Outlet_All: initialize massflow to mext_setup (will be updated each step) === */
  {
    const cs_zone_t *z_outlet_All = cs_boundary_zone_by_name("Outlet_All");
    if (z_outlet_All != NULL) {
      cs_boundary_conditions_open_set_mass_flow_rate_by_value(z_outlet_All, mext_setup);
    }
  }
}

/*----------------------------------------------------------------------------*/
void
cs_user_boundary_conditions(cs_domain_t  *domain,
                            cs_lnum_t           bc_type[])
{
  /* small epsilon */
  const cs_real_t eps = 1e-12;

  /* --- retrieve useful pointers --- */
  const cs_lnum_t *b_face_cells = domain->mesh->b_face_cells;
  const cs_real_3_t *cdgfbo = (const cs_real_3_t *)(domain->mesh_quantities->b_face_cog);
  const cs_real_t *gxyz = cs_glob_physical_constants->gravity;

  const cs_fluid_properties_t *fp = cs_glob_fluid_properties;
  cs_combustion_gas_model_t *cm = cs_glob_combustion_gas_model;

  /* time shifted as you used earlier */
  cs_real_t ttcabs = cs_glob_time_step->t_cur;
  cs_real_t tstart = 100.0;
  cs_real_t tt = ttcabs - tstart;

  /* load pyro mass loss */
  cs_real_t dmdt_pyro = _interpolate(mpyro, n_mpyro, tt);
  bft_printf("DBG: timestep tcur=%g tt=%g dmdt_pyro=%g\n", (cs_real_t)ttcabs, (cs_real_t)tt, (cs_real_t)dmdt_pyro);

  /* update Fire_pan mass flow */
  const cs_zone_t *z_fire = cs_boundary_zone_by_name("Fire_pan");
  if (z_fire != NULL) {
    cs_boundary_conditions_open_set_mass_flow_rate_by_value(z_fire, dmdt_pyro);
  }

  /* get boundary mass flux field (for per-face checks) */
  const cs_lnum_t kbmasf = cs_field_key_id("boundary_mass_flux_id");
  const cs_lnum_t iflmab = cs_field_get_key_int(CS_F_(vel), kbmasf);
  const cs_real_t *bmasfl = cs_field_by_id(iflmab)->val;

  /* set oxidizer composition / enthalpy if needed (unchanged) */
  cs_real_t coefg[CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES];
  {
    cs_lnum_t i;
    for (i = 0; i < CS_COMBUSTION_GAS_MAX_GLOBAL_SPECIES; ++i) coefg[i] = 0.0;
    coefg[1] = 1.0;
  }
  cm->tinoxy = fp->t0;
  cm->hinoxy = cs_gas_combustion_t_to_h(coefg, cm->tinoxy);

  /* --- Zones references --- */
  const cs_zone_t *z_in = cs_boundary_zone_by_name("Inlet_All");
  const cs_zone_t *z_out = cs_boundary_zone_by_name("Outlet_All");

  /* -----------------------------------------------------
     Helper to compute global average of cell rho for a zone:
     - we compute local sum over faces: sum_local (cs_real_t), count_local (int)
     - then use cs_parall_sum to reduce to global
     - return average (or fallback)
     -----------------------------------------------------*/
  {
    /* INLET processing */
    if (z_in != NULL) {
      cs_real_t sum_rho_local = 0.0;
      int count_local = 0;
      cs_lnum_t nf_local = z_in->n_elts;
      /* accumulate local rho from neighbor cells (interior cells) */
      for (cs_lnum_t ii = 0; ii < nf_local; ++ii) {
        cs_lnum_t face_id = z_in->elt_ids[ii];
        bc_type[face_id] = CS_CONVECTIVE_INLET;
        /* safe check */
        if (face_id < 0) continue;
        cs_lnum_t cell_id = b_face_cells[face_id];
        if (cell_id >= 0) {
          sum_rho_local += (cs_real_t) CS_F_(rho)->val[cell_id];
          count_local += 1;
        }
      }

      /* global reduction */
      cs_real_t sum_rho_global = sum_rho_local;
      cs_lnum_t count_global = count_local;
      if (cs_glob_n_ranks > 1) {
        cs_parall_sum(1, CS_REAL_TYPE, &sum_rho_global);
        cs_parall_sum(1, CS_INT_TYPE, &count_global);
      }

      /* compute dp and choose rho for Bernoulli */
      cs_real_t dp_adm = (fp->pther - P0) - Padm; /* dp = (pther-P0) - Padm_act (padm_act loaded from file) */
      cs_real_t rho_used = radm; /* fallback */
      if (count_global > 0) rho_used = (cs_real_t)(sum_rho_global / (cs_real_t)count_global);

      /* convention: if dp_adm > 0 -> domain pushes out -> use local rho (rho_used),
         else dp_adm <=0 -> exterior pushes in -> use imposed radm */
      if (dp_adm <= 0.0) {
        rho_used = radm;
      }

      bft_printf("INLET DBG: dp_adm=%g count_global=%d rho_used(avg)=%g radm=%g Kadm=%g\n",
      (cs_real_t)dp_adm, count_global, (cs_real_t)rho_used, (cs_real_t)radm, (cs_real_t)Kadm);

      /* protect against invalid values */
      if (Kadm <= eps || rho_used <= eps) {
        bft_printf("WARN INLET: skip Bernoulli (Kadm=%g rho_used=%g)\n", (cs_real_t)Kadm, (cs_real_t)rho_used);
      } else {
        cs_real_t vadm = (cs_real_t) sqrt( 2.0 * fabs(dp_adm) / (Kadm * rho_used) );
        cs_real_t signe = (dp_adm > 0.0) ? -1.0 : 1.0; /* -1 if flow exits domain, +1 if enters */
        cs_real_t madm_local = signe * rho_used * vadm * S0;
        bft_printf("INLET RESULT: vadm=%g madm_local=%g (signe=%g)\n",
                   (cs_real_t)vadm, (cs_real_t)madm_local, (cs_real_t)signe);

        /* set mass flow for the whole zone */
        cs_boundary_conditions_open_set_mass_flow_rate_by_value(z_in, madm_local);
        /* optionally: set velocity normal by value instead (uncomment if desired)
           cs_boundary_conditions_open_set_velocity_by_normal_value(z_in, signe * vadm);
         */
      }
    } else {
      bft_printf("INLET_DBG: zone Inlet_All not found on this process.\n");
    }
  }

  /* ----------------------------------------------------- */
  /* OUTLET processing */
  {
      if (z_in != NULL) {
      cs_real_t sum_rho_local = 0.0;
      int count_local = 0;
      cs_lnum_t nf_local = z_out->n_elts;
      /* accumulate local rho from neighbor cells (interior cells) */
      for (cs_lnum_t ii = 0; ii < nf_local; ++ii) {
        cs_lnum_t face_id = z_out->elt_ids[ii];
        bc_type[face_id] = CS_CONVECTIVE_INLET;
        /* safe check */
        if (face_id < 0) continue;
        cs_lnum_t cell_id = b_face_cells[face_id];
        if (cell_id >= 0) {
          sum_rho_local += (cs_real_t) CS_F_(rho)->val[cell_id];
          count_local += 1;
        }
      }

      cs_real_t sum_rho_global = sum_rho_local;
      cs_lnum_t count_global = count_local;
      if (cs_glob_n_ranks > 1) {
        cs_parall_sum(1, CS_GNUM_TYPE, &sum_rho_global);
        cs_parall_sum(1, CS_INT_TYPE, &count_global);
      }

      cs_real_t dp_ext_local = (fp->pther - P0) - Pext;
      cs_real_t rho_used = rext;
      if (count_global > 0) rho_used = (cs_real_t)(sum_rho_global / (cs_real_t)count_global);

      /* convention: if dp_ext > 0 -> domain pushes out -> use local rho (rho_used),
         else dp_ext <=0 -> exterior pushes in -> use imposed rext */
      if (dp_ext_local <= 0.0) {
        rho_used = rext;
      }

      bft_printf("OUTLET DBG: dp_ext=%g count_global=%d rho_used(avg)=%g rext=%g Kext=%g\n",
                 (cs_real_t)dp_ext_local, count_global, (cs_real_t)rho_used, (cs_real_t)rext, (cs_real_t)Kext);

      if (Kext <= eps || rho_used <= eps) {
        bft_printf("WARN OUTLET: skip Bernoulli (Kext=%g rho_used=%g)\n", (cs_real_t)Kext, (cs_real_t)rho_used);
      } else {
        cs_real_t vext = (cs_real_t) sqrt( 2.0 * fabs(dp_ext_local) / (Kext * rho_used) );
        cs_real_t signe_out = (dp_ext_local > 0.0) ? -1.0 : 1.0;
        cs_real_t mext_local = signe_out * rho_used * vext * S0;
        bft_printf("OUTLET RESULT: vext=%g mext_local=%g (signe=%g)\n",
                   (cs_real_t)vext, (cs_real_t)mext_local, (cs_real_t)signe_out);

        cs_boundary_conditions_open_set_mass_flow_rate_by_value(z_out, mext_local);
        /* optionally: set velocity normal by value instead
           cs_boundary_conditions_open_set_velocity_by_normal_value(z_out, signe_out * vext);
         */
      }
    } else {
      bft_printf("OUTLET_DBG: zone Outlet_All not found on this process.\n");
    }
  }

  /* --- Following you can keep per-face BCs (pressure bc, turbulence, enthalpy) ---
     If you need to reapply per-face operations (like in your commented code),
     do it below. For example: apply hydrostatic pressure to Outlet_All faces,
     set turbulence BC on faces where mass flux enters, set enthalpy for oxidizer, etc.
  */

  /* Example: reapply outlet per-face settings (hydrostatic pressure, turbulence, enthalpy) */
  if (z_out != NULL) {
    const cs_lnum_t *elt_ids = z_out->elt_ids;
    cs_lnum_t nf = z_out->n_elts;
    for (cs_lnum_t i = 0; i < nf; ++i) {
      cs_lnum_t face_id = elt_ids[i];

      /* hydrostatic reference pressure */
      cs_real_t pimp = fp->p0 + fp->ro0 * cs_math_3_distance_dot_product(fp->xyzp0, cdgfbo[face_id], gxyz);
      // CS_F_(p)->bc_coeffs->icodcl[face_id] = 1;
      // CS_F_(p)->bc_coeffs->rcodcl1[face_id] = pimp;

      /* if mass flux entering at this face, skip turbulence BC */
      if (bmasfl[face_id] > 0.0) continue;

      /* otherwise set turbulence BC based on local cell velocity */
      cs_lnum_t c_id = b_face_cells[face_id];
      const cs_real_3_t *cvar_vel = (const cs_real_3_t *)CS_F_(vel)->val;
      cs_real_t uref2 = cs_math_3_square_norm(cvar_vel[c_id]);
      uref2 = cs_math_fmax(uref2, cs_math_epzero);

      cs_turbulence_bc_inlet_turb_intensity(face_id, uref2, 0.02, 3.0);

      CS_F_(fm)->bc_coeffs->rcodcl1[face_id]   = 0.0;
      CS_F_(fp2m)->bc_coeffs->rcodcl1[face_id] = 0.0;

      /* enthalpy for incoming oxidizer */
      if (cs_glob_physical_model_flag[CS_COMBUSTION_3PT] == 1)
        CS_F_(h)->bc_coeffs->rcodcl1[face_id] = cm->hinoxy;

      if (cm->isoot == 1) {
        CS_F_(fsm)->bc_coeffs->rcodcl1[face_id] = 0.0;
        CS_F_(npm)->bc_coeffs->rcodcl1[face_id] = 0.0;
      }

      /* set rho_b (boundary density) from ideal gas (for output / debug) */
      CS_F_(rho_b)->val[face_id] = fp->p0 / (cs_physical_constants_r * cm->tinoxy / cm->wmolg[2]);
    }
  }

  /* end of cs_user_boundary_conditions */
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
