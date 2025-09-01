/*============================================================================
 * User definition of physical properties.
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
#include <string.h>

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
 * \file cs_user_physical_properties.cpp
 *
 * \brief User definition of physical properties.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define physical properties.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_physical_properties
void
cs_user_physical_properties(cs_domain_t   *domain)
{
  const cs_lnum_t n_cells = domain->mesh->n_cells;

  /* Fluid properties */

  const cs_fluid_properties_t *fp = cs_glob_fluid_properties;

  /* Variable viscosity as a function of temperature
     ----------------------------------------------- */

  if (fp->ivivar) {

    cs_real_t *cpro_temp = CS_F_(t)->val;
    cs_real_t *cpro_viscl = CS_F_(mu)->val;
    cs_real_t *cpro_yf = cs_field_by_name_try("ym_fuel")->val;
    cs_real_t *cpro_yo = cs_field_by_name_try("ym_oxyd")->val;
    cs_real_t *cpro_yp = cs_field_by_name_try("ym_prod")->val;

    // ivslaw = 10 fit Hanbook for air
    // ivslaw = 11 Sutherland's law
    // ivslaw = 20 mixing law

    int ivslaw = 11;

    /* Fit Handbook for air */

    if (ivslaw == 10) {

      cs_real_t varam = 46.291 * 1.e-7;
      cs_real_t varbm = 0.52546 * 1.e-7;
      cs_real_t varcm = -1.9103e-4 * 1.e-7;
      cs_real_t vardm = 4.4842e-8 * 1.e-7;

      for (cs_lnum_t i = 0; i < n_cells; i++) {
        cs_real_t x_t = cpro_temp[i];
        cs_real_t x_t2 = cs_math_pow2(x_t);
        cs_real_t x_t3 = cs_math_pow3(x_t);
        cpro_viscl[i] = varam + x_t*varbm + varcm*x_t2 + vardm*x_t3;
      }

    }

    /* Sutherland's law */

    else if (ivslaw == 11) {

      cs_real_t tsuthe = 110.4;
      cs_real_t tkelvi = cs_physical_constants_celsius_to_kelvin;
      cs_real_t viscl0 = fp->viscl0;

      for (cs_lnum_t i = 0; i < n_cells; i++) {
        cs_real_t x_t = cpro_temp[i];
        cs_real_t tpkel = x_t/tkelvi;
        tpkel = pow(tpkel,3.0/2.0);
        cpro_viscl[i] = viscl0 *(tkelvi+tsuthe)/(x_t+tsuthe)*tpkel;
      }

    }

    /* Mixing law */

    else if (ivslaw == 20) {

      const int ngaze = cs_glob_combustion_gas_model->n_gas_el_comp;

      cs_real_t  viscoefk[4][CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];
      cs_real_t  visk[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];
      cs_real_t  yi[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];
      cs_real_t  yk[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];
      cs_real_t  xk[CS_COMBUSTION_GAS_MAX_ELEMENTARY_COMPONENTS];

#if 0
      /* c2h4 */
      viscoefk[0][0]= -3.583724530538e-007;
      viscoefk[1][0]=  3.868415673467e-008;
      viscoefk[2][0]= -1.193808567543e-011;
      viscoefk[3][0]=  2.113573044333e-015;
#endif

      /*--- c3h8 ---*/
      viscoefk[0][0]= -1.234349675612e-007;
      viscoefk[1][0]=  3.104728049493e-008;
      viscoefk[2][0]= -9.675403508033e-012;
      viscoefk[3][0]=  1.725116880193e-015;

      /*--- o2 ---*/
      viscoefk[0][1]=  5.488554654900e-006;
      viscoefk[1][1]=  5.699109311534e-008;
      viscoefk[2][1]= -1.793899193773e-011;
      viscoefk[3][1]=  3.451533345962e-015;

      /*--- co2 ---*/
      viscoefk[0][2]=  2.817125838535e-007;
      viscoefk[1][2]=  5.513314232221e-008;
      viscoefk[2][2]=  -1.749297867996e-011;
      viscoefk[3][2]=  3.166176962235e-015;

      /*-- h2o --*/
      viscoefk[0][3]= -2.754263332602e-006;
      viscoefk[1][3]=  5.396477445927e-008;
      viscoefk[2][3]= -9.502751904262e-012;
      viscoefk[3][3]=  9.232335790842e-016;

      /* n2 */
      viscoefk[0][4]= 5.181365609879e-006;
      viscoefk[1][4]= 4.841826612088e-008;
      viscoefk[2][4]= -1.490782620291e-011;
      viscoefk[3][4]=  2.849448433980e-015;

      for (cs_lnum_t i = 0; i < n_cells; i++) {

        cs_real_t sum1 = 0.;
        cs_real_t sum2 = 0.;

        cs_real_t temp = cpro_temp[i];

        yi[0] = cpro_yf[i];
        yi[1] = cpro_yo[i];
        yi[2] = cpro_yp[i];

        /* Calculation of the mole fraction of each species */
        cs_combustion_gas_yg2xye(yi, yk, xk);

        /* Calculation of the viscosity of each species */
        for (int j = 0; j < ngaze; j++) {
          visk[j] = 0.;
          for (int l = 0; l < 4; l++) {
            visk[j] += viscoefk[l][j] * pow(temp, l);
          }
        }

        /* Calculation of the viscosity of the mixture */
        for (int j = 0; j < ngaze; j++) {
          sum1 += visk[j]*xk[j];
          sum2 += xk[j]/visk[j];
        }

        cpro_viscl[i] = 0.5*(sum1 + 1.0/sum2);

      }

    }

  }

  /* Diffusivity as a function of viscosity for scalars
     -------------------------------------------------- */

  /* Mixture fraction */
  cs_field_t *fm = CS_F_(fm);

  const int kivisl = cs_field_key_id("diffusivity_id");
  const int ksigmas = cs_field_key_id("turbulent_schmidt");

  const int ifcvsl = cs_field_get_key_int(fm, kivisl);

  if (ifcvsl >= 0) {

    const int n_fields = cs_field_n_fields();
    const int keysca = cs_field_key_id("scalar_id");
    const int kscavr = cs_field_key_id("first_moment_id");

    for (int f_id = 0; f_id < n_fields; f_id++) {

      cs_field_t *f = cs_field_by_id(f_id);

      /* Here we only handle user or model scalar-type variables
         which are not fluctuations */

      int sc_id = -1;
      if (f->type & CS_FIELD_VARIABLE)
        sc_id = cs_field_get_key_int(f, keysca) - 1;
      if (sc_id < 0)
        continue;

      int variance_id = cs_field_get_key_int(f, kscavr);
      int diffusivity_id = cs_field_get_key_int(f, kivisl);

      if (variance_id > -1 || diffusivity_id < 0)
        continue;

      cs_real_t *cpro_viscl = CS_F_(mu)->val;
      cs_real_t *cpro_visls = cs_field_by_id(diffusivity_id)->val;

      const cs_real_t turb_schmidt = cs_field_get_key_double(f, ksigmas);

      /* Viscosity in kg/(m.s) at cell centers */

      for (cs_lnum_t c_id = 0; c_id < n_cells; c_id++) {
        cpro_visls[c_id] = cpro_viscl[c_id] / turb_schmidt;
      }

    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
