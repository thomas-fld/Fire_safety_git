/*============================================================================
 * Radiation solver operations.
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
#include <string.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*! \file cs_user_radiative_transfer_bcs.cpp */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of radiative transfer boundary conditions.
 *
 * See \ref cs_user_radiative_transfer for examples.
 *
 * \warning the temperature unit here is the Kelvin
 *
 * For each boundary face face_id, a specific output (logging and
 * postprocessing) class id may be assigned. This allows realizing balance
 * sheets by treating them separately for each zone. By default, the
 * output class id is set to the general (input) zone id associated to a face.
 *
 * To access output class ids (both for reading and modifying), use the
 * \ref cs_boundary_zone_face_class_id function.
 * The zone id values are arbitrarily chosen by the user, but must be
 * positive integers; very high numbers may also lead to higher memory
 * consumption.
 *
 * \section cs_user_radiative_transfer_bcs_wall  Wall characteristics
 *
 * The following face characteristics must be set:
 *  - isothp(face_id) boundary face type
 *    * CS_BOUNDARY_RAD_WALL_GRAY:
 *      Gray wall with temperature based on fluid BCs
 *    * CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T:
 *      Gray wall with fixed outside temperature
 *    * CS_BOUNDARY_RAD_WALL_REFL_EXTERIOR_T:
 *      Reflecting wall with fixed outside temperature
 *    * CS_BOUNDARY_RAD_WALL_GRAY_COND_FLUX:
 *      Gray wall with fixed conduction flux
 *    * CS_BOUNDARY_RAD_WALL_REFL_COND_FLUX:
 *      Reflecting wall with fixed conduction flux
 *
 * Depending on the value of isothp, other values may also need to be set:
 *  - rcodcl = conduction flux
 *  - epsp   = emissivity
 *  - xlamp  = conductivity (W/m/K)
 *  - epap   = thickness (m)
 *  - textp  = outside temperature (K)
 *
 * \param[in, out]  domain        pointer to a cs_domain_t structure
 * \param[in]       bc_type       boundary face types
 * \param[in]       isothp        boundary face type for radiative transfer
 * \param[out]      tmin          min allowed value of the wall temperature
 * \param[out]      tmax          max allowed value of the wall temperature
 * \param[in]       tx            relaxation coefficient (0 < tx < 1)
 * \param[in]       dt            time step (per cell)
 * \param[in]       thwall        inside current wall temperature (K)
 * \param[in]       qincid        radiative incident flux  (W/m2)
 * \param[in]       hfcnvp        convective exchange coefficient (W/m2/K)
 * \param[in]       flcnvp        convective flux (W/m2)
 * \param[out]      xlamp         conductivity (W/m/K)
 * \param[out]      epap          thickness (m)
 * \param[out]      epsp          emissivity (>0)
 * \param[out]      textp         outside temperature (K)
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_radiative_transfer_bcs
void
cs_user_radiative_transfer_bcs(cs_domain_t      *domain,
                               const int         bc_type[],
                               int               isothp[],
                               cs_real_t        *tmin,
                               cs_real_t        *tmax,
                               cs_real_t        *tx,
                               const cs_real_t   dt[],
                               const cs_real_t   thwall[],
                               const cs_real_t   qincid[],
                               cs_real_t         hfcnvp[],
                               cs_real_t         flcnvp[],
                               cs_real_t         xlamp[],
                               cs_real_t         epap[],
                               cs_real_t         epsp[],
                               cs_real_t         textp[])
{
  cs_real_t tkelvi = 273.15;
  cs_lnum_t n_b_faces = cs_glob_mesh->n_b_faces;
  int *izfrdp = cs_boundary_zone_face_class_id();

  /* Allocate a temporary array for boundary faces selection */

  cs_lnum_t nlelt;
  cs_lnum_t  *lstelt;
  BFT_MALLOC(lstelt, n_b_faces, cs_lnum_t);

  cs_field_t *fth;

  switch (cs_glob_thermal_model->itherm) {
  case 1:
    fth = CS_F_(t);
    break;
  case 2:
    fth = CS_F_(h);
    break;
  default:
    fth = nullptr;
  }

  const cs_lnum_t ivart
    = cs_field_get_key_int(fth, cs_field_key_id("variable_id")) - 1;

  /* Min and max values for the wall temperatures (clipping otherwise)
   * TMIN and TMAX are given in Kelvin. */

  *tmin = cs_glob_fluid_properties->t0;
  *tmax = 3000;

// #if 0

  const cs_real_3_t *xyz_fbo
    = (const cs_real_3_t *)(cs_glob_mesh_quantities->b_face_cog);

  /* Zone definitions */
  /*------------------*/

  cs_selector_get_b_face_list("Walls_All", &nlelt, lstelt);

  for (cs_lnum_t ilelt = 0; ilelt < nlelt; ilelt++) {

    cs_lnum_t face_id = lstelt[ilelt];

    if (bc_type[face_id] == CS_SMOOTHWALL) {
      /* zone number */
      izfrdp[face_id] = 56;

      /* Type: gray or black wall with fixed outside temperature */
      isothp[face_id] = CS_BOUNDARY_RAD_WALL_GRAY_EXTERIOR_T;

      /* Emissivity */
      epsp[face_id] = 0.7;

      // Ceiling : 50 mm of thermipan (140 kg/m3, 840 J/kg/K, 0.102 W/m/K)
      //        + 300 mm of concrete (2430 kg/m3, 736 J/kg/K, 1.500 W/m/K)
      if (xyz_fbo[face_id][2] >= 3.94) {

        //Conductivity (W/m/K)
        xlamp[face_id] = 0.5;

        // Thickness (m)
        epap[face_id] = 0.35;

        /* Emissivity */
        epsp[face_id] = 0.95;

      }

      // Floor and walls : 300 mm of concrete
      // (2430 kg/m3, 736 J/kg/K, 1.500 W/m/K)
      else if (xyz_fbo[face_id][2] <= 0.02) {

        //Conductivity (W/m/K)
        xlamp[face_id] = 1.5;

        // Thickness (m)
        epap[face_id] = 0.3;
      }

      else {

        //Conductivity (W/m/K)
        xlamp[face_id] = 1.5;

        // Thickness (m)
        epap[face_id] = 0.3;

      }

      // Fixed outside temperature
      textp[face_id] = cs_glob_fluid_properties->t0;
    }

  }

// #endif

  /* Same selection criterion as that used for the 1D wall (uspt1d.f90) */

  /* FIXME je ne connais pas encore le type de face (entree ou paroi) sur ces
     couleurs, determine dans cs_user_boundary_conditions
     cette face sera consideree comme une paroi pour le scalaire enthalpie
     (raycli) mais si la face est donnee comme autre chose dans
     cs_user_boundary_conditions, il y aura incoherence entre le scalaire
     enthalpie et les autres variables */

  /*
  
  cs_selector_get_b_face_list("WALL or CEILING", &nlelt, lstelt);

  for (cs_lnum_t ilelt = 0; ilelt < nlelt; ilelt++) {

    cs_lnum_t face_id = lstelt[ilelt];

    if (bc_type[face_id] == CS_SMOOTHWALL) {
      // zone number 
      izfrdp[face_id] = 56;

      // Type: gray or black wall with fixed temperature (1D law) 
      isothp[face_id] = CS_BOUNDARY_RAD_WALL_GRAY_1D_T;

      // Emissivity 
      epsp[face_id] = 0.9;

    }

  }

  */

  /* Other output zone ids */

  for (cs_lnum_t face_id = 0; face_id < cs_glob_mesh->n_b_faces; face_id++) {

    if (bc_type[face_id] == CS_OUTLET)
      izfrdp[face_id] = 60;
    else if (bc_type[face_id] == CS_FREE_INLET)
      izfrdp[face_id] = 61;
    else if (bc_type[face_id] == CS_INLET)
      izfrdp[face_id] = 62;
    else if (bc_type[face_id] == CS_CONVECTIVE_INLET)
      izfrdp[face_id] = 63;
    else if (bc_type[face_id] == CS_SYMMETRY)
      izfrdp[face_id] = 64;

  }

  /* Deallocate the temporary array */
  BFT_FREE(lstelt);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
