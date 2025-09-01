/*============================================================================
 * Define postprocessing output.
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

#include "stdlib.h"
#include "string.h"

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_geom.h"
#include "cs_interpolate.h"
#include "cs_mesh.h"
#include "cs_selector.h"
#include "cs_parall.h"
#include "cs_post.h"
#include "cs_post_util.h"
#include "cs_probe.h"
#include "cs_time_plot.h"

#include "cs_field_pointer.h"
#include "cs_notebook.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_turbulence_model.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Local (user defined) function definitions
 *============================================================================*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define post-processing writers.
 *
 * The default output format and frequency may be configured, and additional
 * post-processing writers allowing outputs in different formats or with
 * different format options and output frequency than the main writer may
 * be defined.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_postprocess_writers
void
cs_user_postprocess_writers(void)
{
  /* Default writer time dependency */

  fvm_writer_time_dep_t   time_dep = FVM_WRITER_FIXED_MESH;

  /* Default time step or physical time based output frequencies */

  bool       output_at_start = true;
  bool       output_at_end = true;
  int        frequency_n = -1;
  double     frequency_t =  20.0;

  /* Default output format and options */

  const char format_name[] = "EnSight Gold";
  const char format_options[] = "";

  /* Redefine default writer */
  /* ----------------------- */

//  cs_post_define_writer(CS_POST_WRITER_PROBES,   /* writer_id */
//			"",                      /* writer name */
//			"monitoring",
//			"time_plot",             /* format name */
//			"csv",                   /* format options */
//			FVM_WRITER_FIXED_MESH,
//			false,                   /* output at start */
//			false,                   /* output at end */
//			-1,                      /* time step frequency */
//			0.1.);                    /* time value frequency */

  /* Define additional writers */
  /* ------------------------- */

  cs_post_define_writer(1,                       /* writer_id */
			"mid_plane",             /* writer name */
			"postprocessing",        /* directory name */
			format_name,
			format_options,
			time_dep,
			output_at_start,
			output_at_end,
			frequency_n,
			frequency_t);

  cs_post_define_writer(2,                            /* writer_id */
			"ventilation_plane",          /* writer name */
			"postprocessing",             /* directory name */
			format_name,
			format_options,
			time_dep,
			output_at_start,
			output_at_end,
			frequency_n,
			frequency_t);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define post-processing meshes.
 *
 * The main post-processing meshes may be configured, and additional
 * post-processing meshes may be defined as a subset of the main mesh's
 * cells or faces (both interior and boundary).
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_postprocess_meshes
void
cs_user_postprocess_meshes(void)
{
  {
    /* Select interior faces with y = 0.0 */

    const int n_writers = 1;
    const int writer_ids[] = {1};  /* Associate to writers 1 */

    const char *cell_criteria = "plane[-1, 0, 0, 0.0, epsilon = 0.2]";

    cs_post_define_volume_mesh(1,
			       "Median plane",
			       cell_criteria,
			       false,
			       true,
			       n_writers,
			       writer_ids);
  }

  {
    const int n_writers = 1;
    const int writer_ids[] = {2};

    const char *cell_criteria = "plane[0, 0, -1, 3.40, epsilon = 0.2]";

    cs_post_define_volume_mesh(2,
			       "Ventilation plane",
			       cell_criteria,
			       false,
			       true,
			       n_writers,
			       writer_ids);
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define monitoring probes and profiles.
 *
 * Profiles are defined as sets of probes.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_postprocess_probes
void
cs_user_postprocess_probes(void)
{
  /* Define volume (including near-boundary) monitoring probe sets */

  cs_probe_set_t  *pset_v = cs_probe_set_create("probes");

  /* Define boundary monitoring probe sets */

  cs_probe_set_t  *pset_w = cs_probe_set_create("wall_probes");
  cs_probe_set_option(pset_w, "boundary", "true");

  int writer_ids[] = {CS_POST_WRITER_PROBES};

  cs_probe_set_associate_writers(pset_w, 1, writer_ids);

  /* Thermocouples mounted on trees (volume probes) */

  cs_probe_set_t  *pset = cs_probe_set_get("probes");

  int n_thermo = 18;
  int n_probes = 0;

  char label[32];
  double z_gas[18] = {3.90, 3.85, 3.80, 3.55, 3.30, 3.05, 2.80,
		      2.55, 2.30, 2.05, 1.80, 1.55, 1.30, 1.05,
		      0.80, 0.55, 0.30, 0.05};

  /* TG_L2_NE */
  for (int i = 0; i < n_thermo; i++) {
    n_probes += 1;
    cs_real_t x = -1.5, y = 1.25, z = z_gas[i];
    snprintf(label, 31, "%i_TGL2_NE_%g",n_probes,z);
    cs_probe_set_add_probe(pset, x, y, z, label);
  }

  /* TG_L2_SE */
  for (int i = 0; i < n_thermo; i++) {
    n_probes += 1;
    cs_real_t x = 1.5, y = 1.25, z = z_gas[i];
    snprintf(label, 31, "%i_TGL2_SE_%g",n_probes,z);
    cs_probe_set_add_probe(pset, x, y, z, label);
  }

  /* TG_L2_SW */
  for (int i = 0; i < n_thermo; i++) {
    n_probes += 1;
    cs_real_t x = 1.5, y = -1.25, z = z_gas[i];
    snprintf(label, 31, "%i_TGL2_SW_%g",n_probes,z);
    cs_probe_set_add_probe(pset, x, y, z, label);
  }

  /* Gas measurements */

  n_probes += 1;
  snprintf(label, 31, "%i_GAS_L2_FP",n_probes);
  cs_probe_set_add_probe(pset, -0.8,0, 0.35, label);

  n_probes += 1;
  snprintf(label, 31, "%i_GAS_L2_SE_low",n_probes);
  cs_probe_set_add_probe(pset,1.5, 1.25, 0.80, label);

  n_probes += 1;
  snprintf(label, 31, "%i_GAS_L2_SE_high",n_probes);
  cs_probe_set_add_probe(pset,1.5, 1.25, 3.30, label);

  /* Wall probes */

  double z_wall[4] = {0.30, 1.55, 2.60, 3.55};
  double z_wall_c[1] = {2.60};

  n_thermo = 4;
  n_probes = 0;

  pset = cs_probe_set_get("wall_probes");

  for (int i = 0; i < n_thermo; i++) {
    n_probes += 1;
    cs_real_t x = -3.0, y = 0, z = z_wall[i];
    snprintf(label, 31, "%i_TPL2_N_%g",n_probes,z);
    cs_probe_set_add_probe(pset, x, y, z, label);
  }

  for (int i = 0; i < n_thermo; i++) {
    n_probes += 1;
    cs_real_t x = -3.0, y = 1.5, z = z_wall[i];
    snprintf(label, 31, "%i_FLT_L2_N_%g",n_probes,z);
    cs_probe_set_add_probe(pset, x, y, z, label);
  }

  for (int i = 0; i < n_thermo; i++) {
    n_probes += 1;
    cs_real_t x = 3.0, y = 0, z = z_wall_c[i];
    snprintf(label, 31, "%i_TPL2_SC_%g",n_probes,z);
    cs_probe_set_add_probe(pset, x, y, z, label);
  }

  for (int i = 0; i < n_thermo; i++) {
    n_probes += 1;
    cs_real_t x = 0.0, y = 2.50, z = z_wall_c[i];
    snprintf(label, 31, "%i_TPL2_EC_%g",n_probes,z);
    cs_probe_set_add_probe(pset, x, y, z, label);
  }

  for (int i = 0; i < n_thermo; i++) {
    n_probes += 1;
    cs_real_t x = 0.0, y = -2.50, z = z_wall_c[i];
    snprintf(label, 31, "%i_TPL2_WC_%g",n_probes,z);
    cs_probe_set_add_probe(pset, x, y, z, label);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User function for output of values on a post-processing mesh.
 *
 * \param[in]       mesh_name    name of the output mesh for the current call
 * \param[in]       mesh_id      id of the output mesh for the current call
 * \param[in]       cat_id       category id of the output mesh for the
 *                               current call
 * \param[in]       probes       pointer to associated probe set structure if
 *                               the mesh is a probe set, nullptr otherwise
 * \param[in]       n_cells      local number of cells of post_mesh
 * \param[in]       n_i_faces    local number of interior faces of post_mesh
 * \param[in]       n_b_faces    local number of boundary faces of post_mesh
 * \param[in]       n_vertices   local number of vertices faces of post_mesh
 * \param[in]       cell_list    list of cells (0 to n-1) of post-processing
 *                               mesh
 * \param[in]       i_face_list  list of interior faces (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       b_face_list  list of boundary faces (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       vertex_list  list of vertices (0 to n-1) of
 *                               post-processing mesh
 * \param[in]       ts           time step status structure, or nullptr
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_postprocess_values
void
cs_user_postprocess_values(const char            *mesh_name,
                           int                    mesh_id,
                           int                    cat_id,
                           cs_probe_set_t        *probes,
                           cs_lnum_t              n_cells,
                           cs_lnum_t              n_i_faces,
                           cs_lnum_t              n_b_faces,
                           cs_lnum_t              n_vertices,
                           const cs_lnum_t        cell_list[],
                           const cs_lnum_t        i_face_list[],
                           const cs_lnum_t        b_face_list[],
                           const cs_lnum_t        vertex_list[],
                           const cs_time_step_t  *ts)
{
  CS_UNUSED(mesh_name);
  CS_UNUSED(mesh_id);
  CS_UNUSED(cat_id);
  CS_UNUSED(probes);
  CS_UNUSED(n_cells);
  CS_UNUSED(n_i_faces);
  CS_UNUSED(n_b_faces);
  CS_UNUSED(n_vertices);
  CS_UNUSED(cell_list);
  CS_UNUSED(i_face_list);
  CS_UNUSED(b_face_list);
  CS_UNUSED(vertex_list);
  CS_UNUSED(ts);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
