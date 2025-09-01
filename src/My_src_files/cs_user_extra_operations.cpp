/*============================================================================
 * cs_user_extra_operations.cpp - FULL translation from Fortran
 *============================================================================*/

#include "cs_defs.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "cs_headers.h"

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
void
cs_user_extra_operations(cs_domain_t *domain)
{
  const cs_lnum_t nt_cur = cs_glob_time_step->nt_cur;
  const cs_real_t ttcabs = cs_glob_time_step->t_cur;
  const cs_lnum_t nfabor = domain->mesh->n_b_faces;
  const cs_real_t *cpro_wt = CS_F_(t_b)->val;

  const cs_mesh_t *m = cs_glob_mesh;
  const int n_cells = m->n_cells;

  const cs_real_t *rho = CS_F_(rho)->val;

  cs_fluid_properties_t *fp = cs_get_glob_fluid_properties();
  cs_combustion_gas_model_t  *cm = cs_glob_combustion_gas_model;
  const int ngas_g = cm -> n_gas_species;
  const int ngas_e = cm -> n_gas_el_comp;

  const cs_real_t xsoot = cm-> xsoot;
  const cs_real_t rosoot = cm-> rosoot;

  cs_real_t *yg = (cs_real_t *)malloc(ngas_g * sizeof(cs_real_t));
  cs_real_t *ye = (cs_real_t *)malloc(ngas_e * sizeof(cs_real_t));
  cs_real_t *xe = (cs_real_t *)malloc(ngas_e * sizeof(cs_real_t));

  // Variables de champ
  cs_real_t *bmasfl = cs_field_by_name("boundary_mass_flux")->val;
  cs_real_t *brom = cs_field_by_name("boundary_density")->val;
  cs_real_t *bfqinc = cs_field_by_name("rad_incident_flux")->val;
  cs_real_t *bfconv = cs_field_by_name("rad_convective_flux")->val;
  cs_real_t *bhconv = cs_field_by_name("rad_exchange_coefficient")->val;
  cs_real_t *cvar_wt = cs_field_by_name("wall_temp")->val;
  cs_real_t *cvar_qinc = cs_field_by_name("wall_incident_flux")->val;
  cs_real_t *cvar_qconv = cs_field_by_name("wall_convective_flux")->val;
  cs_real_t *cvar_wf = cs_field_by_name("wall_total_flux")->val;
  cs_real_t *cpro_Yf = cs_field_by_name("ym_fuel")->val;
  cs_real_t *cpro_Yo = cs_field_by_name("ym_oxyd")->val;
  cs_real_t *cpro_Yp = cs_field_by_name("ym_prod")->val;
  cs_real_t *bym1 = cs_field_by_name("Xm_O2")->val;
  cs_real_t *bym2 = cs_field_by_name("Xm_CO2")->val;
  cs_real_t *cvar_fm = cs_field_by_name("Cm_Soot")->val;

  if (cm->isoot >= 1) {
    cvar_fm = CS_F_(fsm)->val;
  }

  for (cs_lnum_t cell_id = 0; cell_id < n_cells; cell_id++) {
    yg[0] = cpro_Yf[cell_id];
    yg[1] = cpro_Yo[cell_id];
    yg[2] = cpro_Yp[cell_id];

    cs_combustion_gas_yg2xye(yg,ye,xe);

    bym1[cell_id] = xe[cm -> iio2]; /*xe[2]*/
    bym2[cell_id] = xe[cm -> iico2]; /*xe[3]*/
    /*cvar_fm[cell_id] = cm -> rosoot[cell_id]*/ /*xe[cm -> iic] *rho[cell_id];*/

    cs_real_t xs; /*massic fraction of soot*/
    if (cm -> isoot ==0 && cm -> iic >0) {
      // no soot model but carbon is calculated through iic
      xs = cpro_Yp[cell_id] *cm -> coefeg[2][cm -> iic-1]; /*-1 base 0 ?*/
    }
    else if (cm->isoot == 0) {
      // no soot model, so we calcul it with soot fraction production
      xs = xsoot * cpro_Yp[cell_id];
    }
    else if (cm->isoot >= 1) {
      // soot model already active just copy user field Cm_Soot
      xs = cvar_fm[cell_id];
    }
    else {
      // if nothing exist just keep 0
      xs = 0;
    }
    cvar_fm[cell_id] = xs; /*mass fraction */
    // bft_printf("xs %g",xs);
  }

  free(yg);
  free(ye);
  free(xe);

  cs_real_t mfuel = 0., madm = 0., mext = 0.;
  cs_real_t mfuel_local = 0., madm_local = 0., mext_local = 0.;
  cs_real_t wall_temp = 0., wall_flux = 0.;
  cs_real_t en_masfl = 0., en_o2_masfl = 0.;
  cs_real_t so_masfl = 0., so_o2_masfl = 0.;
  //
  // cs_lnum_t *lstelt = (cs_lnum_t *)malloc(sizeof(cs_lnum_t) * nfabor);
  // cs_lnum_t nlelt;

  const cs_zone_t  *zn = nullptr;

  // Fire_pan
  zn = cs_boundary_zone_by_name("Fire_pan");
  cs_lnum_t n_fpan_local = zn->n_elts;
  for (cs_lnum_t ilelt = 0; ilelt < n_fpan_local; ilelt++) {
    const cs_lnum_t face_id = zn->elt_ids[ilelt];
    mfuel_local -= bmasfl[face_id];
  }

  // Inlets
  zn = cs_boundary_zone_by_name("Inlet_All");
  cs_lnum_t n_inlet_local = zn->n_elts;
  for (cs_lnum_t ilelt = 0; ilelt < n_inlet_local; ilelt++) {
    const cs_lnum_t face_id = zn->elt_ids[ilelt];
    cs_real_t rho = brom[face_id];
    cs_real_t flumab = bmasfl[face_id];
    madm_local -= 3600. * flumab / rho;
    // en_masfl += flumab;
  }

  // Outlets
  zn = cs_boundary_zone_by_name("Outlet_All");
  cs_lnum_t n_outlet_local = zn->n_elts;
  for (cs_lnum_t ilelt = 0; ilelt < n_outlet_local ; ilelt++) {
    const cs_lnum_t face_id = zn->elt_ids[ilelt];
    cs_lnum_t iel = domain->mesh->b_face_cells[face_id];
    cs_real_t rho = brom[face_id];
    cs_real_t flumab = bmasfl[face_id];
    mext_local += 3600. * flumab / rho;
    // so_masfl += flumab;
  }

  // Walls
  zn = cs_boundary_zone_by_name("Walls_All");
  cs_lnum_t n_walls = zn->n_elts;
  for (cs_lnum_t ilelt = 0; ilelt < n_walls; ilelt++) {
    const cs_lnum_t face_id = zn->elt_ids[ilelt];
    cvar_wt[face_id] = cpro_wt[face_id];
    cvar_qinc[face_id] = bfqinc[face_id];
    cvar_qconv[face_id] = bfconv[face_id];
    cvar_wf[face_id] = bfqinc[face_id] + bfconv[face_id];
    wall_temp += cpro_wt[face_id];
    wall_flux += cvar_wf[face_id];
  }

  // /* global reduction */
  cs_real_t mfuel_global = mfuel_local;
  cs_real_t madm_global = madm_local;
  cs_real_t mext_global = mext_local;
  cs_lnum_t n_fpan_global = n_fpan_local;
  cs_lnum_t n_inlet_global = n_inlet_local;
  cs_lnum_t n_outlet_global = n_outlet_local;
  if (cs_glob_n_ranks > 1) {
    cs_parall_sum(1, CS_REAL_TYPE, &mfuel_global);
    cs_parall_sum(1, CS_REAL_TYPE, &madm_global);
    cs_parall_sum(1, CS_REAL_TYPE, &mext_global);
    cs_parall_sum(1, CS_GNUM_TYPE, &n_fpan_global);
    cs_parall_sum(1, CS_GNUM_TYPE, &n_inlet_global);
    cs_parall_sum(1, CS_GNUM_TYPE, &n_outlet_global);
  }

  if (n_fpan_global > 0) {
    mfuel = mfuel_global;
    mfuel /= (n_fpan_global);
    }
  if (n_inlet_global > 0) {
    madm = madm_global;
  }
  if (n_inlet_global > 0) {
    mext = mext_global;
  }
  // bft_printf("mfuel = %e\n", mfuel);
  // bft_printf("n_fpan = %d\n",n_fpan_global);
  // bft_printf("n_inlet = %d\n",n_inlet_global);
  // bft_printf("n_outlet = %d\n",n_outlet_global);
  // bft_printf("n_walls = %d\n",n_walls);

  if (cs_glob_rank_id ==0 ) {
    FILE *fop = fopen("results.dat",nt_cur == cs_glob_time_step->nt_ini ? "w" : "a");
    if (fop != nullptr) {
      if (nt_cur == cs_glob_time_step->nt_ini) {
        fprintf(fop, "#%10s %15s %12s %12s %12s %12s %12s\n","TIME","MFUEL","DP","DADM_MOY","DEXT_MOY","WALL_TEMP","WALL_FLUX");
        fprintf(fop, "#%10s %15s %12s %12s %12s %12s %12s\n","(s)","(kg/s)","(Pa)","(m3/h)","(m3/h)","(deg.C)","(W/m2)");
      }
      fprintf(fop, "#%10.1f %15.5e %12.5f %12.5f %12.5f %12.5f %12.5f\n",ttcabs,mfuel,fp->pther-fp->p0,madm,mext, wall_temp/n_walls,wall_flux/n_walls);
      fclose(fop);
    }
  }
}

END_C_DECLS
