/*============================================================================
 * cs_user_initialization.cpp
 *============================================================================*/

#include "cs_defs.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include <ple_coupling.h>

#include "cs_headers.h"

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * rief User initialization routine (translated from Fortran)
 */
void
cs_user_initialization(cs_domain_t *domain)
{
  cs_real_t *k = nullptr;
  cs_real_t *eps = nullptr;
  cs_real_t *nut = nullptr;

  const cs_lnum_t n_cells = domain->mesh->n_cells;
  const cs_real_t *coor = domain->mesh_quantities->cell_cen;

  // // Field pointers
  // k = CS_F_(k)->val;
  // eps = CS_F_(eps)->val;
  //
  // if (k == nullptr);
  //   bft_error(__FILE__, __LINE__, 0, "Field'k'not available!\n");

  for (cs_lnum_t iel = 0; iel < n_cells; iel++) {

  const cs_real_t x = coor[3 * iel    ];
  const cs_real_t y = coor[3 * iel + 1];
  const cs_real_t z = coor[3 * iel + 2];

  // // Initial values (example: uniform or based on position)
  // k[iel]   = 1.e-4;
  // eps[iel] = 1.e-4;
}
}

END_C_DECLS
