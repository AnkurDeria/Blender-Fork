/*
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * The Original Code is Copyright (C) 2013 Blender Foundation
 * All rights reserved.
 */

/** \file
 * \ingroup admmpd
 */

#include "admmpd_api.h"
#include "admmpd_solver.h"
#include "admmpd_lattice.h"
#include "DNA_object_types.h" // Object
#include "BKE_softbody.h"
#include <iostream>

struct ADMMPD_Data {
  admmpd::ADMMPD_Options *options;
  admmpd::ADMMPD_Data *data;
  admmpd::Lattice *lattice;
};

ADMMPD_Data* admmpd_init(
    BodyPoint *bp,
    int numVerts)
{
  if (!bp)
    return NULL;

  ADMMPD_Data *admmpd_data = new ADMMPD_Data();
  admmpd_data->options = new admmpd::ADMMPD_Options();
  admmpd::ADMMPD_Options *options = admmpd_data->options;
  admmpd_data->data = new admmpd::ADMMPD_Data();
  admmpd::ADMMPD_Data *data = admmpd_data->data;
  admmpd_data->lattice = new admmpd::Lattice();
  admmpd::Lattice *lattice = admmpd_data->lattice;

  // Create initializer
  Eigen::MatrixXd V(numVerts,3);
  V.setZero();
  for (int i = 0; i < numVerts; ++i)
  {
    BodyPoint &bpi = bp[i];
    V(i,0) = bpi.pos[0];
    V(i,1) = bpi.pos[1];
    V(i,2) = bpi.pos[2];
  }

  // Generate a mapping from x -> V
  // Usually, dim(x) < dim(V)
  Eigen::MatrixXd x;
  Eigen::MatrixXi T;
  lattice->generate(V, &x, &T);

  try
  {
    admmpd::Solver().init(x, T, options, data);
  }
  catch(const std::exception &e)
  {
    printf("**ADMMPD Error on init: %s\n", e.what());
    // Should probably return nullptr
  }

  return admmpd_data;
}

void admmpd_cleanup(ADMMPD_Data *admmpd_data)
{
  if (!admmpd_data)
    return;
  
  delete admmpd_data->data;
  delete admmpd_data->options;
  delete admmpd_data->lattice;
  delete admmpd_data;
}

void admmpd_solve(ADMMPD_Data *admmpd_data)
{
  try
  {
    admmpd::Solver().solve(admmpd_data->options,admmpd_data->data);
  }
  catch(const std::exception &e)
  {
    printf("**ADMMPD Error on solve: %s\n", e.what());
  }
}

void admmpd_to_bodypoint(
    ADMMPD_Data* data,
    BodyPoint *bp,
    int numVerts)
{
  if (!data || !bp)
  {
    printf("DATA OR BP NULL!?\n");
    return;
  }

  if (!data->data || !data->lattice)
    return;

  // Map x -> BodyPoint
  // Usually, dim(x) < dim(V)
  for (int i=0; i<numVerts; ++i)
  {
    BodyPoint &bpi = bp[i];
    Eigen::Vector3d xi = data->lattice->get_mapped_vertex(
      i, &data->data->x, &data->data->tets);
    Eigen::Vector3d vi = data->lattice->get_mapped_vertex(
      i, &data->data->v, &data->data->tets);
    bpi.pos[0] = xi[0];
    bpi.pos[1] = xi[1];
    bpi.pos[2] = xi[2];
    bpi.vec[0] = vi[0];
    bpi.vec[1] = vi[1];
    bpi.vec[2] = vi[2];
  }

} // end map to object