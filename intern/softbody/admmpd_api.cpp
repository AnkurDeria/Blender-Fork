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
#include "tetgen_api.h"
#include "DNA_mesh_types.h" // Mesh
#include "DNA_meshdata_types.h" // MVert
#include "BKE_mesh_remesh_voxel.h" // TetGen
#include "BKE_mesh.h" // BKE_mesh_free
#include "BKE_softbody.h" // BodyPoint
#include "MEM_guardedalloc.h" // 

#include <iostream>


struct ADMMPDInternalData {
  admmpd::Options *options;
  admmpd::Data *data;
//  admmpd::Lattice *lattice;
  int in_totverts; // number of input verts
};

void admmpd_dealloc(ADMMPDInterfaceData *iface)
{
  if (iface==NULL)
    return;

  iface->totverts = 0;
  iface->mesh_totverts = 0;
  iface->mesh_totfaces = 0;

  if (iface->data)
  {
    if(iface->data->options)
        delete iface->data->options;
    if(iface->data->data)
        delete iface->data->data;
    delete iface->data;
  }

  iface->data = NULL;
}

int admmpd_init(ADMMPDInterfaceData *iface, float *in_verts, unsigned int *in_faces)
{

  if (iface==NULL)
    return 0;
  if (in_verts==NULL || in_faces==NULL)
    return 0;
  if (iface->mesh_totverts<=0 || iface->mesh_totfaces<=0)
    return 0;

  // Generate tets
  TetGenRemeshData tg;
  init_tetgenremeshdata(&tg);
  tg.in_verts = in_verts;
  tg.in_totverts = iface->mesh_totverts;
  tg.in_faces = in_faces;
  tg.in_totfaces = iface->mesh_totfaces;
  bool success = tetgen_resmesh(&tg);
  if (!success || tg.out_tottets==0)
    return 0;

  // Create initializer for ADMMPD
  iface->totverts = tg.out_totverts;
  int nv = tg.out_totverts;
  int nt = tg.out_tottets;
  Eigen::MatrixXd V(nv,3);
  Eigen::MatrixXi T(nt,4);
  V.setZero();
  for (int i=0; i<nv; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      V(i,j) = tg.out_verts[i*3+j];
    }
  }
  T.setZero();
  for (int i=0; i<nt; ++i)
  {
    T(i,0) = tg.out_tets[i*4+0];
    T(i,1) = tg.out_tets[i*4+1];
    T(i,2) = tg.out_tets[i*4+2];
    T(i,3) = tg.out_tets[i*4+3];
  }

  // Generate solver data
  iface->data = new ADMMPDInternalData();
  iface->data->options = new admmpd::Options();
  admmpd::Options *options = iface->data->options;
  iface->data->data = new admmpd::Data();
  admmpd::Data *data = iface->data->data;

  // Initialize
  bool init_success = false;
  try
  {
    init_success = admmpd::Solver().init(V, T, options, data);
  }
  catch(const std::exception &e)
  {
    printf("**ADMMPD Error on init: %s\n", e.what());
  }

  // Clean up tetgen data
  MEM_freeN(tg.out_tets);
  MEM_freeN(tg.out_facets);
  MEM_freeN(tg.out_verts);
  return int(init_success);
}

void admmpd_copy_from_bodypoint(ADMMPDInterfaceData *iface, const BodyPoint *pts)
{
  if (iface == NULL || pts == NULL)
    return;

  for (int i=0; i<iface->totverts; ++i)
  {
    const BodyPoint *pt = &pts[i];
    for(int j=0; j<3; ++j)
    {
      iface->data->data->x(i,j)=pt->pos[j];
      iface->data->data->v(i,j)=pt->vec[j];
    }
  }
}

void admmpd_copy_to_bodypoint_and_object(ADMMPDInterfaceData *iface, BodyPoint *pts, float (*vertexCos)[3])
{
  if (iface == NULL)
    return;

  for (int i=0; i<iface->totverts; ++i)
  {
    if (pts != NULL)
    {
      BodyPoint *pt = &pts[i];
      for(int j=0; j<3; ++j)
      {
        pt->pos[j] = iface->data->data->x(i,j);
        pt->vec[j] = iface->data->data->v(i,j);
      }
    }

    // TODO eventually replace with mapping. For now
    // we assume TetGen, which the first N vertices are
    // mapped one-to-one.
    if (vertexCos != NULL && i<iface->mesh_totverts)
    {
      vertexCos[i][0] = iface->data->data->x(i,0);
      vertexCos[i][1] = iface->data->data->x(i,1);
      vertexCos[i][2] = iface->data->data->x(i,2);
    }
  }
}

void admmpd_solve(ADMMPDInterfaceData *iface)
{
  
  if (iface == NULL)
    return;

  try
  {
    admmpd::Solver().solve(iface->data->options,iface->data->data);
  }
  catch(const std::exception &e)
  {
    printf("**ADMMPD Error on solve: %s\n", e.what());
  }
}