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
#include "MEM_guardedalloc.h" // 

#include <iostream>


struct ADMMPDInternalData {
  admmpd::Options *options;
  admmpd::Data *data;
//  admmpd::Lattice *lattice;
};


void admmpd_alloc(ADMMPDInterfaceData *iface)
{
  if (iface==NULL)
    return;

  if (iface->in_verts != NULL)
  {
    MEM_freeN(iface->in_verts);
    iface->in_verts = NULL;
  }
  if (iface->in_vel != NULL)
  {
    MEM_freeN(iface->in_vel);
    iface->in_vel = NULL;
  }
  if (iface->in_faces != NULL)
  {
    MEM_freeN(iface->in_faces);
    iface->in_faces = NULL;
  }

  iface->in_verts = (float *)MEM_mallocN(iface->in_totverts*3*sizeof(float), "admmpd_verts");
  iface->in_vel = (float *)MEM_mallocN(iface->in_totverts*3*sizeof(float), "admmpd_vel");
  iface->in_faces = (unsigned int *)MEM_mallocN(iface->in_totfaces*3*sizeof(unsigned int), "admmpd_faces");
}

void admmpd_dealloc(ADMMPDInterfaceData *iface)
{
  if (iface==NULL)
    return;

  iface->in_totverts = 0;
  if (iface->in_verts != NULL)
    MEM_freeN(iface->in_verts);
  if (iface->in_vel != NULL)
    MEM_freeN(iface->in_vel);

  iface->in_totfaces = 0;
  if (iface->in_faces != NULL)
    MEM_freeN(iface->in_faces);

  iface->out_totverts = 0;
  if (iface->out_verts != NULL)
    MEM_freeN(iface->out_verts);
  if (iface->out_vel != NULL)
    MEM_freeN(iface->out_vel);

  if (iface->data)
  {
    if(iface->data->options)
        delete iface->data->options;
    if(iface->data->data)
        delete iface->data->data;
    delete iface->data;
  }

  iface->in_verts = NULL;
  iface->in_vel = NULL;
  iface->in_faces = NULL;
  iface->out_verts = NULL;
  iface->out_vel = NULL;
  iface->data = NULL;
}

int admmpd_init(ADMMPDInterfaceData *iface)
{
  if (iface==NULL)
    return 0;
  if (iface->in_verts==NULL || iface->in_vel==NULL || iface->in_faces==NULL)
    return 0;
  if (iface->in_totverts<=0 || iface->in_totfaces<=0)
    return 0;

  // Generate tets
  TetGenRemeshData tg;
  init_tetgenremeshdata(&tg);
  tg.in_verts = iface->in_verts;
  tg.in_totverts = iface->in_totverts;
  tg.in_faces = iface->in_faces;
  tg.in_totfaces = iface->in_totfaces;
  bool success = tetgen_resmesh(&tg);
  if (!success || tg.out_tottets==0)
    return 0;

  // Resize data
  iface->out_totverts = tg.out_totverts;
  if (iface->out_verts != NULL)
  {
    MEM_freeN(iface->out_verts);
    iface->out_verts = NULL;
  }
  if (iface->out_vel != NULL)
  {
    MEM_freeN(iface->out_vel);
    iface->out_vel = NULL;
  }
  iface->out_verts = (float *)MEM_callocN(
      iface->out_totverts*3*sizeof(float), "ADMMPD_out_verts");
  iface->out_vel = (float *)MEM_callocN(
      iface->out_totverts*3*sizeof(float), "ADMMPD_out_vel");

  // Create initializer for ADMMPD
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
      iface->out_verts[i*3+j] = tg.out_verts[i*3+j];
      iface->out_vel[i*3+j] = 0;
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

void admmpd_solve(ADMMPDInterfaceData *iface)
{
  if (iface == NULL)
    return;

  // Whatever is in out_verts and out_vel needs
  // to be mapped to internal data, as it's used as input
  // when reading from cached data.
  int nv = iface->out_totverts;
  for (int i=0; i<nv; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      iface->data->data->x(i,j) = iface->out_verts[i*3+j];
      iface->data->data->v(i,j) = iface->out_vel[i*3+j];
    }
  }

  try
  {
    admmpd::Solver().solve(iface->data->options,iface->data->data);
    for (int i=0; i<iface->data->data->x.rows(); ++i)
    {
      for (int j=0; j<3; ++j)
      {
        iface->out_verts[i*3+j] = iface->data->data->x(i,j);
        iface->out_vel[i*3+j] = iface->data->data->v(i,j);
      }
    }
  }
  catch(const std::exception &e)
  {
    printf("**ADMMPD Error on solve: %s\n", e.what());
  }
}

void admmpd_get_vertices(ADMMPDInterfaceData *iface, float (*vertexCos)[3], int numVerts)
{
  if (iface == NULL)
    return;

  if (numVerts != iface->in_totverts || numVerts > iface->out_totverts)
  {
    printf("**ADMMPD TODO: PROPER VERTEX MAPPINGS\n");
    return;
  }
  // TODO double check this, but I believe the first n verts
  // created by tetgen are the same as the input. I hope. We'll find out I guess.
  // If not, this function will be a placeholder for the mapping that
  // will have to occur.
  for (int i=0; i<numVerts; ++i)
  {
    vertexCos[i][0] = iface->out_verts[i*3+0];
    vertexCos[i][1] = iface->out_verts[i*3+1];
    vertexCos[i][2] = iface->out_verts[i*3+2];
  }
}