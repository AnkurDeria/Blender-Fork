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
#include "admmpd_types.h"
#include "admmpd_solver.h"
#include "admmpd_embeddedmesh.h"
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
  admmpd::SolverData *data;
  admmpd::EmbeddedMeshData *embmesh;
  int in_totverts; // number of input verts
};

void admmpd_dealloc(ADMMPDInterfaceData *iface)
{
  if (iface==NULL)
    return;

  iface->totverts = 0; // output vertices

  if (iface->idata)
  {
    if(iface->idata->options)
        delete iface->idata->options;
    if(iface->idata->data)
        delete iface->idata->data;
    if(iface->idata->embmesh)
        delete iface->idata->embmesh;
    delete iface->idata;
  }

  iface->idata = NULL;
}

static int admmpd_init_with_tetgen(
  ADMMPDInterfaceData *iface, float *in_verts, unsigned int *in_faces,
  Eigen::MatrixXd *V, Eigen::MatrixXi *T)
{
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
  V->resize(nv,3);
  T->resize(nt,4);
  V->setZero();
  for (int i=0; i<nv; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      V->operator()(i,j) = tg.out_verts[i*3+j];
    }
  }
  T->setZero();
  for (int i=0; i<nt; ++i)
  {
    T->operator()(i,0) = tg.out_tets[i*4+0];
    T->operator()(i,1) = tg.out_tets[i*4+1];
    T->operator()(i,2) = tg.out_tets[i*4+2];
    T->operator()(i,3) = tg.out_tets[i*4+3];
  }
  // Clean up tetgen data
  MEM_freeN(tg.out_tets);
  MEM_freeN(tg.out_facets);
  MEM_freeN(tg.out_verts);
  return 1;
}

static int admmpd_init_with_lattice(
  ADMMPDInterfaceData *iface, float *in_verts, unsigned int *in_faces,
  Eigen::MatrixXd *V, Eigen::MatrixXi *T)
{

  int nv = iface->mesh_totverts;
  Eigen::MatrixXd in_V(nv,3);
  for (int i=0; i<nv; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      in_V(i,j) = in_verts[i*3+j];
    }
  }

  int nf = iface->mesh_totfaces;
  Eigen::MatrixXi in_F(nf,3);
  for (int i=0; i<nf; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      in_F(i,j) = in_faces[i*3+j];
    }
  }

  iface->totverts = 0;
  bool success = admmpd::EmbeddedMesh().generate(in_V,in_F,iface->idata->embmesh,V);
  if (success)
  {
    iface->totverts = V->rows();
    *T = iface->idata->embmesh->tets;
    return 1;
  }

  return 0;
}

int admmpd_init(ADMMPDInterfaceData *iface, float *in_verts, unsigned int *in_faces)
{
  if (iface==NULL)
    return 0;
  if (in_verts==NULL || in_faces==NULL)
    return 0;
  if (iface->mesh_totverts<=0 || iface->mesh_totfaces<=0)
    return 0;

  // Delete any existing data
  admmpd_dealloc(iface);

  // Generate solver data
  iface->idata = new ADMMPDInternalData();
  iface->idata->options = new admmpd::Options();
  admmpd::Options *options = iface->idata->options;
  iface->idata->data = new admmpd::SolverData();
  admmpd::SolverData *data = iface->idata->data;
  iface->idata->embmesh = new admmpd::EmbeddedMeshData();

  // Generate tets and vertices
  Eigen::MatrixXd V;
  Eigen::MatrixXi T;
  int gen_success = 0;
  switch (iface->init_mode)
  {
    default:
    case 0:
      gen_success = admmpd_init_with_tetgen(iface,in_verts,in_faces,&V,&T);
      break;
    case 1:
      gen_success = admmpd_init_with_lattice(iface,in_verts,in_faces,&V,&T);
      break;
  }
  if (!gen_success || iface->totverts==0)
  {
    printf("**ADMMPD Failed to generate tets\n");
    return 0;
  } 

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

  if (!init_success)
  {
    printf("**ADMMPD Failed to initialize\n");
    return 0;
  }

  return 1;
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
      iface->idata->data->x(i,j)=pt->pos[j];
      iface->idata->data->v(i,j)=pt->vec[j];
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
        pt->pos[j] = iface->idata->data->x(i,j);
        pt->vec[j] = iface->idata->data->v(i,j);
      }
    }

    // If we're using TetGen, then we know the first
    // n vertices of the tet mesh are the input surface mesh.
    if (vertexCos != NULL && iface->init_mode==0 && i<iface->mesh_totverts)
    {
      vertexCos[i][0] = iface->idata->data->x(i,0);
      vertexCos[i][1] = iface->idata->data->x(i,1);
      vertexCos[i][2] = iface->idata->data->x(i,2);
    }
  } // end loop all verts

  // If using lattice, get the embedded vertex position
  // from the deformed lattice.
  if (vertexCos != NULL && iface->init_mode==1)
  {
      for (int i=0; i<iface->mesh_totverts; ++i)
      {
        
        Eigen::Vector3d xi = admmpd::EmbeddedMesh().get_mapped_vertex(
          iface->idata->embmesh, &iface->idata->data->x, i);
        vertexCos[i][0] = xi[0];
        vertexCos[i][1] = xi[1];
        vertexCos[i][2] = xi[2];
      }
  }

} // end map ADMMPD to bodypoint and object

void admmpd_solve(ADMMPDInterfaceData *iface)
{
  
  if (iface == NULL)
    return;

  if (iface->idata == NULL || iface->idata->options == NULL || iface->idata->data == NULL)
    return;

  try
  {
    admmpd::Solver().solve(iface->idata->options,iface->idata->data,NULL);
  }
  catch(const std::exception &e)
  {
    printf("**ADMMPD Error on solve: %s\n", e.what());
  }
}