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
#include "admmpd_tetmesh.h"
#include "admmpd_embeddedmesh.h"
#include "admmpd_collision.h"
#include "admmpd_pin.h"

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
  admmpd::TetMeshData *tetmesh; // init_mode=0
  admmpd::EmbeddedMeshData *embmesh; // init_mode=1
  admmpd::Collision *collision;
  admmpd::Pin *pin;
  int in_totverts; // number of input verts
};

void admmpd_dealloc(ADMMPDInterfaceData *iface)
{
  if (iface==NULL)
    return;

  iface->totverts = 0; // output vertices

  if (iface->idata)
  {
    if (iface->idata->options)
      delete iface->idata->options;
    if (iface->idata->data)
      delete iface->idata->data;
    if (iface->idata->tetmesh)
      delete iface->idata->tetmesh;
    if (iface->idata->embmesh)
      delete iface->idata->embmesh;
    if (iface->idata->collision)
      delete iface->idata->collision;
    delete iface->idata;
  }

  iface->idata = NULL;
}

static int admmpd_init_with_tetgen(
  ADMMPDInterfaceData *iface, float *in_verts, unsigned int *in_faces,
  Eigen::MatrixXd *V, Eigen::MatrixXi *T, Eigen::VectorXd *m)
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
  iface->idata->tetmesh->faces.resize(iface->mesh_totfaces,3);
  iface->idata->tetmesh->x_rest.resize(nv,3);
  iface->idata->tetmesh->tets.resize(nt,3);
  for (int i=0; i<iface->mesh_totfaces; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      iface->idata->tetmesh->faces(i,j) = in_faces[i*3+j];
    } 
  }

  for (int i=0; i<nv; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      V->operator()(i,j) = tg.out_verts[i*3+j];
      iface->idata->tetmesh->x_rest(i,j) = tg.out_verts[i*3+j];
    }
  }
  T->setZero();
  for (int i=0; i<nt; ++i)
  {
    T->operator()(i,0) = tg.out_tets[i*4+0];
    T->operator()(i,1) = tg.out_tets[i*4+1];
    T->operator()(i,2) = tg.out_tets[i*4+2];
    T->operator()(i,3) = tg.out_tets[i*4+3];
    iface->idata->tetmesh->tets.row(i) = T->row(i);
  }

  admmpd::TetMesh().compute_masses(
       iface->idata->tetmesh, V, m);

  // Clean up tetgen data
  MEM_freeN(tg.out_tets);
  MEM_freeN(tg.out_facets);
  MEM_freeN(tg.out_verts);
  return 1;
}

static int admmpd_init_with_lattice(
  ADMMPDInterfaceData *iface, float *in_verts, unsigned int *in_faces,
  Eigen::MatrixXd *V, Eigen::MatrixXi *T, Eigen::VectorXd *m)
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
  bool trim_lattice = true;
  bool success = admmpd::EmbeddedMesh().generate(in_V,in_F,iface->idata->embmesh,trim_lattice);
  if (success)
  {
    admmpd::EmbeddedMesh().compute_masses(iface->idata->embmesh, m);
    *T = iface->idata->embmesh->tets;
    *V = iface->idata->embmesh->rest_x;
    iface->totverts = V->rows();
    return 1;
  }

  return 0;
}

int admmpd_init(ADMMPDInterfaceData *iface, ADMMPDInitData *in_mesh)
{
  if (iface==NULL)
    return 0;
  if (in_mesh->verts==NULL || in_mesh->faces==NULL)
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
  iface->idata->tetmesh = new admmpd::TetMeshData();
  iface->idata->embmesh = new admmpd::EmbeddedMeshData();
  iface->idata->collision = NULL;
  iface->idata->pin = NULL;

  // Generate tets and vertices
  Eigen::MatrixXd V; // defo verts
  Eigen::MatrixXi T; // defo tets
  Eigen::VectorXd m; // masses
  int gen_success = 0;
  switch (iface->init_mode)
  {
    default:
    case 0: {
      gen_success = admmpd_init_with_tetgen(iface,in_mesh->verts,in_mesh->faces,&V,&T,&m);
      //iface->idata->collision = new admmpd::TetMeshCollision();
      } break;
    case 1: {
      gen_success = admmpd_init_with_lattice(iface,in_mesh->verts,in_mesh->faces,&V,&T,&m);
      iface->idata->collision = new admmpd::EmbeddedMeshCollision(iface->idata->embmesh);
      iface->idata->pin = new admmpd::EmbeddedMeshPin(iface->idata->embmesh);
    } break;
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
    init_success = admmpd::Solver().init(V, T, m, options, data);
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

void admmpd_update_obstacles(
    ADMMPDInterfaceData *iface,
    float *in_verts_0,
    float *in_verts_1,
    int nv,
    unsigned int *in_faces,
    int nf)
{
    if (iface==NULL || in_verts_0==NULL || in_verts_1==NULL || in_faces==NULL)
      return;
    if (iface->idata==NULL)
      return;
    if (iface->idata->collision==NULL)
      return;

    iface->idata->collision->set_obstacles(
      in_verts_0, in_verts_1, nv, in_faces, nf);
}

void admmpd_update_goals(
    ADMMPDInterfaceData *iface,
    float *goal_k, // goal stiffness, nv
    float *goal_pos, // goal position, nv x 3
    int nv)
{
    if (iface==NULL || goal_k==NULL || goal_pos==NULL)
      return;
    if (iface->idata==NULL)
      return;
    if (iface->idata->pin==NULL)
      return;

    for (int i=0; i<nv; ++i)
    {
      if (goal_k[i] <= 0.f)
        continue;

      Eigen::Vector3d ki = Eigen::Vector3d::Ones() * goal_k[i];
      Eigen::Vector3d qi(goal_pos[i*3+0], goal_pos[i*3+1], goal_pos[i*3+2]);
      iface->idata->pin->set_pin(i,qi,ki);
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
    admmpd::Solver().solve(
        iface->idata->options,
        iface->idata->data,
        iface->idata->collision,
        iface->idata->pin);
  }
  catch(const std::exception &e)
  {
    iface->idata->data->x = iface->idata->data->x_start;
    printf("**ADMMPD Error on solve: %s\n", e.what());
  }
}