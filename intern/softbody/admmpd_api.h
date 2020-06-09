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
 * The Original Code is Copyright (C) 2013 Blender Foundation,
 * All rights reserved.
 */

/** \file
 * \ingroup admmpd
 */

#ifndef __ADMMPD_API_H__
#define __ADMMPD_API_H__

#ifdef __cplusplus
extern "C" {
#endif

//typedef struct Mesh Mesh_;

typedef struct ADMMPDInterfaceData {
    float *in_verts;
    float *in_vel;
    unsigned int *in_faces;
    int in_totfaces;
    int in_totverts;
    // Num output verts might be different than num input verts.
    // This is due to the lattice/tetmesh that is generated
    // in init. They need to be cached by the system
    float *out_verts;
    float *out_vel;
    int out_totverts;
    // Solver data used internally
    struct ADMMPDInternalData *data;
} ADMMPDInterfaceData;

void admmpd_alloc(ADMMPDInterfaceData*, int in_verts, int in_faces);
void admmpd_dealloc(ADMMPDInterfaceData*);
int admmpd_init(ADMMPDInterfaceData*);
int admmpd_cache_valid(ADMMPDInterfaceData*, int numVerts);
void admmpd_solve(ADMMPDInterfaceData*);
void admmpd_map_vertices(ADMMPDInterfaceData*, float (*vertexCos)[3], int numVerts); 

//void admmpd_solve(ADMMPDInterfaceData*);

// Copies the results of the solve (pos, vel) into BodyPoint
//void admmpd_to_bodypoint(
//    ADMMPD_Data *data,
//    BodyPoint *bp,
//    int numVerts);

#ifdef __cplusplus
}
#endif

#endif /* __ADMMPD_API_H__ */
