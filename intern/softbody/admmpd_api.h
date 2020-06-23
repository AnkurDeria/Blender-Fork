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

#ifndef ADMMPD_API_H
#define ADMMPD_API_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ADMMPDInterfaceData {
    // totverts is usually different than mesh_totverts.
    // This is due to the lattice/tetmesh that is generated
    // in init. You can use them as input if reading from cache,
    // as they will be copied to internal solver data before admmpd_solve.
    int totverts; // number of deformable verts (output)
    int mesh_totverts; // number of surface mesh vertices (input)
    int mesh_totfaces; // number of surface mesh faces (input)
    int init_mode; // 0=tetgen, 1=lattice
    // Solver data used internally
    struct ADMMPDInternalData *idata;
} ADMMPDInterfaceData;

// SoftBody bodypoint (contains pos,vec)
typedef struct BodyPoint BodyPoint;

// Clears all solver data and ADMMPDInterfaceData
void admmpd_dealloc(ADMMPDInterfaceData*);

// Initializes solver and allocates internal data
int admmpd_init(ADMMPDInterfaceData*, float *in_verts, unsigned int *in_faces);

// Copies BodyPoint data (from SoftBody)
// to internal vertex position and velocity
void admmpd_copy_from_bodypoint(ADMMPDInterfaceData*, const BodyPoint *pts);

// Copies internal vertex position and velocity data
// to BodyPoints (from SoftBody) AND surface mesh vertices.
// If pts or vertexCos is null, its skipped
void admmpd_copy_to_bodypoint_and_object(
    ADMMPDInterfaceData*,
    BodyPoint *pts,
    float (*vertexCos)[3]);

// Copies out_verts and out_verts to internal data
// Performs solve over the time step
// Copies internal data to out_verts and out_vel
void admmpd_solve(ADMMPDInterfaceData*);

#ifdef __cplusplus
}
#endif

#endif // ADMMPD_API_H
