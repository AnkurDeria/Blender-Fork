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

typedef struct ADMMPD_Data ADMMPD_Data;
typedef struct Object Object;
typedef struct BodyPoint BodyPoint;

ADMMPD_Data* admmpd_init(
    BodyPoint *bp,
    int numVerts);

void admmpd_cleanup(ADMMPD_Data*);

void admmpd_solve(ADMMPD_Data*);

// Copies the results of the solve (pos, vel) into BodyPoint
void admmpd_to_bodypoint(
    ADMMPD_Data *data,
    BodyPoint *bp,
    int numVerts);

#ifdef __cplusplus
}
#endif

#endif /* __ADMMPD_API_H__ */
