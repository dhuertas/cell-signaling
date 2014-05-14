//  Omnet++ project to simulate cell signaling communications 
//  Copyright (C) 2014  Daniel Huertas
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef DEFINES_H
#define DEFINES_H

#define NO_TIME					-1

// Mobility Event Types
#define EV_NONE					0
#define EV_COLLISION			1
#define EV_CHECK				2
#define EV_TRANSFER				3
#define EV_BOUNDARYCOLLISION	4
#define EV_OUTOFNEIGHBORHOOD	5
#define EV_TKENVUPDATE			6
#define EV_STATSUPDATE			7
#define EV_TTLEXPIRE			8
#define EV_EMIT					9
#define EV_RECEIVE				10
#define EV_BROWNIAN				11

// Particle types
#define T_MOLECULE				1
#define T_SIMPLECELL			2
#define T_EMITTER				3
#define T_RECEIVER				4
#define T_EMITTER_RECEIVER		5
#define T_SIGNALING				6

// Collision algorithm modes
#define M_CELLLIST				1	// Only Cell Lists are used
#define M_NNLIST				2	// A combination of Near-Neighbor List and Cell List methods is used

// Boundary modes
#define BM_ELASTIC				1
#define BM_EXPIRE				2
#define BM_PERIODIC				3

// Status
#define STATUS_OFF				0
#define STATUS_ON				1

#define BUFFER_LENGTH 	1048576 // 1024*1024
#define READ			0
#define WRITE			1

#define ISMOBILITY(kind)			\
	kind == EV_COLLISION ||			\
	kind == EV_CHECK ||				\
	kind == EV_TRANSFER ||			\
	kind == EV_BOUNDARYCOLLISION ||	\
	kind == EV_OUTOFNEIGHBORHOOD || \
	kind == EV_BROWNIAN

#define CELLBELONGSTOSIMSPACE(a,b,c,Nx,Ny,Nz)	\
	0 <= a && a < Nx &&							\
	0 <= b && b < Ny &&							\
	0 <= c && c < Nz

typedef struct Point3D {
	double x;
	double y;
	double z;
} point_t;

typedef struct Vector3D {
	double x;
	double y;
	double z;
} vect_t;

typedef struct CellIndex3D {
	unsigned char im; // index mask
		// xxxxxxx. spare
		// .......x whether the particle is outside of sim space (0) or inside (1)
	unsigned int i; // space cell index along x axis
	unsigned int j; // space cell index along y axis
	unsigned int k; // space cell index along z axis
} index_t;

typedef struct Statistics {
	unsigned long long allCollisions;
	unsigned long long particleCollisions;
	unsigned long long wallCollisions;
	unsigned long long transfers;
	unsigned long long expires;
} statistics_t;

// Space cell side points that allows to obtain the side equation and compute 
// the transfer time.
extern char sp[10*6];

#endif
