#ifndef DEFINES_H
#define DEFINES_H

#define NO_TIME             	-1

// Mobility Event Types
#define EV_COLLISION			1
#define EV_CHECK				2
#define EV_TRANSFER				3
#define EV_WALLCOLLISION		4
#define EV_OUTOFNEIGHBORHOOD	5
#define EV_TKENVUPDATE			6

// Collision Dynamics Algorithm
#define CD_CELLLIST				1	// Only Cell Lists are used
#define CD_NNLIST				2 	// A combination of Near-Neighbor List and
									// Cell List methods is used
typedef struct Point {
	double x;
	double y;
	double z;
} point_t;

typedef struct Vector3D {
	double x;
	double y;
	double z;
} vect_t;

// Space cell side points that allows to obtain the side equation and compute 
// the transfer time.
extern char sp[10*6];

#endif
