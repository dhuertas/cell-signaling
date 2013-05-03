#ifndef DEFINES_H
#define DEFINES_H

#define EV_COLLISION    	1
#define EV_TRANSFER     	2
#define EV_CHECK        	3
#define EV_WALLCOLLISION	4
#define EV_TKENVUPDATE		5

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

#endif
