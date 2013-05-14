#ifndef DEFINES_H
#define DEFINES_H

#define NO_TIME             -1

#define EV_COLLISION    	1
#define EV_CHECK        	2
#define EV_TRANSFER         3
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
