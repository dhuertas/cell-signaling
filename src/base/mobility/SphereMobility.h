#ifndef COMPUTESPHERE_H
#define COMPUTESPHERE_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <omnetpp.h>

#include "Mobility.h"
#include "Manager.h"
#include "Sphere.h"

#include "../messages/Transfer_m.h"
#include "../messages/Collision_m.h"
#include "../messages/OutOfNeighborhood_m.h"

class SphereMobility : public Mobility {

	public:

		static double nextCollision(CollisionMessage *,int, Sphere *);
		static double nextWallCollision(CollisionMessage *, Sphere *);
		static double solveCollision(Particle *, Particle *);

};

#endif
