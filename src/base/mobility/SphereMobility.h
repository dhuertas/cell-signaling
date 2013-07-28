#ifndef COMPUTESPHERE_H
#define COMPUTESPHERE_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <omnetpp.h>

#include "Mobility.h"
#include "Manager.h"
#include "Sphere.h"
#include "../messages/Mobility_m.h"

class SphereMobility : public Mobility {

	public:

		static void collisionTime(MobilityMessage *, Sphere *);
		static void wallCollisionTime(MobilityMessage *, Sphere *);
		static double solveCollisionTime(Particle *, Particle *);

};

#endif