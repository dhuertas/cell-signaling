#ifndef COMPUTE_H
#define COMPUTE_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <omnetpp.h>

#include "Manager.h"
#include "Particle.h"
#include "../messages/Mobility_m.h"

class Mobility {

	public:

		static void transferTime(MobilityMessage *, Particle *);
		static void outOfNeighborhoodTime(MobilityMessage *, Particle *);
};

#endif