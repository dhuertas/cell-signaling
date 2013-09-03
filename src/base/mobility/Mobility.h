#ifndef COMPUTE_H
#define COMPUTE_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <omnetpp.h>

#include "Manager.h"
#include "Particle.h"

#include "../messages/Transfer_m.h"
#include "../messages/Collision_m.h"
#include "../messages/OutOfNeighborhood_m.h"

class Mobility {

	public:

		static double nextTransfer(TransferMessage *, Particle *);
		static double outOfNeighborhoodTime(OutOfNeighborhoodMessage *, Particle *);

		static void resetCollisionMessage(CollisionMessage *);

};

#endif
