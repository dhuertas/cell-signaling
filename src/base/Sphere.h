#ifndef SPHERE_H
#define SPHERE_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <coutvector.h>
#include <omnetpp.h>

#include "Manager.h"
#include "Circle.h"
#include "../messages/Mobility_m.h"

class Sphere : public Circle, public cSimpleModule {

	private:

		Manager *manager;

		MobilityMessage *transferMsg;
		MobilityMessage *collisionMsg;
		MobilityMessage *wallCollisionMsg;
		MobilityMessage *outOfNeighborhoodMsg;

	protected:

		cOutVector collisionVector;

	public:

		Sphere() : Circle() {};
		Sphere(double, double, double, double , double, double, double, double);

		// Initialize the event queue
		void firstEventTime(void);

		// Find the next event time for a given particle
		void nextEventTime(void);

		double scheduledCollisionTime(void);

		void handleMobilityMessage(MobilityMessage *);

		void updateStateAfterTransfer(MobilityMessage *);
		void updateStateAfterCollision(MobilityMessage *);
		void updateStateAfterWallCollision(MobilityMessage *);

		// Near-Neighbor List methods
		void createNearNeighborList(void);
		void updateNearNeighborList(void);

		// tk Environment related methods
		void tkEnvDrawShape(void);
		void tkEnvUpdatePosition(void);
		void tkEnvUpdatePosition(double);

		// Gets and sets
		Manager *getManager(void) { return manager; }
		MobilityMessage *getScheduledMobilityMessage(void);

		void setManager(std::string);
};

#endif
