#ifndef SPHERE_H
#define SPHERE_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <coutvector.h>
#include <omnetpp.h>

#include "Manager.h"
#include "Circle.h"

#include "messages/Transfer_m.h"
#include "messages/Collision_m.h"
#include "messages/OutOfNeighborhood_m.h"

class Sphere : public Circle, public cSimpleModule {

	private:

		Manager *manager;

// Self messages
		TransferMessage *transferMsg;
		CollisionMessage *collisionMsg;
		OutOfNeighborhoodMessage *outOfNeighborhoodMsg;

	protected:

		cOutVector collisionVector;

	public:

		Sphere() : Circle() {};
		Sphere(double, double, double, double , double, double, double, double);

		void initMessages(void);

// Initialize the event queue
		void initEvents(void);

		void handleMobilityMessage(cMessage *);

		void handleTransfer(TransferMessage *);
		void handleCollision(CollisionMessage *);
		void handleWallCollision(CollisionMessage *);
		void handleOutOfNeighborhood(void);

		void adjustCollision(double, Particle *);

// Near-Neighbor List methods
		void createNearNeighborList(void);
		void updateNearNeighborList(void);

// tk Environment related methods
		void tkEnvDrawShape(void);
		void tkEnvUpdatePosition(void);
		void tkEnvUpdatePosition(double);

// Gets and sets
		Manager * getManager(void) { return manager; }
		TransferMessage * getTransferMessage(void);
		CollisionMessage * getCollisionMessage(void);

		void setManager(std::string);

};

#endif
