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
#include "messages/BrownianMotion_m.h"

class Sphere : public Circle, public cSimpleModule {

	private:

		// Self messages
		TransferMessage *transferMsg;
		CollisionMessage *collisionMsg;
		OutOfNeighborhoodMessage *outOfNeighborhoodMsg;
		BrownianMotionMessage *brownianMotionMsg;

	protected:

		Manager *manager;

		// Log collisions
		int logCollisions;

		// Log position
		int logPosition;

		// Track time between collisions
		cOutVector *collisionTimeVector;

		// Vectors to get the mean free path
		cOutVector *xCollisionPositionVector;
		cOutVector *yCollisionPositionVector;
		cOutVector *zCollisionPositionVector;

	public:

		Sphere() : Circle() {};
		Sphere(double, double, double, double , double, double, double, double);

		void initMobilityMessages(void);
		void deleteMobilityMessages(void);

		// Initialize the event queue
		void initializeMobility(void);
		void finishMobility(void);
		void finishMobility(Particle *from);

		void handleMobilityMessage(cMessage *);

		void handleTransfer(TransferMessage *);
		void handleCollision(CollisionMessage *);
		void handleBoundaryCollision(CollisionMessage *);
		void handleWallCollision(CollisionMessage *);
		void handleBrownianMotion(BrownianMotionMessage *);
		void handleOutOfNeighborhood(void);

		void adjustCollision(double, Particle *);

		// Near-Neighbor List methods
		void createNearNeighborList(void);
		void updateNearNeighborList(void);

		// cOutVector methods
		void logCollisionTime(double stime) {
		    double st = simTime().dbl();
			if (logCollisions && collisionTimeVector != NULL) {
				collisionTimeVector->recordWithTimestamp(st, stime);
			}

			if (logCollisions && 
				xCollisionPositionVector != NULL &&
				yCollisionPositionVector != NULL &&
				zCollisionPositionVector != NULL) {

				xCollisionPositionVector->recordWithTimestamp(st, position.x);
				yCollisionPositionVector->recordWithTimestamp(st, position.y);
				zCollisionPositionVector->recordWithTimestamp(st, position.z);
			}
		}

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
