#ifndef SPHERE_H
#define SPHERE_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <omnetpp.h>

#include "Manager.h"
#include "Circle.h"

class Sphere : public Circle, public cSimpleModule {

	private:

		Manager *manager;

	protected:

	public:

		Sphere() : Circle() {};
		Sphere(double, double, double, double , double, double, double, double);

		// Find the next event time for a given particle
		void nextEventTime();

		cMessage * computeTransferTime(void);
		cMessage * computeCollisionTime(void);
		cMessage * computeWallCollisionTime(void);
		double solveCollisionTime(Particle *);

		void updateStateAfterCollision(cMessage *);
		void updateStateAfterWallCollision(cMessage *);

		// tk Environment related methods
		void tkEnvDrawShape(void);
		void tkEnvUpdatePosition(void);
		void tkEnvUpdatePosition(double);

		// Gets and sets
		Manager *getManager(void) { return manager; }
		void setManager(std::string);

};

#endif
