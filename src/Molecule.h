#ifndef MOLECULE_H
#define MOLECULE_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <coutvector.h>
#include <omnetpp.h>
#include <string.h>

#include <iostream>

#include "base/Sphere.h"
#include "messages/TimeToLive_m.h"

class Molecule : public Sphere {

	private:

		// Amount of time that will pass before it expires.
		double timeToLive;

		// Update the cOutVectors periodically.
		double statsRefreshRate;

		// Self messages
		TimeToLiveMessage *timeToLiveMsg;

	protected:

		// Vectors to track particle position over time.
		cOutVector *xPositionVector;
		cOutVector *yPositionVector;
		cOutVector *zPositionVector;

	public:

		~Molecule();

		// cSimpleModule inheritance
		virtual void initialize(int);
		virtual int numInitStages(void) const;
		virtual void handleMessage(cMessage *);
		virtual void finish();

};

#endif
