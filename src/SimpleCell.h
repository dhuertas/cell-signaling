#ifndef CELL_H
#define CELL_H

#include <string.h>
#include <iostream>

#include "base/Sphere.h"
#include "receiver/MoleculeReceiver.h"
#include "messages/TimeToLive_m.h"
#include "Molecule.h"

class SimpleCell : public Sphere {

	private:

		// Amount of time that will pass before it expires.
		double timeToLive;

		// Update the cOutVectors periodically.
		double statsRefreshRate;

		// Self messages
		TimeToLiveMessage *timeToLiveMsg;

	protected:

	public:

		~SimpleCell();

		void expire();

		void scheduleExpire(double time);

		bool isSignaling(cMessage *msg);

		void handleSignaling(cMessage *msg);

		//
		// cSimpleModule inheritance
		//
		virtual void initialize(int);

		virtual int numInitStages() const;

		virtual void handleMessage(cMessage *msg);

		virtual void finish();

};

#endif
