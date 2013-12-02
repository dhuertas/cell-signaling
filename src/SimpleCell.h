#ifndef CELL_H
#define CELL_H

#include <string.h>
#include <iostream>

#include "base/Sphere.h"
#include "emitter/MoleculeEmitter.h"
#include "receiver/MoleculeReceiver.h"

#include "messages/TimeToLive_m.h"

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

		// cSimpleModule inheritance
		virtual void initialize(int);
		virtual int numInitStages() const;
		virtual void handleMessage(cMessage *);
		virtual void finish();

	friend class MoleculeEmitter;
	friend class MoleculeReceiver;

};

#endif
