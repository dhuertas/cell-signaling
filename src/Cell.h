#ifndef CELL_H
#define CELL_H

#include <string.h>
#include <iostream>

#include "base/Sphere.h"
#include "messages/TimeToLive_m.h"
class Cell : public Sphere {

	private:

		double emitEvery;   // Dumb var to test simulation
		int emitCount;

		// Amount of time that will pass before it expires.
		double timeToLive;

		// Update the cOutVectors periodically.
		double statsRefreshRate;

		// Self messages
		TimeToLiveMessage *timeToLiveMsg;

	protected:

	public:

		~Cell();

		// cSimpleModule inheritance
		virtual void initialize(int);
		virtual int numInitStages() const;
		virtual void handleMessage(cMessage *);
		virtual void finish();
		
};

#endif
