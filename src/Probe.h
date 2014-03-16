#ifndef PROBE_H
#define PROBE_H

#include <cmessage.h>
#include <csimplemodule.h>
#include <cqueue.h>
#include <omnetpp.h>

#include <cmath>
#include <vector>
#include <list>

#include "base/Defines.h"

class Probe : public cSimpleModule {

	private:

	protected:

		point_t position;

		double radius;

	public:

		~Probe();

		//
		// cSimpleModule inheritance
		//
		virtual void initialize(int stage);

		virtual int numInitStages() const;

		virtual void handleMessage(cMessage *);

		virtual void finish();
};

#endif