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

class Molecule : public Sphere {

	private:

		Manager *manager;

	protected:

	public:

		~Molecule();

		// cSimpleModule inheritance
		virtual void initialize(int);
		virtual int numInitStages(void) const;
		virtual void handleMessage(cMessage *);
		virtual void finish();

};

#endif
