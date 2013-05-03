#ifndef MOLECULE_H
#define MOLECULE_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <omnetpp.h>
#include <string.h>
#include <iostream>

#include "Manager.h"
#include "base/Sphere.h"

class Molecule : public Sphere {

	private:

		int identifier;
		Manager *manager;

	protected:

	public:

	~Molecule();
	int getIdentifier(void) { return identifier; };
	void setIdentifer(int id) { identifier = id; };

	// cSimpleModule inheritance
	virtual void initialize(int);
	virtual int numInitStages(void) const;
	virtual void handleMessage(cMessage *);
	virtual void finish();

};

#endif
