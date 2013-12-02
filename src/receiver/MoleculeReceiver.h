#ifndef MOLECULERECEIVER_H
#define MOLECULERECEIVER_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <coutvector.h>
#include <omnetpp.h>

#include "Receiver.h"

class MoleculeReceiver : public Receiver, public cSimpleModule {

	private:
	
	protected:
	
	public:

		MoleculeReceiver();

		~MoleculeReceiver();

		virtual void initialize();

		virtual void handleMessage(cMessage *);

		virtual void finish();

};

#endif