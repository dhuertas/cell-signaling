#ifndef MOLECULEEMITTER_H
#define MOLECULEEMITTER_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <coutvector.h>
#include <omnetpp.h>

#include "Emitter.h"

class MoleculeEmitter : public Emitter, public cSimpleModule {

	private:

	protected:

		double emissionStart;

		double emissionDuration;

		double emissionRate;

		uint8_t emissionFunction;

	public:

		MoleculeEmitter();

		~MoleculeEmitter();

		virtual void initialize();

		virtual void handleMessage(cMessage *);

		virtual void finish();

};

#endif