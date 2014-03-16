#ifndef MOLECULEEMITTER_H
#define MOLECULEEMITTER_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <coutvector.h>
#include <omnetpp.h>
#include <random.h>

#include <list>
#include <vector>

#include "Emitter.h"
#include "SimpleCell.h"
#include "Manager.h"

class MoleculeEmitter : public Emitter, public cSimpleModule {

	private:

		Manager *manager;

		SimpleCell *mobility;

		bool checkOverlap(point_t, double);

	protected:

		double emissionStart;

		double emissionDuration;

		double emissionRate;

		double emissionParticleRadius;

		double emissionParticleMass;

		double emissionTimeToLive;

		uint8_t emissionFunction;

		int emissionBoundariesMode;

		double emissionVelocity;

		double emissionListRadius;

		double emissionRefreshListRadius;

	public:

		MoleculeEmitter();
		~MoleculeEmitter();

		virtual void initialize(int);
		virtual int numInitStages() const;
		virtual void handleMessage(cMessage *);
		virtual void finish();

		Molecule * createMolecule();

		void setManager(std::string);
};

#endif
