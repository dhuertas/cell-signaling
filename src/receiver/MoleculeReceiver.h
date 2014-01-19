#ifndef MOLECULERECEIVER_H
#define MOLECULERECEIVER_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <coutvector.h>
#include <omnetpp.h>

#include "Receiver.h"
#include "../SimpleCell.h"
#include "../base/Manager.h"

// Forward declarations
class SimpleCell;

class MoleculeReceiver : public Receiver, public cSimpleModule {

	private:

		Manager *manager;

		SimpleCell *mobility;

		// Update the cOutVectors periodically.
		double statsRefreshRate;

	protected:
	
		unsigned long long received;

		cOutVector particlesReceivedVector; 

	public:

		MoleculeReceiver();

		~MoleculeReceiver();

		void setManager(std::string);

		void registerReception(int particleType);

		//
		// cSimpleModule inheritance
		//
		virtual void initialize(int);

		virtual int numInitStages() const;

		virtual void handleMessage(cMessage *);

		virtual void finish();
};

#endif
