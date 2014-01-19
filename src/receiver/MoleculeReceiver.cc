#include "MoleculeReceiver.h"

Define_Module(MoleculeReceiver);

MoleculeReceiver::MoleculeReceiver() : Receiver() {

}

MoleculeReceiver::~MoleculeReceiver() {

}

/*
 *
 * @param {integer} stage
 */
void MoleculeReceiver::initialize(int stage) {

	if (stage == 0) {
		// First stage manager initialization
		enabled = par("enabled");

	} else if (stage == 1) {
		// Particles initialization
	} else if (stage == 2) {
		// Second stage manager initialization
	} else if (stage == 3) {
		// Particle emitters and receivers initialization

		mobility = (SimpleCell *)getParentModule()->getSubmodule("mobility");

		if (enabled && getParentModule()->getSubmodule("emitter")->par("enabled")) {
			mobility->setParticleType(T_EMITTER_RECEIVER);
		} else if (enabled) {
			mobility->setParticleType(T_RECEIVER);
		}

		particlesReceivedVector.setName("particlesReceived");

		statsRefreshRate = par("statsRefreshRate");

		// Subscribe to manager
		setManager("manager");

		if (statsRefreshRate > 0) {
			scheduleAt(simTime() + statsRefreshRate/1000,
				new cMessage("refresh", EV_STATSUPDATE));
		}
	}

}

/*
 * Returns the number of initialization stages.
 */
int MoleculeReceiver::numInitStages() const {
	return 4;
}

/*
 *
 * @param {cMessage *} msg
 */
void MoleculeReceiver::handleMessage(cMessage *msg) {

	int kind = msg->getKind();
	double st;

	if (kind == EV_RECEIVE) {

		received++;		

	} else if (kind == EV_STATSUPDATE) {

		st = simTime().dbl();
		particlesReceivedVector.recordWithTimestamp(st, received);

		received = 0;

		if (statsRefreshRate > 0) {
			scheduleAt(st + statsRefreshRate/1000, msg);
		}

	}
}

/*
 *
 */
void MoleculeReceiver::finish() {

}

/*
 * Sets the manager attribute.
 * 
 * @param {string} param: the name of the manager module
 */
void MoleculeReceiver::setManager(std::string param) {

	try {
		manager = (Manager *)simulation.
			getSystemModule()->getSubmodule(param.c_str());
	} catch (cException *e) {
		EV << "setManager error" << "\n";
	}

}

/*
 *
 */
void MoleculeReceiver::registerReception(int particleType) {
	received++;
	EV << "Molecule received with type " << particleType << "\n";
	// TODO log the received particle type
}