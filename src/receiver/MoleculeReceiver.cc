//  This file is part of the Cell-Signaling project. Cell-Signaling is an
//  Omnet++ project to simulate cell signaling communications.
//  Copyright (C) 2014  Daniel Huertas
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "MoleculeReceiver.h"

Define_Module(MoleculeReceiver);

MoleculeReceiver::MoleculeReceiver() : Receiver() {

	received = 0;

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

		// Initialize statistics
		received = 0;

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

		statsRefreshRate_ = par("statsRefreshRate");

		// Set the reference to the manager module
		setManager("manager");

		if (statsRefreshRate_ > 0) {
			scheduleAt(simTime() + statsRefreshRate_,
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

		if (statsRefreshRate_ > 0) {
			scheduleAt(st + statsRefreshRate_, msg);
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
		manager_ = (Manager *)simulation.
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