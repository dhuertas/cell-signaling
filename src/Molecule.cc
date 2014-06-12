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

#include "Molecule.h"

using namespace std;

Define_Module(Molecule);

Molecule::~Molecule() {}

/*
 *
 */
void Molecule::initialize(int stage) {

	std::stringstream buffer;

	if (stage == 0) {
		// Manager module initializes during this stage
	} else if (stage == 1) {

		active = true;

		setParticleType(T_MOLECULE);

		// Initial position
		setX(par("xpos").doubleValue());
		setY(par("ypos").doubleValue());
		setZ(par("zpos").doubleValue());

		// Velocity
		setVx(par("vx").doubleValue());
		setVy(par("vy").doubleValue());
		setVz(par("vz").doubleValue());

		// Cell radius
		setRadius(par("radius").doubleValue());
		setMass(par("mass").doubleValue());


		// Subscribe to manager
		setManager("manager");
		getManager()->subscribe(this);

		// Brownian Motion parameters
		setDiffusion(par("diffusion").doubleValue());
		setInertia(par("inertia").doubleValue());
		setViscosity(par("viscosity").doubleValue());

		// Compute Brownian Motion Standard Deviation
		//        _____________
		//  \    /             |
		//   \  /  4*M_PI*D*dt
		//    \/
		setBMStdDev(sqrt(4*M_PI*par("diffusion").doubleValue()*(manager->getDeltaTime())));

		// Near-Neighbor List radius
		setListRadius(par("listRadius").doubleValue());

		// Near Neighbor List refresh list radius
		setRefreshListRadius(par("refreshListRadius").doubleValue());

		boundariesMode = par("boundariesMode");

		timeToLive = par("timeToLive");

		if (timeToLive > 0) {
			timeToLiveMsg = new TimeToLiveMessage("expire", EV_TTLEXPIRE);
			scheduleAt(simTime() + timeToLive, timeToLiveMsg);
		}

		logCollisions = par("logCollisions");

		if (logCollisions > 0) {
			collisionTimeVector = new cOutVector("collisionTime");
			xCollisionPositionVector = new cOutVector("xCollisionPosition");
			yCollisionPositionVector = new cOutVector("yCollisionPosition");
			zCollisionPositionVector = new cOutVector("zCollisionPosition");
		}

		logPosition = par("logPosition");

		if (logPosition > 0) {
			xPositionVector = new cOutVector("xPosition");
			yPositionVector = new cOutVector("yPosition");
			zPositionVector = new cOutVector("zPosition");
		}

		statsRefreshRate = par("statsRefreshRate");

		if (statsRefreshRate > 0) {
			scheduleAt(simTime() + statsRefreshRate/1000,
				new cMessage("refresh", EV_STATSUPDATE));
		}

		// update Molecule position in the tk environment
		tkEnvUpdatePosition();

		// draw module shape in the tk environment
		tkEnvDrawShape();

	} else {

	}

}

/*
 * The molecules must be initialized during the second stage.
 */
int Molecule::numInitStages() const {

	return 2;

}

/*
 * Handles every message that the molecule receives.
 *
 * @param msg pointer to a cMessage object
 */
void Molecule::handleMessage(cMessage *msg) {

	int kind = msg->getKind();

	if (isSignaling(msg)) {

		handleSignaling(msg);

	} else if (ISMOBILITY(kind)) {

		handleMobilityMessage(msg);

	} else if (kind == EV_STATSUPDATE) {

		double st = simTime().dbl();
		double dt = st - lastCollisionTime;

		// Put the statistics logged so far to cout vectors
		xPositionVector->recordWithTimestamp(st, position.x + velocity.x * dt);
		yPositionVector->recordWithTimestamp(st, position.y + velocity.y * dt);
		zPositionVector->recordWithTimestamp(st, position.z + velocity.z * dt);

		if (statsRefreshRate > 0) {
			scheduleAt(st + statsRefreshRate/1000, msg);
		}

	} else if (kind == EV_TTLEXPIRE) {
		expire();
	}

}

/*
 * Clean and close everything.
 */
void Molecule::finish() {

	// Unsubscribe from the manager
	getManager()->unsubscribe(this);

	// Delete mobility messages
	deleteMobilityMessages();

	if (timeToLive > 0) {
		cancelAndDelete(timeToLiveMsg);
	}

}

/*
 * Molecule must expire. This function gets called when the timeToLive 
 * parameter is set and the event EV_TTLEXPIRE arrives.
 */
void Molecule::expire() {
	// Methods called from other modules must have this macro
	Enter_Method_Silent();

	active = false;

	manager->registerExpire();

	// Unsubscribe from the manager
	getManager()->unsubscribe(this);

	// Get out of the simulation space gracefully
	this->finishMobility();

	// Delete mobility messages
	deleteMobilityMessages();

	if (timeToLive > 0) {
		cancelAndDelete(timeToLiveMsg);
	}

	this->deleteModule();

}

/*
 *
 * @param {double} time: simulation time to schedule the expire message
 */
void Molecule::scheduleExpire(double time) {
	// Methods called from other modules must have this macro
	Enter_Method_Silent();

	if (timeToLive > 0) {
		if (timeToLiveMsg->isScheduled()) {
			cancelEvent(timeToLiveMsg);
		}
	} else {
		timeToLiveMsg = new TimeToLiveMessage("expire", EV_TTLEXPIRE);
	}

	scheduleAt(time, timeToLiveMsg);

}

/*
 * Tells whether the molecule has to process a collision as a signaling.
 *
 * @param {cMessage *} msg
 * @return {bool}
 */
bool Molecule::isSignaling(cMessage *msg) {

	int partnerParticleType;

	int kind = msg->getKind();
	CollisionMessage *cmsg;

	cmsg = NULL;
	partnerParticleType = 0;

	if (kind == EV_COLLISION) {

		cmsg = (CollisionMessage *)msg;

		if (particleType == T_SIGNALING) {

			partnerParticleType = cmsg->getPartner()->getParticleType();

			if (partnerParticleType == T_RECEIVER ||
				partnerParticleType == T_EMITTER_RECEIVER) {
				return true;
			}
		}
	}

	return false;
}

/*
 * Handle the signaling process. This is a molecule and the other is a receiver.
 *
 * @param {cMessage *} msg
 */
void Molecule::handleSignaling(cMessage *msg) {

	MoleculeReceiver *receiver;
	Particle *p;

	CollisionMessage *cmsg;

	cmsg = (CollisionMessage *)msg;
	p = cmsg->getPartner();

	receiver = (MoleculeReceiver *)((SimpleCell *)p)->getParentModule()
		->getSubmodule("receiver");

	receiver->registerReception(getParticleType());

	expire();
}
