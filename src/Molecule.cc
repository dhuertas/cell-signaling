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

		// Near-Neighbor List radius
		setListRadius(par("listRadius").doubleValue());

		// Near Neighbor List refresh list radius
		setRefreshListRadius(par("refreshListRadius").doubleValue());

		// Subscribe to manager
		setManager("manager");
		getManager()->subscribe(this);

		boundariesMode = par("boundariesMode");

		timeToLive = par("timeToLive");

		if (timeToLive > 0) {

			timeToLiveMsg = new TimeToLiveMessage("expire", EV_TTLEXPIRE);

			scheduleAt(simTime() + timeToLive, timeToLiveMsg);

		}

		statsRefreshRate = par("statsRefreshRate");

		if (statsRefreshRate > 0) {
			xPositionVector = new cOutVector("x position");
			yPositionVector = new cOutVector("y position");
			zPositionVector = new cOutVector("z position");
		}

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