#include "SimpleCell.h"

using namespace std;

Define_Module(SimpleCell);

SimpleCell::~SimpleCell() {

}

/*
 * Cell initialization
 */
void SimpleCell::initialize(int stage) {

	std::stringstream buffer;

	if (stage == 0) {
		// Manager module initializes during this stage
	} else if (stage == 1) {

		setParticleType(T_SIMPLECELL);

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
		setMass(par("radius").doubleValue());

		// Near-Neighbor List radius
		setListRadius(par("listRadius").doubleValue());

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
			scheduleAt(simTime() + statsRefreshRate/1000,
				new cMessage("refresh", EV_STATSUPDATE));
		}

		// update Cell position in the tk environment
		tkEnvUpdatePosition();

		// draw module shape in the tk environment
		tkEnvDrawShape();

	}

}

/*
 * Returns the number of initialization stages.
 */
int SimpleCell::numInitStages() const {
	return 2;
}

/*
 * Handles every message that the cell receives.
 *
 * @param msg pointer to a cMessage object
 */
void SimpleCell::handleMessage(cMessage *msg) {

	int kind = msg->getKind();

	if (isSignaling(msg)) {

		handleSignaling(msg);

	} else if (ISMOBILITY(kind)) {

		handleMobilityMessage(msg);

	} else if (kind == EV_STATSUPDATE) {

		if (statsRefreshRate > 0) {
			scheduleAt(simTime() + statsRefreshRate/1000, msg);
		}

	} else if (kind == EV_TTLEXPIRE) {

		expire();

	}

}

/*
 * Clean and close everything.
 */
void SimpleCell::finish() {

	// Unsubscribe from the manager
	getManager()->unsubscribe(this);

	// Delete mobility messages
	deleteMobilityMessages();

	if (timeToLive > 0) {
		cancelAndDelete(timeToLiveMsg);
	}

}

/*
 * Cell must expire. This function gets called when the timeToLive 
 * parameter is set and the event EV_TTLEXPIRE arrives.
 */
void SimpleCell::expire() {

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
 */
void SimpleCell::scheduleExpire(double time) {
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
 *
 */
bool SimpleCell::isSignaling(cMessage *msg) {

	int kind = msg->getKind();
	CollisionMessage *cmsg;

	cmsg = NULL;

	if (kind == EV_COLLISION) {

		cmsg = (CollisionMessage *)msg;

		if (particleType == T_RECEIVER ||
			particleType == T_EMITTER_RECEIVER) {

			if (cmsg->getPartner()->getParticleType() == T_SIGNALING) {
				return true;
			}
		}
	}

	return false;
}

/*
 * Handle the signaling process. This is a receiver and the other is a molecule.
 *
 * @param {cMessage *} msg
 */
void SimpleCell::handleSignaling(cMessage *msg) {

	MoleculeReceiver *receiver;
	Particle *p;

	CollisionMessage *cmsg;

	cmsg = (CollisionMessage *)msg;
	p = cmsg->getPartner();

	receiver = ((MoleculeReceiver *)getParentModule()->getSubmodule("receiver"));
	receiver->registerReception(p->getParticleType());

	// TODO change the following line, perhaps using gates
	((Molecule *)p)->scheduleExpire(cmsg->getCollisionTime());

	// Change the event type to CHECK so we can find the next event
	// for the receiver
	msg->setKind(EV_CHECK);
	handleMobilityMessage(msg);
}
