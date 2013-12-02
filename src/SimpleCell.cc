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
 *
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

    if (ISMOBILITY(kind)) {

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