#include "Cell.h"

using namespace std;

Define_Module(Cell);

Cell::~Cell() {

}

/*
 *
 */
void Cell::initialize(int stage) {

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

		lastCollisionTime = 0;

		// Subscribe to manager		
		setManager("managerName");
		getManager()->subscribe(this);
	
		// update Cell position in the tk environment
		tkEnvUpdatePosition();

		// draw module shape in the tk environment
		tkEnvDrawShape();


	}

	// Get the time of the first collision

}

/*
 *
 */
int Cell::numInitStages() const {
	return 2;
}

/*
 * Handles every message that the cell receives.
 *
 * @param msg pointer to a cMessage object
 */
void Cell::handleMessage(cMessage *msg) {

    int kind = msg->getKind();

// Step 2. Handle the event
    if (kind == EV_TRANSFER) {
        // Update the molecule space cell
        updateStateAfterTransfer((MobilityMessage *)msg);
        nextEventTime();

    } else if (kind == EV_WALLCOLLISION) {
        // Update the molecule data
        updateStateAfterWallCollision((MobilityMessage *)msg);
        nextEventTime();

    } else if (kind == EV_COLLISION) {

        updateStateAfterCollision((MobilityMessage *)msg);
        nextEventTime();

    } else if (kind == EV_CHECK) {
// TODO our collision time has turned invalid. We must check again for the next
// event.

// Two cases are possible:
// 1. We obtained an expected collision time, but the partner has a smaller
// collision time. Therefore we can calculate next collision at the expected
// collision time.
//
// 2. Someone else has forced us to recompute our collision time since it
// expects a collision. If we have a scheduled collision event, we must cancel
// it telling the partner to check again for its next collision time (thus
// going to case 1).
        nextEventTime();

    } else {

    }

}

/*
 * Clean and close everything.
 */
void Cell::finish() {

	// Unsubscribe from the manager
	getManager()->unsubscribe(this);

	// All events related to this cell should be discarded
}
