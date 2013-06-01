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

    if (strcmp(msg->getName(), "mobility") == 0) {

        handleMobilityMessage((MobilityMessage *)msg);

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
