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

		lastCollisionTime = 0;

		// Subscribe to manager
		setManager("managerName");
		getManager()->subscribe(this);

		// update Molecule position in the tk environment
		tkEnvUpdatePosition();

		// draw module shape in the tk environment
		tkEnvDrawShape();
	} else {

	}

}

/*
 *
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

	if (kind == EV_CHECK) {

	} else if (kind == EV_WALLCOLLISION) {
		// Update the molecule data
		updateStateAfterWallCollision(msg);
		delete msg;

	}

	nextEventTime();
}

/*
 * Clean and close everything.
 */
void Molecule::finish() {
	
	// Unsubscribe from the manager
	getManager()->unsubscribe(this);

}
