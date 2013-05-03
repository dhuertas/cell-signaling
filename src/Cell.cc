#include "Cell.h"

using namespace std;

Define_Module(Cell);

Cell::~Cell() {

}

void Cell::initialize(int stage) {

	cPar *managerName;
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
		managerName = & simulation.getSystemModule()->par("managerName");
		manager = (Manager *)simulation.getSystemModule()
				->getSubmodule(managerName->stringValue());
		
		if (manager != NULL) {
			manager->subscribe(this);
		} else {
			// TODO stop simulation
		}
	
		// update Cell position in the tk environment
		tkEnvUpdatePosition();

		// draw module shape in the tk environment
		tkEnvDrawShape();


	}

	// Get the time of the first collision

}

int Cell::numInitStages() const {
	return 2;
}

void Cell::handleMessage(cMessage *msg) {

	if (msg->getKind() == 1 && emitEvery > 0) {  // kind 1: emit molecule

		emitCount++;

		if (emitCount == emitEvery) {
			emitCount = 0;

			EV << "Molecule emitted";
		}

		delete msg;

		cMessage * internalMsg = new cMessage("emit", 1);

		scheduleAt(simTime()+10, internalMsg);
	}

}

void Cell::finish() {

	// Unsubscribe from the manager
	manager->unsubscribe(this);

	// All events related to this cell should be discarded
}
