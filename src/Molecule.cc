#include "Molecule.h"

using namespace std;

Define_Module(Molecule);

Molecule::~Molecule() {}

void Molecule::initialize(int stage) {

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

		// update Molecule position in the tk environment
		tkEnvUpdatePosition();

		// draw module shape in the tk environment
		tkEnvDrawShape();
	} else {

	}

}

int Molecule::numInitStages() const {

	return 2;

}

void Molecule::handleMessage(cMessage *msg) {

}

void Molecule::finish() {
	
	// Unsubscribe from the manager
	manager->unsubscribe(this);

}
