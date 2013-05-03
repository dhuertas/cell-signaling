#include <cmath>

#include "Manager.h"

using namespace std;

Define_Module(Manager);

Manager::~Manager() {}

void Manager::subscribe(Particle * p) {

	// Add the particle pointer to the particles vector
	particles.push_back(p);

	// Also add the particle to its corresponding space cell
	if (spaceCellSize > 0) {
		attachParticleToSpaceCell(p);
	}

}

void Manager::unsubscribe(Particle * p) {

	// Remove the particle pointer from the space cell structure
	detachParticleFromSpaceCell(p);

	// Remove the particle pointer from the particles vector
	particles.remove(p);

}

void Manager::initialize(int stage) {

	int N;
    double diameter;

	cPar *managerName;
	cModule *module;
	std::list<Particle *>::const_iterator p;
	std::stringstream buffer;

	// The manager node should be the first module to be initialized
	if (stage == 0) {

	    managerName = & simulation.getSystemModule()->par("managerName");
		
		// Set the name of the manager so we can have later access from
		// the other nodes.
		this->setName(managerName->stringValue());

		// Get the simulation space size
		spaceSizeX = simulation.getSystemModule()->par("spaceSizeX");
		spaceSizeY = simulation.getSystemModule()->par("spaceSizeY");
		spaceSizeZ = simulation.getSystemModule()->par("spaceSizeZ");

		tkEnvRefreshRate = simulation.getSystemModule()->par("refreshRate");

		// Set network size for tkenv
	    buffer << spaceSizeY;

	    module = simulation.getSystemModule();
	    module->getDisplayString().setTagArg("bgb", 0, buffer.str().c_str());
		buffer.str(std::string()); // clear buffer

		buffer << spaceSizeX;
	    module->getDisplayString().setTagArg("bgb", 1, buffer.str().c_str());
		buffer.str(std::string()); // clear buffer

		// We must wait to determine the space cell size till all the initial 
		// particles have been subscribed
		spaceCellSize = 0;

	} else if (stage == 1) {
		// the rest of the modules are being initialized ...
	} else if (stage == 2) {

		// All the particles are in the simulation space now. We can determine
		// the spaceCellSize
		for (p = particles.begin(); p != particles.end(); ++p) {
		    diameter = 2*(*p)->getRadius();

			if (spaceCellSize < diameter) {
				spaceCellSize = diameter;
			}

		}

		// Put every subscribed particle in its corresponding space cell
		Nx = ceil(spaceSizeX/spaceCellSize);
		Ny = ceil(spaceSizeY/spaceCellSize);
		Nz = ceil(spaceSizeZ/spaceCellSize);

		N = Nx*Ny*Nz;

		spaceCells.resize(N);

		for (p = particles.begin(); p != particles.end(); ++p) {
			attachParticleToSpaceCell(*p);
		}

		// Self message to refresh the tk environment
		scheduleAt(simTime()+tkEnvRefreshRate, new cMessage("refresh", EV_TKENVUPDATE));

	}

}

int Manager::numInitStages() const {

	return 3;

}

void Manager::attachParticleToSpaceCell(Particle *p) {

	int n, i, j, k;

	// The cell position for a given particle is a function of the form
	// n = f(x, y, z, spaceCellSize)
	i = floor(p->getX()/spaceCellSize);
	j = floor(p->getY()/spaceCellSize);
	k = floor(p->getZ()/spaceCellSize);

	n = i*(Nz*Ny) + j*(Nz) + k;

	p->setSpaceCell(n);
	spaceCells.at(n).push_back(p);

}

void Manager::detachParticleFromSpaceCell(Particle *p) {

	spaceCells.at(p->getSpaceCell()).remove(p);
	p->setSpaceCell(-1);

}

void Manager::transferParticle(Particle *p) {

	// Detach particle from its space cell
	detachParticleFromSpaceCell(p);

	// Obtain the particle position after the transfer event occured

	// Attach the particle to its new space cell
	attachParticleToSpaceCell(p);

}

void Manager::handleMessage(cMessage *msg) {

	if (msg->getKind() == EV_TKENVUPDATE) {

		tkEnvUpdateNetwork();

		// Self-message
		scheduleAt(simTime().dbl()+tkEnvRefreshRate, msg);

	}
}

void Manager::finish() {

}

void Manager::tkEnvUpdateNetwork() {

	std::list<Particle *>::const_iterator p;

	// Update particle positions
	for (p = particles.begin(); p != particles.end(); ++p) {

	    (*p)->tkEnvUpdatePosition(simTime().dbl());

	}

}
