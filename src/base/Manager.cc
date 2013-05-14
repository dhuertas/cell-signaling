#include <cmath>

#include "Manager.h"

using namespace std;

Define_Module(Manager);

Manager::~Manager() {}

/*
 * Subscribe particle to the simulation space
 *
 * @param {Particle *} p 
 */
void Manager::subscribe(Particle * p) {

	// Add the particle pointer to the particles vector
	particles.push_back(p);

	// Also add the particle to its corresponding space cell
	if (spaceCellSize > 0) {
		attachParticleToSpaceCell(p, -1);
	}

}

/*
 * Unsubscribe particle from the simulation space
 *
 * @param {Particle *} p
 */
void Manager::unsubscribe(Particle * p) {

	// Remove the particle pointer from the space cell structure
	detachParticleFromSpaceCell(p, -1);

	// Remove the particle pointer from the particles vector
	particles.remove(p);

}

/*
 * Initialize the simulation space and its modules.
 *
 * @param {Integer} stage
 */
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
		setNumberOfSpaceCellsX(ceil(spaceSizeX/spaceCellSize));
		setNumberOfSpaceCellsY(ceil(spaceSizeY/spaceCellSize));
		setNumberOfSpaceCellsZ(ceil(spaceSizeZ/spaceCellSize));

		N = Nx*Ny*Nz;

		spaceCells.resize(N);

		for (p = particles.begin(); p != particles.end(); ++p) {
			attachParticleToSpaceCell(*p, -1);
		}

		// Self message to refresh the tk environment
		scheduleAt(simTime() + tkEnvRefreshRate, 
			new cMessage("refresh", EV_TKENVUPDATE));

		// Make that every subscribed particle compute its next event time
		for (p = particles.begin(); p != particles.end(); ++p) {

			(*p)->firstEventTime();
		}
	}

}

/*
 * The manager must be initialized and act during the first and third stages.
 *
 * @return {const integer}
 */
int Manager::numInitStages() const {

	return 3;

}

/*
 * Attach a particle to a space cell. If the value of the second argument is 
 * negative means that the particle has been just subscribed and its space cell
 * must be calculated. 
 *
 * @param {Particle *} p
 * @param {integer} to
 */
void Manager::attachParticleToSpaceCell(Particle *p, int to) {

	int n, i, j, k;

	if (to == -1) {

		i = floor(p->getX()/spaceCellSize);
		j = floor(p->getY()/spaceCellSize);
		k = floor(p->getZ()/spaceCellSize);

		n = i*Ny*Nz + j*Nz + k;

		p->setSpaceCell(n);
		spaceCells.at(n).push_back(p);

	} else {

		spaceCells.at(to).push_back(p);

	}

}

/*
 * Detach a particle from a space cell.
 *
 * @param {Particle *} p
 * @param {integer} from
 */
void Manager::detachParticleFromSpaceCell(Particle *p, int from) {

	if (from == -1) {

		spaceCells.at(p->getSpaceCell()).remove(p);

	} else {

		spaceCells.at(from).remove(p);

	}

}

/*
 * Transfer a particle from one space cell to another.
 *
 * @param {Particle *} p
 * @param {integer} from
 * @param {integer} to
 */
void Manager::transferParticle(Particle *p, int from, int to) {

	// Detach particle from its space cell
	detachParticleFromSpaceCell(p, from);

	// Obtain the particle position after the transfer event occured

	// Attach the particle to its new space cell
	attachParticleToSpaceCell(p, to);

}

/*
 * Return a list containing the pointers to the particles in a given
 * space cell.
 *
 * @param {Integer} n space cell index
 * @return list of particle pointers
 */
std::list<Particle *> Manager::getSpaceCellParticles(int n) {

	return spaceCells.at(n);

}

/*
 * Handles every message that the manager module receives.
 *
 * @param {cMessage *} msg
 */
void Manager::handleMessage(cMessage *msg) {

	int kind = msg->getKind();

	if (kind == EV_TKENVUPDATE) {

		tkEnvUpdateNetwork();

		// Self-message
		scheduleAt(simTime().dbl() + tkEnvRefreshRate, msg);

	}

}

/*
 * Clean and close everything.
 */
void Manager::finish() {

}

/*
 * Updates the position of all the particles in the tk environment during
 * the simulation.
 */
void Manager::tkEnvUpdateNetwork() {

	std::list<Particle *>::const_iterator p;

	// Update particle positions
	for (p = particles.begin(); p != particles.end(); ++p) {

		(*p)->tkEnvUpdatePosition(simTime().dbl());

	}

}
