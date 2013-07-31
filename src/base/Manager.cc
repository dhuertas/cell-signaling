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

	int N, m;
	int i, j, k;

    unsigned int count;

	double diameter;

	// cPar *managerName;
	cModule *module;
	std::list<Particle *>::const_iterator p;
	std::stringstream buffer;

// Initialize variables
	N = 0;
	m = 0;
	count = 0;

	diameter = 0.0;

// The manager node should be the first module to be initialized
	if (stage == 0) {

// Set whether it is a performance simulation
		setPerformanceSimulation(simulation.getSystemModule()->par("performanceSimulation"));

// Set the mode of operation for the molecule dynamics 
		mode = par("mode");

// Set the name of the manager so we can have later access from
// the other nodes.
		this->setName(par("name").stringValue());

// Get the simulation space size
		setSpaceSizeX(simulation.getSystemModule()->par("spaceSizeX"));
		setSpaceSizeY(simulation.getSystemModule()->par("spaceSizeY"));
		setSpaceSizeZ(simulation.getSystemModule()->par("spaceSizeZ"));

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

// All the particles are in the simulation space now. We can determine the 
// spaceCellSize
		count = 0;

		for (p = particles.begin(); p != particles.end(); ++p) {

			diameter = 2*(*p)->getRadius();

			if (spaceCellSize < diameter) {
				spaceCellSize = diameter;
			}

// Also set the simulation mode (Cell List or Near-Neighbor List)
			(*p)->setMode(mode);
			(*p)->setIdentifier(count);

			count++;

		}

// Put every subscribed particle in its corresponding space cell
		setNumberOfSpaceCellsX(ceil(spaceSizeX/spaceCellSize));
		setNumberOfSpaceCellsY(ceil(spaceSizeY/spaceCellSize));
		setNumberOfSpaceCellsZ(ceil(spaceSizeZ/spaceCellSize));

		N = Nx*Ny*Nz;

		spaceCells.resize(N);

// Piece of code for the perfomance tests. Set each particle position.
		count = 0;
		m = (int)round(pow(particles.size(), 1/3.0));
		i = 0; j = 0; k = 0;

		if (performanceSimulation) {
			
			for (p = particles.begin(); p != particles.end(); ++p) {

				(*p)->setX((0.5 + i)*spaceSizeX/m);
				(*p)->setY((0.5 + j)*spaceSizeY/m);
				(*p)->setZ((0.5 + k)*spaceSizeZ/m);

				k++;

				if (k%m == 0) {

					k = 0; j++;

					if (j%m == 0) {

						j = 0; i++;

						if (i%m == 0) {
							i = 0;
						}

					}

				}

			}

		}

		switch (mode) {

			case M_NNLIST:
// To start with, set the listRadius equal to the spaceCellSize
				for (p = particles.begin(); p != particles.end(); ++p) {
					(*p)->setListRadius(spaceCellSize);
					attachParticleToSpaceCell(*p, -1);
				}

				break;

			case M_CELLLIST:
			default:

				for (p = particles.begin(); p != particles.end(); ++p) {
					attachParticleToSpaceCell(*p, -1);
				}

				break;
		}

		if (mode == M_NNLIST) {
			for (p = particles.begin(); p != particles.end(); ++p) {
				(*p)->createNearNeighborList();
			}
		}

// Self message to refresh the tk environment
		scheduleAt(simTime() + tkEnvRefreshRate, 
			new cMessage("refresh", EV_TKENVUPDATE));

// Make that every subscribed particle compute its next event time
		for (p = particles.begin(); p != particles.end(); ++p) {
			(*p)->firstEventTime();
		}

		tkEnvUpdateNetwork();
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
		scheduleAt(simTime() + tkEnvRefreshRate, msg);

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
