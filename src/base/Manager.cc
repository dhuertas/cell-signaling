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

    unsigned int count;

	double diameter;

	cModule *module;
	std::list<Particle *>::const_iterator p;
	std::string particleDistribution;
	std::stringstream buffer;

	// Initialize variables
	N = 0;
	count = 0;
	diameter = 0.0;

	// The manager node should be the first module to be initialized
	if (stage == 0) {

		// Set the mode of operation for the molecule dynamics 
		setMode(par("mode"));

		// Set the name of the manager so we can have later access from
		// the other nodes.
		setName(par("name").stringValue());

		EV << "particle distribution: " << particleDistribution << "\n";

		// Get the simulation space size
		setSpaceSizeX(simulation.getSystemModule()->par("spaceSizeX"));
		setSpaceSizeY(simulation.getSystemModule()->par("spaceSizeY"));
		setSpaceSizeZ(simulation.getSystemModule()->par("spaceSizeZ"));

		EV << "space size: ";
		EV << "X=" << spaceSize.x << ", ";
		EV << "Y=" << spaceSize.y << ", ";
		EV << "Z=" << spaceSize.z << "\n";

		tkEnvRefreshRate = par("refreshRate");
		enableWebServer = par("enableWebServer");

		// Set network size for tkenv
		buffer << spaceSize.y;

		module = simulation.getSystemModule();
		module->getDisplayString().setTagArg("bgb", 0, buffer.str().c_str());
		buffer.str(std::string()); // clear buffer

		buffer << spaceSize.x;
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
		count = 0;

		for (p = particles.begin(); p != particles.end(); ++p) {

			diameter = 2*(*p)->getRadius();

			if (spaceCellSize < diameter) {
				spaceCellSize = diameter;
			}

			// Also set the simulation mode (Cell List or Near-Neighbor List)
			(*p)->setMode(mode);

			// Initialize particle attributes
			(*p)->setParticleId(count);
			(*p)->setLastCollisionTime(0);
			(*p)->initMessages();

			count++;

		}

		// Put every subscribed particle in its corresponding space cell
		setNumberOfSpaceCellsX(ceil(spaceSize.x/spaceCellSize));
		setNumberOfSpaceCellsY(ceil(spaceSize.y/spaceCellSize));
		setNumberOfSpaceCellsZ(ceil(spaceSize.z/spaceCellSize));

		N = Nx*Ny*Nz;

		spaceCells.resize(N);

		particleDistribution = par("particleDistribution").stringValue();

		if (particleDistribution.compare("undefined") != 0) {

			if (particleDistribution.compare("random") == 0) {
				randomDistribution(spaceSize, &particles);
			}

			if (particleDistribution.compare("cube") == 0) {
				cubeDistribution(spaceSize, &particles);
			}

			if (particleDistribution.compare("sphere") == 0) {
				//sphereDistribution(spaceSize, &particles, NULL, 0);
			}

			if (particleDistribution.compare("twosided") == 0) {
				//twoSidedDistribution(spaceSize, &particles, NULL, NULL);
			}

			if (particleDistribution.compare("highdensity") == 0) {
				//highDensityDistribution(spaceSize, &particles, NULL);
			}

		}

		for (p = particles.begin(); p != particles.end(); ++p) {
			if (mode == M_NNLIST) {
				(*p)->setListRadius(spaceCellSize);
			}
			attachParticleToSpaceCell(*p, -1);
		}

		if (mode == M_NNLIST) {
			for (p = particles.begin(); p != particles.end(); ++p) {
				(*p)->createNearNeighborList();
			}
		}

		// Self message to refresh the tk environment
		if (tkEnvRefreshRate > 0) {
			scheduleAt(simTime() + tkEnvRefreshRate, 
				new cMessage("refresh", EV_TKENVUPDATE));
		}

		// Make that every subscribed particle compute its next event time
		for (p = particles.begin(); p != particles.end(); ++p) {
			(*p)->initEvents();
		}

		tkEnvUpdateNetwork();

		// Start the web server
		if (enableWebServer == 1) {
			startWebServerThread();
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

		// Self message to refresh the tk environment
		if (tkEnvRefreshRate > 0) {
			scheduleAt(simTime() + tkEnvRefreshRate, msg);
		}

	}

}

/*
 * Clean and close everything.
 */
void Manager::finish() {

	if (enableWebServer == 1) {
		stopWebServerThread();
	}

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

/*
 * Spawns a new thread that starts a Web Server
 */
void Manager::startWebServerThread() {

	int res;

	pipe(quitFd);

	settings_t settings;

	// Gather simulation settings
	settings.numberOfParticles = particles.size();
	settings.simSpaceSize = spaceSize;

	webServerArgs.quitFd = quitFd[READ];
	webServerArgs.settings = settings;
	webServerArgs.particles = &particles;

	res = pthread_create(&webServerThread, NULL, startServerThread, &webServerArgs);

	if (res == 0) {
	    EV << "Web server started" << endl;
	} else {
	    EV << "Error starting the web server" << endl;
	}

}

/*
 * Terminates the Web Server thread
 */
void Manager::stopWebServerThread() {

	endServerThread(quitFd[WRITE]);

	pthread_join(webServerThread, NULL);

	close(quitFd[READ]);
	close(quitFd[WRITE]);

}
