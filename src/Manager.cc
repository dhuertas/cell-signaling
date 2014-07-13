//  This file is part of the Cell-Signaling project. Cell-Signaling is an
//  Omnet++ project to simulate cell signaling communications.
//  Copyright (C) 2014  Daniel Huertas
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
//	if (spaceCellSize > 0) {
//		attachParticleToSpaceCell(p, -1);
//	}

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

	double diameter, tempSpaceCellSize;

	cModule *module;
	std::list<Particle *>::const_iterator p;
	std::string particleDistribution;
	std::stringstream buffer;

	// Initialize variables
	N = 0;
	diameter = 0.0;
	tempSpaceCellSize = 0;

	// The manager node should be the first module to be initialized
	if (stage == 0) {

		// Initialize the statistics data structure
		clearStatistics();

		allCollisionsVector.setName("allCollisions");
		particleCollisionsVector.setName("particleCollisions");
		wallCollisionsVector.setName("wallCollisions");
		transfersVector.setName("transfers");
		expiresVector.setName("expires");

		// Set the mode of operation for the molecule dynamics 
		setMode(par("mode"));

		// Set the name of the manager so we can have later access from
		// the other nodes.
		setName(par("name").stringValue());

		// The global list radius. When set, overwrites the list radius
		// of the particles (default value is 0).
		setListRadius(par("listRadius"));

		EV << "particle distribution: " << particleDistribution << "\n";

		setDeltaTime(par("deltaTime").doubleValue());

		// Get the simulation space size
		setSpaceSizeX(simulation.getSystemModule()->par("spaceSizeX"));
		setSpaceSizeY(simulation.getSystemModule()->par("spaceSizeY"));
		setSpaceSizeZ(simulation.getSystemModule()->par("spaceSizeZ"));

		// Get the space cell size. If set to 0 we must wait untill all the 
		// initial particles have been subscribed.
		setSpaceCellSize(par("spaceCellSize"));

		EV << "space size: ";
		EV << "X=" << spaceSize.x << ", ";
		EV << "Y=" << spaceSize.y << ", ";
		EV << "Z=" << spaceSize.z << "\n";

		EV << "Space cell size: ";

		if (spaceCellSize == 0) {
			EV << "auto" << "\n";
		} else {
			EV << spaceCellSize << "\n";
		}

		tkEnvRefreshRate = par("tkRefreshRate");
		statsRefreshRate_ = par("statsRefreshRate");
		enableWebServer = par("enableWebServer");

		// Set network size for tkenv
		buffer << spaceSize.y;

		module = simulation.getSystemModule();
		module->getDisplayString().setTagArg("bgb", 0, buffer.str().c_str());
		buffer.str(std::string()); // clear buffer

		buffer << spaceSize.x;
		module->getDisplayString().setTagArg("bgb", 1, buffer.str().c_str());
		buffer.str(std::string()); // clear buffer

	} else if (stage == 1) {
		// the rest of the modules are being initialized ...
	} else if (stage == 2) {

		// All the particles are in the simulation space now. We can determine 
		// the space cell size in case it has been set to 0 (default).
		lastParticleId = 0;

		for (p = particles.begin(); p != particles.end(); ++p) {

			diameter = 2*(*p)->getRadius();

			if (tempSpaceCellSize < diameter) {
				tempSpaceCellSize = diameter;
			}

			// Initialize the particle attributes
			(*p)->setParticleId(lastParticleId);

			(*p)->setLastCollisionTime(0);

			(*p)->initMobilityMessages();

			lastParticleId++;
		}

		if (spaceCellSize == 0) {
			spaceCellSize = tempSpaceCellSize;
		}

		// Put every subscribed particle in its corresponding space cell
		setNumberOfSpaceCellsX(ceil(spaceSize.x/spaceCellSize));
		setNumberOfSpaceCellsY(ceil(spaceSize.y/spaceCellSize));
		setNumberOfSpaceCellsZ(ceil(spaceSize.z/spaceCellSize));

		// N is the total number of space cells that the simulation space has.
		N = Nx*Ny*Nz;

		spaceCells.resize(N);

		particleDistribution = par("particleDistribution").stringValue();

		if (particleDistribution.compare("undefined") != 0) {

			if (particleDistribution.compare("uniform") == 0) {

				uniformDistribution3(spaceSize, &particles);

			} else if (particleDistribution.compare("cube") == 0) {

				cubeDistribution(spaceSize, &particles);

			} else if (particleDistribution.compare("sphere") == 0) {

			    point_t c = {spaceSize.x/2, spaceSize.y/2, spaceSize.z/2};
				sphereDistribution(spaceSize, &particles, c, 0);

			} else if (particleDistribution.compare("highdensity") == 0) {

			    point_t c = {spaceSize.x/2, spaceSize.y/2, spaceSize.z/2};
				highDensityDistribution(spaceSize, &particles, c);

			} else if (particleDistribution.compare("densepacked") == 0) {

				point_t c = {spaceSize.x/2, spaceSize.y/2, spaceSize.z/2};
				densepacked(spaceSize, &particles, c);

			}
		}

		for (p = particles.begin(); p != particles.end(); ++p) {
			// If the manager has set a list radius, overwrite the list
			// radius of the particles
			if (mode == M_NNLIST && listRadius_ > 0) {
				(*p)->setListRadius(listRadius_);
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
			scheduleAt(simTime() + tkEnvRefreshRate/1000, 
				new cMessage("refresh", EV_TKENVUPDATE));
		}

		if (statsRefreshRate_ > 0) {
			scheduleAt(simTime() + statsRefreshRate_/1000,
				new cMessage("refresh", EV_STATSUPDATE));
		}

		// Make that every subscribed particle compute its next event time
		for (p = particles.begin(); p != particles.end(); ++p) {
			(*p)->initializeMobility();
		}

		if (ev.isGUI()) {
			tkEnvUpdateNetwork();
		}

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

	int i, j, k, n;
	point_t *pos = NULL;

	if (to < 0) {

		pos = p->getPosition();

		i = floor(pos->x/spaceCellSize);
		j = floor(pos->y/spaceCellSize);
		k = floor(pos->z/spaceCellSize);

		n = i*Ny*Nz + j*Nz + k;

		p->setSpaceCell(n);
		p->setPrevSpaceCell(-1);

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

	if (from < 0) {
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

	// Attach the particle to its new space cell
	attachParticleToSpaceCell(p, to);

}

/*
 * Return a list containing the pointers to the particles in a given
 * space cell.
 *
 * @param {Integer} n n-th space cell index
 * @return {std::list<Particle *> *} a pointer to a list of particle pointers
 */
std::list<Particle *> *Manager::getSpaceCellParticles(int n) {

	return &spaceCells.at(n);

}

int Manager::getNextParticleId() {

	int result = lastParticleId;

	lastParticleId++;

	return result;
}
/*
 * Handles every message that the manager module receives.
 *
 * @param {cMessage *} msg
 */
void Manager::handleMessage(cMessage *msg) {

	int kind = msg->getKind();

	simtime_t st = simTime();

	if (kind == EV_TKENVUPDATE) {

		tkEnvUpdateNetwork();

		// Self message to refresh the tk environment
		if (tkEnvRefreshRate > 0) {
			scheduleAt(st + tkEnvRefreshRate/1000, msg);
		}

	} else if (kind == EV_STATSUPDATE) {

		// Put the statistics logged so far to cout vectors
		allCollisionsVector.recordWithTimestamp(st, stats.allCollisions);
		particleCollisionsVector.recordWithTimestamp(st, stats.particleCollisions);
		wallCollisionsVector.recordWithTimestamp(st, stats.wallCollisions);
		transfersVector.recordWithTimestamp(st, stats.transfers);
		expiresVector.recordWithTimestamp(st, stats.expires);

		// Clear the stats data structure
		clearStatistics();

		if (statsRefreshRate_ > 0) {
			scheduleAt(st + statsRefreshRate_/1000, msg);
		}
	}

}

/*
 * Clean and close everything.
 */
void Manager::finish() {

	if (enableWebServer == 1) {
		// stopWebServerThread(); // TODO the server hangs the tk environment
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

	if(pipe(quitFd)) {
	    EV << "pipe failed\n";
	    return;
	}

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

/*
 * Sets the statistics data structure to zero.
 */
void Manager::clearStatistics() {
	memset(&stats, 0, sizeof(stats));
}

/*
 * Increment the collisions counter
 */
void Manager::registerCollision() {
	stats.particleCollisions++;
	stats.allCollisions++;
}

/*
 * Increment the wall collisions counter
 */
void Manager::registerWallCollision() {
	stats.wallCollisions++;
	stats.allCollisions++;
}

/*
 * Increment the transfer counter
 */
void Manager::registerTransfer() {
	stats.transfers++;
}

/*
 * Increment the expires counter
 */
void Manager::registerExpire() {
	stats.expires++;
}
