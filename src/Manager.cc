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

  // Add the particle pointer to the particles list
  particles_.push_back(p);

  // Also add the particle to its corresponding space cell
  space_.subscribe(p);

}

/*
 * Unsubscribe particle from the simulation space
 *
 * @param {Particle *} p
 */
void Manager::unsubscribe(Particle * p) {

  // Remove the particle pointer from the space cell structure
  space_.unsubscribe(p);

  // Remove the particle pointer from the particles list
  particles_.remove(p);

}

/*
 * Initialize the simulation space and its modules.
 *
 * @param {Integer} stage
 */
void Manager::initialize(int stage) {

  int N;

  double diameter, tempSpaceCellSize;

  vect_t spaceSize;

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

    allCollisionsVector_.setName("allCollisions");
    particleCollisionsVector_.setName("particleCollisions");
    wallCollisionsVector_.setName("wallCollisions");
    transfersVector_.setName("transfers");
    expiresVector_.setName("expires");

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
    spaceSize.x = simulation.getSystemModule()->par("spaceSizeX");
    spaceSize.y = simulation.getSystemModule()->par("spaceSizeY");
    spaceSize.z = simulation.getSystemModule()->par("spaceSizeZ");

    space_.setSpaceSize(spaceSize);

    // Get the space cell size. If set to 0 we must wait until all the
    // initial particles have been subscribed.
    space_.setMaxSpaceCellSize(par("maxSpaceCellSize"));

    tkEnvRefreshRate_ = par("tkRefreshRate");
    statsRefreshRate_ = par("statsRefreshRate");
    enableWebServer_ = par("enableWebServer");

    // Set network size for tkenv
    buffer << space_.getSpaceSize()->y;

    module = simulation.getSystemModule();
    module->getDisplayString().setTagArg("bgb", 0, buffer.str().c_str());
    buffer.str(std::string()); // clear buffer

    buffer << space_.getSpaceSize()->x;
    module->getDisplayString().setTagArg("bgb", 1, buffer.str().c_str());
    buffer.str(std::string()); // clear buffer

  } else if (stage == 1) {
    // the rest of the modules are being initialized ...
  } else if (stage == 2) {

    // All the particles are in the simulation space now. We can determine 
    // the space cell size in case it has been set to 0 (default).
    lastParticleId_ = 0;

    for (p = particles_.begin(); p != particles_.end(); ++p) {

      diameter = 2*(*p)->getRadius();

      if (tempSpaceCellSize < diameter) {
        tempSpaceCellSize = diameter;
      }

      // Initialize the particle attributes
      (*p)->setParticleId(lastParticleId_);

      (*p)->setLastCollisionTime(0);

      (*p)->initMobilityMessages();

      lastParticleId_++;
    }

    if (space_.getMaxSpaceCellSize() == 0) {
      space_.setMaxSpaceCellSize(tempSpaceCellSize);
    }

    particleDistribution = par("particleDistribution").stringValue();

    if (particleDistribution.compare("undefined") != 0) {

      if (particleDistribution.compare("uniform") == 0) {

        uniformDistribution3(*space_.getSpaceSize(), &particles_);

      } else if (particleDistribution.compare("cube") == 0) {

        cubeDistribution(*space_.getSpaceSize(), &particles_);

      } else if (particleDistribution.compare("sphere") == 0) {

        point_t c = {
          space_.getSpaceSize()->x/2, 
          space_.getSpaceSize()->y/2, 
          space_.getSpaceSize()->z/2};

        sphereDistribution(*space_.getSpaceSize(), &particles_, c, 0);

      } else if (particleDistribution.compare("highdensity") == 0) {

        point_t c = {
          space_.getSpaceSize()->x/2, 
          space_.getSpaceSize()->y/2, 
          space_.getSpaceSize()->z/2};

        highDensityDistribution(*space_.getSpaceSize(), &particles_, c);

      } else if (particleDistribution.compare("densepacked") == 0) {

        point_t c = {
          space_.getSpaceSize()->x/2, 
          space_.getSpaceSize()->y/2, 
          space_.getSpaceSize()->z/2};
        densepacked(*space_.getSpaceSize(), &particles_, c);

      }
    }

    for (p = particles_.begin(); p != particles_.end(); ++p) {
      // If the manager has set a list radius, overwrite the list
      // radius of the particles
      if (mode_ == M_NNLIST && listRadius_ > 0) {
        (*p)->setListRadius(listRadius_);
      }

      space_.attachParticleToSpaceCell(*p, IDX_NULL);
    }

    if (mode_ == M_NNLIST) {
      for (p = particles_.begin(); p != particles_.end(); ++p) {
        (*p)->createNearNeighborList();
      }
    }

    // Self message to refresh the tk environment
    if (tkEnvRefreshRate_ > 0) {
      scheduleAt(simTime() + tkEnvRefreshRate_/1000, 
        new cMessage("refresh", EV_TKENVUPDATE));
    }

    if (statsRefreshRate_ > 0) {
      scheduleAt(simTime() + statsRefreshRate_/1000,
        new cMessage("refresh", EV_STATSUPDATE));
    }

    // Make that every subscribed particle compute its next event time
    for (p = particles_.begin(); p != particles_.end(); ++p) {
      (*p)->initializeMobility();
    }

    if (ev.isGUI()) {
      tkEnvUpdateNetwork();
    }

    // Start the web server
    if (enableWebServer_ == 1) {
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
 * @param {index_t} to
 */
void Manager::attachParticleToSpaceCell(Particle *p, index_t to) {

  space_.attachParticleToSpaceCell(p, to);

}

/*
 * Detach a particle from a space cell.
 *
 * @param {Particle *} p
 * @param {index_t} from
 */
void Manager::detachParticleFromSpaceCell(Particle *p, index_t from) {

  space_.detachParticleFromSpaceCell(p, from);

}

/*
 * Transfer a particle from one space cell to another.
 *
 * @param {Particle *} p
 * @param {index_t} from
 * @param {index_t} to
 */
void Manager::transferParticle(Particle *p, index_t from, index_t to) {

  space_.transferParticle(p, from, to);

}

/*
 *
 */
std::vector<Particle *> *Manager::getNeighborParticles(index_t idx, std::vector<Particle *> *list) {

  return space_.getNeighborParticles(idx, list);
}

/*
 *
 */
int Manager::getNextParticleId() {

  int result = lastParticleId_;

  lastParticleId_++;

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
    if (tkEnvRefreshRate_ > 0) {
      scheduleAt(st + tkEnvRefreshRate_/1000, msg);
    }

  } else if (kind == EV_STATSUPDATE) {

    // Put the statistics logged so far to cout vectors
    allCollisionsVector_.recordWithTimestamp(st, stats_.allCollisions);
    particleCollisionsVector_.recordWithTimestamp(st, stats_.particleCollisions);
    wallCollisionsVector_.recordWithTimestamp(st, stats_.wallCollisions);
    transfersVector_.recordWithTimestamp(st, stats_.transfers);
    expiresVector_.recordWithTimestamp(st, stats_.expires);

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

  if (enableWebServer_ == 1) {
    // stopWebServerThread(); // TODO the server hangs the tk environment here
  }

}

/*
 * Updates the position of all the particles in the tk environment during
 * the simulation.
 */
void Manager::tkEnvUpdateNetwork() {

  std::list<Particle *>::const_iterator p;

  // Update particle positions
  for (p = particles_.begin(); p != particles_.end(); ++p) {
    (*p)->tkEnvUpdatePosition(simTime().dbl());
  }

}

/*
 * Spawns a new thread that starts a Web Server
 */
void Manager::startWebServerThread() {

  if(pipe(quitFd_)) {
      EV << "pipe failed\n";
      return;
  }

  settings_t settings;

  // Gather simulation settings
  settings.numberOfParticles = particles_.size();
  settings.simSpaceSize = *(space_.getSpaceSize());

  webServerArgs_.quitFd = quitFd_[READ];
  webServerArgs_.settings = settings;
  webServerArgs_.particles = &particles_;

  if (pthread_create(&webServerThread_, NULL, startServerThread, &webServerArgs_) == 0) {
      EV << "Web server started" << endl;
  } else {
      EV << "Error starting the web server" << endl;
  }

}

/*
 * Terminates the Web Server thread
 */
void Manager::stopWebServerThread() {

  endServerThread(quitFd_[WRITE]);

  pthread_join(webServerThread_, NULL);

  close(quitFd_[READ]);
  close(quitFd_[WRITE]);

}

/*
 * Sets the statistics data structure to zero.
 */
void Manager::clearStatistics() {
  memset(&stats_, 0, sizeof(stats_));
}

/*
 * Increment the collisions counter
 */
void Manager::registerCollision() {
  stats_.particleCollisions++;
  stats_.allCollisions++;
}

/*
 * Increment the wall collisions counter
 */
void Manager::registerWallCollision() {
  stats_.wallCollisions++;
  stats_.allCollisions++;
}

/*
 * Increment the transfer counter
 */
void Manager::registerTransfer() {
  stats_.transfers++;
}

/*
 * Increment the expires counter
 */
void Manager::registerExpire() {
  stats_.expires++;
}
