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

/*
 * Constructor
 */
Manager::Manager() : 
  depth_(0),
  deltaTime_(0.0),
  space_(NULL),
  it_(NULL),
  maxSpaceSize_(0),
  nextParticleId_(0) {

}

/*
 * Destructor
 */
Manager::~Manager() {

  if (space_ != NULL) {
    delete space_;
  }

  if (it_ != NULL) {
    delete it_;
  }
}

/*
 * Subscribe particle to the simulation space
 *
 * @param {Particle *} p 
 */
void Manager::subscribe(Particle *p) {

  index3_t *idx = p->getSpaceCellIdx();

  if (idx->depth == 0) {
    point3_t *pos = p->getPosition();
    double radius = p->getRadius();

    while (2*radius / (maxSpaceSize_ / (1 << (idx->depth+1))) <= 1 && idx->depth+1 <= depth_) {
      idx->depth++;
    }

    idx->i = floor(pos->x / (maxSpaceSize_ / (1 << idx->depth)));
    idx->j = floor(pos->y / (maxSpaceSize_ / (1 << idx->depth)));
    idx->k = floor(pos->z / (maxSpaceSize_ / (1 << idx->depth)));
  }

  space_->insertAt(*idx, (cObject *)p);

  p->setParticleId(nextParticleId_);
  nextParticleId_++;
}

/*
 * Unsubscribe particle from the simulation space
 *
 * @param {Particle *} p
 */
void Manager::unsubscribe(Particle *p) {

  index3_t *idx = p->getSpaceCellIdx();
  space_->removeFrom(*idx, (cObject *)p);
}

/*
 * Initialize the simulation space and its modules.
 *
 * @param {Integer} stage
 */
void Manager::initialize(int stage) {

  cModule *module;
  std::stringstream buffer;

  // The manager node should be the first module to be initialized
  if (stage == 0) {

    nextParticleId_ = 0;

    depth_ = par("depth");

    deltaTime_ = par("deltaTime").doubleValue();

    space_ = new OctreeNode(NULL, NULL, 0, depth_, OCTREE);

    it_ = space_->getList()->iterator();

    // Initialize the statistics data structure
    clearStatistics();

    allCollisionsVector_.setName("all-collisions");
    particleCollisionsVector_.setName("particle-collisions");
    wallCollisionsVector_.setName("wall-collisions");
    transfersVector_.setName("transfers");
    expiresVector_.setName("expires");

    // Set the name of the manager so we can have later access from
    // the other nodes.
    setName(par("name").stringValue());

    // Get the simulation space size
    spaceSize_.x = simulation.getSystemModule()->par("spaceSizeX");
    spaceSize_.y = simulation.getSystemModule()->par("spaceSizeY");
    spaceSize_.z = simulation.getSystemModule()->par("spaceSizeZ");

    maxSpaceSize_ = max(maxSpaceSize_, spaceSize_.x);
    maxSpaceSize_ = max(maxSpaceSize_, spaceSize_.y);
    maxSpaceSize_ = max(maxSpaceSize_, spaceSize_.z);

    tkEnvRefreshRate_ = par("tkRefreshRate");
    statsRefreshRate_ = par("statsRefreshRate");

    // Set network size for tkenv
    buffer << spaceSize_.y;

    module = simulation.getSystemModule();
    module->getDisplayString().setTagArg("bgb", 0, buffer.str().c_str());
    buffer.str(std::string()); // clear buffer

    buffer << spaceSize_.x;
    module->getDisplayString().setTagArg("bgb", 1, buffer.str().c_str());
    buffer.str(std::string()); // clear buffer

  } else if (stage == 1) {
    // the rest of the modules are being initialized ...
  } else if (stage == 2) {

    it_->reset();
    while (it_->hasNext()) {
      // Step 1 - initialize mobility messages
      Particle *p = (Particle *)it_->next();
      p->setLastCollisionTime(0);
      p->initializeMobilityMessages();
      p->initializeMobility();
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

    if (ev.isGUI()) {
      tkEnvUpdateNetwork();
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
 * Transfer a particle from one space cell to another.
 *
 * @param {Particle *} p
 * @param {index3_t *} from
 * @param {index3_t *} to
 */
void Manager::transferParticle(Particle *p, index3_t *from, index3_t *to) {

  space_->removeFrom(*from, (cObject *)p);
  space_->insertAt(*to, (cObject *)p);
}

/*
 *
 * @param {index3_t *} idx
 * @param {vector<Particle *> *} container
 */
void Manager::getNeighborParticles(index3_t *idx, std::vector<Particle *> *container) {

  index3_t neighborIdx;
  index3_t currentIdx, tmpIdx;
  index3_t maxIdx;

  unsigned int imageIdx;

  maxIdx.i = floor(spaceSize_.x / (maxSpaceSize_ / (1 << idx->depth)));
  maxIdx.j = floor(spaceSize_.y / (maxSpaceSize_ / (1 << idx->depth)));
  maxIdx.k = floor(spaceSize_.z / (maxSpaceSize_ / (1 << idx->depth)));

  // Starting at current layer, go up in the octree to find neighbors 
  // in upper layers. Only get particles from lower layers when looking 
  // for neighbors in the starting layer. Once at upper layers only get
  // neighbors at the same layer (or otherwise all particles would be 
  // retrieved and duplicated!)
  for (unsigned int layer = idx->depth; layer != UINT_MAX; layer--)  {

    if (layer < idx->depth) {
      tmpIdx = currentIdx;
      currentIdx = getParentIdx(&tmpIdx);
    } else {
      currentIdx = *idx;
    }

    // Loop through the current 27 neighbor space cells
    for (int8_t i = -1; i <= 1; i++)
    for (int8_t j = -1; j <= 1; j++)
    for (int8_t k = -1; k <= 1; k++) {

      neighborIdx.depth = currentIdx.depth;

      // If a space cell idx belongs to an image (out of the simulation domain),
      // set the corresponding image idx to each particle

      // Images are identified as 3 pairs of 2 bits
      // --xx---- x axis
      // ----xx-- y axis
      // ------xx z axis
      // 00 means real domain
      // 10 means particles belong to the left image domain for that axis
      // 01 means particles belong to the right image domain for that axis
      imageIdx = 0;

      neighborIdx.i = currentIdx.i + i;
      neighborIdx.j = currentIdx.j + j;
      neighborIdx.k = currentIdx.k + k;

      // 0U - 1U == UINT_MAX
      if (maxIdx.k < neighborIdx.k && neighborIdx.k < UINT_MAX) {
          imageIdx |= 0x01;
          neighborIdx.k = 0;
      } else if (neighborIdx.k == UINT_MAX) {
          imageIdx |= 0x02;
          neighborIdx.k = maxIdx.k;
      }

      if (maxIdx.j < neighborIdx.j && neighborIdx.j < UINT_MAX) {
          imageIdx |= 0x04;
          neighborIdx.j = 0;
      } else if (neighborIdx.j == UINT_MAX) {
          imageIdx |= 0x08;
          neighborIdx.j = maxIdx.j;
      }

      if (maxIdx.i < neighborIdx.i && neighborIdx.i < UINT_MAX) {
          imageIdx |= 0x10;
          neighborIdx.i = 0;
      } else if (neighborIdx.i == UINT_MAX) {
          imageIdx |= 0x20;
          neighborIdx.i = maxIdx.i;
      }

      OctreeNode *node = space_->findOctreeNode(neighborIdx);

      if (node != NULL) {
        if (layer == idx->depth) {
          node->getAllElements(container, imageIdx);
        } else {
          node->getElements(container, imageIdx);
        }
      }
    }
  }
}

index3_t Manager::getParentIdx(index3_t *idx) {
  
  index3_t parentIdx;

  if (idx->depth == 0) {
    memcpy(&parentIdx, idx, sizeof(index3_t));
    return parentIdx;
  }

  parentIdx.i = idx->i/2;
  parentIdx.j = idx->j/2;
  parentIdx.k = idx->k/2;
  parentIdx.depth = idx->depth-1;

  return parentIdx;
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

}

/*
 * Updates the position of all the particles in the tk environment during
 * the simulation.
 */
void Manager::tkEnvUpdateNetwork() {

  double currentTime = simTime().dbl();

  it_->reset();
  while (it_->hasNext()) {
    Particle *p = (Particle *)it_->next();
    p->tkEnvUpdatePosition(currentTime);
  }
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
