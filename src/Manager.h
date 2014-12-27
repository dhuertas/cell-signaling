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

#ifndef MANAGER_H
#define MANAGER_H

#include <cmessage.h>
#include <csimplemodule.h>
#include <cqueue.h>
#include <omnetpp.h>

#include <cmath>
#include <vector>
#include <list>

#include "Particle.h"
#include "ParticleDistributionHelper.h"
#include "OctreeNode.h"

class Manager : public cSimpleModule {

 protected:

  unsigned int depth_;

  // Delta time
  double deltaTime_;

  vector3_t spaceSize_;

  // The simulation space
  OctreeNode *space_;

  LinkedListIterator *it_;

  double maxSpaceSize_;

  // Contains the particle Id from the last added particle to the domain.
  long long unsigned nextParticleId_;

  // File Descriptors to communicate with the Web Server thread
  int quitFd_[2];

  // TK environment refresh rate
  // It allows the manager module to send self-messages in order to
  // update the position of each particle.
  double tkEnvRefreshRate_;

  // Update the cOutVectors periodically.
  double statsRefreshRate_;

  // Manager name
  std::string name_;

  // Statistics & vectors
  statistics_t stats_;

  cOutVector allCollisionsVector_;

  cOutVector particleCollisionsVector_;

  cOutVector wallCollisionsVector_;

  cOutVector transfersVector_;

  cOutVector expiresVector_;

 public:

  Manager();

  ~Manager();

  //
  // cSimpleModule inheritance
  //
  virtual void initialize(int stage);

  virtual int numInitStages() const;

  virtual void handleMessage(cMessage *);

  virtual void finish();

  //
  // Handle particles
  //
  void subscribe(Particle *);

  void unsubscribe(Particle *);

  void transferParticle(Particle *, index3_t *from, index3_t *to);

  //
  // Update the tk environment
  //
  void tkEnvUpdateNetwork(void);

  //
  // Gets and sets
  //
  vector3_t *getSpaceSize(void) { return &spaceSize_; }

  unsigned int getDepth() { return depth_; }

  double getMaxSpaceSize() { return maxSpaceSize_; }

  double getDeltaTime(void) { return deltaTime_; };

  void getNeighborParticles(index3_t *idx, std::vector<Particle *> *container);

  long long unsigned getNextParticleId(void) { return nextParticleId_; }

  long long unsigned getNextIncParticleId(void) { 
    long long unsigned n = nextParticleId_; 
    nextParticleId_++;
    return n;
  }

  void setMaxSpaceSize(double spaceSize) { maxSpaceSize_ = spaceSize; }

  void setDeltaTime(double dt) { deltaTime_ = dt; }

  //
  // Statistics
  //
  void clearStatistics(void);

  void registerCollision(void);

  void registerWallCollision(void);

  void registerTransfer(void);

  void registerExpire(void);

};

#endif
