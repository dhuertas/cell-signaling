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

#ifndef PARTICLE_H
#define PARTICLE_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <coutvector.h>
#include <omnetpp.h>

#include "Common.h"
#include <vector>

class Manager;

class Particle : public cSimpleModule {

 private:

 protected:

  Manager *manager_;

  bool active_;

  point3_t position_;

  vector3_t velocity_;

  double mass_;

  double radius_;

  double lastCollisionTime_;

  unsigned int particleId_;

  unsigned int particleType_;

  index3_t spaceCellIdx_;

  // Whether the particle must bounce, expire, etc.
  unsigned int boundariesMode_;

  unsigned int imageIdx_;

  //
  // Brownian motion parameters
  // 
  double diffusion_;

  double brownianMotionStdDev_;

 public:

  Particle();

  //
  // cSimpleModule inheritance
  //
  virtual void initialize(int stage);

  virtual int numInitStages() const;

  virtual void handleMessage(cMessage *);

  virtual void finish();

  //
  // Gets and sets
  //
  Manager *getManager()  { return manager_; }

  point3_t *getPosition(void) { return &position_; }

  vector3_t *getVelocity(void) { return &velocity_; }

  double getMass(void) { return mass_; }

  double getRadius(void) { return radius_; }

  double getLastCollisionTime(void) { return lastCollisionTime_; }

  unsigned int getParticleId(void) { return particleId_; }

  unsigned int getParticleType(void) { return particleType_; }

  index3_t *getSpaceCellIdx(void) { return &spaceCellIdx_; }

  unsigned int getBoundariesMode(void) const { return boundariesMode_; }

  unsigned int getImageIdx() { return imageIdx_; }

  double getDiffusion(void) { return diffusion_; }

  double getBrownianMotionStdDev(void) { return brownianMotionStdDev_; }

  void setManager(Manager *manager) { manager_ = manager; }

  void setPosition(point3_t p) { position_ = p; }

  void setVelocity(vector3_t v) { velocity_ = v; }

  void setMass(double m) { mass_ = m; }

  void setRadius(double r) { radius_ = r; }

  void setLastCollisionTime(double t) { lastCollisionTime_ = t; }

  void setParticleId(unsigned int id) { particleId_ = id; }

  void setParticleType(unsigned int t) { particleType_ = t; }

  void setSpaceCellIdx(index3_t idx) { spaceCellIdx_ = idx; }

  void setBoundariesMode(int mode) { boundariesMode_ = mode; }

  void setImageIdx(unsigned int imageIdx) { imageIdx_ = imageIdx; }

  void setDiffusion(double d) { diffusion_ = d; }

  void setBrownianMotionStdDev(double sd) { brownianMotionStdDev_ = sd; };

  void setActive(bool a) { active_ = a; };

  void setManager(std::string name);

  //
  // Boolean methods
  //
  bool isActive (void) { return active_; };

  //
  // tk Environment methods
  //
  virtual void tkEnvDrawShape(void) = 0;

  virtual void tkEnvUpdatePosition(void) = 0;

  virtual void tkEnvUpdatePosition(double) = 0;

  //
  // Event related methods
  //
  virtual void initializeMobilityMessages() = 0;

  virtual void initializeMobility() = 0;

  virtual void finishMobility() = 0;

  virtual void expire() = 0;
};

#endif
