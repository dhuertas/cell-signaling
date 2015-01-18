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

#include "Molecule.h"

using namespace std;

Define_Module(Molecule);

Molecule::~Molecule() {}

/*
 *
 */
void Molecule::initialize(int stage) {

  std::stringstream buffer;

  if (stage == 0) {
    // Manager module initializes during this stage
  } else if (stage == 1) {

    active_ = true;

    setParticleType(T_MOLECULE);

    // Initial position
    position_.x = par("xpos").doubleValue();
    position_.y = par("ypos").doubleValue();
    position_.z = par("zpos").doubleValue();

    //setX(par("xpos").doubleValue());
    //setY(par("ypos").doubleValue());
    //setZ(par("zpos").doubleValue());

    // Velocity
    velocity_.x = par("vx").doubleValue();
    velocity_.y = par("vy").doubleValue();
    velocity_.z = par("vz").doubleValue();

    //setVx(par("vx").doubleValue());
    //setVy(par("vy").doubleValue());
    //setVz(par("vz").doubleValue());

    // Cell radius
    setRadius(par("radius").doubleValue());
    setMass(par("mass").doubleValue());


    // Subscribe to manager
    setManager("manager");
    getManager()->subscribe(this);

    // Brownian Motion parameters
    setDiffusion(par("diffusion").doubleValue());

    // Compute Brownian Motion Standard Deviation
    //        _____________
    //  \    /             |
    //   \  /  4*M_PI*D*dt
    //    \/
    double diffusion = par("diffusion").doubleValue();
    double dt = manager_->getDeltaTime();
    setBrownianMotionStdDev(sqrt(4*M_PI*diffusion*dt));

    setBoundariesMode(par("boundariesMode"));

    timeToLive_ = par("timeToLive");

    if (timeToLive_ > 0) {
      timeToLiveMsg_ = new TimeToLiveMessage("expire", EV_TTLEXPIRE);
      scheduleAt(simTime() + timeToLive_, timeToLiveMsg_);
    }

    logCollisions_ = par("logCollisions").boolValue();

    if (logCollisions_ > 0) {
      collisionTimeVector_ = new cOutVector("collisionTime");
      xCollisionPositionVector_ = new cOutVector("xCollisionPosition");
      yCollisionPositionVector_ = new cOutVector("yCollisionPosition");
      zCollisionPositionVector_ = new cOutVector("zCollisionPosition");
    }

    logPosition_ = par("logPosition").boolValue();

    if (logPosition_ > 0) {
      xPositionVector_ = new cOutVector("xPosition");
      yPositionVector_ = new cOutVector("yPosition");
      zPositionVector_ = new cOutVector("zPosition");
    }

    statsRefreshRate_ = par("statsRefreshRate");

    if (statsRefreshRate_ > 0) {
      scheduleAt(simTime() + statsRefreshRate_/1000,
        new cMessage("refresh", EV_STATSUPDATE));
    }

    // update Molecule position in the tk environment
    tkEnvUpdatePosition();

    // draw module shape in the tk environment
    tkEnvDrawShape();

  } else {

  }

}

/*
 * The molecules must be initialized during the second stage.
 */
int Molecule::numInitStages() const {

  return 2;

}

/*
 * Handles every message that the molecule receives.
 *
 * @param msg pointer to a cMessage object
 */
void Molecule::handleMessage(cMessage *msg) {

  int kind = msg->getKind();

  if (isSignaling(msg)) {

    handleSignaling(msg);

  } else if (ISMOBILITY(kind)) {

    handleMobilityMessage(msg);

  } else if (kind == EV_STATSUPDATE) {

    double st = simTime().dbl();
    double dt = st - lastCollisionTime_;

    // Put the statistics logged so far to cout vectors
    xPositionVector_->recordWithTimestamp(st, position_.x + velocity_.x * dt);
    yPositionVector_->recordWithTimestamp(st, position_.y + velocity_.y * dt);
    zPositionVector_->recordWithTimestamp(st, position_.z + velocity_.z * dt);

    if (statsRefreshRate_ > 0) {
      scheduleAt(st + statsRefreshRate_/1000, msg);
    }

  } else if (kind == EV_TTLEXPIRE) {
    expire();
  }

}

/*
 * Clean and close everything.
 */
void Molecule::finish() {

  // Unsubscribe from the manager
  getManager()->unsubscribe(this);

  // Delete mobility messages
  deleteMobilityMessages();

  if (timeToLive_ > 0) {
    cancelAndDelete(timeToLiveMsg_);
  }

}

/*
 * Molecule must expire. This function gets called when the timeToLive 
 * parameter is set and the event EV_TTLEXPIRE arrives.
 */
void Molecule::expire() {
  // Methods called from other modules must have this macro
  Enter_Method_Silent();

  active_ = false;

  manager_->registerExpire();

  // Unsubscribe from the manager
  manager_->unsubscribe(this);

  // Get out of the simulation space gracefully
  this->finishMobility();

  // Delete mobility messages
  deleteMobilityMessages();

  if (timeToLive_ > 0) {
    cancelAndDelete(timeToLiveMsg_);
  }

  this->deleteModule();

}

/*
 *
 * @param {double} time: simulation time to schedule the expire message
 */
void Molecule::scheduleExpire(double time) {
  // Methods called from other modules must have this macro
  Enter_Method_Silent();

  if (timeToLive_ > 0) {
    if (timeToLiveMsg_->isScheduled()) {
      cancelEvent(timeToLiveMsg_);
    }
  } else {
    timeToLiveMsg_ = new TimeToLiveMessage("expire", EV_TTLEXPIRE);
  }

  scheduleAt(time, timeToLiveMsg_);

}

/*
 * Tells whether the molecule has to process a collision as a signaling.
 *
 * @param {cMessage *} msg
 * @return {bool}
 */
bool Molecule::isSignaling(cMessage *msg) {

  int partnerParticleType;

  int kind = msg->getKind();
  CollisionMessage *cmsg;

  cmsg = NULL;
  partnerParticleType = 0;

  if (kind == EV_COLLISION) {

    cmsg = (CollisionMessage *)msg;

    if (particleType_ == T_SIGNALING) {

      partnerParticleType = cmsg->getPartner()->getParticleType();

      if (partnerParticleType == T_RECEIVER ||
        partnerParticleType == T_EMITTER_RECEIVER) {
        return true;
      }
    }
  }

  return false;
}

/*
 * Handle the signaling process. This is a molecule and the other is a receiver.
 *
 * @param {cMessage *} msg
 */
void Molecule::handleSignaling(cMessage *msg) {

  MoleculeReceiver *receiver;
  Particle *p;

  CollisionMessage *cmsg;

  cmsg = (CollisionMessage *)msg;
  p = cmsg->getPartner();

  receiver = (MoleculeReceiver *)((SimpleCell *)p)->getParentModule()
    ->getSubmodule("receiver");

  receiver->registerReception(getParticleType());

  expire();
}
