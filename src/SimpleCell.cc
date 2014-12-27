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

#include "SimpleCell.h"

using namespace std;

Define_Module(SimpleCell);

SimpleCell::~SimpleCell() {

}

/*
 * Cell initialization.
 *
 * @param {Integer} stage
 */
void SimpleCell::initialize(int stage) {

  std::stringstream buffer;

  if (stage == 0) {
    // Manager module initializes during this stage
  } else if (stage == 1) {

    setParticleType(T_SIMPLECELL);

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
    setMass(par("radius").doubleValue());

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

    statsRefreshRate_ = par("statsRefreshRate");

    if (statsRefreshRate_ > 0) {
      scheduleAt(simTime() + statsRefreshRate_/1000,
        new cMessage("refresh", EV_STATSUPDATE));
    }

    logCollisions_ = par("logCollisions");

    if (logCollisions_ > 0) {
      collisionTimeVector_ = new cOutVector("collisionTime");
    }

    // update Cell position in the tk environment
    tkEnvUpdatePosition();

    // draw module shape in the tk environment
    tkEnvDrawShape();
  }
}

/*
 * Returns the number of initialization stages.
 */
int SimpleCell::numInitStages() const {
  return 2;
}

/*
 * Handles every message that the cell receives.
 *
 * @param msg pointer to a cMessage object
 */
void SimpleCell::handleMessage(cMessage *msg) {

  int kind = msg->getKind();

  if (isSignaling(msg)) {

    handleSignaling(msg);

  } else if (ISMOBILITY(kind)) {

    handleMobilityMessage(msg);

  } else if (kind == EV_STATSUPDATE) {

    if (statsRefreshRate_ > 0) {
      scheduleAt(simTime() + statsRefreshRate_/1000, msg);
    }

  } else if (kind == EV_TTLEXPIRE) {

    expire();

  }

}

/*
 * Clean and close everything.
 */
void SimpleCell::finish() {

  // Unsubscribe from the manager
  getManager()->unsubscribe(this);

  // Delete mobility messages
  deleteMobilityMessages();

  if (timeToLive_ > 0) {
    cancelAndDelete(timeToLiveMsg_);
  }

}

/*
 * Cell must expire. This function gets called when the timeToLive 
 * parameter is set and the event EV_TTLEXPIRE arrives.
 */
void SimpleCell::expire() {

  manager_->registerExpire();

  // Unsubscribe from the manager
  manager_->unsubscribe(this);

  // Get out of the simulation space gracefully
  finishMobility();

  // Delete mobility messages
  deleteMobilityMessages();

  if (timeToLive_ > 0) {
    cancelAndDelete(timeToLiveMsg_);
  }

  deleteModule();
}

/*
 * Schedules a TTL message in the FES for the current cell.
 *
 * @param {double} time
 */
void SimpleCell::scheduleExpire(double time) {
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
 * Returns whether a message comes from a signaling process or not.
 *
 * @param {cMessage *} msg
 * @return {boolean} True if it is a signaling message
 */
bool SimpleCell::isSignaling(cMessage *msg) {

  int kind = msg->getKind();
  CollisionMessage *cmsg;

  cmsg = NULL;

  if (kind == EV_COLLISION) {

    cmsg = (CollisionMessage *)msg;

    if (particleType_ == T_RECEIVER ||
      particleType_ == T_EMITTER_RECEIVER) {

      if (cmsg->getPartner()->getParticleType() == T_SIGNALING) {
        return true;
      }
    }
  }

  return false;
}

/*
 * Handle the signaling process. This is a receiver and the other is a molecule.
 *
 * @param {cMessage *} msg
 */
void SimpleCell::handleSignaling(cMessage *msg) {

  MoleculeReceiver *receiver;
  Particle *p;

  CollisionMessage *cmsg;

  cmsg = (CollisionMessage *)msg;
  p = cmsg->getPartner();

  receiver = ((MoleculeReceiver *)getParentModule()->getSubmodule("receiver"));
  receiver->registerReception(p->getParticleType());

  // TODO change the following line, perhaps using gates
  ((Molecule *)p)->scheduleExpire(cmsg->getCollisionTime());

  // Change the event type to CHECK so we can find the next event
  // for the receiver
  msg->setKind(EV_CHECK);
  handleMobilityMessage(msg);
}
