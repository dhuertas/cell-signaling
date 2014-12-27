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

#include "Sphere.h"
#include "mobility/SphereMobility.h"
#include "../Molecule.h"
#include "../SimpleCell.h"
#include "../receiver/MoleculeReceiver.h"

using namespace std;

/*
 * Constructor
 */
Sphere::Sphere() :
  transferMsg_(NULL),
  collisionMsg_(NULL),
  manager_(NULL),
  logCollisions_(false),
  logPosition_(false),
  collisionTimeVector_(NULL),
  xCollisionPositionVector_(NULL),
  yCollisionPositionVector_(NULL),
  zCollisionPositionVector_(NULL) {

}

/*
 *
 */
Sphere::~Sphere() {

}

/*
 * Initialize the self messages. This method is called by the manager module
 * during network initialization.
 */
void Sphere::initializeMobilityMessages() {
  // Methods called from other modules must have this macro
  Enter_Method_Silent();

  transferMsg_ = new TransferMessage("mobility", EV_TRANSFER);
  transferMsg_->setManager(manager_);

  collisionMsg_ = new CollisionMessage("mobility", EV_NONE);
  collisionMsg_->setManager(manager_);

  if (manager_->getDeltaTime() > 0) {
    brownianMotionMsg_ = new BrownianMotionMessage("mobility", EV_BROWNIAN);
    brownianMotionMsg_->setManager(manager_);
  }

  SphereMobility::resetCollisionMessage(collisionMsg_);

}

/*
 *
 */
void Sphere::deleteMobilityMessages() {

  if (transferMsg_ != NULL) {
    cancelAndDelete(transferMsg_);
    transferMsg_ = NULL;
  }

  if (collisionMsg_ != NULL) {
    cancelAndDelete(collisionMsg_);
    collisionMsg_ = NULL;
  }

  if (brownianMotionMsg_ != NULL) {
    cancelAndDelete(brownianMotionMsg_);
    brownianMotionMsg_ = NULL;
  }

}

/*
 * Initialize the event queue by computing the first event for the sphere. This
 * method is called by the manager in order to initialize the event queue.
 */
void Sphere::initializeMobility() {
  // Methods called from other modules must have this macro
  Enter_Method_Silent();

  double transferTime = NO_TIME;
  double scheduledCollisionTime = NO_TIME;

  double collisionTime = NO_TIME;
  double boundaryCollisionTime = NO_TIME;
  double brownianMotionTime = NO_TIME;

  double minTime = NO_TIME;

  vector<double> times;
  vector<double>::const_iterator t;

  // Compute the first collision and the first transfer
  transferTime = SphereMobility::nextTransferTime(transferMsg_, this);

  if (transferTime != NO_TIME) {
    scheduleAt(transferTime, transferMsg_);
  }

  collisionTime = SphereMobility::nextCollisionTime(collisionMsg_, this);
  boundaryCollisionTime = SphereMobility::nextBoundaryCollisionTime(collisionMsg_, this);

  if (manager_->getDeltaTime() > 0) {
    brownianMotionTime = SphereMobility::brownianMotion(brownianMotionMsg_, this);
  }

  if (collisionMsg_->isScheduled()) {
    scheduledCollisionTime = collisionMsg_->getCollisionTime();
    times.push_back(scheduledCollisionTime);
  }

  if (collisionTime != NO_TIME) {
    times.push_back(collisionTime);
  }

  if (boundaryCollisionTime != NO_TIME) {
    times.push_back(boundaryCollisionTime);
  }

  if (brownianMotionTime != NO_TIME) {
    times.push_back(brownianMotionTime);
  }

  minTime = boundaryCollisionTime;

  for (t = times.begin(); t!= times.end(); ++t) {
    if ((*t) < minTime) minTime = (*t);
  }

  if (minTime == collisionTime && collisionTime != NO_TIME) {

    if (collisionMsg_->isScheduled()) {
      ((Sphere *)collisionMsg_->getPrevPartner())->getCollisionMessage()->setKind(EV_CHECK);
      ((Sphere *)collisionMsg_->getPrevPartner())->getCollisionMessage()->setPartner(NULL);
    }

    cancelEvent(collisionMsg_);

    collisionMsg_->setKind(EV_COLLISION);
    collisionMsg_->setCollisionTime(collisionTime);

    scheduleAt(collisionTime, collisionMsg_);

    ((Sphere *)collisionMsg_->getPartner())->adjustCollision(collisionTime, this);

  } else if (minTime == boundaryCollisionTime && boundaryCollisionTime != NO_TIME) {

    if (collisionMsg_->isScheduled()) {
      ((Sphere *)collisionMsg_->getPrevPartner())->getCollisionMessage()->setKind(EV_CHECK);
      ((Sphere *)collisionMsg_->getPrevPartner())->getCollisionMessage()->setPartner(NULL);
    }

    cancelEvent(collisionMsg_);

    collisionMsg_->setKind(EV_BOUNDARYCOLLISION);
    collisionMsg_->setCollisionTime(boundaryCollisionTime);

    scheduleAt(boundaryCollisionTime, collisionMsg_);

  } else if (minTime == brownianMotionTime && brownianMotionTime != NO_TIME) {

    if (collisionMsg_->isScheduled()) {
      ((Sphere *)collisionMsg_->getPrevPartner())->getCollisionMessage()->setKind(EV_CHECK);
      ((Sphere *)collisionMsg_->getPrevPartner())->getCollisionMessage()->setPartner(NULL);
    }

    cancelEvent(collisionMsg_);

    scheduleAt(brownianMotionTime, brownianMotionMsg_);

  } else if (minTime == scheduledCollisionTime) {
    // Leave it as it is scheduled
  } else {
    // minTime equals NO_TIME;
  }

}

/*
 * This method gets called when the sphere must leave the simulation space
 * (either it has expired, crossed a boundary, etc).
 */
void Sphere::finishMobility() {

  active_ = false;

  if (transferMsg_->isScheduled()) {
    cancelEvent(transferMsg_);
  }

  if (brownianMotionMsg_ != NULL && brownianMotionMsg_->isScheduled()) {
    cancelEvent(brownianMotionMsg_);
  }

  if (collisionMsg_->isScheduled()) {

    cancelEvent(collisionMsg_);

    // Change the event type of the third party sphere to EV_CHECK
    if (collisionMsg_->getPartner() != NULL) {
      ((Sphere *)collisionMsg_->getPartner())->getCollisionMessage()->setKind(EV_CHECK);
      ((Sphere *)collisionMsg_->getPartner())->getCollisionMessage()->setPartner(NULL);
    }
  }
}

/*
 * Method overwrite. Finish mobility when called from another sphere (it has
 * been absorved by a receiver, combined with another molecule, etc).
 *
 * @param {Particle *} from
 */
void Sphere::finishMobility(Particle *from) {
  // Methods called from other modules must have this macro
  Enter_Method_Silent();

  active_ = false;

  if (transferMsg_->isScheduled()) {
    cancelEvent(transferMsg_);
  }

  if (brownianMotionMsg_ != NULL && brownianMotionMsg_->isScheduled()) {
    cancelEvent(brownianMotionMsg_);
  }

  if (collisionMsg_->isScheduled()) {

    cancelEvent(collisionMsg_);

    // Change the event type of the third party sphere to EV_CHECK
    if (collisionMsg_->getPartner() != NULL) {
      ((Sphere *)collisionMsg_->getPartner())->getCollisionMessage()->setKind(EV_CHECK);
      ((Sphere *)collisionMsg_->getPartner())->getCollisionMessage()->setPartner(NULL);
    }
  }
}

/*
 * Handles the mobility of the sphere
 *
 * @param {cMessage *} msg
 */
void Sphere::handleMobilityMessage(cMessage *msg) {

  double transferTime = NO_TIME;
  double collisionTime = NO_TIME;
  double boundaryCollisionTime = NO_TIME;
  double brownianMotionTime = NO_TIME;

  // Step 1. Find the next event in the queue.

  int kind = msg->getKind();

  // Step 2. Handle the event.

  if (kind == EV_TRANSFER) {

    SphereMobility::handleTransfer((TransferMessage *)msg, this);

  } else if (kind == EV_COLLISION) {

    SphereMobility::handleCollision((CollisionMessage *)msg, this);
    SphereMobility::resetCollisionMessage(collisionMsg_);

    // Reset brownian motion
    if (manager_->getDeltaTime() > 0) {
      if (brownianMotionMsg_->isScheduled()) {
        cancelEvent(brownianMotionMsg_);
      }
      // Compute the next velocity with brownian motion
      brownianMotionTime = SphereMobility::brownianMotion(brownianMotionMsg_, this);
      scheduleAt(brownianMotionTime, brownianMotionMsg_);
    }

  } else if (kind == EV_BOUNDARYCOLLISION) {

    SphereMobility::handleBoundaryCollision((CollisionMessage *)msg, this);
    SphereMobility::resetCollisionMessage(collisionMsg_);

    // Reset brownian motion
    if (manager_->getDeltaTime() > 0) {
      if (brownianMotionMsg_->isScheduled()) {
        cancelEvent(brownianMotionMsg_);
      }
      // Compute the next velocity with brownian motion
      brownianMotionTime = SphereMobility::brownianMotion(brownianMotionMsg_, this);
      scheduleAt(brownianMotionTime, brownianMotionMsg_);
    }

  } else if (kind == EV_BROWNIAN) {

    // Update position, velocity and lastCollisionTime
    SphereMobility::handleBrownianMotion((BrownianMotionMessage *)msg, this);

    // Compute the next velocity with brownian motion
    brownianMotionTime = SphereMobility::brownianMotion(brownianMotionMsg_, this);
    scheduleAt(brownianMotionTime, brownianMotionMsg_);

    // Collisions must be rescheduled
    if (collisionMsg_->isScheduled()) cancelEvent(collisionMsg_);
    SphereMobility::resetCollisionMessage(collisionMsg_);

  } else if (kind == EV_CHECK) {

    SphereMobility::resetCollisionMessage(collisionMsg_);

    // Reset brownian motion
    if (manager_->getDeltaTime() > 0) {
      if (brownianMotionMsg_->isScheduled()) {
        cancelEvent(brownianMotionMsg_);
      }
      // Compute the next velocity with brownian motion
      brownianMotionTime = SphereMobility::brownianMotion(brownianMotionMsg_, this);
      scheduleAt(brownianMotionTime, brownianMotionMsg_);
    }
  }

  // Step 3. Compute the next transfer time for the particle corresponding to the
  // event.

  if (transferMsg_->isScheduled()) cancelEvent(transferMsg_);

  transferTime = SphereMobility::nextTransferTime(transferMsg_, this);

  if (transferTime != NO_TIME) {
    scheduleAt(transferTime, transferMsg_);
  }

  // Step 4. Compute the next collision time with particles in appropriate 
  // neighboring cells.

  collisionTime = SphereMobility::nextCollisionTime(collisionMsg_, this);
  boundaryCollisionTime = SphereMobility::nextBoundaryCollisionTime(collisionMsg_, this);

  // Step 5. Adjust the position of the event and its new partner’s event in the 
  // event queue. Since a wall collision changes the path of a particle, we only 
  // keep either a particle collision or a wall collision for each particle.

  if (collisionMsg_->isScheduled()) {

    if (collisionTime < collisionMsg_->getCollisionTime() && collisionTime != NO_TIME) {

      if (collisionMsg_->getPrevPartner() != NULL) {
        ((Sphere *)collisionMsg_->getPrevPartner())->getCollisionMessage()->setKind(EV_CHECK);
        ((Sphere *)collisionMsg_->getPrevPartner())->getCollisionMessage()->setPartner(NULL);
      }

      cancelEvent(collisionMsg_);

      collisionMsg_->setKind(EV_COLLISION);
      collisionMsg_->setCollisionTime(collisionTime);

      scheduleAt(collisionTime, collisionMsg_);

      if (collisionMsg_->getPartner() != NULL) {
        ((Sphere *)collisionMsg_->getPartner())->adjustCollision(collisionTime, this);
      }
    }

  } else {

    if (collisionTime < boundaryCollisionTime && collisionTime != NO_TIME) {

      collisionMsg_->setKind(EV_COLLISION);
      collisionMsg_->setCollisionTime(collisionTime);

      scheduleAt(collisionTime, collisionMsg_);

      if (collisionMsg_->getPartner() != NULL) {
        ((Sphere *)collisionMsg_->getPartner())->adjustCollision(collisionTime, this);
      }

    } else {

      collisionMsg_->setKind(EV_BOUNDARYCOLLISION);
      collisionMsg_->setCollisionTime(boundaryCollisionTime);

      collisionMsg_->setPartner(NULL);
      collisionMsg_->setPrevPartner(NULL);

      if (boundaryCollisionTime != NO_TIME) {
        scheduleAt(boundaryCollisionTime, collisionMsg_);
      }
    }
  }

  // Step 6. Return to Step 1.
}

/*
 * This method gets called when the sphere is a partner in a collision event.
 *
 * @param {double} newTime: the collision event time
 * @param {Particle *} from: the particle who is handling the collision
 */
void Sphere::adjustCollision(double newTime, Sphere *from) {

  // Methods called from other modules must have this macro
  Enter_Method_Silent();

  Particle* partner = NULL;
  CollisionMessage* partnerCollisionMsg = NULL;

  if (collisionMsg_->isScheduled()) cancelEvent(collisionMsg_);

  partner = dynamic_cast<Particle*>(collisionMsg_->getPartner());

  // Change the event type of the third party sphere to EV_CHECK
  if (partner != NULL && from->getParticleId() != partner->getParticleId()) {

    partnerCollisionMsg = (((Sphere *)partner)->getCollisionMessage());

    if (partnerCollisionMsg != NULL) {
      partnerCollisionMsg->setKind(EV_CHECK);
      partnerCollisionMsg->setPartner(NULL);
    }
  }

  SphereMobility::resetCollisionMessage(collisionMsg_);

  collisionMsg_->setKind(EV_CHECK);
  collisionMsg_->setCollisionTime(newTime);
  collisionMsg_->setPartner(from);

  scheduleAt(newTime, collisionMsg_);

}

/*
 *
 */
void Sphere::logCollisionTime(double stime) {

  double st = simTime().dbl();
  if (logCollisions_ && collisionTimeVector_ != NULL) {
    collisionTimeVector_->recordWithTimestamp(st, stime);
  }

  if (logCollisions_ && 
    xCollisionPositionVector_ != NULL &&
    yCollisionPositionVector_ != NULL &&
    zCollisionPositionVector_ != NULL) {

    xCollisionPositionVector_->recordWithTimestamp(st, position_.x);
    yCollisionPositionVector_->recordWithTimestamp(st, position_.y);
    zCollisionPositionVector_->recordWithTimestamp(st, position_.z);
  }
}

/*
 * Draws the shape of the cell in the tk environment.
 */
void Sphere::tkEnvDrawShape() {

  std::stringstream buffer;
  cModule *parent = getParentModule();

  // We will use the shape drawing tool to draw a circle around the particle
  // center
  buffer << 2*radius_;

  getDisplayString().setTagArg("b", 0, buffer.str().c_str());
  getDisplayString().setTagArg("b", 1, buffer.str().c_str());

  getDisplayString().setTagArg("b", 2, "oval");
  getDisplayString().setTagArg("b", 3, "white");
  getDisplayString().setTagArg("b", 4, "black");
  getDisplayString().setTagArg("b", 5, 1);

  if (parent != NULL && strcmp(parent->getName(), "cell") == 0) {
    parent->getDisplayString().setTagArg("b", 0, buffer.str().c_str());
    parent->getDisplayString().setTagArg("b", 1, buffer.str().c_str());

    parent->getDisplayString().setTagArg("b", 2, "oval");
    parent->getDisplayString().setTagArg("b", 3, "white");
    parent->getDisplayString().setTagArg("b", 4, "black");
    parent->getDisplayString().setTagArg("b", 5, 1);
  }
}

/*
 * Update the module position in the tk environment. This method is used when
 * the particle is initialized.
 */
void Sphere::tkEnvUpdatePosition() {

  std::stringstream buffer;
  cModule *parent = getParentModule();

  buffer << position_.y;

  // Set position string for tkenv
  if (strcmp(getName(), "molecule") == 0) {
    getDisplayString().setTagArg("p", 0, buffer.str().c_str());
  }

  if (parent != NULL && strcmp(parent->getName(), "cell") == 0) {
    parent->getDisplayString().setTagArg("p", 0, buffer.str().c_str());
  }

  buffer.str(std::string()); // clear buffer

  buffer << position_.x;

  if (strcmp(getName(), "molecule") == 0) {
    getDisplayString().setTagArg("p", 1, buffer.str().c_str());
  }

  if (parent != NULL && strcmp(parent->getName(), "cell") == 0) {
    parent->getDisplayString().setTagArg("p", 1, buffer.str().c_str());
  }

  buffer.str(std::string());
}

/*
 * Update the module position in the tk environment. This method is used 
 * during the simulation.
 *
 * @param {Double} t Time
 */
void Sphere::tkEnvUpdatePosition(double t) {

  std::stringstream buffer;
  cModule *parent = getParentModule();

  double lc = getLastCollisionTime();

  // Set position string for tkenv
  buffer << position_.y + velocity_.y*(t-lc);

  if (strcmp(getName(), "molecule") == 0) {
    getDisplayString().setTagArg("p", 0, buffer.str().c_str());
  }

  // Also move the parent's shape
  if (parent != NULL && strcmp(parent->getName(), "cell") == 0) {
    parent->getDisplayString().setTagArg("p", 0, buffer.str().c_str());
  }

  buffer.str(std::string()); // clear buffer

  buffer << position_.x + velocity_.x*(t-lc);

  if (strcmp(getName(), "molecule") == 0) {
    getDisplayString().setTagArg("p", 1, buffer.str().c_str());
  }

  // Also move the parent's shape
  if (parent != NULL && strcmp(parent->getName(), "cell") == 0) {
    parent->getDisplayString().setTagArg("p", 1, buffer.str().c_str());
  }

  buffer.str(std::string());
}