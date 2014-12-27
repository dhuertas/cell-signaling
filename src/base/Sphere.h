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

#ifndef SPHERE_H
#define SPHERE_H

#include "Manager.h"
#include "Particle.h"

#include "messages/Transfer_m.h"
#include "messages/Collision_m.h"
#include "messages/OutOfNeighborhood_m.h"
#include "messages/BrownianMotion_m.h"

class Sphere : public Particle {

 private:

  // Self messages
  TransferMessage *transferMsg_;
  CollisionMessage *collisionMsg_;

  BrownianMotionMessage *brownianMotionMsg_;

 protected:

  Manager *manager_;

  // Log collisions
  bool logCollisions_;

  // Log position
  bool logPosition_;

  // Track time between collisions
  cOutVector *collisionTimeVector_;

  // Vectors to get the mean free path
  cOutVector *xCollisionPositionVector_;
  cOutVector *yCollisionPositionVector_;
  cOutVector *zCollisionPositionVector_;

 public:

  Sphere();

  ~Sphere();

  void initializeMobilityMessages(void);

  void initializeMobility(void);

  void deleteMobilityMessages(void);

  void finishMobility(void);

  void finishMobility(Particle *from);

  void handleMobilityMessage(cMessage *);

  void adjustCollision(double time, Sphere *from);

  void logCollisionTime(double time);

  //
  // tk Environment related methods
  //
  void tkEnvDrawShape(void);

  void tkEnvUpdatePosition(void);

  void tkEnvUpdatePosition(double);

  //
  // Gets and sets
  //
  Manager * getManager(void) { return manager_; }
  
  TransferMessage * getTransferMessage(void) { return transferMsg_; }
  
  CollisionMessage * getCollisionMessage(void) { return collisionMsg_; }

};

#endif
