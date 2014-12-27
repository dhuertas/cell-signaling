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

#ifndef SPHEREMOBILITY_H
#define SPHEREMOBILITY_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <omnetpp.h>

#include "Mobility.h"
#include "Manager.h"
#include "Sphere.h"

#include "../messages/Transfer_m.h"
#include "../messages/Collision_m.h"
#include "../messages/OutOfNeighborhood_m.h"
#include "../messages/BrownianMotion_m.h"

class SphereMobility : public Mobility {

 public:

  static double nextCollisionTime(CollisionMessage *msg, Sphere *s);

  static double nextBoundaryCollisionTime(CollisionMessage *msg, Sphere *s);

  static double leaveBoundedSpace(CollisionMessage *msg, Sphere *s);

  static double solveCollision(Sphere *sa, Sphere *sb);

  static double brownianMotion(BrownianMotionMessage *msg, Sphere *s);

  static void handleCollision(CollisionMessage *msg, Sphere *s);

  static void handleBoundaryCollision(CollisionMessage *msg, Sphere *s);

  static void handleWallCollision(CollisionMessage *msg, Sphere *s);

  static void handlePeriodicBoundary(CollisionMessage *msg, Sphere *s);

  static void handleBrownianMotion(BrownianMotionMessage *msg, Sphere *s);

};

#endif
