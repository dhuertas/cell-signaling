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

#ifndef COMPUTE_H
#define COMPUTE_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <omnetpp.h>

#include "Manager.h"
#include "Particle.h"

#include "../messages/Transfer_m.h"
#include "../messages/Collision_m.h"
#include "../messages/OutOfNeighborhood_m.h"

class Mobility {

 public:

  static uint8_t sides_[];

  static double nextTransferTime(TransferMessage *msg, Particle *p);

  static void handleTransfer(TransferMessage *msg, Particle *p);

  static void resetCollisionMessage(CollisionMessage *msg);

};

#endif
