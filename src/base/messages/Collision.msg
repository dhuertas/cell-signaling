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

cplusplus {{
  #include "Common.h"
  #include "Manager.h"
  #include "Particle.h"

  typedef Manager *ManagerPtr;
  typedef Particle *ParticlePtr;

  typedef point3_t Point;
  typedef vector3_t Vector;

  //typedef uint8_t Hits;
}}

class noncobject ManagerPtr;
class noncobject ParticlePtr;
class noncobject Point;
class noncobject Vector;
//class noncobject Hits;

message CollisionMessage {

  ManagerPtr manager;

  Point position;

  Vector velocity;

  double collisionTime;

  int hits;

  ParticlePtr partner;
  ParticlePtr prevPartner;
}