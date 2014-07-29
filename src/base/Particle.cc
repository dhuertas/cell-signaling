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

#include "Particle.h"

using namespace std;

/*
 * Constructor.
 */
Particle::Particle() {}

/*
 * Constructor overload.
 */
Particle::Particle(
  double x,
  double y,
  double vx,
  double vy,
  double m) {

  position_.x = x;
  position_.y = y;
  position_.z = 0;

  velocity_.x = vx;
  velocity_.y = vy;
  velocity_.z = 0;

  spaceCell_ = -1;

  spaceCellIdx_.flags = 0x00;
  spaceCellIdx_.i = 0;
  spaceCellIdx_.j = 0;
  spaceCellIdx_.k = 0;
  spaceCellIdx_.layer = 0;

  prevSpaceCellIdx_.flags = 0x00;
  prevSpaceCellIdx_.i = 0;
  prevSpaceCellIdx_.j = 0;
  prevSpaceCellIdx_.k = 0;
  spaceCellIdx_.layer = 0;


  mass_ = m;
  lastCollisionTime_ = 0;

  listRadius_ = 1;

}