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
Particle::Particle() :
  manager_(NULL),
  active_(false),
  mass_(0.0),
  radius_(0.0),
  lastCollisionTime_(0.0),
  particleId_(0),
  particleType_(0),
  boundariesMode_(0),
  imageIdx_(0),
  diffusion_(0),
  brownianMotionStdDev_(0.0) {

  memset(&position_, 0, sizeof(point3_t));
  memset(&velocity_, 0, sizeof(vector3_t));
  memset(&spaceCellIdx_, 0, sizeof(index3_t));
}

/*
 *
 */
void Particle::initialize(int stage) {

}

/*
 *
 */
int Particle::numInitStages() const {
  return 0;
}

/*
 *
 */
void Particle::handleMessage(cMessage *) {

}

/*
 *
 */
void Particle::finish() {

}

/*
 * Sets the manager attribute.
 * 
 * @param {string} param: the name of the manager module
 */
void Particle::setManager(std::string name) {

  try {
    manager_ = (Manager *)simulation.
      getSystemModule()->getSubmodule(name.c_str());
  } catch (cException *e) {
    EV << "setManager error" << "\n";
  }
}