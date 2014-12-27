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

#include "Probe.h"

using namespace std;

Define_Module(Probe);

Probe::~Probe() {}

/*
 *
 */
void Probe::initialize(int stage) {
  
  if (stage == 0) {
    // First stage manager initialization
  } else if (stage == 1) {
    // Particles initialization
    // Probes initialization
    
    name_ = par("name").stringValue();
    
    position_.x = par("xpos").doubleValue();
    position_.y = par("ypos").doubleValue();
    position_.z = par("zpos").doubleValue();

    radius_ = par("radius").doubleValue();

    statsRefreshRate_ = par("statsRefreshRate");

    type_ = par("type");

    moleculeDensityVector.setName(name_.c_str());

    setManager("manager");

    if (statsRefreshRate_ > 0) {
      scheduleAt(simTime() + statsRefreshRate_/1000,
        new cMessage("refresh", EV_STATSUPDATE));
    }   

  }
}

/*
 *
 */
int Probe::numInitStages() const {

  return 2;

}

/*
 *
 */
void Probe::handleMessage(cMessage *msg) {

  int kind = msg->getKind();

  simtime_t st = simTime();

  if (kind == EV_STATSUPDATE) {

    moleculeDensityVector.recordWithTimestamp(st, getMoleculeDensity());

    if (statsRefreshRate_ > 0) {
      scheduleAt(st + statsRefreshRate_/1000, msg);
    }
  }
}

/*
 *
 */
void Probe::finish() {

}

/*
 *
 */
double Probe::getMoleculeDensity() {
  // Find cell coordinates that fall within the given sphere (with position 
  // pos and radius r)
  //int a, b, c, l;
  //int n;

  //index3_t idx;

  //int Nx, Ny, Nz;

  //int count = 0;

  // TODO rewrite this

  // return the quotient of the counter and the volume of the probe
  return 0;
  //count/((4.0/3.0)*M_PI*radius_*radius_*radius_);

}

/*
 *
 */
void Probe::setManager(std::string param) {

  try {
    manager_ = (Manager *)simulation.
      getSystemModule()->getSubmodule(param.c_str());
  } catch (cException *e) {
    EV << "setManager error" << "\n";
  }

}
