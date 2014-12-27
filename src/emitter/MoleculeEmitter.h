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

#ifndef MOLECULEEMITTER_H
#define MOLECULEEMITTER_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <coutvector.h>
#include <omnetpp.h>
#include <random.h>
#include <math.h>

#include <list>
#include <vector>


#include "Emitter.h"
#include "SimpleCell.h"
#include "Manager.h"

class MoleculeEmitter : public Emitter, public cSimpleModule {

 private:

  Manager *manager_;

  SimpleCell *mobility_;

  bool checkOverlap(point3_t, double);

  std::list<Molecule *> preloadedMolecules_;

 protected:

  double emissionStart_;

  double emissionDuration_;

  double emissionRate_;

  double emissionParticleRadius_;

  double emissionParticleMass_;

  double emissionTimeToLive_;

  uint8_t emissionFunction_;

  int emissionBoundariesMode_;

  double emissionVelocity_;

  double emissionListRadius_;

  double emissionRefreshListRadius_;

  double emissionDiffusion_;

  bool preloadMolecules_;

 public:

  MoleculeEmitter();

  ~MoleculeEmitter();

  virtual void initialize(int);

  virtual int numInitStages() const;

  virtual void handleMessage(cMessage *);

  virtual void finish();

  Molecule *createMolecule();

  void setManager(std::string);

};

#endif
