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

#ifndef OCTREE_H
#define OCTREE_H

#include <cmessage.h>
#include <csimplemodule.h>
#include <cqueue.h>
#include <omnetpp.h>

#include <vector>
#include <list>

#include "Particle.h"

class Octree {

  typedef std::vector<std::list<Particle*> > layer_t;

 private:

 protected:

  unsigned int numLayers_;

  double maxSpaceCellSize_;

  vect_t spaceSize_;

  std::vector<layer_t> spaceCellsTree_;

 public:

  Octree();

  Octree(unsigned int nl, double msl);

  ~Octree();

  std::list<Particle *> *getSpaceCellParticles(index_t idx);

  std::list<Particle *> getNeighborParticles(index_t idx);

  void transferParticle(Particle *p, index_t from, index_t to);

  void subscribe(Particle *p);

  void unsubscribe(Particle *p);

  void attachParticleToSpaceCell(Particle *, index_t idx);

  void detachParticleFromSpaceCell(Particle *, index_t idx);

  //
  // gets and sets
  //
  inline unsigned int getNumLayers(void) { return numLayers_; };

  inline double getMaxSpaceCellSize(void) { return maxSpaceCellSize_; };

  inline vect_t *getSpaceSize(void) { return &spaceSize_; };

  // Returns the side length of a given cell and layer
  // D_cell = D_max/(2^layer)
  inline double getSpaceCellSideLength(index_t idx) {
    return maxSpaceCellSize_/(1 << idx.layer);
  }

  inline void setNumLayers(unsigned int nl) { numLayers_ = nl; };

  inline void setMaxSpaceCellSize(double sl) { maxSpaceCellSize_ = sl; };

  inline void setSpaceSize(vect_t ss) { spaceSize_ = ss; }; 

};

#endif