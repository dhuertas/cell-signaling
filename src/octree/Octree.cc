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

#include "Octree.h"

#include "Defines.h"
#include <math.h>

/*
 *
 */
Octree::Octree() {

}

/*
 *
 */
Octree::Octree(unsigned int numLayers, double maxSideLength) {

  numLayers_ = numLayers;

  maxSpaceCellSize_ = maxSideLength;

  spaceCellsTree_.resize(numLayers_);

}

/*
 *
 */
Octree::~Octree() {

}

/*
 * Returns the list of particles from a given space cell 
 */
std::list<Particle *> *Octree::getSpaceCellParticles(index_t idx) {

  //unsigned int Nx = ceil(spaceSize_.x/(maxSpaceCellSize_/(1 << idx.l)));
  unsigned int Ny = ceil(spaceSize_.y/(maxSpaceCellSize_/(1 << idx.layer)));
  unsigned int Nz = ceil(spaceSize_.z/(maxSpaceCellSize_/(1 << idx.layer)));

  return &spaceCellsTree_[idx.layer][(Ny*Nz)*idx.i + Nz*idx.j + idx.k];
}

/*
 *
 */
std::vector<Particle *> *Octree::getNeighborParticles(index_t idx, std::vector<Particle *> *list) {

  // Get the particles from each of the listed space cells
  //for (sc = spaceCells.begin(); sc != spaceCells.end(); ++sc) {
  //
  //  particleList = manager->getSpaceCellParticles(*sc);
  //  particles.insert(particles.end(), particleList->begin(), particleList->end());
  //}

  return list;
}

/*
 *
 */
void Octree::transferParticle(Particle *p, index_t from, index_t to) {

  // Detach particle from its space cell
  detachParticleFromSpaceCell(p, from);

  // Attach the particle to its new space cell
  attachParticleToSpaceCell(p, to);
}

/*
 *
 */
void Octree::subscribe(Particle *p) {

}

/*
 *
 */
void Octree::unsubscribe(Particle *p) {

}

/*
 *
 */
void Octree::attachParticleToSpaceCell(Particle *p, index_t idx) {

  point_t *pos = NULL;

  unsigned int Nx, Ny, Nz; 
  unsigned long n;

  Nx = Ny = Nz = 0;
  n = 0;

  if (IDX_ISNULL(idx)) {

    // Set enable
    idx.flags |= IDX_ENABLED;

    idx.layer = floor(log(maxSpaceCellSize_/(2*p->getRadius()))/LOG2);

    Nx = ceil(spaceSize_.x/(maxSpaceCellSize_/(1 << idx.layer)));
    Ny = ceil(spaceSize_.y/(maxSpaceCellSize_/(1 << idx.layer)));
    Nz = ceil(spaceSize_.z/(maxSpaceCellSize_/(1 << idx.layer)));

    pos = p->getPosition();

    // Note that the integer already performs the floor() operation
    //idx.i = floor(pos->x/(maxSpaceCellSize_/(1 << idx.layer)));
    //idx.j = floor(pos->y/(maxSpaceCellSize_/(1 << idx.layer)));
    //idx.k = floor(pos->z/(maxSpaceCellSize_/(1 << idx.layer)));

    idx.i = pos->x/(maxSpaceCellSize_/(1 << idx.layer));
    idx.j = pos->y/(maxSpaceCellSize_/(1 << idx.layer));
    idx.k = pos->z/(maxSpaceCellSize_/(1 << idx.layer));

    n = Ny*Nz*idx.i + Nz*idx.j + idx.k;

    // Initialize particle space cell index
    p->setSpaceCellIdx(idx);
    p->setPrevSpaceCellIdx(IDX_NULL);

  } else {

    Nx = ceil(spaceSize_.x/(maxSpaceCellSize_/(1 << idx.layer)));
    Ny = ceil(spaceSize_.y/(maxSpaceCellSize_/(1 << idx.layer)));
    Nz = ceil(spaceSize_.z/(maxSpaceCellSize_/(1 << idx.layer)));

    n = Ny*Nz*idx.i + Nz*idx.j + idx.k;
  }

  spaceCellsTree_[idx.layer].at(n).push_back(p);

}

/*
 *
 */
void Octree::detachParticleFromSpaceCell(Particle *p, index_t idx) {

  unsigned int Nx, Ny, Nz; 
  unsigned long n; 

  Nx = Ny = Nz = 0;
  n = 0;

  if (IDX_ISNULL(idx)) {
    idx = p->getSpaceCellIdx();
  }

  Nx = ceil(spaceSize_.x/(maxSpaceCellSize_/(1 << idx.layer)));
  Ny = ceil(spaceSize_.y/(maxSpaceCellSize_/(1 << idx.layer)));
  Nz = ceil(spaceSize_.z/(maxSpaceCellSize_/(1 << idx.layer)));

  n = Ny*Nz*idx.i + Nz*idx.j + idx.k;

  spaceCellsTree_[idx.layer].at(n).remove(p);
}