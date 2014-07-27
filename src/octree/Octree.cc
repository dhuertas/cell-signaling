#include "Octree.h"

#include <math.h>

Octree::Octree() {

}

Octree::Octree(unsigned int numLayers, double maxSideLength) {

  numLayers_ = numLayers;

  maxSideLength_ = maxSideLength;

  spaceCellsTree_.resize(numLayers_);

}

Octree::~Octree() {

}

/*
 * Returns the side length of a given cell and layer
 * D_cell = D_max/(2^layer)
 */
double Octree::getSpaceCellSideLength(index_t idx) {

  return maxSideLength_/(1 << idx.l); 
}

/*
 * Returns the list of particles from a given space cell 
 */
std::list<Particle *> *Octree::getSpaceCellParticles(index_t idx) {

  //unsigned int Nx = ceil(spaceSize_.x/(maxSideLength_/(1 << idx.l)));
  unsigned int Ny = ceil(spaceSize_.y/(maxSideLength_/(1 << idx.l)));
  unsigned int Nz = ceil(spaceSize_.z/(maxSideLength_/(1 << idx.l)));

  // TODO finish computing list index
  return &spaceCellsTree_[idx.l][(Ny*Nz)*idx.i + Nz*idx.j + idx.k];
}

void Octree::transferParticle(Particle *p, index_t from, index_t to) {

}

void Octree::subscribe(Particle *p) {

}

void Octree::unsubscribe(Particle *p) {

}

void Octree::attachParticleToSpaceCell(Particle *, index_t idx) {

}

void Octree::detachParticleFromSpaceCell(Particle *, index_t idx) {

}