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

#include "MoleculeEmitter.h"

Define_Module(MoleculeEmitter);

MoleculeEmitter::MoleculeEmitter() : Emitter() {

}

MoleculeEmitter::~MoleculeEmitter() {

}

/*
 * Initialize the molecule emitter
 */
void MoleculeEmitter::initialize(int stage) {

  if (stage == 0) {
    // First stage manager initialization
    enabled_ = par("enabled");

  } else if (stage == 1) {
    // Particles initialization
  } else if (stage == 2) {
    // Second stage manager initialization
  } else if (stage == 3) {
    // Particle emitters and receivers initialization

    // Emission variables
    emissionStart_ = par("emissionStart");
    emissionDuration_ = par("emissionDuration");
    emissionRate_ = par("emissionRate");
    emissionParticleRadius_ = par("emissionParticleRadius");
    emissionParticleMass_ = par("emissionParticleMass");
    emissionTimeToLive_ = par("emissionTimeToLive");
    emissionBoundariesMode_ = par("emissionBoundariesMode");
    emissionVelocity_ = par("emissionVelocity");
    emissionDiffusion_ = par("emissionDiffusion");

    preloadMolecules_ = par("preloadMolecules");

    mobility_ = (SimpleCell *)getParentModule()->getSubmodule("mobility");

    if (enabled_ && getParentModule()->getSubmodule("receiver")->par("enabled")) {
      mobility_->setParticleType(T_EMITTER_RECEIVER);
    } else if (enabled_) {
      mobility_->setParticleType(T_EMITTER);
    }

    if (preloadMolecules_) {
      // Do not create more modules during the initialization process,
      // since omnet will automatically initialize them...
      scheduleAt(simTime(), new cMessage("preload", EV_PRELOAD));
    }

    if (emissionStart_ > 0) {
      scheduleAt(simTime() + emissionStart_, new cMessage("emit", EV_EMIT));
    }

    // Subscribe to manager
    setManager("manager");
  }
}

/*
 * Returns the number of initialization stages.
 */
int MoleculeEmitter::numInitStages() const {
  return 4;
}

/*
 *
 */
void MoleculeEmitter::handleMessage(cMessage *msg) {

  int kind = msg->getKind();
  double st = NO_TIME;

  Molecule *molecule;

  if (kind == EV_EMIT) {

    // TODO improve this
    st = simTime().dbl();

    molecule = createMolecule();
    molecule->callInitialize();

    molecule->setParticleType(T_SIGNALING);
    molecule->setLastCollisionTime(st);

    molecule->initializeMobilityMessages();

    molecule->initializeMobility();

    if (st < emissionStart_ + emissionDuration_) {
      scheduleAt(st + 1/emissionRate_, msg);
    }

  } else if (kind == EV_PRELOAD) {

    // Compute the number of molecules to be released
    uint32_t numberOfMolecules = ceil(emissionRate_*emissionDuration_);

    // Create molecules and save them for later use
    for (uint32_t i = 0; i < numberOfMolecules; i++) {

      cModuleType *moduleType = cModuleType::get("cellsignaling.src.Molecule");
      Molecule *m = (Molecule *)moduleType->create("molecule", simulation.getSystemModule());

      preloadedMolecules_.push_back(m);
    }
  }
}

/*
 *
 */
Molecule * MoleculeEmitter::createMolecule() {

  bool overlap;

  double dt;
  double theta, phi;
  double r, e;
  double epr;
  double vm;

  point3_t pos, c;
  vector3_t *ss;
  vector3_t v;

  point3_t *mpos = NULL;
  vector3_t *mvel = NULL;

  // force enter the first while loop
  overlap = true;

  epr = emissionParticleRadius_;
  e = 0.01*epr;
  r = mobility_->getRadius() + epr + e;

  ss = manager_->getSpaceSize();

  mpos = mobility_->getPosition();
  mvel = mobility_->getVelocity();

  // Create molecule or get a preloaded one
  Molecule *m = NULL;

  if (preloadedMolecules_.size() > 0) {
    m = preloadedMolecules_.front();
    preloadedMolecules_.pop_front(); // remove first molecule from the list
  } else {
    cModuleType *moduleType = cModuleType::get("cellsignaling.src.Molecule");
    m = (Molecule *)moduleType->create("molecule", simulation.getSystemModule()); 
  }

  dt = simTime().dbl() - mobility_->getLastCollisionTime();

  // set up parameters and gate sizes before we set up its submodules
  theta = dblrand()*M_PI;
  phi = dblrand()*M_PI*2;

  // force enter the first while loop
  pos.x = 0;
  pos.y = 0;
  pos.z = 0;

  c.x = mpos->x + mvel->x*dt;
  c.y = mpos->y + mvel->y*dt;
  c.z = mpos->z + mvel->z*dt;

  while (pos.x - epr <= 0 || pos.x + epr >= ss->x ||
    pos.y - epr <= 0 || pos.y + epr >= ss->y ||
    pos.z - epr <= 0 || pos.z + epr >= ss->z ||
    overlap) {

    theta = dblrand()*M_PI;
    phi = dblrand()*M_PI*2;

    pos.x = c.x + r*sin(theta)*cos(phi);
    pos.y = c.y + r*sin(theta)*sin(phi);
    pos.z = c.z + r*cos(theta);

    // Check whether the emitted molecule overlaps with surrounding
    // particles
    // TODO rethink this part since it is too costly
    overlap = checkOverlap(pos, emissionParticleRadius_);

  }

  m->par("xpos") = pos.x;
  m->par("ypos") = pos.y;
  m->par("zpos") = pos.z;

  v.x = pos.x - c.x;
  v.y = pos.y - c.y;
  v.z = pos.z - c.z;

  vm = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);

  v.x = emissionVelocity_*v.x/vm;
  v.y = emissionVelocity_*v.y/vm;
  v.z = emissionVelocity_*v.z/vm;

  // Particle parameters
  m->par("vx") = v.x;
  m->par("vy") = v.y;
  m->par("vz") = v.z;

  m->par("mass") = emissionParticleMass_;
  m->par("radius") = emissionParticleRadius_;

  m->par("timeToLive") = emissionTimeToLive_;

  m->par("boundariesMode") = emissionBoundariesMode_;
  // m->par("statsRefreshRate");

  m->par("diffusion") = emissionDiffusion_;

  m->finalizeParameters();

  // m->setGateSize("in", 1);
  // m->setGateSize("out", 1);

  // create internals, and schedule it
  m->buildInside();
  m->scheduleStart(simTime());

  return m;

}

/*
 *
 */
void MoleculeEmitter::finish() {

  while (preloadedMolecules_.size() > 0) {
    preloadedMolecules_.front()->deleteModule();
    preloadedMolecules_.pop_front();
  }
}

/*
 * Sets the manager attribute.
 * 
 * @param {string} param: the name of the manager module
 */
void MoleculeEmitter::setManager(std::string param) {

  try {
    manager_ = (Manager *)simulation.
      getSystemModule()->getSubmodule(param.c_str());
  } catch (cException *e) {
    EV << "setManager error" << "\n";
  }

}

/*
 * Function to check whether an emitted particle overlaps.
 *
 * @return {bool}
 */
bool MoleculeEmitter::checkOverlap(point3_t pos, double radius) {

  double dx, dy, dz;
  double rb;

  point3_t *cb = NULL;

  index3_t idx = {0, 0, 0, 0, 0};

  std::vector<Particle *> particles;
  std::vector<Particle *>::iterator p;

  double maxSpaceSize = manager_->getMaxSpaceSize();
  unsigned int maxDepth = manager_->getDepth();

  while (2*radius / (maxSpaceSize / (1 << (idx.depth+1))) <= 1 && idx.depth+1 <= maxDepth) {
    idx.depth++;
  }

  idx.i = floor(pos.x / (maxSpaceSize / (1 << idx.depth)));
  idx.j = floor(pos.y / (maxSpaceSize / (1 << idx.depth)));
  idx.k = floor(pos.z / (maxSpaceSize / (1 << idx.depth)));

  manager_->getNeighborParticles(&idx, &particles);

  for (p = particles.begin(); p != particles.end(); ++p) {

    cb = (*p)->getPosition();
    rb = (*p)->getRadius();

    dx = pos.x - cb->x;
    dy = pos.y - cb->y;
    dz = pos.z - cb->z;

    if (sqrt(dx*dx + dy*dy + dz*dz) <= radius + rb) {
      return true;
    }
  }

  // TODO check for wall overlap
  return false;
}
