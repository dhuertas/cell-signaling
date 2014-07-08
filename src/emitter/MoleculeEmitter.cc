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
		enabled = par("enabled");

	} else if (stage == 1) {
		// Particles initialization
	} else if (stage == 2) {
		// Second stage manager initialization
	} else if (stage == 3) {
		// Particle emitters and receivers initialization

		// Emission variables
		emissionStart = par("emissionStart");
		emissionDuration = par("emissionDuration");
		emissionRate = par("emissionRate");
		emissionParticleRadius = par("emissionParticleRadius");
		emissionParticleMass = par("emissionParticleMass");
		emissionTimeToLive = par("emissionTimeToLive");
		emissionBoundariesMode = par("emissionBoundariesMode");
		emissionVelocity = par("emissionVelocity");
		emissionListRadius = par("emissionListRadius");
		emissionRefreshListRadius = par("emissionRefreshListRadius");
		emissionDiffusion = par("emissionDiffusion");

		preloadMolecules = par("preloadMolecules");

		mobility = (SimpleCell *)getParentModule()->getSubmodule("mobility");

		if (enabled && getParentModule()->getSubmodule("receiver")->par("enabled")) {
			mobility->setParticleType(T_EMITTER_RECEIVER);
		} else if (enabled) {
			mobility->setParticleType(T_EMITTER);
		}

		if (preloadMolecules) {
			// Compute the number of molecules to be released
			uint32_t numberOfMolecules = ceil(emissionRate*emissionDuration);

			// Create molecules and save them for later use
			for (uint32_t i = 0; i < numberOfMolecules; i++) {

				cModuleType *moduleType = cModuleType::get("cellsignaling.src.Molecule");
				Molecule *m = (Molecule *)moduleType->create("molecule", simulation.getSystemModule());

				preloadedMolecules.push_back(m);
			}
		}

		if (emissionStart > 0) {
			scheduleAt(simTime() + emissionStart, 
			new cMessage("emit", EV_EMIT));
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
	double st;

	Molecule *molecule;

	if (kind == EV_EMIT) {

		st = simTime().dbl();

		molecule = createMolecule();
		molecule->callInitialize();
		
		molecule->setParticleId(manager->getNextParticleId());
		molecule->setParticleType(T_SIGNALING);
		molecule->setLastCollisionTime(st);

		molecule->initMobilityMessages();

		manager->attachParticleToSpaceCell(molecule, -1);

		if (manager->getMode() == M_NNLIST) {
			molecule->setListRadius(2*emissionParticleRadius);
			molecule->createNearNeighborList();
		}

		molecule->initializeMobility();

		if (st < emissionStart + emissionDuration) {
			scheduleAt(st + 1/emissionRate, msg);
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

	point_t pos, c;
	vect_t *ss;
	vect_t v;

	point_t *mpos = NULL;
	vect_t *mvel = NULL;

	// force enter the first while loop
	overlap = true;

	epr = emissionParticleRadius;
	e = 0.01*epr;
	r = mobility->getRadius() + epr + e;

	ss = manager->getSpaceSize();

	mpos = mobility->getPosition();
	mvel = mobility->getVelocity();

	// Create molecule or get a preloaded one
	Molecule *m = NULL;
	if (preloadedMolecules.size() > 0) {
		m = preloadedMolecules.front();
		preloadedMolecules.pop_front(); // remove first molecule from the list
	} else {
		cModuleType *moduleType = cModuleType::get("cellsignaling.src.Molecule");
		m = (Molecule *)moduleType->create("molecule", simulation.getSystemModule());	
	}

	dt = simTime().dbl() - mobility->getLastCollisionTime();

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
		overlap = checkOverlap(pos, emissionParticleRadius);

	}

	m->par("xpos") = pos.x;
	m->par("ypos") = pos.y;
	m->par("zpos") = pos.z;

	v.x = pos.x - c.x;
	v.y = pos.y - c.y;
	v.z = pos.z - c.z;

	vm = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);

	v.x = emissionVelocity*v.x/vm;
	v.y = emissionVelocity*v.y/vm;
	v.z = emissionVelocity*v.z/vm;

	// Particle parameters
	m->par("vx") = v.x;
	m->par("vy") = v.y;
	m->par("vz") = v.z;

	m->par("mass") = emissionParticleMass;
	m->par("radius") = emissionParticleRadius;

	m->par("timeToLive") = emissionTimeToLive;

	// Verlet List parameters
	m->par("listRadius") = emissionListRadius;
	m->par("refreshListRadius") = emissionRefreshListRadius;

	m->par("boundariesMode") = emissionBoundariesMode;
	// m->par("statsRefreshRate");

	m->par("diffusion") = emissionDiffusion;

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

	while (preloadedMolecules.size() > 0) {
		preloadedMolecules.front()->deleteModule();
		preloadedMolecules.pop_front();
	}
}

/*
 * Sets the manager attribute.
 * 
 * @param {string} param: the name of the manager module
 */
void MoleculeEmitter::setManager(std::string param) {

	try {
		manager = (Manager *)simulation.
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
bool MoleculeEmitter::checkOverlap(point_t ca, double ra) {

	int a, b, c;			// Nested "for" loops
	int i, j, k;			// Indexes to access the current space cell
	int *Nx, *Ny, *Nz;		// Number of space cells (or divisions) in each axis

	double dx, dy, dz;
	double rb;
	double spaceCellSize;

	point_t *cb = NULL;

	std::vector<Particle *> particles;
	std::vector<Particle *>::iterator p;

	std::vector<int> spaceCells;
	std::vector<int>::iterator sc;

	std::list<Particle *> *particleList;
	std::list<Particle *>::const_iterator pl;

	spaceCellSize = manager->getSpaceCellSize();

	Nx = manager->getNumberOfSpaceCellsX();
	Ny = manager->getNumberOfSpaceCellsY();
	Nz = manager->getNumberOfSpaceCellsZ();

	i = floor(ca.x/spaceCellSize);
	j = floor(ca.y/spaceCellSize);
	k = floor(ca.z/spaceCellSize);

	if (manager->getMode() == M_NNLIST) {

		particles = mobility->getNeighborParticles();

	} else {

		for (a = -1; a <= 1; a++)
		for (b = -1; b <= 1; b++)
		for (c = -1; c <= 1; c++) {

			if (CELLBELONGSTOSIMSPACE(i+a, j+b, k+c, *Nx, *Ny, *Nz)) {
				spaceCells.push_back((i+a)*(*Ny)*(*Nz)+(j+b)*(*Nz)+(k+c));
			}

		}

		// Get the particles from each of the listed space cells
		for (sc = spaceCells.begin(); sc != spaceCells.end(); ++sc) {

			particleList = manager->getSpaceCellParticles(*sc);

			for (pl = particleList->begin(); pl != particleList->end(); ++pl) {
				particles.push_back(*pl);
			}

		}

		for (p = particles.begin(); p != particles.end(); ++p) {

			cb = (*p)->getPosition();
			rb = (*p)->getRadius();

			dx = ca.x - cb->x;
			dy = ca.y - cb->y;
			dz = ca.z - cb->z;

			if (sqrt(dx*dx + dy*dy + dz*dz) <= ra + rb) {
				return true;
			}

		}

	}

	return false;

}
