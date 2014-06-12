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
		
		name = par("name").stringValue();
		
		position.x = par("xpos").doubleValue();
		position.y = par("ypos").doubleValue();
		position.z = par("zpos").doubleValue();

		radius = par("radius").doubleValue();

		statsRefreshRate = par("statsRefreshRate");

		type = par("type");

		moleculeDensityVector.setName(name.c_str());

		setManager("manager");

		if (statsRefreshRate > 0) {
			scheduleAt(simTime() + statsRefreshRate/1000,
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

		if (statsRefreshRate > 0) {
			scheduleAt(st + statsRefreshRate/1000, msg);
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
	int a, b, c, l;
	int i, j, k, n;

	int *Nx, *Ny, *Nz;

	int count = 0;

	double spaceCellSize = 0.0;

	point_t *ppos = NULL;

	list<Particle *> *particleList;
	list<Particle *>::const_iterator p;

	Nx = manager->getNumberOfSpaceCellsX();
	Ny = manager->getNumberOfSpaceCellsY();
	Nz = manager->getNumberOfSpaceCellsZ();

	spaceCellSize = manager->getSpaceCellSize();

	i = floor(position.x/spaceCellSize);
	j = floor(position.y/spaceCellSize);
	k = floor(position.z/spaceCellSize);

	n = i*(*Ny)*(*Nz) + j*(*Nz) + k;

	l = ceil(radius/spaceCellSize);

	// Loop each cell list looking for the particles that fall in the sphere
	for (a = -l; a <= l; a++)
	for (b = -l; b <= l; b++)
	for (c = -l; c <= l; c++) {
		if (CELLBELONGSTOSIMSPACE(i+a, j+b, k+c, *Nx, *Ny, *Nz)) {

			n = (i+a)*(*Ny)*(*Nz) + (j+b)*(*Nz) + (k+c);
			particleList = manager->getSpaceCellParticles(n);

			for (p = particleList->begin(); p != particleList->end(); ++p) {

				ppos = (*p)->getPosition();

				if ((*p)->getParticleType() == type && 
					((position.x-ppos->x)*(position.x-ppos->x) + 
					(position.y-ppos->y)*(position.y-ppos->y) +
					(position.z-ppos->z)*(position.z-ppos->z) <= radius*radius)) {

					count++;
				}
			}
		}
	}

	// return the quotient of the counter and the volume of the probe
	return count/((4.0/3.0)*M_PI*radius*radius*radius);

}

/*
 *
 */
void Probe::setManager(std::string param) {

	try {
		manager = (Manager *)simulation.
			getSystemModule()->getSubmodule(param.c_str());
	} catch (cException *e) {
		EV << "setManager error" << "\n";
	}

}
