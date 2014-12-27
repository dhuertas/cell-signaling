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

#ifndef PARTICLEDISTRIBUTION_H
#define PARTICLEDISTRIBUTION_H

#include <stdlib.h>
#include <time.h>
#include <list>
#include <cmath>
#include <algorithm>

// Omnet++
#include <random.h>

// project
#include "Common.h"
#include "Particle.h"


// Functions that enable the Manager module to place particles following 
// different distributions.

// Place each particle at a random position.
void uniformDistribution(vector3_t, std::list<Particle *> *);
void uniformDistribution2(vector3_t, std::list<Particle *> *);
void uniformDistribution3(vector3_t, std::list<Particle *> *);

// Place each particle following a cube pattern.
void cubeDistribution(vector3_t, std::list<Particle *> *);

// Place each particle randomly on a sphere surface.
void sphereDistribution(vector3_t, std::list<Particle *> *, point3_t, double);

// Place each particle equally distributed on a sphere surface
void sphereEquallyDistributed(vector3_t, std::list<Particle *> *, point3_t, double);

// Put all the particles in the center of the simulation space one near the other.
void highDensityDistribution(vector3_t, std::list<Particle *> *, point3_t);

void densepacked(vector3_t, std::list<Particle *> *, point3_t);

// Detect whether two sphere particles are overlaping or not.
bool checkOverlap(point3_t*, double, point3_t*, double);

double dblRandNormal(double, double);

#endif