#ifndef PARTICLEDISTRIBUTION_H
#define PARTICLEDISTRIBUTION_H

#include <stdlib.h>
#include <time.h>
#include <list>
#include <cmath>

#include "Defines.h"
#include "Particle.h"


// Functions that enable the Manager module to place particles following 
// different distributions.

// Place each particle at a random position.
void randomDistribution(vect_t, std::list<Particle *> *);

// Place each particle following a cube pattern.
void cubeDistribution(vect_t, std::list<Particle *> *);

// Place each particle following a sphere surface.
void sphereDistribution(vect_t, std::list<Particle *> *, point_t, double);

// Create two distant groups of particles. The particles in each group are
// placed randomly.
void twoSidedDistribution(vect_t, std::list<Particle *> *, point_t, point_t);

// Put all the particles in the center of the simulation space one near the other.
void highDensityDistribution(vect_t, std::list<Particle *> *, point_t);

// Detect whether two sphere particles are overlaping or not.
bool checkOverlap(point_t, double, point_t, double);

#endif