#include "ParticleDistribution.h"

using namespace std;
/*
 * Place each particle at a random position.
 *
 * @param {vect_t} spaceSize: the simulation space (x, y and z)
 * @param {std::list<Particle *> *} particles
 */
void randomDistribution(vect_t spaceSize, std::list<Particle *> *particles) {

	bool overlap;

	point_t pos;
	vect_t range;

	std::list<Particle *>::iterator p, q;

	p = particles->begin();

	srand(time(NULL));

	while (p != particles->end()) {

		overlap = false;
		
		range.x = spaceSize.x - 2*(*p)->getRadius();
		range.y = spaceSize.y - 2*(*p)->getRadius();
		range.z = spaceSize.z - 2*(*p)->getRadius();

		pos.x = (*p)->getRadius() + range.x*rand()/RAND_MAX;
		pos.y = (*p)->getRadius() + range.y*rand()/RAND_MAX;
		pos.z = (*p)->getRadius() + range.z*rand()/RAND_MAX;

		for (q = particles->begin(); q != p; ++q) {

			if (checkOverlap(pos, (*p)->getRadius(), (*q)->getPosition(), (*q)->getRadius())) {
				overlap = true;
				break;
			}

		}

		if ( ! overlap) {
			(*p)->setPosition(pos);
			++p;
		}

	}

}

/*
 * Place each particle following a cube pattern.
 *
 * @param {vect_t} spaceSize: the simulation space (x, y and z)
 * @param {std::list<Particle *> *} particles
 */
void cubeDistribution(vect_t spaceSize, std::list<Particle *> *particles) {

	int m;
	int i, j, k;

	std::list<Particle *>::iterator p;

	m = (int)round(pow(particles->size(), 1/3.0));

	i = 0;
	j = 0;
	k = 0;

	for (p = particles->begin(); p != particles->end(); ++p) {

		(*p)->setX((0.5 + i)*spaceSize.x/m);
		(*p)->setY((0.5 + j)*spaceSize.y/m);
		(*p)->setZ((0.5 + k)*spaceSize.z/m);

		k++;

		if (k%m == 0) {

			k = 0; j++;

			if (j%m == 0) {

				j = 0; i++;

				if (i%m == 0) {
					i = 0;
				}

			}

		}

	}

}

/* 
 * Place each particle following a sphere surface.
 *
 * @param {vect_t} spaceSize: the simulation space (x, y and z)
 * @param {std::list<Particle *> *} particles
 * @param {point_t} c: the center of the sphere
 * @param {double} radius: the sphere radius containing the center of the particles
 */
void sphereDistribution(vect_t spaceSize, std::list<Particle *> *particles, point_t c, double radius) {

	std::list<Particle *>::iterator p;

	for (p = particles->begin(); p != particles->end(); ++p) {

	}

}

/*
 * Create two distant groups of particles. The particles in each group are
 * placed randomly.
 *
 * @param {vect_t} spaceSize: the simulation space (x, y and z)
 * @param {std::list<Particle *> *} particles
 * @param {point_t} ca: the center of the first group
 * @param {point_t} cb: the center of the second group
 */
void twoSidedDistribution(vect_t spaceSize, std::list<Particle *> *particles, point_t ca, point_t cb) {

	std::list<Particle *>::iterator p;

	for (p = particles->begin(); p != particles->end(); ++p) {

	}

}

/*
 * Put all the particles one near the other at a certain point.
 *
 * @param {vect_t} spaceSize: the simulation space (x, y and z)
 * @param {std::list<Particle *> *} particles
 * @param {point_t} c: center
 */
void highDensityDistribution(vect_t spaceSize, std::list<Particle *> *particles, point_t c) {

	std::list<Particle *>::iterator p;

	for (p = particles->begin(); p != particles->end(); ++p) {

	}

}

/*
 * Detect whether two sphere particles are overlaping or not.
 *
 * @param {point_t} ca: the center coordinates of the first particle
 * @param {double} ra: the radius of the first particle
 * @param {point_t} cb: the center coordinates of the second particle
 * @param {double} rb: the radius of the second particle
 * @return {bool}
 */
bool checkOverlap(point_t ca, double ra, point_t cb, double rb) {

	double dx = (ca.x-cb.x);
	double dy = (ca.y-cb.y);
	double dz = (ca.z-cb.z);

	return  sqrt(dx*dx + dy*dy + dz*dz) < ra + rb ? true : false;

}
