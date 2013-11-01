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

	while (p != particles->end()) {

		overlap = false;
		
		range.x = spaceSize.x - 2*(*p)->getRadius();
		range.y = spaceSize.y - 2*(*p)->getRadius();
		range.z = spaceSize.z - 2*(*p)->getRadius();

		// We are using the Omnet internal random number generator dblrand()
		pos.x = (*p)->getRadius() + range.x*dblrand();
		pos.y = (*p)->getRadius() + range.y*dblrand();
		pos.z = (*p)->getRadius() + range.z*dblrand();

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
 * Place each particle following a sphere surface. The particles are uniformly 
 * distributed over a sphere surface. More info here:
 *
 * - http://mathworld.wolfram.com/SpherePointPicking.html
 *
 * @param {vect_t} spaceSize: the simulation space (x, y and z)
 * @param {std::list<Particle *> *} particles
 * @param {point_t} c: the center of the sphere
 * @param {double} radius: the sphere radius containing the center of the particles
 */
void sphereDistribution(vect_t spaceSize, std::list<Particle *> *particles, point_t c, double radius) {

	bool overlap;

	double u, v;
	double theta, phi;

	double minx, miny, minz;
	double maxRadius;

	point_t pos;

	std::list<Particle *>::iterator p, q;

	u = 0;
	v = 0;

	maxRadius = 0;

	if (c.x == 0) c.x = spaceSize.x/2;
	if (c.y == 0) c.y = spaceSize.y/2;
	if (c.z == 0) c.z = spaceSize.z/2;

	if (radius == 0) {
		// Use the maximum radius possible given the space size
		minx = std::min(spaceSize.x - c.x, c.x - 0);
		miny = std::min(spaceSize.y - c.y, c.y - 0);
		minz = std::min(spaceSize.z - c.z, c.z - 0);

		radius = std::min(std::min(minx, miny), minz);

		// Find the max radius from the particles and subtract half of it (and
		// a little more).
		for (p = particles->begin(); p != particles->end(); ++p) {
			maxRadius = std::max((*p)->getRadius(), maxRadius);
		}

		radius -= maxRadius*(1.0 + 0.0001);

	} else if (radius < 0) {
		radius = -radius;
	}

	p = particles->begin();

	// WARNING: if the space size is not big enough this function may hang the
	// simulation
	while (p != particles->end()) {

		overlap = false;

		u = 0;
		v = 0;

		// dblrand() returns values in [0,1) while u and v must be in (0, 1)
		while (u == 0) u = dblrand();
		while (v == 0) v = dblrand();

		theta = 2*M_PI*u;
		phi = acos(2*v - 1);

		pos.x = c.x + radius*sin(theta)*cos(phi);
		pos.y = c.y + radius*sin(theta)*sin(phi);
		pos.z = c.z + radius*cos(theta);

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
 * Put all the particles one near the other at a certain point.
 *
 * @param {vect_t} spaceSize: the simulation space (x, y and z)
 * @param {std::list<Particle *> *} particles
 * @param {point_t} c: center
 */
void highDensityDistribution(vect_t spaceSize, std::list<Particle *> *particles, point_t c) {

	bool overlap;

	uint8_t overlapCount;
	double movingVariance;

	point_t pos;

	std::list<Particle *>::iterator p, q;

	overlapCount = 100;

	p = particles->begin();

	while (p != particles->end()) {

		overlap = false;

		pos.x = dblRandNormal(c.x, movingVariance);
		pos.y = dblRandNormal(c.y, movingVariance);
		pos.z = dblRandNormal(c.z, movingVariance);

		if (pos.x - (*p)->getRadius() < 0 || pos.x + (*p)->getRadius() > spaceSize.x ||
			pos.y - (*p)->getRadius() < 0 || pos.y + (*p)->getRadius() > spaceSize.y ||
			pos.z - (*p)->getRadius() < 0 || pos.z + (*p)->getRadius() > spaceSize.z) {
			continue;	
		}

		for (q = particles->begin(); q != p; ++q) {

			if (checkOverlap(pos, (*p)->getRadius(), (*q)->getPosition(), (*q)->getRadius())) {
				overlap = true;
				break;
			}

		}

		if ( ! overlap) {
			(*p)->setPosition(pos);
			++p;
		} else {
			overlapCount--;
			if (overlapCount == 0) {
				overlapCount = 100;
			}
			movingVariance++;
		}

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

/*
 * Returns a random number following a Normal distribution N(mean, var)
 * 
 * @param {double} mean
 * @param {double} vari
 */
double dblRandNormal(double mean, double var) {

	double x1, x2, w, y1;
	// double y2;

	do {
		x1 = 2.0*dblrand() - 1.0;
		x2 = 2.0*dblrand() - 1.0;
		w = x1*x1 + x2*x2;
	} while (w >= 1.0);

	w = sqrt((-2.0*log(w))/w);
	y1 = x1*w;
	//y2 = x2*w;

	return mean + y1*sqrt(var);
}
