#include "Particle.h"

using namespace std;

/*
 * Constructor.
 */
Particle::Particle() {}

/*
 * Constructor overload.
 */
Particle::Particle(
	double x,
	double y,
	double vx,
	double vy,
	double m) {

	position.x = x;
	position.y = y;
	position.z = 0;

	velocity.x = vx;
	velocity.y = vy;
	velocity.z = 0;

	spaceCell = -1;
	mass = m;
	lastCollisionTime = 0;

	listRadius = 1;

}

/*
 * Create and populate the Near-Neighbor list
 */
void Particle::createNNList(std::list<Particle*> *l) {

	double dx, dy, dz;
	double lrs; // listRadiusSquared

	std::list<Particle *>::const_iterator pa;

	for (pa = l->begin(); pa != l->end(); ++pa) {

		if ((*pa) == this) continue;

// Two particles are said to be neighbor when the sum of their listRadius is greater
// than the distance between their centroids.
		dx = this->getX() - (*pa)->getX();
		dy = this->getY() - (*pa)->getY();
		dz = this->getZ() - (*pa)->getZ();

		lrs = listRadius+(*pa)->getListRadius();
		lrs *= lrs;

		if (lrs > dx*dx+dy*dy+dz*dz) {
			this->addParticleToNNList((*pa));
		}

	}

}

/*
 * Add a particle to our Verlet list
 *
 * @param {Particle *} p
 */
void Particle::addParticleToNNList(Particle *p) {

	neighbourParticles.push_back(p);

}

/*
 * Add a particle to our Verlet list
 *
 * @param {Particle *} p
 */
void Particle::updateNNList(std::list<Particle*> *l) {

// Empty previous list
	neighbourParticles.clear();

// Look for the particles closer than listRadius
	this->createNNList(l);
}