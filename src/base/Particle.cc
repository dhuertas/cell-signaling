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

}

/*
 * Create and populate the Verlet list
 */
void Particle::createVerletList(std::list<Particle*> *l) {

    bool found;

    double dx, dy, dz;

    std::list<Particle *>::const_iterator pa, pn;

    for (pa = l->begin(); pa != l->end(); ++pa) {

        if ((*pa) == this) continue;

// Check first if the particle has already been added to our list
        found = false;

        for (pn = neighbourParticles.begin();
            pn != neighbourParticles.end(); ++pn) {

            if ((*pa) == (*pn)) {
// Particle pn is already in our list
                found = true;
                break;

            }
        }

        if ( ! found) {
// Compute distance and, if it is closer than or equal to listRadius, add the
// particle to our list and add ourselves to its list
            dx = this->getX() - (*pa)->getX();
            dy = this->getY() - (*pa)->getY();
            dz = this->getZ() - (*pa)->getZ();

            if (listRadius*listRadius <= dx*dx + dy*dy + dz*dz) {
                l->push_back((*pa));
                (*pa)->addParticleToVerletList(this);
            }
        }
    }

// Eventually We will have to update our list.
// TODO to be continued ...
}

/*
 * Add a particle to our Verlet list
 *
 * @param {Particle *} p
 */
void Particle::addParticleToVerletList(Particle *p) {

    neighbourParticles.push_back(p);

}
