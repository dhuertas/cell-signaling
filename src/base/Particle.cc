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