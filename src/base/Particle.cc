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

	spaceCellIdx.im = 0x00;
	spaceCellIdx.i = 0;
	spaceCellIdx.j = 0;
	spaceCellIdx.k = 0;

	prevSpaceCellIdx.im = 0x00;
	prevSpaceCellIdx.i = 0;
	prevSpaceCellIdx.j = 0;
	prevSpaceCellIdx.k = 0;

	mass = m;
	lastCollisionTime = 0;

	listRadius = 1;

}