//  Omnet++ project to simulate cell signaling communications 
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

#include "SphereMobility.h"

using namespace std;

/*
 * Computes the collision time with the neighboring particles and stores the
 * smallest one.
 *
 * @param {CollisionMessage *} msg
 * @param {integer} kind message kind
 * @param {Sphere *} s
 * @return {double} the smallest computed collision time
 */
double SphereMobility::nextCollision(CollisionMessage *msg, int kind, Sphere *s) {

	int a, b, c;			// Nested "for" loops
	int i, j, k, n;			// Indexes to access the current space cell
	int *Nx, *Ny, *Nz;		// Number of space cells (or divisions) in each axis
	int collisionCounter;

	double collisionTime, sTime, temp;
	double prevCollisionTime;
	double partnerCollisionTime;

	vector<int> spaceCells;
	vector<int>::iterator sc;

	vector<Particle *> particles;
	vector<Particle *>::iterator p;

	list<Particle *> *particleList;
	list<Particle *>::const_iterator pl;

	Manager *manager;
	Particle *partner, *prevPartner; // The collision partner

	CollisionMessage *partnerCollisionMsg;

	manager = msg->getManager();
	sTime = simTime().dbl();

	Nx = manager->getNumberOfSpaceCellsX();
	Ny = manager->getNumberOfSpaceCellsY();
	Nz = manager->getNumberOfSpaceCellsZ();

	n = s->getSpaceCell();

	collisionCounter = 0;

	// i, j and k are the indexes of the space cell for each axis
	i =  n/((*Nz)*(*Ny));
	j = (n%((*Nz)*(*Ny)))/(*Nz);
	k = (n%((*Nz)*(*Ny)))%(*Nz);

	collisionTime = NO_TIME;
	prevCollisionTime = NO_TIME;
	partnerCollisionTime = NO_TIME;

	partner = NULL;
	prevPartner = NULL;

	if (manager->getMode() == M_NNLIST) {

		particles = s->getNeighborParticles();

	} else {
		// M_CELLLIST / default
		// Get the list of space cell indexes to ask the manager for the 
		// particles

		// It should be faster to directly compute this here rather than
		// calling a manager method, so we are not passing vectors around
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
			particles.insert(particles.end(), particleList->begin(), particleList->end());
		}

	}

	// Initialize the collision time. If the collision event is still scheduled
	// it should return the current scheduled collision time. Otherwise it should
	// return NO_TIME.
	prevCollisionTime = msg->getCollisionTime();

	if (prevCollisionTime != NO_TIME && prevCollisionTime > sTime) {

		partner = msg->getPartner();
		prevPartner = msg->getPartner();

		collisionTime = prevCollisionTime;
		collisionCounter++;

	}

	// Loop through the retrieved particles
	for (p = particles.begin(); p != particles.end(); ++p) {

		if (*p == s) continue;

		// Solve particle to particle collision
		temp = solveCollision(s, *p);

		// Only keep the collision time if it's smaller than the rest of computed 
		// collision times so far and is also smaller than the partner collision time.
		partnerCollisionMsg = ((Sphere *)(*p))->getCollisionMessage();

		if (partnerCollisionMsg != NULL) {
		    partnerCollisionTime = partnerCollisionMsg->getCollisionTime();
		} else {
		    partnerCollisionTime = NO_TIME;
		}

		if (temp != NO_TIME && sTime < temp) {
			// Collision found!
			if (collisionCounter == 0) {

				if (temp < partnerCollisionTime || partnerCollisionTime == NO_TIME) {

					collisionTime = temp;
					partner = *p;
					collisionCounter++;

				}

			} else {

				if (temp < collisionTime) {

					if (temp < partnerCollisionTime || partnerCollisionTime == NO_TIME) {

						collisionTime = temp;
						partner = *p;
						collisionCounter++;

					}

				}

			}

		}

	}

	msg->setPartner(partner);
	msg->setPrevPartner(prevPartner);

	return collisionTime;

}

/*
 *
 */
double SphereMobility::nextBoundaryCollision(CollisionMessage *msg, Sphere *s) {

	int boundariesMode = s->getBoundariesMode();

	if (boundariesMode == BM_ELASTIC) {

		return SphereMobility::nextWallCollision(msg, s);

	} else if (boundariesMode == BM_EXPIRE) {

		return leaveBoundedSpace(msg, s);

	} else if (boundariesMode == BM_PERIODIC) {

	}

	return NO_TIME;

}

/*
 * Computes the collision time with the simulation space walls and stores the
 * smallest one.
 *
 * @param {CollisionMessage *} event
 * @param {Sphere *} s
 * @return {double} the smallest computed time
 */
double SphereMobility::nextWallCollision(CollisionMessage *msg, Sphere *s) {

	int collisionCounter;

	double rad;
	double boundaryCollisionTime, lastCollisionTime, temp;

	point_t P, Q, R;
	vect_t V, W, N;
	vect_t *S;

	point_t *pos = NULL;
	vect_t *vel = NULL;

	vector<int> sides, hits;
	vector<int>::const_iterator side, hit;

	Manager *manager;

	manager = s->getManager();

	boundaryCollisionTime = NO_TIME;
	lastCollisionTime = s->getLastCollisionTime();

	S = manager->getSpaceSize();

	rad = s->getRadius();

	pos = s->getPosition();
	vel = s->getVelocity();

	if (vel->x == 0 && vel->y == 0 && vel->z == 0) return NO_TIME;

	// In the simulation space we have 6 possible sides but since we know the 
	// direction of the particle we need to check only 3 (at most).
	if (vel->x > 0) sides.push_back(1);
	else if (vel->x < 0) sides.push_back(4);

	if (vel->y > 0) sides.push_back(3);
	else if (vel->y < 0) sides.push_back(2);

	if (vel->z > 0) sides.push_back(0);
	else if (vel->z < 0) sides.push_back(5);

	collisionCounter = 0;

	for (side = sides.begin(); side != sides.end(); ++side) {
		// Select 3 points for each side following the previous rules
		P.x = sp[0*6 + *side]*S->x;
		P.y = sp[1*6 + *side]*S->y;
		P.z = sp[2*6 + *side]*S->z;

		Q.x = sp[3*6 + *side]*S->x;
		Q.y = sp[4*6 + *side]*S->y;
		Q.z = sp[5*6 + *side]*S->z;

		R.x = sp[6*6 + *side]*S->x;
		R.y = sp[7*6 + *side]*S->y;
		R.z = sp[8*6 + *side]*S->z;

		V.x = Q.x - P.x; V.y = Q.y - P.y; V.z = Q.z - P.z;
		W.x = R.x - P.x; W.y = R.y - P.y; W.z = R.z - P.z;

		// Cross product to find the Normal vector
		N.x = V.y * W.z - V.z * W.y;
		N.y = V.z * W.x - V.x * W.z;
		N.z = V.x * W.y - V.y * W.x;

		// This is the plane equation N.x*(x-P.x)+N.y*(y-P.y)+N.z*(z-P.z) = 0
		// Substitute for:
		//     x = xi + vx*t + R*vnx, vnx = vx/sqrt(vx2+vy2+vz2), R = radius
		//     y = yi + vy*t + R*vny, vny = vy/sqrt(vx2+vy2+vz2), "
		//     z = zi + vz*t + R*vnz, vnz = vz/sqrt(vx2+vy2+vz2), "
		// and find for t.
		temp = (N.x*vel->x + N.y*vel->y + N.z*vel->z);

		if (temp != 0) {

			if (*side == 0)
				temp = (N.x*(P.x - pos->x) + N.y*(P.y - pos->y) + N.z*(P.z - (pos->z + rad)))/temp;
			else if (*side == 1)
				temp = (N.x*(P.x - (pos->x + rad)) + N.y*(P.y - pos->y) + N.z*(P.z - pos->z))/temp;
			else if (*side == 2)
				temp = (N.x*(P.x - pos->x) + N.y*(P.y - (pos->y - rad)) + N.z*(P.z - pos->z))/temp;
			else if (*side == 3)
				temp = (N.x*(P.x - pos->x) + N.y*(P.y - (pos->y + rad)) + N.z*(P.z - pos->z))/temp;
			else if (*side == 4)
				temp = (N.x*(P.x - (pos->x - rad)) + N.y*(P.y - pos->y) + N.z*(P.z - pos->z))/temp;
			else
				temp = (N.x*(P.x - pos->x) + N.y*(P.y - pos->y) + N.z*(P.z - (pos->z - rad)))/temp;

			// Solution found. "temp" is the amount of time the particle takes 
			// to go from its last event point to the simulation space side.
			if (collisionCounter == 0) {

				boundaryCollisionTime = temp;
				hits.push_back(*side);
				collisionCounter++;

			} else {

				if (0 < temp && temp < boundaryCollisionTime) {

					boundaryCollisionTime = temp;

					hits.clear();
					hits.push_back(*side);
					collisionCounter++;

				} else if (0 < temp && temp == boundaryCollisionTime) {
					// The particle hit two or more sides at the same time
					boundaryCollisionTime = temp;

					hits.push_back(*side);
					collisionCounter++;

				} else {
					// It's greater
				}

			}

		} else {
			// Line and plane do not intersect
		}

	}

	// The future particle position
	msg->setX(pos->x + vel->x*boundaryCollisionTime);
	msg->setY(pos->y + vel->y*boundaryCollisionTime);
	msg->setZ(pos->z + vel->z*boundaryCollisionTime);

	// Compute the future time value
	boundaryCollisionTime += lastCollisionTime;

	// Compute the future velocity vector
	for (hit = hits.begin(); hit != hits.end(); ++hit) {
		msg->setVx((*hit == 1 || *hit == 4) ? -vel->x : vel->x);
		msg->setVy((*hit == 2 || *hit == 3) ? -vel->y : vel->y);
		msg->setVz((*hit == 0 || *hit == 5) ? -vel->z : vel->z);
	}

	return boundaryCollisionTime;

}

/*
 * Still using the time when the sphere surface collides with the wall.
 * @param {CollisionMessage *} msg
 * @param {Sphere *} s
 * @return {double} the computed time
 */
double SphereMobility::leaveBoundedSpace(CollisionMessage *msg, Sphere *s) {

	int collisionCounter;

	double rad;

	double boundaryCollisionTime, lastCollisionTime, temp;

	point_t P, Q, R;
	vect_t V, W, N;
	vect_t *S;

	point_t *pos = NULL;
	vect_t *vel = NULL;

	vector<int> sides, hits;
	vector<int>::const_iterator side, hit;

	Manager *manager;

	manager = s->getManager();

	boundaryCollisionTime = NO_TIME;
	lastCollisionTime = s->getLastCollisionTime();

	S = manager->getSpaceSize();

	rad = s->getRadius();

	pos = s->getPosition();
	vel = s->getVelocity();

	if (vel->x == 0 && vel->y == 0 && vel->z == 0) return NO_TIME;

	// In the simulation space we have 6 possible sides but since we know the 
	// direction of the particle we need to check only 3 (at most).
	if (vel->x > 0) sides.push_back(1);
	else if (vel->x < 0) sides.push_back(4);

	if (vel->y > 0) sides.push_back(3);
	else if (vel->y < 0) sides.push_back(2);

	if (vel->z > 0) sides.push_back(0);
	else if (vel->z < 0) sides.push_back(5);

	collisionCounter = 0;

	for (side = sides.begin(); side != sides.end(); ++side) {
		// Select 3 points for each side following the previous rules
		P.x = sp[0*6 + *side]*S->x;
		P.y = sp[1*6 + *side]*S->y;
		P.z = sp[2*6 + *side]*S->z;

		Q.x = sp[3*6 + *side]*S->x;
		Q.y = sp[4*6 + *side]*S->y;
		Q.z = sp[5*6 + *side]*S->z;

		R.x = sp[6*6 + *side]*S->x;
		R.y = sp[7*6 + *side]*S->y;
		R.z = sp[8*6 + *side]*S->z;

		V.x = Q.x - P.x; V.y = Q.y - P.y; V.z = Q.z - P.z;
		W.x = R.x - P.x; W.y = R.y - P.y; W.z = R.z - P.z;

		// Cross product to find the Normal vector
		N.x = V.y * W.z - V.z * W.y;
		N.y = V.z * W.x - V.x * W.z;
		N.z = V.x * W.y - V.y * W.x;

		// This is the plane equation N.x*(x-P.x)+N.y*(y-P.y)+N.z*(z-P.z) = 0
		// Substitute for:
		//     x = xi + vx*t + R*vnx, vnx = vx/sqrt(vx2+vy2+vz2), R = radius
		//     y = yi + vy*t + R*vny, vny = vy/sqrt(vx2+vy2+vz2), "
		//     z = zi + vz*t + R*vnz, vnz = vz/sqrt(vx2+vy2+vz2), "
		// and find for t.
		temp = (N.x*vel->x + N.y*vel->y + N.z*vel->z);

		if (temp != 0) {

			if (*side == 0)
				temp = (N.x*(P.x - pos->x) + N.y*(P.y - pos->y) + N.z*(P.z - (pos->z + rad)))/temp;
			else if (*side == 1)
				temp = (N.x*(P.x - (pos->x + rad)) + N.y*(P.y - pos->y) + N.z*(P.z - pos->z))/temp;
			else if (*side == 2)
				temp = (N.x*(P.x - pos->x) + N.y*(P.y - (pos->y - rad)) + N.z*(P.z - pos->z))/temp;
			else if (*side == 3)
				temp = (N.x*(P.x - pos->x) + N.y*(P.y - (pos->y + rad)) + N.z*(P.z - pos->z))/temp;
			else if (*side == 4)
				temp = (N.x*(P.x - (pos->x - rad)) + N.y*(P.y - pos->y) + N.z*(P.z - pos->z))/temp;
			else
				temp = (N.x*(P.x - pos->x) + N.y*(P.y - pos->y) + N.z*(P.z - (pos->z - rad)))/temp;

			// Solution found. "temp" is the amount of time the particle takes 
			// to go from its last event point to the simulation space side.
			if (collisionCounter == 0) {

				boundaryCollisionTime = temp;
				hits.push_back(*side);
				collisionCounter++;

			} else {

				if (0 < temp && temp < boundaryCollisionTime) {

					boundaryCollisionTime = temp;

					hits.clear();
					hits.push_back(*side);
					collisionCounter++;

				} else if (0 < temp && temp == boundaryCollisionTime) {
					// The particle hit two or more sides at the same time
					boundaryCollisionTime = temp;

					hits.push_back(*side);
					collisionCounter++;

				} else {
					// It's greater
				}

			}

		} else {
			// Line and plane do not intersect
		}

	}

	// The future particle position
	msg->setX(pos->x + vel->x*boundaryCollisionTime);
	msg->setY(pos->y + vel->y*boundaryCollisionTime);
	msg->setZ(pos->z + vel->z*boundaryCollisionTime);

	// Compute the future time value
	boundaryCollisionTime += lastCollisionTime;

	// Compute the future velocity vector
	for (hit = hits.begin(); hit != hits.end(); ++hit) {
		msg->setVx((*hit == 1 || *hit == 4) ? -vel->x : vel->x);
		msg->setVy((*hit == 2 || *hit == 3) ? -vel->y : vel->y);
		msg->setVz((*hit == 0 || *hit == 5) ? -vel->z : vel->z);
	}

	return boundaryCollisionTime;

}

/*
 * Solves the Sphere-Sphere collision problem and returns the collision time
 * 
 * @param {Particle *} pa
 * @param {Particle *} pb
 * @return {double} the collision time
 */
double SphereMobility::solveCollision(Particle * pa, Particle * pb) {
	// Distance between centers A and B when t = tc (time of collision):
	//                 ______________________________________________
	//           \    / ( Ax + Avx*(tc-ta) - (Bx + Bvx*(tc-tb) )² +
	//  D(A, B) = \  /  ( Ay + Avy*(tc-ta) - (By + Bvy*(tc-tb) )² +   = Ra + Rb
	//             \/   ( Az + Avz*(tc-ta) - (Bz + Bvz*(tc-tb) )²
	//
	// ta: time when the previous collision took place for particle A
	// tb: same for particle B

	// (dxi + dvx*tc)² + (dyi + dvy*tc)² + (dyi + dvy*tc)²= (A.r + B.r)²
	double ta = pa->getLastCollisionTime();
	double tb = pb->getLastCollisionTime();

	point_t *posa = NULL;
	point_t *posb = NULL;

	vect_t *vela = NULL;
	vect_t *velb = NULL;

	posa = pa->getPosition();
	posb = pb->getPosition();

	vela = pa->getVelocity();
	velb = pb->getVelocity();

	double dxi = posa->x - posb->x - vela->x*ta + velb->x*tb;
	double dyi = posa->y - posb->y - vela->y*ta + velb->y*tb;
	double dzi = posa->z - posb->z - vela->z*ta + velb->z*tb;

	double dvx = vela->x - velb->x;
	double dvy = vela->y - velb->y;
	double dvz = vela->z - velb->z;

	double radd = pb->getRadius() + pa->getRadius();

	// a*t² + b*t + c = 0
	double a = dvx*dvx + dvy*dvy + dvz*dvz;
	double b = 2*(dxi*dvx + dyi*dvy + dzi*dvz);
	double c = dxi*dxi + dyi*dyi + dzi*dzi - radd*radd;

	if (b*b >= 4*a*c) {
	// Return the smaller solution
		return (-b - sqrt(b*b - 4*a*c))/(2*a);
	} else {
		return NO_TIME;
	}

}

/*
 * Computes the next position of the sphere after a delta time following
 * the Boltzmann-Maxwell distribution, with parameters:
 * - D: diffusion
 * - dt: delta time
 * Note: delta time is a fixed value for now
 * @param {Particle *}
 * @return {double} next delta simulation time
 */
double SphereMobility::brownianMotion(BrownianMotionMessage *msg, Particle *p) {

	double nextDeltaTime = NO_TIME;
	double *BMStdDev = NULL;
	double dt;

	point_t *pos = NULL;
	point_t nextPos;

	dt = msg->getManager()->getDeltaTime();

	if (dt > 0) {

		nextDeltaTime = simTime().dbl() + dt;

		pos = p->getPosition();
		BMStdDev = p->getBMStdDev();

		nextPos.x = normal(pos->x, *BMStdDev);
		nextPos.y = normal(pos->y, *BMStdDev);
		nextPos.z = normal(pos->z, *BMStdDev);

		msg->setVx((nextPos.x - pos->x)/dt);
		msg->setVy((nextPos.y - pos->y)/dt);
		msg->setVz((nextPos.z - pos->z)/dt);

		msg->setBrownianMotionTime(nextDeltaTime);

	}

	return nextDeltaTime;
}
