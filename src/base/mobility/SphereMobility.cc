#include "SphereMobility.h"

using namespace std;

/*
 * Computes the collision time with the neighboring particles and returns the
 * smallest one.
 */
void SphereMobility::collisionTime(MobilityMessage *collisionMsg, Sphere *s) {

	int a, b, c;			// Nested "for" loops

	int collisionCounter;

	int i, j, k, N, n;		// Indexes to access the current space cell
	int Nx, Ny, Nz;			// Number of space cells (or divisions) in each axis

	double collisionTime, sTime, temp;

	list<Particle *> particles;
	list<Particle *>::const_iterator p;

	Manager * manager;
	Particle * partner; // The collision partner

	manager = collisionMsg->getManager();

	collisionTime = NO_TIME;
	sTime = simTime().dbl();

	Nx = manager->getNumberOfSpaceCellsX();
	Ny = manager->getNumberOfSpaceCellsY();
	Nz = manager->getNumberOfSpaceCellsZ();

	n = s->getSpaceCell();

// i, j and k are the indexes of the space cell for each axis
	i = n / (Nz*Ny);
	j = (n % (Nz*Ny)) / Nz;
	k = (n % (Nz*Ny)) % Nz;

	collisionCounter = 0;

	switch (manager->getMode()) {

		case M_NNLIST:

			// TODO better return a pointer to the list instead of replicating it
			particles = s->getNeighborParticles();

			for (p = particles.begin(); p != particles.end(); ++p) {

				if (*p == s) continue;

// Solve particle to particle collision
				temp = solveCollisionTime(s, *p);

				if (collisionCounter == 0) {

					if (sTime <= temp) {
						collisionTime = temp;
						partner = (*p);
						collisionCounter++;
					}

				} else {

					if (sTime <= temp  && temp < collisionTime) {
						collisionTime = temp;
						partner = (*p);
						collisionCounter++;
					}

				}
			}

			break;

		case M_CELLLIST:
		default:

			for (a = -1; a <= 1; a++)
			for (b = -1; b <= 1; b++)
			for (c = -1; c <= 1; c++) {

// The neighbor cell must be contained in the simulation space
				if (CELLBELONGSTOSIMSPACE(i+a, j+b, k+c, Nx, Ny, Nz)) {

					N = (i+a)*Ny*Nz + (j+b)*Nz + (k+c);
					// TODO better return a pointer to the list instead of replicating it
					particles = manager->getSpaceCellParticles(N);

					for (p = particles.begin(); p != particles.end(); ++p) {

						if (*p == s) {
							continue;
						}

// Solve particle to particle collision
						temp = solveCollisionTime(s, *p);

						if (collisionCounter == 0) {

							if (sTime <= temp) {
								collisionTime = temp;
								partner = (*p);
								collisionCounter++;
							}

						} else {

							if (sTime <= temp  && temp < collisionTime) {
								collisionTime = temp;
								partner = (*p);
								collisionCounter++;
							}

						}

					}
				}
			}

			break;
	}

	collisionMsg->setKind(EV_COLLISION);
	collisionMsg->setEventTime(collisionTime);
	collisionMsg->setPartner(partner);

}

void SphereMobility::wallCollisionTime(MobilityMessage *wallCollisionMsg, Sphere *s) {
	
	int collisionCounter;

	double Sx, Sy, Sz;
	double rad, x, y, z, vx, vy, vz; // Particle radius, position and velocity
	double wallCollisionTime, lastCollisionTime, temp;

	point_t P, Q, R;
	vect_t V, W, N;

	vector<int> sides, hits;
	vector<int>::const_iterator side, hit;

	Manager * manager;

	manager = wallCollisionMsg->getManager();

	wallCollisionTime = NO_TIME;
	lastCollisionTime = s->getLastCollisionTime();

	Sx = manager->getSpaceSizeX();
	Sy = manager->getSpaceSizeY();
	Sz = manager->getSpaceSizeZ();

	rad = s->getRadius();

	x = s->getX();
	y = s->getY();
	z = s->getZ();

	vx = s->getVx();
	vy = s->getVy();
	vz = s->getVz();

// In the simulation space we have 6 possible sides but since we know the 
// direction of the particle we need to check only 3 (at most).
	if (vx > 0) sides.push_back(1);
	else if (vx < 0) sides.push_back(4);

	if (vy > 0) sides.push_back(3);
	else if (vy < 0) sides.push_back(2);

	if (vz > 0) sides.push_back(0);
	else if (vz < 0) sides.push_back(5);

	collisionCounter = 0;

	for (side = sides.begin(); side != sides.end(); ++side) {
// Select 3 points for each side following the previous rules
		P.x = sp[0*6 + *side]*Sx;
		P.y = sp[1*6 + *side]*Sy;
		P.z = sp[2*6 + *side]*Sz;

		Q.x = sp[3*6 + *side]*Sx;
		Q.y = sp[4*6 + *side]*Sy;
		Q.z = sp[5*6 + *side]*Sz;

		R.x = sp[6*6 + *side]*Sx;
		R.y = sp[7*6 + *side]*Sy;
		R.z = sp[8*6 + *side]*Sz;

		V.x = Q.x - P.x; V.y = Q.y - P.y; V.z = Q.z - P.z;
		W.x = R.x - P.x; W.y = R.y - P.y; W.z = R.z - P.z;

// Cross product to find the Normal vector
		N.x = V.y * W.z - V.z * W.y;
		N.y = V.z * W.x - V.x * W.z;
		N.z = V.x * W.y - V.y * W.x;

// This is our plane equation N.x*(x-P.x)+N.y*(y-P.y)+N.z*(z-P.z) = 0
// Substitute for:
//     x = xi + vx*t + R*vnx, vnx = vx/sqrt(vx2+vy2+vz2), R = radius
//     y = yi + vy*t + R*vny, vny = vy/sqrt(vx2+vy2+vz2), "
//     z = zi + vz*t + R*vnz, vnz = vz/sqrt(vx2+vy2+vz2), "
// and find for t.
		temp = (N.x*vx + N.y*vy + N.z*vz);

		if (temp != 0) {

			if (*side == 0)
				temp = (N.x*(P.x - x) + N.y*(P.y - y) + N.z*(P.z - (z + rad)))/temp;
			else if (*side == 1)
				temp = (N.x*(P.x - (x + rad)) + N.y*(P.y - y) + N.z*(P.z - z))/temp;
			else if (*side == 2)
				temp = (N.x*(P.x - x) + N.y*(P.y - (y - rad)) + N.z*(P.z - z))/temp;
			else if (*side == 3)
				temp = (N.x*(P.x - x) + N.y*(P.y - (y + rad)) + N.z*(P.z - z))/temp;
			else if (*side == 4)
				temp = (N.x*(P.x - (x - rad)) + N.y*(P.y - y) + N.z*(P.z - z))/temp;
			else
				temp = (N.x*(P.x - x) + N.y*(P.y - y) + N.z*(P.z - (z - rad)))/temp;

// Solution found. "temp" is the amount of time the particle takes to go from its last 
// event point to the simulation space side.
			if (collisionCounter == 0) {

				wallCollisionTime = temp;
				hits.push_back(*side);
				collisionCounter++;

			} else {

				if (0 < temp && temp < wallCollisionTime) {

					wallCollisionTime = temp;

					hits.clear();
					hits.push_back(*side);
					collisionCounter = 1;

				} else if (0 < temp && temp == wallCollisionTime) {
// The particle hit two or more sides at the same time
					hits.push_back(*side);
					collisionCounter++;

				}

			}

		} else {
// Line and plane doesn't intersect
		}

	}

	wallCollisionMsg->setKind(EV_WALLCOLLISION);

// Compute the future values
	wallCollisionMsg->setEventTime(wallCollisionTime + lastCollisionTime);

// The future particle position
	wallCollisionMsg->setX(x + vx*wallCollisionTime);
	wallCollisionMsg->setY(y + vy*wallCollisionTime);
	wallCollisionMsg->setZ(z + vz*wallCollisionTime);

// The future velocity vector
	for (hit = hits.begin(); hit != hits.end(); ++hit) {

		wallCollisionMsg->setVx((*hit == 1 || *hit == 4) ? -vx : vx);
		wallCollisionMsg->setVy((*hit == 2 || *hit == 3) ? -vy : vy);
		wallCollisionMsg->setVz((*hit == 0 || *hit == 5) ? -vz : vz);

	}

}

double SphereMobility::solveCollisionTime(Particle * pa, Particle * pb) {
// Distance between centers A and B when t = tc (time of collision):
//                 ______________________________________________
//           \    / ( Ax + Avx*(tc-ta) - (Bx + Bvx*(tc-tb) )² + |
//  D(A, B) = \  /  ( Ay + Avy*(tc-ta) - (By + Bvy*(tc-tb) )² +   = Ra + Rb
//             \/   ( Az + Avz*(tc-ta) - (Bz + Bvz*(tc-tb) )²
//
// ta: when the previous collision take place for particle A
// tb: same for particle B

// (dxi + dvx*tc)² + (dyi + dvy*tc)² + (dyi + dvy*tc)²= (A.r + B.r)²

	double ta = pa->getLastCollisionTime();
	double tb = pb->getLastCollisionTime();

	double dxi = pa->getX() - pb->getX() - pa->getVx()*ta + pb->getVx()*tb;
	double dyi = pa->getY() - pb->getY() - pa->getVy()*ta + pb->getVy()*tb;
	double dzi = pa->getZ() - pb->getZ() - pa->getVz()*ta + pb->getVz()*tb;

	double dvx = pa->getVx() - pb->getVx();
	double dvy = pa->getVy() - pb->getVy();
	double dvz = pa->getVz() - pb->getVz();

	double radd = pb->getRadius() + pa->getRadius();

	// a*t² + b*t + c = 0
	double a = dvx*dvx + dvy*dvy + dvz*dvz;
	double b = 2*(dxi*dvx + dyi*dvy + dzi*dvz);
	double c = dxi*dxi + dyi*dyi + dzi*dzi - radd*radd;

	double result = -1;

	if (b*b >= 4*a*c) {
		
		double dtba = (-b + sqrt(b*b - 4*a*c))/(2*a);
		double dtbb = (-b - sqrt(b*b - 4*a*c))/(2*a);

		result = dtbb < dtba ? dtbb : dtba;

	}

	return result;

}