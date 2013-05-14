#include "Sphere.h"

using namespace std;

/*
 * Constructor.
 *
 * @param {double} x
 * @param {double} y
 * @param {double} z
 * @param {double} vx
 * @param {double} vy
 * @param {double} vz
 * @param {double} radius
 * @param {double} mass
 */
Sphere::Sphere(
	double x, 
	double y, 
	double z, 
	double vx, 
	double vy, 
	double vz, 
	double radius, 
	double mass) : Circle(x, y, vx, vy, radius, mass) {

	position.z = z;
	velocity.z = vz;

}

/*
 * Sets the manager attribute.
 * 
 * @param {string} param: the name of the manager module.
 */
void Sphere::setManager(string param) {

	cPar *managerName;

	try {

		managerName = & simulation.getSystemModule()->par(param.c_str());
		manager = (Manager *)simulation.getSystemModule()
			->getSubmodule(managerName->stringValue());

	} catch (cException *e) {
// TODO stop simulation
	}

}

/*
 * Draws the shape of the cell in the tk environment.
 */
void Sphere::tkEnvDrawShape() {

	std::stringstream buffer;

// We will use the shape drawing tool to draw a circle around the particle
// center instead of using
	buffer << 2*getRadius();
	getDisplayString().setTagArg("b", 0, buffer.str().c_str());
	getDisplayString().setTagArg("b", 1, buffer.str().c_str());

	getDisplayString().setTagArg("b", 2, "oval");
	getDisplayString().setTagArg("b", 3, "white");
	getDisplayString().setTagArg("b", 4, "black");
	getDisplayString().setTagArg("b", 5, 1);
}

/*
 * Update the module position in the tk environment. This method is used when
 * the particle is initialized.
 */
void Sphere::tkEnvUpdatePosition() {

	std::stringstream buffer;

// Set position string for tkenv
	buffer << getY();

	getDisplayString().setTagArg("p", 0, buffer.str().c_str());
	buffer.str(std::string()); // clear buffer

	buffer << getX();

	getDisplayString().setTagArg("p", 1, buffer.str().c_str());
	buffer.str(std::string());

}

/*
 * Update the module position in the tk environment. This method is used 
 * during the simulation.
 *
 * @param {Double} t
 */
void Sphere::tkEnvUpdatePosition(double t) {

	std::stringstream buffer;

// Set position string for tkenv
	buffer << getY() + getVy()*(t-getLastCollisionTime());

	getDisplayString().setTagArg("p", 0, buffer.str().c_str());
	buffer.str(std::string()); // clear buffer

	buffer << getX() + getVx()*(t-getLastCollisionTime());

	getDisplayString().setTagArg("p", 1, buffer.str().c_str());
	buffer.str(std::string());

}

/*
 * Initialize the event queue by computing the first event for the particle.
 * This function is called by the manager in order to initialize the event
 * queue.
 */
void Sphere::firstEventTime() {

// Methods called from other modules must have this macro
	Enter_Method_Silent();

	transferMsg = new MobilityMessage("transfer", EV_TRANSFER);
	collisionMsg = new MobilityMessage("collision", EV_COLLISION);
	wallCollisionMsg =  new MobilityMessage("collision", EV_WALLCOLLISION);

// Compute the first collision and the first transfer
	computeTransferTime();

	if (transferMsg->getEventTime() > 0) {
		scheduleAt(transferMsg->getEventTime(), transferMsg);
	}

// Schedule only the first event
	computeWallCollisionTime();
	computeCollisionTime();

	if (0 < collisionMsg->getEventTime() &&
		collisionMsg->getEventTime() < wallCollisionMsg->getEventTime()) {

		scheduleAt(collisionMsg->getEventTime(), collisionMsg);

	} else if (0 < wallCollisionMsg->getEventTime() &&
		wallCollisionMsg->getEventTime() < collisionMsg->getEventTime()) {

		scheduleAt(wallCollisionMsg->getEventTime(), wallCollisionMsg);

	} else {
// TODO Something that I didn't think of just happened
	}

}

/*
 * Generates a cMessage object for the next event. That event can be a transfer
 * event, a particle-particle collision event or a wall bounce event.
 */
void Sphere::nextEventTime() {
// Methods called from other modules must have this macro
    Enter_Method_Silent();

// Steps 3, 4 and 5.
    double transferTime,
        collisionTime,
        wallCollisionTime,
        partnerCollisionTime;

// Step 3. Get the next space cell transfer time
    if (transferMsg->isScheduled()) cancelEvent(transferMsg);
    if (collisionMsg->isScheduled()) cancelEvent(collisionMsg);
    if (wallCollisionMsg->isScheduled()) cancelEvent(wallCollisionMsg);

    computeTransferTime();
    transferTime = transferMsg->getEventTime();

    computeCollisionTime();
    collisionTime = collisionMsg->getEventTime();

    computeWallCollisionTime();
    wallCollisionTime = wallCollisionMsg->getEventTime();

    if (transferTime == NO_TIME) {
// Particle is not moving ... but we may have detected a collision
        if (collisionTime == NO_TIME) {
// Nope
        } else {
// Some other particle expects to collide with us
            collisionMsg->setKind(EV_CHECK);
            scheduleAt(collisionTime, collisionMsg);

        }

    } else {

        if (collisionTime == NO_TIME) {
// No collision time with other particle has been found
            if (wallCollisionTime <= transferTime) {

                scheduleAt(wallCollisionTime, wallCollisionMsg);

            } else {

                scheduleAt(transferTime, transferMsg);

            }

        } else {

            if (collisionTime <= transferTime &&
                collisionTime <= wallCollisionTime) {

                partnerCollisionTime = collisionMsg->getPartner()
                    ->scheduledCollisionTime();

                if (partnerCollisionTime == NO_TIME) {
// Partner particle does not have any scheduled collision event
                    scheduleAt(collisionTime, collisionMsg);
                    collisionMsg->getPartner()->nextEventTime();

                } else {

                    if (collisionTime < partnerCollisionTime) {

                        scheduleAt(collisionTime, collisionMsg);
                        collisionMsg->getPartner()->nextEventTime();

                    } else {

                        collisionMsg->setKind(EV_CHECK);
                        scheduleAt(collisionTime, collisionMsg);

                    }

                }


            } else if (wallCollisionTime < collisionTime &&
                wallCollisionTime <= transferTime) {

                scheduleAt(wallCollisionTime, wallCollisionMsg);

            } else {

                scheduleAt(transferTime, transferMsg);

            }

        }

    }

}

/*
 * Computes the transfer time of the particle with the sides of its space cell
 * and returns the smallest one.
 */
void Sphere::computeTransferTime() {

	int i, j, k, n;			// Indexes to access the current space cell
	int Nx, Ny, Nz;			// Number of space cells (or divisions) in each axis

	int counter;

	int sp[] = { // side points
	//  0, 1, 2, 3, 4, 5 <- dice side
		1, 1, 0, 0, 0, 0,
		0, 0, 0, 1, 0, 0,
		1, 0, 0, 0, 0, 0,
		1, 1, 1, 1, 0, 1,
		1, 1, 0, 1, 0, 0,
		1, 0, 0, 0, 1, 0,
		0, 1, 0, 0, 0, 0,
		0, 0, 0, 1, 1, 1,
		1, 1, 1, 1, 0, 0
	};

	double x, y, z, vx, vy, vz;	// Particle position and velocity

	double spaceCellSize,
		transferTime,
		sTime,
		collisionTime,
		temp;

	point_t P, Q, R;
	vect_t V, W, N;

	vector<int> sides, hits;
	vector<int>::const_iterator side, hit;

	transferTime = NO_TIME;

	sTime = simTime().dbl();
	collisionTime = getLastCollisionTime();

	spaceCellSize = manager->getSpaceCellSize();

	Nx = manager->getNumberOfSpaceCellsX();
	Ny = manager->getNumberOfSpaceCellsY();
	Nz = manager->getNumberOfSpaceCellsZ();

	n = getSpaceCell();

	x = getX(); y = getY(); z = getZ();
	vx = getVx(); vy = getVy(); vz = getVz();

	x += vx*(sTime - collisionTime);
	y += vy*(sTime - collisionTime);
	z += vz*(sTime - collisionTime);

// i, j and k are the indexes of the space cell for each axis
	i = n / (Nz*Ny);
	j = (n % (Nz*Ny)) / Nz;
	k = (n % (Nz*Ny)) % Nz;

// In a space cell we have 6 possible sides but since we know the direction of
// the particle we need to check only 3 (at most).
	if (vx > 0) sides.push_back(1);
	else if (vx < 0) sides.push_back(4);

	if (vy > 0) sides.push_back(3);
	else if (vy < 0) sides.push_back(2);

	if (vz > 0) sides.push_back(0);
	else if (vz < 0) sides.push_back(5);

	counter = 0;

	for (side = sides.begin(); side != sides.end(); ++side) {
// Select 3 points for each side following the previous rules (sp[])
		P.x = (i + sp[0*6 + *side])*spaceCellSize;
		P.y = (j + sp[1*6 + *side])*spaceCellSize;
		P.z = (k + sp[2*6 + *side])*spaceCellSize;

		Q.x = (i + sp[3*6 + *side])*spaceCellSize;
		Q.y = (j + sp[4*6 + *side])*spaceCellSize;
		Q.z = (k + sp[5*6 + *side])*spaceCellSize;

		R.x = (i + sp[6*6 + *side])*spaceCellSize;
		R.y = (j + sp[7*6 + *side])*spaceCellSize;
		R.z = (k + sp[8*6 + *side])*spaceCellSize;

		V.x = Q.x - P.x; V.y = Q.y - P.y; V.z = Q.z - P.z;
		W.x = R.x - P.x; W.y = R.y - P.y; W.z = R.z - P.z;

// Cross product to find the Normal vector
		N.x = V.y * W.z - V.z * W.y;
		N.y = V.z * W.x - V.x * W.z;
		N.z = V.x * W.y - V.y * W.x;

// This is our plane equation N.x*(x-P.x)+N.y*(y-P.y)+N.z*(z-P.z) = 0
// Substitute for: x = xi + vx*t, y = yi + vy*t, z = zi + vz*t and solve for t.
		temp = (N.x*vx + N.y*vy + N.z*vz);

// Check that temp really is greater than 0 ...
		if (temp*temp > 0) {

			temp = (N.x*(P.x-x) + N.y*(P.y-y) + N.z*(P.z-z))/temp;

// We found a solution :)
// 
// "temp" is the amount of time the particle takes to go from its last
// event point to the space cell side where it is bounded now.
			if (counter == 0) {

				transferTime = temp;
				hits.push_back(*side);

				if (temp == 0) break;

			} else {

				if (temp == 0) {
// The particle is on a space cell side and needs to update its space cell
// value and recalculate the transfer time.
					transferTime = temp;
					hits.push_back(*side);
					break;

				} else if (0 < temp && temp < transferTime) {

					transferTime = temp;

					hits.clear();
					hits.push_back(*side);

				} else if (0 < temp && temp == transferTime) {
// The particle hit two or more sides at the same time
					hits.push_back(*side);

				}

			}

			counter++;

		} else {
// Line and plane doesn't intersect
		}

	}

	if (transferTime > 0) {
		transferTime += sTime;
	}

	transferMsg->setEventTime(transferTime);
	transferMsg->setPrevSpaceCell(getSpaceCell());

	for (hit = hits.begin(); hit != hits.end(); ++hit) {

		if (*hit == 0) k = (k < Nz-1) ? k+1 : k;
		else if (*hit == 1) i = (i < Nx-1) ? i+1 : i;
		else if (*hit == 2) j = (j > 0) ? j-1 : j;
		else if (*hit == 3) j = (j < Ny-1) ? j+1 : j;
		else if (*hit == 4) i = (i > 0) ? i-1 : i;
		else k = (k > 0) ? k-1 : k;

	}

	transferMsg->setNextSpaceCell(i*Nz*Ny + j*Nz + k);

}

/*
 * Computes the collision time with the neighboring particles and returns the
 * smallest one.
 */
void Sphere::computeCollisionTime() {

	int a, b, c;			// Nested "for" loops

	int counter;

	int i, j, k, N, n;		// Indexes to access the current space cell
	int Nx, Ny, Nz;			// Number of space cells (or divisions) in each axis

	double collisionTime, sTime, temp;

	list<Particle *> particles;
	list<Particle *>::const_iterator p;

	Particle * partner; // The collision partner

	collisionTime = NO_TIME;
	sTime = simTime().dbl();

	Nx = manager->getNumberOfSpaceCellsX();
	Ny = manager->getNumberOfSpaceCellsY();
	Nz = manager->getNumberOfSpaceCellsZ();

	n = getSpaceCell();

// i, j and k are the indexes of the space cell for each axis
	i = n / (Nz*Ny);
	j = (n % (Nz*Ny)) / Nz;
	k = (n % (Nz*Ny)) % Nz;

	counter = 0;

	for (a = -1; a <= 1; a++) {
		for (b = -1; b <= 1; b++) {
			for (c = -1; c <= 1; c++) {

// The neighbor cell must be contained in the simulation space
				if (0 <= i+a && i+a < Nx && 
					0 <= j+b && j+b < Ny && 
					0 <= k+c && k+c < Nz) {

					N = (i+a)*Ny*Nz + (j+b)*Nz + (k+c);
					particles = manager->getSpaceCellParticles(N);

					for (p = particles.begin(); p != particles.end(); ++p) {

						if (*p == this) {
							continue;
						}

// Solve particle to particle collision
						temp = solveCollisionTime(*p);

						if (counter == 0) {

							if (sTime <= temp) {
								collisionTime = temp;
								partner = (*p);
							}

						} else {

							if (sTime <= temp  && temp < collisionTime) {
								collisionTime = temp;
								partner = (*p);
							}

						}
						
						counter++;

					}
				}

			}
		}
	}

	collisionMsg->setKind(EV_COLLISION);
	collisionMsg->setEventTime(collisionTime);
	collisionMsg->setPartner(partner);

}

/*
 *
 */
void Sphere::computeWallCollisionTime() {

	int counter;

	int sp[] = { // side points
	//  0, 1, 2, 3, 4, 5 <- dice side
		1, 1, 0, 0, 0, 0,
		0, 0, 0, 1, 0, 0,
		1, 0, 0, 0, 0, 0,
		1, 1, 1, 1, 0, 1,
		1, 1, 0, 1, 0, 0,
		1, 0, 0, 0, 1, 0,
		0, 1, 0, 0, 0, 0,
		0, 0, 0, 1, 1, 1,
		1, 1, 1, 1, 0, 0
	};

	double Sx, Sy, Sz;
	double rad, x, y, z, vx, vy, vz, vm; // Particle radius, position and velocity
	double wallCollisionTime,
		lastCollisionTime,
		temp;

	point_t P, Q, R;
	vect_t V, W, N;

	vector<int> sides, hits;
	vector<int>::const_iterator side, hit;

	wallCollisionTime = NO_TIME;
	lastCollisionTime = getLastCollisionTime();

	Sx = manager->getSpaceSizeX();
	Sy = manager->getSpaceSizeY();
	Sz = manager->getSpaceSizeZ();

	rad = getRadius();

	x = getX(); y = getY(); z = getZ();
	vx = getVx(); vy = getVy(); vz = getVz();

	vm = sqrt(vx*vx + vy*vy + vz*vz);

// In the simulation space we have 6 possible sides but since we know the 
// direction of the particle we need to check only 3 (at most).
	if (vx > 0) sides.push_back(1);
	else if (vx < 0) sides.push_back(4);

	if (vy > 0) sides.push_back(3);
	else if (vy < 0) sides.push_back(2);

	if (vz > 0) sides.push_back(0);
	else if (vz < 0) sides.push_back(5);

	counter = 0;

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

// Check that temp really is greater than 0 ...
		if (temp*temp > 0) {

			temp = (N.x*(P.x - x - rad*vx/vm) +
				N.y*(P.y - y - rad*vy/vm) +
				N.z*(P.z - z - rad*vz/vm))/temp;

// We found a solution :)
// 
// "temp" is the amount of time the particle takes to go from its last 
// event point to the space cell side where it is bounded now.
			if (counter == 0) {
				// First
				wallCollisionTime = temp;
				hits.push_back(*side);

			} else {

				if (temp > 0 && temp < wallCollisionTime) {

					wallCollisionTime = temp;

					hits.clear();
					hits.push_back(*side);

				} else if (temp > 0 && temp == wallCollisionTime) {
// The particle hit two or more sides at the same time
					hits.push_back(*side);

				}

			}

		} else {
// Line and plane doesn't intersect
		}

		counter++;

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

/*
 * Solves the particle to particle collision problem. When no collision occurs
 * it returns a negative number (-1).
 *
 * @param {Particle *} p
 * @return {double} the collision time
 */
double Sphere::solveCollisionTime(Particle *p) {

// Distance between centers A and B when t = tc (time of collision):
//                 ______________________________________________
//           \    / ( Ax + Avx*(tc-ta) - (Bx + Bvx*(tc-tb) )² + |
//  D(A, B) = \  /  ( Ay + Avy*(tc-ta) - (By + Bvy*(tc-tb) )² +   = Ra + Rb
//             \/   ( Az + Avz*(tc-ta) - (Bz + Bvz*(tc-tb) )²
//
// ta: when the previous collision take place for particle A
// tb: same for particle B

// (dxi + dvx*tc)² + (dyi + dvy*tc)² + (dyi + dvy*tc)²= (A.r + B.r)²

	double ta = getLastCollisionTime();
	double tb = p->getLastCollisionTime();

	double dxi = getX() - p->getX() - getVx()*ta + p->getVx()*tb;
	double dyi = getY() - p->getY() - getVy()*ta + p->getVy()*tb;
	double dzi = getZ() - p->getZ() - getVz()*ta + p->getVz()*tb;

	double dvx = getVx() - p->getVx();
	double dvy = getVy() - p->getVy();
	double dvz = getVz() - p->getVz();

	double radd = p->getRadius() + getRadius();

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

/*
 * Returns the eventTime from the scheduled collision event.
 *
 * @return {double} the event time 
 */
double Sphere::scheduledCollisionTime() {

	if (collisionMsg->isScheduled()) {

		return collisionMsg->getEventTime();

	} else if (wallCollisionMsg->isScheduled()) {

	    return wallCollisionMsg->getEventTime();

	} else {

	    return NO_TIME;

	}

}

/*
 * Tells the manager to update the space cell position for the particle.
 *
 * @param {MobilityMessage *} msg
 */
void Sphere::updateStateAfterTransfer(MobilityMessage *msg) {

	int prevSpaceCell, nextSpaceCell;

	prevSpaceCell = msg->getPrevSpaceCell();
	nextSpaceCell = msg->getNextSpaceCell();

	manager->transferParticle(this, prevSpaceCell, nextSpaceCell);

	setSpaceCell(nextSpaceCell);

}

/*
 * Updates the particle position, velocity and the last collision time values
 * and does the same for the partner particle.
 *
 * @param {MobilityMessage *} msg
 */
void Sphere::updateStateAfterCollision(MobilityMessage *msg) {

	double tc;
	double v1m, v2m;
	double e1m, e2m;
	double m1, m2;
	double th, temp;

	point_t C1, C2, P, Q;

	vect_t v1, v2;      // Velocity vectors for particle 1 and 2
	vect_t e1, e2, e3;  // Orthonormal basis
	vect_t v1e, v2e;    // Velocity of particles 1 and 2 for the basis found

	Particle * p;

	tc = msg->getEventTime();
	p = msg->getPartner();

	m1 = getMass();
	m2 = p->getMass();

	v1.x = getVx();
	v1.y = getVy();
	v1.z = getVz();

	v1m = sqrt(v1.x*v1.x + v1.y*v1.y + v1.z*v1.z);

	v2.x = p->getVx();
	v2.y = p->getVy();
	v2.z = p->getVz();

	v2m = sqrt(v2.x*v2.x + v2.y*v2.y + v2.z*v2.z);

// Compute the particle center at the collision time.
	C1.x = getX() + v1.x*(tc-getLastCollisionTime());
	C1.y = getY() + v1.y*(tc-getLastCollisionTime());
	C1.z = getZ() + v1.z*(tc-getLastCollisionTime());

// Compute the partner center at the collision time.
	C2.x = p->getX() + v2.x*(tc - p->getLastCollisionTime());
	C2.y = p->getY() + v2.y*(tc - p->getLastCollisionTime());
	C2.z = p->getZ() + v2.z*(tc - p->getLastCollisionTime());

// Find an orthonormal basis
	e1.x = C2.x - C1.x;
	e1.y = C2.y - C1.y;
	e1.z = C2.z - C1.z;

	e1m = sqrt(e1.x*e1.x + e1.y*e1.y + e1.z*e1.z);
	e1.x /= e1m;
	e1.y /= e1m;
	e1.z /= e1m;

// Find the collision point
	P.x = C1.x + getRadius()*e1.x;
	P.y = C1.y + getRadius()*e1.y;
	P.z = C1.z + getRadius()*e1.z;

// Now that we have the collision plane: A(x-x0) + B(y-y0) + C(z-z0) = 0
// e1.x*(x-P.x) + e1.y*(y-P.y) + e1.z*(z-P.z) = 0, we can find one point
// contained in that plane and form an orthonormal basis.

	if (e1.x != 0 && (e1.y != 0 || e1.z != 0)) {

		Q.x = (e1.y*P.y + e1.z*P.z)/e1.x + P.x;
		Q.y = 0;
		Q.z = 0;

	} else if (e1.y != 0 && (e1.x != 0 || e1.z != 0)) {

		Q.x = 0;
		Q.y = (e1.x*P.x + e1.z*P.z)/e1.y + P.y;
		Q.z = 0;

	} else if (e1.z != 0 && (e1.x != 0 || e1.y != 0)) {

		Q.x = 0;
		Q.y = 0;
		Q.z = (e1.x*P.x + e1.y*P.y)/e1.z + P.z;

	} else {
// TODO something went wrong ... particles are at the same point
	}

// Find the PQ vector
	e2.x = Q.x - P.x;
	e2.y = Q.y - P.y;
	e2.z = Q.z - P.z;

	e2m = sqrt(e2.x*e2.x + e2.y*e2.y + e2.z*e2.z);

	e2.x /= e2m;
	e2.y /= e2m;
	e2.z /= e2m;

	e3.x = e1.y*e2.z - e1.z*e2.y;
	e3.y = e1.z*e2.x - e1.x*e2.z;
	e3.z = e1.x*e2.y - e1.y*e2.x;

// Find the velocity vectors in the new basis
	if (v1m > 0) {

		th = acos((v1.x*e1.x + v1.y*e1.y + v1.z*e1.z)/v1m);
		v1e.x = v1m*cos(th);

		th = acos((v1.x*e2.x + v1.y*e2.y + v1.z*e2.z)/v1m);
		v1e.y = v1m*cos(th);

		th = acos((v1.x*e3.x + v1.y*e3.y + v1.z*e3.z)/v1m);
		v1e.z = v1m*cos(th);

	} else {

		v1e.x = 0;
		v1e.y = 0;
		v1e.z = 0;

	}

	if (v2m > 0) {

		th = acos((v2.x*e1.x + v2.y*e1.y + v2.z*e1.z)/v2m);
		v2e.x = v2m*cos(th);

		th = acos((v2.x*e2.x + v2.y*e2.y + v2.z*e2.z)/v2m);
		v2e.y = v2m*cos(th);

		th = acos((v2.x*e3.x + v2.y*e3.y + v2.z*e3.z)/v2m);
		v2e.z = v2m*cos(th);

	} else {

		v2e.x = 0;
		v2e.y = 0;
		v2e.z = 0;

	}

// Solve the exchange of momentum between particles. Note that only the normal
// velocity component changes its value.
	temp = ((m1-m2)*v1e.x + 2*m2*v2e.x)/(m1 + m2);
	v2e.x = ((m2-m1)*v2e.x + 2*m1*v1e.x)/(m1 + m2);
	v1e.x = temp;

// Revert the velocity vectors and set the new velocity
	setVx(v1e.x*e1.x + v1e.y*e2.x + v1e.z*e3.x);
	setVy(v1e.x*e1.y + v1e.y*e2.y + v1e.z*e3.y);
	setVz(v1e.x*e1.z + v1e.y*e2.z + v1e.z*e3.z);

	p->setVx(v2e.x*e1.x + v2e.y*e2.x + v2e.z*e3.x);
	p->setVy(v2e.x*e1.y + v2e.y*e2.y + v2e.z*e3.y);
	p->setVz(v2e.x*e1.z + v2e.y*e2.z + v2e.z*e3.z);

// Update the particles position
	setX(C1.x);
	setY(C1.y);
	setZ(C1.z);

	p->setX(C2.x);
	p->setY(C2.y);
	p->setZ(C2.z);

// Update the last collision times
	setLastCollisionTime(tc);
	p->setLastCollisionTime(tc);

}

/*
 * Updates the position, velocity and the last collision time when the particle
 * collides with a wall.
 *
 * @param {MobilityMessage *} msg
 */
void Sphere::updateStateAfterWallCollision(MobilityMessage *msg) {

	setX(msg->getX());
	setY(msg->getY());
	setZ(msg->getZ());

	setVx(msg->getVx());
	setVy(msg->getVy());
	setVz(msg->getVz());

	setLastCollisionTime(msg->getEventTime());

}
