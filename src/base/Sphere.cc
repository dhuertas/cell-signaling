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

    EV << "x: " << getX()*getVx() << ", y: " << getY() << ", z: " << getZ() << "\n";

}

/*
 * Update the module position in the tk environment. This method is used 
 * during the simulation.
 *
 * @param {Double} t
 */
void Sphere::tkEnvUpdatePosition(double t) {

	std::stringstream buffer;

	double lc = getLastCollisionTime();

// Set position string for tkenv
	buffer << getY() + getVy()*(t-getLastCollisionTime());

	getDisplayString().setTagArg("p", 0, buffer.str().c_str());
	buffer.str(std::string()); // clear buffer

	buffer << getX() + getVx()*(t-getLastCollisionTime());

	getDisplayString().setTagArg("p", 1, buffer.str().c_str());
	buffer.str(std::string());

	EV << "x: " << getX() + getVx()*(t - lc) <<
        ", y: " << getY() + getVy()*(t - lc) <<
        ", z: " << getZ() + getVz()*(t - lc) << "\n";

}

/*
 * Initialize the event queue by computing the first event for the particle.
 * This function is called by the manager in order to initialize the event
 * queue.
 */
void Sphere::firstEventTime() {

// Methods called from other modules must have this macro
	Enter_Method_Silent();

	double transferTime, collisionTime, wallCollisionTime;

	transferMsg = new MobilityMessage("mobility", EV_TRANSFER);
	collisionMsg = new MobilityMessage("mobility", EV_COLLISION);
	wallCollisionMsg =  new MobilityMessage("mobility", EV_WALLCOLLISION);

// Compute the first collision and the first transfer
	computeTransferTime();
    computeWallCollisionTime();
    computeCollisionTime();

    transferTime = transferMsg->getEventTime();
    collisionTime = collisionMsg->getEventTime();
    wallCollisionTime = wallCollisionMsg->getEventTime();

	if (transferTime > 0) {
		scheduleAt(transferTime, transferMsg);
	}

// Schedule only the first event

	if (0 < collisionTime && collisionTime < wallCollisionTime) {

		scheduleAt(collisionTime, collisionMsg);

	} else if (0 < wallCollisionTime &&
	        (wallCollisionTime < collisionTime || collisionTime == NO_TIME)) {

		scheduleAt(wallCollisionTime, wallCollisionMsg);

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
// Must check that the partner is solving the collision (its next event is not
// of type EV_CHECK)
                        //collisionMsg->getPartner()->
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

	int collisionCounter;

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

	collisionCounter = 0;

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

	int collisionCounter;

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
	double rad, x, y, z, vx, vy, vz; // Particle radius, position and velocity
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

	// vm = sqrt(vx*vx + vy*vy + vz*vz);

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


// We found a solution :)
// 
// "temp" is the amount of time the particle takes to go from its last 
// event point to the simulation space side.
			if (collisionCounter == 0) {
				// First
				wallCollisionTime = temp;
				hits.push_back(*side);
				collisionCounter++;

			} else {

				if (temp > 0 && temp < wallCollisionTime) {

					wallCollisionTime = temp;

					hits.clear();
					hits.push_back(*side);
					collisionCounter = 1;

				} else if (temp > 0 && temp == wallCollisionTime) {
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
 * Handles the mobility of the sphere
 *
 * @param {MobilityMessage *} msg
 */
void Sphere::handleMobilityMessage(MobilityMessage *msg) {

    int kind = msg->getKind();

// Step 2. Handle the event
    if (kind == EV_TRANSFER) {
        // Update the molecule space cell
        updateStateAfterTransfer((MobilityMessage *)msg);
        nextEventTime();

    } else if (kind == EV_WALLCOLLISION) {
        // Update the molecule data
        updateStateAfterWallCollision((MobilityMessage *)msg);
        nextEventTime();

    } else if (kind == EV_COLLISION) {

        updateStateAfterCollision((MobilityMessage *)msg);
        nextEventTime();

    } else if (kind == EV_CHECK) {
// TODO our collision time has turned invalid. We must check again for the next
// event.

// Two cases are possible:
// 1. We obtained an expected collision time, but the partner has a smaller
// collision time. Therefore we can calculate next collision at the expected
// collision time.
//
// 2. Someone else has forced us to recompute our collision time since it
// expects a collision. If we have a scheduled collision event, we must cancel
// it telling the partner to check again for its next collision time (thus
// going to case 1).
        nextEventTime();

    } else {

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

    double tc, m1, m2, tmp;
    double v1n, v1e1, v1e2, v2n, v2e1, v2e2;

    point_t c1, c2;
    vect_t v1, n, e1, e2;

    Particle *p;

    tc = msg->getEventTime();
    p = msg->getPartner();

// Find the center position of the spheres
    c1.x = this->getX() + this->getVx()*(tc - this->getLastCollisionTime());
    c1.y = this->getY() + this->getVy()*(tc - this->getLastCollisionTime());
    c1.z = this->getZ() + this->getVz()*(tc - this->getLastCollisionTime());

    c2.x = p->getX() + p->getVx()*(tc - p->getLastCollisionTime());
    c2.y = p->getY() + p->getVy()*(tc - p->getLastCollisionTime());
    c2.z = p->getZ() + p->getVz()*(tc - p->getLastCollisionTime());

    m1 = this->getMass();
    m2 = p->getMass();

// Change frame of reference of the system to one of the spheres
    v1.x = this->getVx() - p->getVx();
    v1.y = this->getVy() - p->getVy();
    v1.z = this->getVz() - p->getVz();

// Find the normal vector of the plane of collision
    n.x = c2.x - c1.x;
    n.y = c2.y - c1.y;
    n.z = c2.z - c1.z;

    tmp = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);

    n.x /= tmp;
    n.y /= tmp;
    n.z /= tmp;

// Find e1 as the perpendicular vector to both n and v, and then e2 as the one
// perpendicular to n and e1
    e1.x = n.y*v1.z - n.z*v1.y;
    e1.y = n.z*v1.x - n.x*v1.z;
    e1.z = n.x*v1.y - n.y*v1.x;

// Normalize the vector found, e1
    tmp = sqrt(e1.x*e1.x + e1.y*e1.y + e1.z*e1.z);

    e1.x /= tmp;
    e1.y /= tmp;
    e1.z /= tmp;

// Find the velocity vectors in the new basis and if ...
    v1n  = v1.x*n.x  + v1.y*n.y  + v1.z*n.z;

    if (e1.x == 0.0 && e1.y == 0.0 && e1.z == 0.0) {
// n and v are parallel, we can solve directly
        tmp = (m1 - m2)*v1n/(m1 + m2);
        v2n = 2*m1*v1n/(m1 + m2);
        v1n = tmp;

// Revert the frame of reference, the velocity vectors and set the new velocity
        setVx(v1n*n.x + p->getVx());
        setVy(v1n*n.y + p->getVy());
        setVz(v1n*n.z + p->getVz());

        p->setVx(v2n*n.x + p->getVx());
        p->setVy(v2n*n.y + p->getVy());
        p->setVz(v2n*n.z + p->getVz());

    } else {

        e2.x = e1.y*n.z - e1.z*n.y;
        e2.y = e1.z*n.x - e1.x*n.z;
        e2.z = e1.x*n.y - e1.y*n.x;

// Find the rest of the components
        v1e1 = v1.x*e1.x + v1.y*e1.y + v1.z*e1.z;
        v1e2 = v1.x*e2.x + v1.y*e2.y + v1.z*e2.z;

        v2n  = 0.0;
        v2e1 = 0.0;
        v2e2 = 0.0;

// Find the new velocity in the normal component (remember that v2n initially
// is 0.0)
        tmp = (m1 - m2)*v1n/(m1 + m2);
        v2n = 2*m1*v1n/(m1 + m2);
        v1n = tmp;

// Revert the frame of reference, the velocity vectors and set the new velocity
        setVx(v1n*n.x + v1e1*e1.x + v1e2*e2.x + p->getVx());
        setVy(v1n*n.y + v1e1*e1.y + v1e2*e2.y + p->getVy());
        setVz(v1n*n.z + v1e1*e1.z + v1e2*e2.z + p->getVz());

        p->setVx(v2n*n.x + v2e1*e1.x + v2e2*e2.x + p->getVx());
        p->setVy(v2n*n.y + v2e1*e1.y + v2e2*e2.y + p->getVy());
        p->setVz(v2n*n.z + v2e1*e1.z + v2e2*e2.z + p->getVz());

    }

// Update the particles position
    setX(c1.x);
    setY(c1.y);
    setZ(c1.z);

    p->setX(c2.x);
    p->setY(c2.y);
    p->setZ(c2.z);

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
