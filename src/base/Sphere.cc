#include "Sphere.h"
#include "../messages/WallCollision_m.h"
#include "../messages/Transfer_m.h"

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
 */
void Sphere::tkEnvUpdatePosition(double t) {

	std::stringstream buffer;

	// Set position string for tkenv
	buffer << getY() + getVy()*(t-getLastCollisionTime());
	EV << "y: " << buffer.str().c_str() << "\n";

	getDisplayString().setTagArg("p", 0, buffer.str().c_str());
	buffer.str(std::string()); // clear buffer

	buffer << getX() + getVx()*(t-getLastCollisionTime());
	EV << "x: " << buffer.str().c_str() << "\n";
	getDisplayString().setTagArg("p", 1, buffer.str().c_str());
	buffer.str(std::string());

}

/*
 * Generates a cMessage object for the next event. That event can be a transfer
 * event, a particle-particle collision event or a wall bounce event.
 *
 * @param {cMessage *} msg
 */
void Sphere::nextEventTime() {  // TODO rethink this function

    Enter_Method_Silent(); // Methods called from other modules must have this macro

    double transferTime, 
    	wallCollisionTime;

	cMessage *transferMsg, 
		*collisionMsg, 
		*wallCollisionMsg;

	// Step 1. Get the next space cell transfer time
	transferMsg = computeTransferTime();
	transferTime = ((TransferMessage *)transferMsg)->getTransferTime();

	// Step 2. Get the next collision time with the surrounding particles or 
	// walls. We only consider collisions with the particles in the 
	// neighbouring space cells (9+9+(8+1) cells)
	// collisionMsg = computeCollisionTime();

	// Step 3. Get the next wall collision time
	wallCollisionMsg = computeWallCollisionTime();
	wallCollisionTime = ((WallCollisionMessage *)wallCollisionMsg)
		->getCollisionTime();

	if (wallCollisionTime > 0) {

		scheduleAt(((WallCollisionMessage *)wallCollisionMsg)
			->getCollisionTime(), wallCollisionMsg);

	} else delete wallCollisionMsg;

}

/*
 * Computes the transfer time of the particle with the sides of its space cell
 * and returns the smallest one.
 * 
 * @return {double}: transfer time
 */
cMessage * Sphere::computeTransferTime() {

	int i, j, k, n, r, t;		// Indexes to access the current space cell
	int Ny, Nz;					// Number of space cells (or divisions) in each axis

	int transferSpaceCell, counter;

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
		transferTime;

	point_t P, Q, R;
	vect_t V, W, N;

	vector<int> sides;
	vector<int>::const_iterator side;

	cMessage * msg = new TransferMessage("TransferMessage", EV_TRANSFER);

	transferSpaceCell = -1;
	transferTime = -1;

	spaceCellSize = manager->getSpaceCellSize();

	Ny = manager->getNumberOfSpaceCellsY();
	Nz = manager->getNumberOfSpaceCellsZ();

	n = getSpaceCell();

	x = getX();
	y = getY();
	z = getZ();
	vx = getVx();
	vy = getVy();
	vz = getVz();

	// i, j and k are the indexes of the space cell for each axis
	i = n / (Nz*Ny);
	r = n % (Nz*Ny);
	j = r / Nz;
	k = r % Nz;

	// In a space cell we have 6 possible sides but since we know the direction
	// of the particle we need to check only 3 (at most).
	if (vx > 0) {
		sides.push_back(1);
	} else if (vx < 0) {
		sides.push_back(4);
	}

	if (vy > 0) {
		sides.push_back(3);
	} else if (vy < 0) {
		sides.push_back(2);
	}

	if (vz > 0) {
		sides.push_back(0);
	} else if (vz < 0) {
		sides.push_back(5);
	}

	counter = 0;

	for (side = sides.begin(); side != sides.end(); ++side) {
		// Select 3 points for each side following the previous rules
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

		// Cross product
		N.x = V.y * W.z - V.z * W.y;
		N.y = V.z * W.x - V.x * W.z;
		N.z = V.x * W.y - V.y * W.x;

		// This is our plane equation N.x*(x-P.x)+N.y*(y-P.y)+N.z*(z-P.z) = 0
		// Substitute for: x = xi + vx*t, y = yi + vy*t, z = zi + vz*t and find
		// for t.
		t = (N.x*vx + N.y*vy + N.z*vz);
		
		if (t != 0) {
			t = (N.x*(P.x-x) + N.y*(P.y-y) + N.z*(P.z-z))/t;
			// We found a solution :)
			// 
			// "t" is the amount of time the particle takes to go from its
			// last collision point to the space cell side where it is bounded 
			// now.

			// WARNING! The particle could hit the three sides at the same time.
			// TODO solve for multiple sides
			if (counter == 0) {

				transferTime = t + getLastCollisionTime();
				// transferSpaceCell = 0;

			} else {

				if (t > 0 && t + getLastCollisionTime() < transferTime) {
					transferTime = t + getLastCollisionTime();
					// transferSpaceCell = 0;
				}

			}
			

		} else {
			// Line and plane doesn't intersect
		}

		counter++;

	}

	((TransferMessage *)msg)->setTransferTime(transferTime);
	((TransferMessage *)msg)->setPrevSpaceCell(getSpaceCell());
	((TransferMessage *)msg)->setNextSpaceCell(transferSpaceCell);

	return msg;

}

/*
 * Computes the collision time with the neighbouring particles and returns the
 * smallest one.
 * 
 * @return {double}: collision time
 */
cMessage * Sphere::computeCollisionTime() {

	int a, b, c;				// Nested "for" loops

	int i, j, k, N, n, r;		// Indexes to access the current space cell
	int Nx, Ny, Nz;				// Number of space cells (or divisions) in each axis

	double collisionTime, t;

	list<Particle *> particles;
	list<Particle *>::const_iterator p;

	cMessage * msg;

	collisionTime = -1;

	Nx = manager->getNumberOfSpaceCellsX();
	Ny = manager->getNumberOfSpaceCellsY();
	Nz = manager->getNumberOfSpaceCellsZ();

	n = getSpaceCell();

	// i, j and k are the indexes of the space cell for each axis
	i = n / (Nz*Ny);
	r = n % (Nz*Ny);
	j = r / Nz;
	k = r % Nz;

	for (a = -1; a <= 1; a++) {
		for (b = -1; b <= 1; b++) {
			for (c = -1; c <= 1; c++) {

				// The neighbour cell must be contained in the simulation space
				if (0 <= i+a && i+a < Nx && 
					0 <= j+b && j+b < Ny && 
					0 <= k+c && k+c < Nz) {

					N = (i+a)*Ny*Nz + (j+b)*Nz + (k+c);
					particles = manager->getSpaceCellParticles(N);

					for (p = particles.begin(); p != particles.end(); ++p) {
						// Solve particle to particle problem
						t = solveCollisionTime(*p);

						if (collisionTime > t) {
							collisionTime = t;
							// candidate = (*p);
						}
					}
				}

			}
		}
	}

	return msg;

}

/*
 *
 */
cMessage * Sphere::computeWallCollisionTime() {

	int counter, wallCollision;

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

	double Sx, Sy, Sz, t;
	double x, y, z, vx, vy, vz;	// Particle position and velocity
	double wallCollisionTime;

	point_t P, Q, R;
	vect_t V, W, N;

	vector<int> sides;
	vector<int>::const_iterator side;

	cMessage * msg = new WallCollisionMessage("WallCollision", EV_WALLCOLLISION);

	wallCollision = -1;
	wallCollisionTime = -1;

	Sx = manager->getSpaceSizeX();
	Sy = manager->getSpaceSizeY();
	Sz = manager->getSpaceSizeZ();

	x = getX();
	y = getY();
	z = getZ();
	vx = getVx();
	vy = getVy();
	vz = getVz();

	// In the simulation space we have 6 possible sides but since we know the 
	// direction of the particle we need to check only 3 (at most).
	if (vx > 0) {
		sides.push_back(1);
	} else if (vx < 0) {
		sides.push_back(4);
	}

	if (vy > 0) {
		sides.push_back(3);
	} else if (vy < 0) {
		sides.push_back(2);
	}

	if (vz > 0) {
		sides.push_back(0);
	} else if (vz < 0) {
		sides.push_back(5);
	}

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

		// Cross product
		N.x = V.y * W.z - V.z * W.y;
		N.y = V.z * W.x - V.x * W.z;
		N.z = V.x * W.y - V.y * W.x;

		// This is our plane equation N.x*(x-P.x)+N.y*(y-P.y)+N.z*(z-P.z) = 0
		// Substitute for: x = xi + vx*t, y = yi + vy*t, z = zi + vz*t and find
		// for t.
		t = (N.x*vx + N.y*vy + N.z*vz);

		if (t != 0) {
			t = (N.x*(P.x-x) + N.y*(P.y-y) + N.z*(P.z-z))/t;
			// We found a solution :)
			// 
			// "t" is the amount of time the particle takes to go from its
			// last collision point to the space cell side where it is bounded 
			// now.

			// WARNING! The particle could hit two or the three sides of the space
			// cell.
			// TODO solve for multiple sides
			if (counter == 0) {
				// First
				wallCollisionTime = t + getLastCollisionTime();
				wallCollision = *side;

			} else {

				if (t > 0 && t + getLastCollisionTime() < wallCollisionTime) {
					wallCollisionTime = t + getLastCollisionTime();
					wallCollision = *side;
				}

			}

		} else {
			// Line and plane doesn't intersect
		}

		counter++;

	}

	// Compute the future values
	((WallCollisionMessage *)msg)->setCollisionTime(wallCollisionTime);

	((WallCollisionMessage *)msg)->setX(x + vx*(wallCollisionTime - getLastCollisionTime()));
	((WallCollisionMessage *)msg)->setY(y + vy*(wallCollisionTime - getLastCollisionTime()));
	((WallCollisionMessage *)msg)->setZ(z + vz*(wallCollisionTime - getLastCollisionTime()));

    ((WallCollisionMessage *)msg)->setVx((wallCollision == 1 || wallCollision == 4) ? -vx : vx);
    ((WallCollisionMessage *)msg)->setVy((wallCollision == 2 || wallCollision == 3) ? -vy : vy);
    ((WallCollisionMessage *)msg)->setVz((wallCollision == 0 || wallCollision == 5) ? -vz : vz);

	return msg;

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
 *
 */
void Sphere::updateStateAfterCollision(cMessage *msg) {

}

/*
 *
 */
void Sphere::updateStateAfterWallCollision(cMessage *msg) {

	setX(((WallCollisionMessage *)msg)->getX());
	setY(((WallCollisionMessage *)msg)->getY());
	setZ(((WallCollisionMessage *)msg)->getZ());

	setVx(((WallCollisionMessage *)msg)->getVx());
	setVy(((WallCollisionMessage *)msg)->getVy());
	setVz(((WallCollisionMessage *)msg)->getVz());

	setLastCollisionTime(((WallCollisionMessage *)msg)->getCollisionTime());

}
