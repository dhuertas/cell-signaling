#include "Mobility.h"

using namespace std;

// Space cell side points that allows to obtain the side equation and compute 
// the transfer time. Defined in src/base/Defines.h file
char sp[10*6] = { 
	//  0, 1, 2, 3, 4, 5 <- dice side
		1, 1, 0, 0, 0, 0, // Px
		0, 0, 0, 1, 0, 0, // Py
		1, 0, 0, 0, 0, 0, // Pz
		1, 1, 1, 1, 0, 1, // Qx
		1, 1, 0, 1, 0, 0, // Qy
		1, 0, 0, 0, 1, 0, // Qz
		0, 1, 0, 0, 0, 0, // Rx
		0, 0, 0, 1, 1, 1, // Ry
		1, 1, 1, 1, 0, 0  // Rz
	};

/*
 * Computes the transfer time of the particle with the sides of its space cell
 * and returns the smallest one.
 *
 * @param {TransferMessage *} msg
 * @param {Particle *} p
 */
double Mobility::nextTransfer(TransferMessage *msg, Particle *p) {

	int i, j, k, n;			// Indexes to access the current space cell
	int Nx, Ny, Nz;			// Number of space cells (or divisions) in each axis

	int counter;

	double x, y, z, vx, vy, vz;	// Particle position and velocity

	double spaceCellSize, transferTime, sTime, lastCollisionTime, temp;

	point_t P, Q, R;
	vect_t V, W, N;

	vector<int> sides, hits;
	vector<int>::const_iterator side, hit;

	Manager *manager;

	manager = msg->getManager();

	spaceCellSize = manager->getSpaceCellSize();

	transferTime = NO_TIME;

	sTime = simTime().dbl();
	lastCollisionTime = p->getLastCollisionTime();

	temp = NO_TIME;

	Nx = manager->getNumberOfSpaceCellsX();
	Ny = manager->getNumberOfSpaceCellsY();
	Nz = manager->getNumberOfSpaceCellsZ();

	n = p->getSpaceCell();

	x = p->getX();
	y = p->getY();
	z = p->getZ();
	
	vx = p->getVx();
	vy = p->getVy();
	vz = p->getVz();

	x += vx*(sTime - lastCollisionTime);
	y += vy*(sTime - lastCollisionTime);
	z += vz*(sTime - lastCollisionTime);

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

		if (temp != 0) {

			temp = (N.x*(P.x-x) + N.y*(P.y-y) + N.z*(P.z-z))/temp;

// Solution found. "temp" is the amount of time the centroid of the particle 
// takes to go from its current position to the space cell side where it is 
// bounded.
			if (counter == 0) {

				transferTime = temp;
				hits.push_back(*side);

				if (temp == 0) {
// The centroid of the particle is on a space cell side and needs to update its
// space cell value and recalculate the transfer time.
					break;
				}

			} else {

				if (temp == 0) {
// The centroid of the particle is on a space cell side and needs to update its
// space cell value and recalculate the transfer time.
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

	} else if (transferTime == 0) {
// If transferTime equals the simuation time sTime it means that temp 
// equals = 0, thus the centroid of the particle belongs to the plane it is 
// crossing. Set an event transfer at the same simulation time so it will 
// update the NextSpaceCell value and compute again the next transfer time.
	} else {
// transfer time not found (NO_TIME)
	}

	msg->setTransferTime(transferTime);
	msg->setPrevSpaceCell(p->getSpaceCell());

	for (hit = hits.begin(); hit != hits.end(); ++hit) {

		if (*hit == 0) k = (k < Nz-1) ? k+1 : k;
		else if (*hit == 1) i = (i < Nx-1) ? i+1 : i;
		else if (*hit == 2) j = (j > 0) ? j-1 : j;
		else if (*hit == 3) j = (j < Ny-1) ? j+1 : j;
		else if (*hit == 4) i = (i > 0) ? i-1 : i;
		else k = (k > 0) ? k-1 : k;

	}

	msg->setNextSpaceCell(i*Nz*Ny + j*Nz + k);

	return transferTime;

}

/*
 * Computes the time when the particle leaves its neighborhood area and its
 * Near-Neighbor List must be updated.
 *
 * @param {OutOfNeighborhoodMessage *} msg
 * @param {Particle *} p
 * @return {double} the computed time
 */
double Mobility::outOfNeighborhoodTime(OutOfNeighborhoodMessage *msg, Particle *p) {

	double vx, vy, vz, vm, lr;

	double outOfNeighborhoodTime;
	double sTime;

	vx = p->getVx();
	vy = p->getVy(); 
	vz = p->getVz();

	sTime = simTime().dbl();
	outOfNeighborhoodTime = NO_TIME;

	lr = p->getListRadius();

	vm = sqrt(vx*vx + vy*vy + vz*vz);

	if (vm != 0) {

		outOfNeighborhoodTime = lr/vm + sTime;

	}

    return outOfNeighborhoodTime;

}

/*
 * Set default values for the collision message.
 *
 * @param {CollisionMessage *} msg
 */
void Mobility::resetCollisionMessage(CollisionMessage *msg) {

	msg->setCollisionTime(NO_TIME);

	msg->setX(-1);
	msg->setY(-1);
	msg->setZ(-1);

	msg->setVx(0);
	msg->setVy(0);
	msg->setVz(0);

	msg->setPartner(NULL);
	msg->setPrevPartner(NULL);

}
