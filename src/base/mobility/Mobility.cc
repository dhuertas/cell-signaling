//  This file is part of the Cell-Signaling project. Cell-Signaling is an
//  Omnet++ project to simulate cell signaling communications.
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
	int *Nx, *Ny, *Nz;		// Number of space cells (or divisions) in each axis

	int counter;

	double x, y, z;

	point_t *pos = NULL;
	vect_t *vel = NULL;

	double spaceCellSize;
	double temp;
	double transferTime;
	double sTime;
	double lastCollisionTime;

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

	pos = p->getPosition();
	vel = p->getVelocity();

	x = pos->x + vel->x*(sTime - lastCollisionTime);
	y = pos->y + vel->y*(sTime - lastCollisionTime);
	z = pos->z + vel->z*(sTime - lastCollisionTime);

	// i, j and k are the indexes of the space cell for each axis
	i =  n/((*Nz)*(*Ny));
	j = (n%((*Nz)*(*Ny)))/(*Nz);
	k = (n%((*Nz)*(*Ny)))%(*Nz);

	// In a space cell we have 6 possible sides but since we know the direction 
	// of the particle we need to check only 3 (at most).
	if (vel->x > 0) sides.push_back(1);
	else if (vel->x < 0) sides.push_back(4);

	if (vel->y > 0) sides.push_back(3);
	else if (vel->y < 0) sides.push_back(2);

	if (vel->z > 0) sides.push_back(0);
	else if (vel->z < 0) sides.push_back(5);

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
		// Substitute for: x = xi + vx*t, y = yi + vy*t, z = zi + vz*t and 
		// solve for t.
		temp = (N.x*vel->x + N.y*vel->y + N.z*vel->z);

		if (temp != 0) {

			temp = (N.x*(P.x-x) + N.y*(P.y-y) + N.z*(P.z-z))/temp;

			// Solution found. "temp" is the amount of time the centroid of the
			// particle takes to go from its current position to the space cell 
			// side where it is bounded.
			if (counter == 0) {

				if (temp == 0) {
					// The centroid of the particle is on a space cell side and
					// needs to update its space cell value and recalculate the
					// transfer time.
				    transferTime = temp;
				    hits.push_back(*side);
				    counter++;
					break;
				} else if (0 < temp) {
				    transferTime = temp;
				    hits.push_back(*side);
				    counter++;
				} else {
				    // We don't want it
				}

			} else {

				if (temp == 0) {
					// The centroid of the particle is on a space cell side and
					// needs to update its space cell value and recalculate the 
					// transfer time.
					transferTime = temp;
					hits.push_back(*side);
					counter++;
					break;

				} else if (0 < temp && temp < transferTime) {

					transferTime = temp;

					hits.clear();
					hits.push_back(*side);
					counter++;

				} else if (0 < temp && temp == transferTime) {
					// The particle hit two or more sides at the same time
					hits.push_back(*side);
					counter++;

				}

			}

		} else {
			// Line and plane doesn't intersect
		}

	}

	if (transferTime > 0) {

		transferTime += sTime;

	} else if (transferTime == 0) {
		// If transferTime equals the simuation time sTime means that temp 
		// equals 0, thus the centroid of the particle belongs to the plane it 
		// is crossing. Set an event transfer at the same simulation time so it
		// will update the NextSpaceCell value and compute again the next 
		// transfer time.
		transferTime = sTime;
	} else {
		// transfer time not found (NO_TIME)
	}

	msg->setTransferTime(transferTime);
	msg->setPrevSpaceCell(p->getSpaceCell());

	for (hit = hits.begin(); hit != hits.end(); ++hit) {

		if (*hit == 0) k = (k < *Nz-1) ? k+1 : k;
		else if (*hit == 1) i = (i < *Nx-1) ? i+1 : i;
		else if (*hit == 2) j = (j > 0) ? j-1 : j;
		else if (*hit == 3) j = (j < *Ny-1) ? j+1 : j;
		else if (*hit == 4) i = (i > 0) ? i-1 : i;
		else k = (k > 0) ? k-1 : k;

	}

	msg->setNextSpaceCell(i*(*Nz)*(*Ny)+j*(*Nz)+k);

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

//	double vx, vy, vz, vm, lr, rlr;
	double vm, lr, rlr;
	double outOfNeighborhoodTime;
	double sTime;

	vect_t *vel = NULL;

	vel = p->getVelocity();

	lr = 0;
	rlr = 0;

	sTime = simTime().dbl();
	outOfNeighborhoodTime = NO_TIME;

	lr = p->getListRadius();
	rlr = p->getRefreshListRadius();

	if (rlr > 0) lr = rlr;

	vm = sqrt(vel->x*vel->x + vel->y*vel->y + vel->z*vel->z);

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
