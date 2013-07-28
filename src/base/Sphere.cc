#include "Sphere.h"
#include "mobility/SphereMobility.h"

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

	if (mode == M_NNLIST) {
		buffer << getListRadius();

		getDisplayString().setTagArg("r",0, buffer.str().c_str());
		buffer.str(std::string());
	}

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

	double transferTime;
	double collisionTime;
	double wallCollisionTime;
	double outOfNeighborhoodTime;
	double smallestTime;

	vector<double> times;
	vector<double>::const_iterator t;

	transferMsg = new MobilityMessage("mobility", EV_TRANSFER);
	collisionMsg = new MobilityMessage("mobility", EV_COLLISION);
	wallCollisionMsg =  new MobilityMessage("mobility", EV_WALLCOLLISION);
	outOfNeighborhoodMsg = new MobilityMessage("mobility", EV_OUTOFNEIGHBORHOOD);

	transferMsg->setManager(manager);
	collisionMsg->setManager(manager);
	wallCollisionMsg->setManager(manager);
	outOfNeighborhoodMsg->setManager(manager);

// Compute the first collision and the first transfer
	SphereMobility::transferTime(transferMsg, this);
	SphereMobility::wallCollisionTime(wallCollisionMsg, this);
	SphereMobility::collisionTime(collisionMsg, this);

	transferTime = transferMsg->getEventTime();
	collisionTime = collisionMsg->getEventTime();
	wallCollisionTime = wallCollisionMsg->getEventTime();

	times.push_back(transferTime);
	times.push_back(collisionTime);
	times.push_back(wallCollisionTime);

	switch (mode) {
		
		case M_NNLIST:

			SphereMobility::outOfNeighborhoodTime(outOfNeighborhoodMsg, this);

			outOfNeighborhoodTime = outOfNeighborhoodMsg->getEventTime();

			times.push_back(outOfNeighborhoodTime);

// Initialize the smallestTime so it is not a negative value (NO_TIME)
			if (transferTime != NO_TIME) {
				smallestTime = transferTime;
			} else if (collisionTime != NO_TIME) {
				smallestTime = collisionTime;
			} else {
// No event will be scheduled for this sphere for now
				smallestTime = 0;
			}

			for (t = times.begin(); t != times.end(); ++t) {
				if (0 < (*t) && (*t) < smallestTime) smallestTime = (*t);
			}

			if (smallestTime == transferTime) {
				scheduleAt(transferTime, transferMsg);
			} else if (smallestTime == collisionTime) {
				scheduleAt(collisionTime, collisionMsg);
			} else if (smallestTime == wallCollisionTime) {
				scheduleAt(wallCollisionTime, wallCollisionMsg);
			} else if (smallestTime == outOfNeighborhoodTime) {
				scheduleAt(outOfNeighborhoodTime, outOfNeighborhoodMsg);
			} else {
// Nothing can be scheduled for this sphere
			}

			break;

		case M_CELLLIST:
		default:

// Initialize the smallestTime so it is not a negative value (NO_TIME)
			if (transferTime != NO_TIME) {
				smallestTime = transferTime;
			} else if (collisionTime != NO_TIME) {
				smallestTime = collisionTime;
			} else {
// No event will be scheduled for this sphere for now
				smallestTime = 0;
			}

			for (t = times.begin(); t != times.end(); ++t) {
				if (0 < (*t) && (*t) < smallestTime) smallestTime = (*t);
			}

			if (smallestTime == transferTime) {
				scheduleAt(transferTime, transferMsg);
			} else if (smallestTime == collisionTime) {
				scheduleAt(collisionTime, collisionMsg);
			} else if (smallestTime == wallCollisionTime) {
				scheduleAt(wallCollisionTime, wallCollisionMsg);
			} else {
// Nothing can be scheduled for this sphere
			}

			break;
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
	double transferTime;
	double collisionTime;
	double wallCollisionTime;
	double outOfNeighborhoodTime;

	double partnerCollisionTime;

	double smallestTime;

	vector<double> times;
	vector<double>::const_iterator t;

	MobilityMessage *msg;

// Step 3. Get the next space cell transfer time
	if (transferMsg->isScheduled()) cancelEvent(transferMsg);
	if (collisionMsg->isScheduled()) cancelEvent(collisionMsg);
	if (wallCollisionMsg->isScheduled()) cancelEvent(wallCollisionMsg);

	SphereMobility::transferTime(transferMsg, this);
	SphereMobility::collisionTime(collisionMsg, this);
	SphereMobility::wallCollisionTime(wallCollisionMsg, this);

	transferTime = transferMsg->getEventTime();
	collisionTime = collisionMsg->getEventTime();
	wallCollisionTime = wallCollisionMsg->getEventTime();

	times.push_back(transferTime);
	times.push_back(collisionTime);
	times.push_back(wallCollisionTime);

	switch (mode) {

		case M_NNLIST:

			if (outOfNeighborhoodMsg->isScheduled()) cancelEvent(outOfNeighborhoodMsg);

			SphereMobility::outOfNeighborhoodTime(outOfNeighborhoodMsg, this);

			outOfNeighborhoodTime = outOfNeighborhoodMsg->getEventTime();

			times.push_back(outOfNeighborhoodTime);

// Initialize the smallestTime so it is not a negative value (NO_TIME)
			if (collisionTime != NO_TIME) {
				smallestTime = collisionTime;
			} else if (wallCollisionTime != NO_TIME) {
				smallestTime = wallCollisionTime;
			} else if (transferTime != NO_TIME) {
				smallestTime = transferTime;
			} else if (outOfNeighborhoodTime != NO_TIME) {
				smallestTime = outOfNeighborhoodTime;
			} else {
// No event will be scheduled for this sphere for now
				smallestTime = 0;
			}

			for (t = times.begin(); t != times.end(); ++t) {
				if (0 < (*t) && (*t) < smallestTime) smallestTime = (*t);
			}

			if (smallestTime == outOfNeighborhoodTime) {

				scheduleAt(outOfNeighborhoodTime, outOfNeighborhoodMsg);

			} else if (smallestTime == transferTime) {

				scheduleAt(transferTime, transferMsg);

			} else if (smallestTime == wallCollisionTime) {

				scheduleAt(wallCollisionTime, wallCollisionMsg);

			} else if (smallestTime == collisionTime) {

				partnerCollisionTime = collisionMsg->getPartner()->scheduledCollisionTime();
				
				if (partnerCollisionTime == NO_TIME) {
// Partner particle does not have any scheduled collision event
					scheduleAt(collisionTime, collisionMsg);
					collisionMsg->getPartner()->updateNearNeighborList();
					collisionMsg->getPartner()->nextEventTime();

				} else {

					if (collisionTime < partnerCollisionTime) {

						scheduleAt(collisionTime, collisionMsg);
						collisionMsg->getPartner()->updateNearNeighborList();
						collisionMsg->getPartner()->nextEventTime();

					} else {
// Must check that the partner is solving the collision (its next event is not of type EV_CHECK)
						msg = ((Sphere *) collisionMsg->getPartner())->getScheduledMobilityMessage();

						if (msg->getKind() == EV_CHECK) {
// Partner expects that we are solving the collision
							collisionMsg->setKind(EV_COLLISION);
							scheduleAt(collisionTime, collisionMsg);
						} else {
// Partner is solving the collision event
							collisionMsg->setKind(EV_CHECK);
							scheduleAt(collisionTime, collisionMsg);

						}

					}

				}

			} else {
// Do nothing
			}

			break;

		case M_CELLLIST:
		default:

// Initialize the smallestTime so it is not a negative value (NO_TIME)
			if (transferTime != NO_TIME) {
				smallestTime = transferTime;
			} else if (collisionTime != NO_TIME) {
				smallestTime = collisionTime;
			} else if (wallCollisionTime != NO_TIME) {
				smallestTime = wallCollisionTime;
			} else {
// No event will be scheduled for this sphere for now
				smallestTime = 0;
			}

			for (t = times.begin(); t != times.end(); ++t) {
				if (0 < (*t) && (*t) < smallestTime) smallestTime = (*t);
			}

			if (smallestTime == transferTime) {

				scheduleAt(transferTime, transferMsg);

			} else if (smallestTime == wallCollisionTime) {

				scheduleAt(wallCollisionTime, wallCollisionMsg);

			} else if (smallestTime == collisionTime) {
				
				partnerCollisionTime = collisionMsg->getPartner()->scheduledCollisionTime();
				
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
						msg = ((Sphere *) collisionMsg->getPartner())->getScheduledMobilityMessage();

						if (msg->getKind() == EV_CHECK) {
// Partner expects that we are solving the collision
							collisionMsg->setKind(EV_COLLISION);
							scheduleAt(collisionTime, collisionMsg);

						} else {
// Partner is solving the collision event
							collisionMsg->setKind(EV_CHECK);
							scheduleAt(collisionTime, collisionMsg);

						}

					}

				}

			} else {
// Do nothing
			}

			break;
	}

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

	switch (mode) {

		case M_NNLIST:

// Step 2. Handle the event
			if (kind == EV_TRANSFER) {
				// Update the molecule space cell
				updateStateAfterTransfer((MobilityMessage *)msg);
				updateNearNeighborList();
				
				nextEventTime();

			} else if (kind == EV_WALLCOLLISION) {
				// Update the molecule data
				updateStateAfterWallCollision((MobilityMessage *)msg);
				updateNearNeighborList();
				
				nextEventTime();

			} else if (kind == EV_COLLISION) {

				updateStateAfterCollision((MobilityMessage *)msg);
				updateNearNeighborList();
				
				nextEventTime();

			} else if (kind == EV_OUTOFNEIGHBORHOOD) {

				updateStateAfterTransfer((MobilityMessage *)msg);
				updateNearNeighborList();
				
				nextEventTime();

			} else if (kind == EV_CHECK) {

				updateNearNeighborList();
				nextEventTime();

			} else {

			}

			break;

		case M_CELLLIST:
		default:

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

				nextEventTime();

			} else {

			}

			break;
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

// Refresh partner Near-NeighborList
	p->updateNearNeighborList();

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

/*
 * Create and populate the Near-Neighbor list
 */
void Sphere::createNearNeighborList() {

	int a, b, c;
	int i, j, k, n;
	int Nx, Ny, Nz, N; // Number of space cells (or divisions) in each axis

	double dx, dy, dz;
	double lrs; // listRadiusSquared

	double sTime;

	std::list<Particle *> particles;
	std::list<Particle *>::const_iterator pa;

	Nx = manager->getNumberOfSpaceCellsX();
	Ny = manager->getNumberOfSpaceCellsY();
	Nz = manager->getNumberOfSpaceCellsZ();

	n = getSpaceCell();

	sTime = simTime().dbl();

// i, j and k are the indexes of the space cell for each axis
	i = n / (Nz*Ny);
	j = (n % (Nz*Ny)) / Nz;
	k = (n % (Nz*Ny)) % Nz;

	for (a = -1; a <= 1; a++)
	for (b = -1; b <= 1; b++)
	for (c = -1; c <= 1; c++) {

// The neighbor cell must be contained in the simulation space
		if (CELLBELONGSTOSIMSPACE(i+a, j+b, k+c, Nx, Ny, Nz)) {

			N = (i+a)*Ny*Nz + (j+b)*Nz + (k+c);
			particles = manager->getSpaceCellParticles(N);

			for (pa = particles.begin(); pa != particles.end(); ++pa) {

				if ((*pa) == this) continue;

// Two particles are said to be neighbor when the sum of their listRadius is 
// greater than the distance between their centroids.
				dx = getX() + getVx()*(sTime - lastCollisionTime) - 
					((*pa)->getX() + (*pa)->getVx()*(sTime - (*pa)->getLastCollisionTime()));
				dy = getY() + getVy()*(sTime - lastCollisionTime) - 
					((*pa)->getY() + (*pa)->getVy()*(sTime - (*pa)->getLastCollisionTime()));
				dz = getZ() + getVz()*(sTime - lastCollisionTime) - 
					((*pa)->getZ() + (*pa)->getVz()*(sTime - (*pa)->getLastCollisionTime()));

				lrs = listRadius + (*pa)->getListRadius();
				lrs *= lrs;

				if (dx*dx+dy*dy+dz*dz < lrs) {
					this->neighborParticles.push_back(*pa);
				}

			}
		}
	}

}

/*
 * Add a particle to our Verlet list
 */
void Sphere::updateNearNeighborList() {

// Empty previous list
	neighborParticles.clear();

// Look for the particles closer than listRadius
	this->createNearNeighborList();

}

/*
 * Returns the scheduled mobility message
 * 
 * @return {MobilityMessage *}
 */
MobilityMessage *Sphere::getScheduledMobilityMessage() {

	MobilityMessage *msg;

	switch (mode) {

		case M_NNLIST:

			if (collisionMsg->isScheduled()) {
				msg = collisionMsg;
			} else if (wallCollisionMsg->isScheduled()) {
				msg = wallCollisionMsg;
			} else if (outOfNeighborhoodMsg->isScheduled()) {
				msg = outOfNeighborhoodMsg;
			} else {
				msg = transferMsg;
			}

			break;

		case M_CELLLIST:
		default:

			if (collisionMsg->isScheduled()) {
				msg = collisionMsg;
			} else if (wallCollisionMsg->isScheduled()) {
				msg = wallCollisionMsg;
			} else {
				msg = transferMsg;
			}

			break;
	}

	return msg;
	
}
