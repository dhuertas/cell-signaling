#include "Sphere.h"
#include "mobility/SphereMobility.h"
#include "../Molecule.h"
#include "../SimpleCell.h"
#include "../receiver/MoleculeReceiver.h"
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
 * Initialize the self messages. This method is called by the manager module
 * during network initialization.
 */
void Sphere::initMobilityMessages() {
	// Methods called from other modules must have this macro
	Enter_Method_Silent();

	transferMsg = new TransferMessage("mobility", EV_TRANSFER);
	transferMsg->setManager(manager);

	collisionMsg = new CollisionMessage("mobility", EV_NONE);
	collisionMsg->setManager(manager);

	if (manager->getMode() == M_NNLIST) {
		outOfNeighborhoodMsg = new OutOfNeighborhoodMessage("mobility", EV_OUTOFNEIGHBORHOOD);	
	}

	if (manager->getDeltaTime() > 0) {
		brownianMotionMsg = new BrownianMotionMessage("mobility", EV_BROWNIAN);
		brownianMotionMsg->setManager(manager);
	}

	SphereMobility::resetCollisionMessage(collisionMsg);

}

/*
 *
 */
void Sphere::deleteMobilityMessages() {

	cancelAndDelete(transferMsg);
	cancelAndDelete(collisionMsg);

	transferMsg = NULL;
	collisionMsg = NULL;

	if (manager->getMode() == M_NNLIST) {
		cancelAndDelete(outOfNeighborhoodMsg);
		outOfNeighborhoodMsg = NULL;
	}

	if (manager->getDeltaTime() > 0) {
		cancelAndDelete(brownianMotionMsg);
		brownianMotionMsg = NULL;
	}

}

/*
 * Initialize the event queue by computing the first event for the sphere. This
 * method is called by the manager in order to initialize the event queue.
 */
void Sphere::initializeMobility() {
	// Methods called from other modules must have this macro
	Enter_Method_Silent();

	double transferTime = NO_TIME;
	double collisionTime = NO_TIME;
	double boundaryCollisionTime = NO_TIME;
	double brownianMotionTime = NO_TIME;
	double scheduledCollisionTime = NO_TIME;

	double minTime = NO_TIME;

	vector<double> times;
	vector<double>::const_iterator t;

	// Compute the first collision and the first transfer
	transferTime = SphereMobility::nextTransfer(transferMsg, this);

	if (transferTime != NO_TIME) {
		scheduleAt(transferTime, transferMsg);
	}

	if (manager->getMode() == M_NNLIST) {
		this->handleOutOfNeighborhood();
	}

	collisionTime = SphereMobility::nextCollision(collisionMsg, 0, this);
	boundaryCollisionTime = SphereMobility::nextBoundaryCollision(collisionMsg, this);
	brownianMotionTime = SphereMobility::brownianMotion(brownianMotionMsg, this);

	if (collisionMsg->isScheduled()) {
		scheduledCollisionTime = collisionMsg->getCollisionTime();
		times.push_back(scheduledCollisionTime);
	}

	if (collisionTime != NO_TIME) {
		times.push_back(collisionTime);
	}

	if (boundaryCollisionTime != NO_TIME) {
		times.push_back(boundaryCollisionTime);
	}

	if (brownianMotionTime != NO_TIME) {
		times.push_back(brownianMotionTime);
	}

	minTime = boundaryCollisionTime;

	for (t = times.begin(); t!= times.end(); ++t) {
		if ((*t) < minTime) minTime = (*t);
	}

	if (minTime == collisionTime && collisionTime != NO_TIME) {

		if (collisionMsg->isScheduled()) {
			((Sphere *)collisionMsg->getPrevPartner())->getCollisionMessage()->setKind(EV_CHECK);
			((Sphere *)collisionMsg->getPrevPartner())->getCollisionMessage()->setPartner(NULL);
		}

		cancelEvent(collisionMsg);

		collisionMsg->setKind(EV_COLLISION);
		collisionMsg->setCollisionTime(collisionTime);

		scheduleAt(collisionTime, collisionMsg);

		((Sphere *)collisionMsg->getPartner())->adjustCollision(collisionTime, this);

	} else if (minTime == boundaryCollisionTime && boundaryCollisionTime != NO_TIME) {

		if (collisionMsg->isScheduled()) {
			((Sphere *)collisionMsg->getPrevPartner())->getCollisionMessage()->setKind(EV_CHECK);
			((Sphere *)collisionMsg->getPrevPartner())->getCollisionMessage()->setPartner(NULL);
		}

		cancelEvent(collisionMsg);

		collisionMsg->setKind(EV_BOUNDARYCOLLISION);
		collisionMsg->setCollisionTime(boundaryCollisionTime);

		scheduleAt(boundaryCollisionTime, collisionMsg);

	} else if (minTime == brownianMotionTime) {

		if (collisionMsg->isScheduled()) {
			((Sphere *)collisionMsg->getPrevPartner())->getCollisionMessage()->setKind(EV_CHECK);
			((Sphere *)collisionMsg->getPrevPartner())->getCollisionMessage()->setPartner(NULL);
		}

		cancelEvent(collisionMsg);

		scheduleAt(brownianMotionTime, brownianMotionMsg);

	} else if (minTime == scheduledCollisionTime) {
		// Leave it as it is scheduled
	} else {
		// minTime equals NO_TIME;
	}

}

/*
 * This method gets called when the sphere must leave the simulation space
 * (either it has expired, crossed a boundary, etc).
 */
void Sphere::finishMobility() {

	if (collisionMsg->isScheduled()) {

		cancelEvent(collisionMsg);

		// Change the event type of the third party sphere to EV_CHECK
		if (collisionMsg->getPartner() != NULL) {
			((Sphere *)collisionMsg->getPartner())->getCollisionMessage()
				->setKind(EV_CHECK);
			((Sphere *)collisionMsg->getPartner())->getCollisionMessage()
				->setPartner(NULL);
		}

	}

	if (transferMsg->isScheduled()) {
		cancelEvent(transferMsg);
	}

	if (outOfNeighborhoodMsg != NULL) {
		if (outOfNeighborhoodMsg->isScheduled()) {
			cancelEvent(outOfNeighborhoodMsg);
		}
	}

	if (brownianMotionMsg != NULL) {
		if (brownianMotionMsg->isScheduled()) {
			cancelEvent(brownianMotionMsg);
		}
	}

	active = false;

}

/*
 * Method overwrite. Finish mobility when called from another sphere (it has
 * been absorved by a receiver, combined with another molecule, etc).
 *
 * @param {Particle *} from
 */
void Sphere::finishMobility(Particle *from) {
	// Methods called from other modules must have this macro
	Enter_Method_Silent();

	if (collisionMsg->isScheduled()) {

		cancelEvent(collisionMsg);

		// Change the event type of the third party sphere to EV_CHECK
		if (collisionMsg->getPartner() != NULL && 
		from->getParticleId() != collisionMsg->getPartner()->getParticleId()) {

			((Sphere *)collisionMsg->getPartner())->getCollisionMessage()
				->setKind(EV_CHECK);
			((Sphere *)collisionMsg->getPartner())->getCollisionMessage()
				->setPartner(NULL);
		}

	}

	if (transferMsg->isScheduled()) {
		cancelEvent(transferMsg);
	}

	if (outOfNeighborhoodMsg != NULL) {
		if (outOfNeighborhoodMsg->isScheduled()) {
			cancelEvent(outOfNeighborhoodMsg);
		}
	}

	if (brownianMotionMsg != NULL) {
		if (brownianMotionMsg->isScheduled()) {
			cancelEvent(brownianMotionMsg);
		}
	}

	active = false;

}

/*
 * This method gets called when the sphere is a partner in a collision event.
 *
 * @param {double} newTime: the collision event time
 * @param {Particle *} from: the particle who is handling the collision
 */
void Sphere::adjustCollision(double newTime, Particle *from) {

	// Methods called from other modules must have this macro
	Enter_Method_Silent();

	Particle* partner = NULL;
	CollisionMessage* partnerCollisionMsg = NULL;

	if (collisionMsg->isScheduled()) cancelEvent(collisionMsg);

	partner = dynamic_cast<Particle*>(collisionMsg->getPartner());

	// Change the event type of the third party sphere to EV_CHECK
	if (partner != NULL && from->getParticleId() != partner->getParticleId()) {

		partnerCollisionMsg = dynamic_cast<CollisionMessage*>(((Sphere *)partner)->getCollisionMessage());

		if (partnerCollisionMsg != NULL && 
			partnerCollisionMsg->getManager() == manager) {

			partnerCollisionMsg->setKind(EV_CHECK);
			partnerCollisionMsg->setPartner(NULL);
		}

	}

	SphereMobility::resetCollisionMessage(collisionMsg);

	collisionMsg->setKind(EV_CHECK);
	collisionMsg->setCollisionTime(newTime);
	collisionMsg->setPartner(from);

	scheduleAt(newTime, collisionMsg);

}

/*
 * Handles the mobility of the sphere
 *
 * @param {cMessage *} msg
 */
void Sphere::handleMobilityMessage(cMessage *msg) {

	double transferTime = NO_TIME;
	double collisionTime = NO_TIME;
	double boundaryCollisionTime = NO_TIME;
	double brownianMotionTime = NO_TIME;

	// Step 1. Find the next event in the queue.

	int kind = msg->getKind();

	// Step 2. Handle the event.

	if (kind == EV_TRANSFER) {

		this->handleTransfer((TransferMessage *)msg);

	} else if (kind == EV_OUTOFNEIGHBORHOOD) {

		// The list will be updated anyway ... 
		// TODO remove this "else if" section?

	} else if (kind == EV_COLLISION) {

		this->handleCollision((CollisionMessage *)msg);
		SphereMobility::resetCollisionMessage(collisionMsg);

	} else if (kind == EV_BOUNDARYCOLLISION) {

		this->handleBoundaryCollision((CollisionMessage *)msg);
		SphereMobility::resetCollisionMessage(collisionMsg);

	} else if (kind == EV_BROWNIAN) {

		// Update position, velocity and lastCollisionTime
		this->handleBrownianMotion((BrownianMotionMessage *)msg);

		//  Compute the next velocity with brownian motion
		brownianMotionTime = SphereMobility::brownianMotion(brownianMotionMsg, this);
		scheduleAt(brownianMotionTime, brownianMotionMsg);

		// Collisions must be rescheduled
		if (collisionMsg->isScheduled()) cancelEvent(collisionMsg);
		SphereMobility::resetCollisionMessage(collisionMsg);

	} else if (kind == EV_CHECK) {

		SphereMobility::resetCollisionMessage(collisionMsg);

	}

	if (manager->getMode() == M_NNLIST) {
		this->handleOutOfNeighborhood();
	}

	// Step 3. Compute the next transfer time for the particle corresponding to the
	// event.

	if (transferMsg->isScheduled()) cancelEvent(transferMsg);

	transferTime = SphereMobility::nextTransfer(transferMsg, this);

	if (transferTime != NO_TIME) {
		scheduleAt(transferTime, transferMsg);
	}

	// Step 4. Compute the next collision time with particles in appropriate 
	// neighboring cells.

	collisionTime = SphereMobility::nextCollision(collisionMsg, kind, this);
	boundaryCollisionTime = SphereMobility::nextBoundaryCollision(collisionMsg, this);

	// Step 5. Adjust the position of the event and its new partner’s event in the 
	// event queue. Since a wall collision changes the path of a particle, we only 
	// keep either a particle collision or a wall collision for each particle.

	if (collisionMsg->isScheduled()) {

		if (collisionTime < collisionMsg->getCollisionTime() && collisionTime != NO_TIME) {

			if (collisionMsg->getPrevPartner() != NULL) {
				((Sphere *)collisionMsg->getPrevPartner())->getCollisionMessage()->setKind(EV_CHECK);
				((Sphere *)collisionMsg->getPrevPartner())->getCollisionMessage()->setPartner(NULL);
			}

			cancelEvent(collisionMsg);

			collisionMsg->setKind(EV_COLLISION);
			collisionMsg->setCollisionTime(collisionTime);

			scheduleAt(collisionTime, collisionMsg);

			if (collisionMsg->getPartner() != NULL) {
				((Sphere *)collisionMsg->getPartner())->adjustCollision(collisionTime, this);
			}

		}

	} else {

		if (collisionTime < boundaryCollisionTime && collisionTime != NO_TIME) {

			collisionMsg->setKind(EV_COLLISION);
			collisionMsg->setCollisionTime(collisionTime);

			scheduleAt(collisionTime, collisionMsg);

			if (collisionMsg->getPartner() != NULL) {
				((Sphere *)collisionMsg->getPartner())->adjustCollision(collisionTime, this);
			}

		} else {

			collisionMsg->setKind(EV_BOUNDARYCOLLISION);
			collisionMsg->setCollisionTime(boundaryCollisionTime);

			collisionMsg->setPartner(NULL);
			collisionMsg->setPrevPartner(NULL);

			if (boundaryCollisionTime != NO_TIME) {
				scheduleAt(boundaryCollisionTime, collisionMsg);
			}

		}

	}

	// Step 6. Return to Step 1.

}

/*
 * Tells the manager to update the space cell position.
 *
 * @param {TransferMessage *} msg
 */
void Sphere::handleTransfer(TransferMessage *msg) {

	manager->transferParticle(this, 
		msg->getPrevSpaceCell(), 
		msg->getNextSpaceCell());

	setPrevSpaceCell(msg->getPrevSpaceCell());
	setSpaceCell(msg->getNextSpaceCell());

	// Statistics
	manager->registerTransfer();

}

/*
 * Updates the particle position, velocity and the last collision time values
 * and does the same for the partner sphere.
 *
 * @param {CollisionMessage *} msg
 */
void Sphere::handleCollision(CollisionMessage *msg) {

	double tc, m1, m2, tmp;
	double v1n, v1e1, v1e2, v2n, v2e1, v2e2;

	point_t *ppos = NULL; // partner position pointer
	vect_t *pvel = NULL; // partner velocity pointer

	point_t c1, c2;
	vect_t v1, n, e1, e2;

	Particle *p;

	p = msg->getPartner();
	tc = msg->getCollisionTime();

	ppos = p->getPosition();
	pvel = p->getVelocity();

	// Find the center position of the spheres
	c1.x = position.x + velocity.x*(tc - lastCollisionTime);
	c1.y = position.y + velocity.y*(tc - lastCollisionTime);
	c1.z = position.z + velocity.z*(tc - lastCollisionTime);

	c2.x = ppos->x + pvel->x*(tc - p->getLastCollisionTime());
	c2.y = ppos->y + pvel->y*(tc - p->getLastCollisionTime());
	c2.z = ppos->z + pvel->z*(tc - p->getLastCollisionTime());

	m1 = this->getMass();
	m2 = p->getMass();

	// Change frame of reference of the system to one of the spheres
	v1.x = velocity.x - pvel->x;
	v1.y = velocity.y - pvel->y;
	v1.z = velocity.z - pvel->z;

	// Find the normal vector of the plane of collision
	n.x = c2.x - c1.x;
	n.y = c2.y - c1.y;
	n.z = c2.z - c1.z;

	tmp = sqrt(n.x*n.x + n.y*n.y + n.z*n.z);

	n.x /= tmp;
	n.y /= tmp;
	n.z /= tmp;

	// Find e1 as the perpendicular vector to both n and v, and then e2 as the 
	// one perpendicular to n and e1
	e1.x = n.y*v1.z - n.z*v1.y;
	e1.y = n.z*v1.x - n.x*v1.z;
	e1.z = n.x*v1.y - n.y*v1.x;

	// Normalize the vector found, e1
	tmp = sqrt(e1.x*e1.x + e1.y*e1.y + e1.z*e1.z);

	if (tmp > 0) {
		e1.x /= tmp;
		e1.y /= tmp;
		e1.z /= tmp;
	}

	// Find the velocity vectors in the new basis and if ...
	v1n  = v1.x*n.x  + v1.y*n.y  + v1.z*n.z;

	if (e1.x == 0.0 && e1.y == 0.0 && e1.z == 0.0) {
		// n and v are parallel, we can solve directly
		tmp = (m1 - m2)*v1n/(m1 + m2);
		v2n = 2*m1*v1n/(m1 + m2);
		v1n = tmp;

		// Revert the frame of reference, the velocity vectors and set the new
		// velocity
		velocity.x = (v1n*n.x + pvel->x);
		velocity.y = (v1n*n.y + pvel->y);
		velocity.z = (v1n*n.z + pvel->z);

		p->setVx(v2n*n.x + pvel->x);
		p->setVy(v2n*n.y + pvel->y);
		p->setVz(v2n*n.z + pvel->z);

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

		// Find the new velocity in the normal component (remember that v2n 
		// initially is 0.0)
		tmp = (m1 - m2)*v1n/(m1 + m2);
		v2n = 2*m1*v1n/(m1 + m2);
		v1n = tmp;

		// Revert the frame of reference, the velocity vectors and set the new 
		// velocity
		velocity.x = (v1n*n.x + v1e1*e1.x + v1e2*e2.x + pvel->x);
		velocity.y = (v1n*n.y + v1e1*e1.y + v1e2*e2.y + pvel->y);
		velocity.z = (v1n*n.z + v1e1*e1.z + v1e2*e2.z + pvel->z);

		p->setVx(v2n*n.x + v2e1*e1.x + v2e2*e2.x + pvel->x);
		p->setVy(v2n*n.y + v2e1*e1.y + v2e2*e2.y + pvel->y);
		p->setVz(v2n*n.z + v2e1*e1.z + v2e2*e2.z + pvel->z);

	}

	// Update the particles position
	setPosition(c1);
	p->setPosition(c2);

	// Update the last collision times
	setLastCollisionTime(tc);
	p->setLastCollisionTime(tc);

	logCollisionTime(tc);
	((Sphere *)p)->logCollisionTime(tc);

	// Statistics
	manager->registerCollision();

}

void Sphere::handleBoundaryCollision(CollisionMessage *msg) {

	if (boundariesMode == BM_ELASTIC) {

		handleWallCollision(msg);

	} else if (boundariesMode == BM_EXPIRE) {

		expire();

	} else if (boundariesMode == BM_PERIODIC) {
		// TODO complete this
	}

}

/*
 * Updates the position, velocity and the last collision time when the particle
 * collides with a wall.
 *
 * @param {CollisionMessage *} msg
 */
void Sphere::handleWallCollision(CollisionMessage *msg) {

	position.x = msg->getX();
	position.y = msg->getY();
	position.z = msg->getZ();

	velocity.x = msg->getVx();
	velocity.y = msg->getVy();
	velocity.z = msg->getVz();

	lastCollisionTime = msg->getCollisionTime();

	// Statistics
	manager->registerWallCollision();
}

/*
 *
 */
void Sphere::handleBrownianMotion(BrownianMotionMessage *msg) {

	double dt = manager->getDeltaTime();

	point_t temp = position;

	position.x = velocity.x*dt + temp.x;
	position.y = velocity.y*dt + temp.y;
	position.z = velocity.z*dt + temp.z;

	// New velocity for the next brownian motion step
	velocity.x = msg->getVx();
	velocity.y = msg->getVy();
	velocity.z = msg->getVz();

	lastCollisionTime = msg->getBrownianMotionTime();

}

/*
 * Handler for the out-of-neighborhood event
 */
void Sphere::handleOutOfNeighborhood() {

	double outOfNeighborhoodTime;

	outOfNeighborhoodTime = NO_TIME;

	if (outOfNeighborhoodMsg->isScheduled()) {
		cancelEvent(outOfNeighborhoodMsg);
	}

	this->updateNearNeighborList();

	outOfNeighborhoodTime = SphereMobility::outOfNeighborhoodTime(outOfNeighborhoodMsg, this);

	outOfNeighborhoodMsg->setOutOfNeighborhoodTime(outOfNeighborhoodTime);

	if (outOfNeighborhoodTime != NO_TIME) {
		scheduleAt(outOfNeighborhoodTime, outOfNeighborhoodMsg);
	}

}

/*
 * Create and populate the Near-Neighbor list
 */
void Sphere::createNearNeighborList() {

	int a, b, c;
	int i, j, k, n;
	// Number of space cells (or divisions) in each axis
	int *Nx, *Ny, *Nz, N;

	double dt, dx, dy, dz;
	double lrs; // listRadiusSquared

	double sTime;

	std::list<Particle *> *particles;
	std::list<Particle *>::const_iterator pa;

	point_t *ppos = NULL;
	vect_t *pvel = NULL;

	Nx = manager->getNumberOfSpaceCellsX();
	Ny = manager->getNumberOfSpaceCellsY();
	Nz = manager->getNumberOfSpaceCellsZ();

	n = getSpaceCell();

	sTime = simTime().dbl();

	// i, j and k are the indexes of the space cell for each axis
	i =  n/((*Nz)*(*Ny));
	j = (n%((*Nz)*(*Ny)))/(*Nz);
	k = (n%((*Nz)*(*Ny)))%(*Nz);

	for (a = -1; a <= 1; a++)
	for (b = -1; b <= 1; b++)
	for (c = -1; c <= 1; c++) {

		// The neighbor cell must be contained in the simulation space
		if (CELLBELONGSTOSIMSPACE(i+a, j+b, k+c, *Nx, *Ny, *Nz)) {

			N = (i+a)*(*Ny)*(*Nz)+(j+b)*(*Nz)+(k+c);
			particles = manager->getSpaceCellParticles(N);

			for (pa = particles->begin(); pa != particles->end(); ++pa) {

				if ((*pa) == this) continue;

				dt = sTime - (*pa)->getLastCollisionTime();
				// Two particles are said to be neighbor when the sum of their listRadius is 
				// greater than the distance between their centroids.
				ppos = (*pa)->getPosition();
				pvel = (*pa)->getVelocity();

				dx = position.x + velocity.x*(sTime - lastCollisionTime) - 
						(ppos->x + pvel->x*dt);
				dy = position.y + velocity.y*(sTime - lastCollisionTime) - 
						(ppos->y + pvel->y*dt);
				dz = position.z + velocity.z*(sTime - lastCollisionTime) - 
						(ppos->z + pvel->z*dt);

				lrs = listRadius + (*pa)->getListRadius();

				if (dx*dx + dy*dy + dz*dz < lrs*lrs) {
						neighborParticles.push_back(*pa);
				}

			}
		}
	}
}

/*
 * Add nearby spheres to the list.
 *
 * TODO need to find a way to remove as few particles as possible
 * instead of rebuilding the neighbor list every time.
 */
void Sphere::updateNearNeighborList() {

	// Empty previous list
	neighborParticles.clear();

	// Look for the particles closer than listRadius
	this->createNearNeighborList();
}

/*
 * Returns the scheduled transfer mobility message.
 * 
 * @return {TransferMessage *}
 */
TransferMessage * Sphere::getTransferMessage() {

	if (active && transferMsg != NULL) {
		return this->transferMsg;
	}

	return NULL;
}

/*
 * Returns the scheduled collision mobility message.
 * 
 * @return {CollisionMessage *}
 */
CollisionMessage * Sphere::getCollisionMessage() {

	if (active && this->collisionMsg != NULL) {
		return this->collisionMsg;
	}

	return NULL;
}

/*
 * Draws the shape of the cell in the tk environment.
 */
void Sphere::tkEnvDrawShape() {

	std::stringstream buffer;
	cModule *parent = getParentModule();

	// We will use the shape drawing tool to draw a circle around the particle
	// center instead of using
	buffer << 2*radius;

	getDisplayString().setTagArg("b", 0, buffer.str().c_str());
	getDisplayString().setTagArg("b", 1, buffer.str().c_str());

	getDisplayString().setTagArg("b", 2, "oval");
	getDisplayString().setTagArg("b", 3, "white");
	getDisplayString().setTagArg("b", 4, "black");
	getDisplayString().setTagArg("b", 5, 1);

	if (parent != NULL && strcmp(parent->getName(), "cell") == 0) {
		parent->getDisplayString().setTagArg("b", 0, buffer.str().c_str());
		parent->getDisplayString().setTagArg("b", 1, buffer.str().c_str());

		parent->getDisplayString().setTagArg("b", 2, "oval");
		parent->getDisplayString().setTagArg("b", 3, "white");
		parent->getDisplayString().setTagArg("b", 4, "black");
		parent->getDisplayString().setTagArg("b", 5, 1);
	}
}

/*
 * Update the module position in the tk environment. This method is used when
 * the particle is initialized.
 */
void Sphere::tkEnvUpdatePosition() {

	std::stringstream buffer;
	cModule *parent = getParentModule();

	buffer << position.y;

	// Set position string for tkenv
	if (strcmp(getName(), "molecule") == 0) {
		getDisplayString().setTagArg("p", 0, buffer.str().c_str());
	}

	if (parent != NULL && strcmp(parent->getName(), "cell") == 0) {
		parent->getDisplayString().setTagArg("p", 0, buffer.str().c_str());
	}

	buffer.str(std::string()); // clear buffer

	buffer << position.x;

	if (strcmp(getName(), "molecule") == 0) {
		getDisplayString().setTagArg("p", 1, buffer.str().c_str());
	}

	if (parent != NULL && strcmp(parent->getName(), "cell") == 0) {
		parent->getDisplayString().setTagArg("p", 1, buffer.str().c_str());
	}

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
	cModule *parent = getParentModule();

	double lc = getLastCollisionTime();

	// Set position string for tkenv
	buffer << position.y + velocity.y*(t-lc);

	if (strcmp(getName(), "molecule") == 0) {
		getDisplayString().setTagArg("p", 0, buffer.str().c_str());
	}

	// Also move the parent's shape
	if (parent != NULL && strcmp(parent->getName(), "cell") == 0) {
		parent->getDisplayString().setTagArg("p", 0, buffer.str().c_str());
	}

	buffer.str(std::string()); // clear buffer

	buffer << position.x + velocity.x*(t-lc);

	if (strcmp(getName(), "molecule") == 0) {
		getDisplayString().setTagArg("p", 1, buffer.str().c_str());
	}

	// Also move the parent's shape
	if (parent != NULL && strcmp(parent->getName(), "cell") == 0) {
		parent->getDisplayString().setTagArg("p", 1, buffer.str().c_str());
	}

	buffer.str(std::string());
}

/*
 * Sets the manager attribute.
 * 
 * @param {string} param: the name of the manager module
 */
void Sphere::setManager(std::string param) {

	try {
		manager = (Manager *)simulation.
			getSystemModule()->getSubmodule(param.c_str());
	} catch (cException *e) {
		EV << "setManager error" << "\n";
	}
}
