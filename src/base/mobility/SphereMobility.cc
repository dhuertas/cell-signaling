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

#include "SphereMobility.h"

using namespace std;

/*
 * Computes the collision time with the neighboring particles and stores the
 * smallest one.
 *
 * @param {CollisionMessage *} msg
 * @param {Sphere *} s
 * @return {double} the smallest computed collision time
 */
double SphereMobility::nextCollisionTime(CollisionMessage *msg, Sphere *s) {

  unsigned int collisionCounter = 0;

  double temp = NO_TIME;

  double prevCollisionTime = NO_TIME;
  double partnerCollisionTime = NO_TIME;
  double collisionTime = NO_TIME;

  // TODO Improve this
  double currentTime = simTime().dbl();

  Manager *manager = s->getManager();
  std::vector<Particle *> particles;

  index3_t *idx = s->getSpaceCellIdx();

  Sphere *partner = NULL;
  Sphere *prevPartner = NULL;

  CollisionMessage *partnerCollisionMsg = NULL;

  manager->getNeighborParticles(idx, &particles);

  prevCollisionTime = msg->getCollisionTime();

  if (prevCollisionTime != NO_TIME && prevCollisionTime > currentTime) {

    partner = (Sphere *)msg->getPartner();
    prevPartner = (Sphere *)msg->getPartner();

    collisionTime = prevCollisionTime;
    collisionCounter++;
  }

  std::vector<Particle *>::iterator it;
  for (it = particles.begin(); it != particles.end(); ++it) {

    Sphere *candidate = (Sphere *)(*it);
    if (candidate == s) continue;

    // Solve particle to particle collision
    temp = solveCollision(s, candidate);

    // Only keep the collision time if it's smaller than the rest of computed
    // collision times so far and is also smaller than the partner collision time.
    partnerCollisionMsg = candidate->getCollisionMessage();

    if (partnerCollisionMsg != NULL) {
      partnerCollisionTime = partnerCollisionMsg->getCollisionTime();
    } else {
      partnerCollisionTime = NO_TIME;
    }

    if (temp != NO_TIME && currentTime < temp) {
      // Collision found!
      if (
        (collisionCounter == 0 && (partnerCollisionTime > temp || partnerCollisionTime == NO_TIME)) ||
        (collisionCounter  > 0 && (partnerCollisionTime > temp || partnerCollisionTime == NO_TIME) && temp < collisionTime)) {
        
        collisionTime = temp;
        partner = candidate;
        collisionCounter++;
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
double SphereMobility::nextBoundaryCollisionTime(CollisionMessage *msg, Sphere *s) {

  double radius = s->getRadius();

  point3_t *pos = s->getPosition();
  vector3_t *vel = s->getVelocity();

  point3_t P, Q, R;
  vector3_t V, W, N;

  point3_t nextPos;
  vector3_t nextVel;

  Manager *manager = s->getManager();

  vector3_t *spaceSize = manager->getSpaceSize();

  double lastCollisionTime = s->getLastCollisionTime();

  unsigned int boundariesMode = s->getBoundariesMode();

  uint8_t side, sides, count;
  uint8_t hits;

  sides = hits = count = 0;

  double temp = 0;
  double collisionTime = 0;

  if (vel->x == 0 && vel->y == 0 && vel->z == 0) {
    return NO_TIME;
  }

  // In a cubic space we have 6 possible sides but since we know the direction 
  // of the particle we only need to check 3 at most.
  // sides order
  // --x----- side 5 0x20
  // ---x---- side 4 0x10
  // ----x--- side 3 0x08
  // -----x-- side 2 0x04
  // ------x- side 1 0x02
  // -------x side 0 0x01
  if (vel->x > 0) sides |= 0x02;
  if (vel->x < 0) sides |= 0x10;
  if (vel->y > 0) sides |= 0x08;
  if (vel->y < 0) sides |= 0x04;
  if (vel->z > 0) sides |= 0x01;
  if (vel->z < 0) sides |= 0x20;

  side = 0;

  while (sides > 0) {
    if ( ! (sides & 0x01)) {
      side++;
      sides >>= 1;
      continue;
    }

    P.x = sides_[0*6 + side]*spaceSize->x;
    P.y = sides_[1*6 + side]*spaceSize->y;
    P.z = sides_[2*6 + side]*spaceSize->z;

    Q.x = sides_[3*6 + side]*spaceSize->x;
    Q.y = sides_[4*6 + side]*spaceSize->y;
    Q.z = sides_[5*6 + side]*spaceSize->z;

    R.x = sides_[6*6 + side]*spaceSize->x;
    R.y = sides_[7*6 + side]*spaceSize->y;
    R.z = sides_[8*6 + side]*spaceSize->z;

    V.x = Q.x - P.x; V.y = Q.y - P.y; V.z = Q.z - P.z;
    W.x = R.x - P.x; W.y = R.y - P.y; W.z = R.z - P.z;

    // Cross product to find the Normal vector
    N.x = V.y * W.z - V.z * W.y;
    N.y = V.z * W.x - V.x * W.z;
    N.z = V.x * W.y - V.y * W.x;

    // This is the plane equation N.x*(x-P.x)+N.y*(y-P.y)+N.z*(z-P.z) = 0
    // Replace by:
    //     x = xi + vx*t + R*vnx, vnx = vx/sqrt(vx2+vy2+vz2), R = radius
    //     y = yi + vy*t + R*vny, vny = vy/sqrt(vx2+vy2+vz2), "
    //     z = zi + vz*t + R*vnz, vnz = vz/sqrt(vx2+vy2+vz2), "
    // and find for t.
    temp = (N.x*vel->x + N.y*vel->y + N.z*vel->z);

    if (temp != 0) {

      if (boundariesMode == BM_PERIODIC  || boundariesMode == BM_EXPIRE) {
        // Use sphere center
        temp = (N.x*(P.x - pos->x) + N.y*(P.y - pos->y) + N.z*(P.z - pos->z)) / temp;
      } else {
        // Use sphere shell
        switch (side) {
          case 0: temp = (N.x*(P.x - pos->x) + N.y*(P.y - pos->y) + N.z*(P.z - (pos->z + radius))) / temp; break;
          case 1: temp = (N.x*(P.x - (pos->x + radius)) + N.y*(P.y - pos->y) + N.z*(P.z - pos->z)) / temp; break;
          case 2: temp = (N.x*(P.x - pos->x) + N.y*(P.y - (pos->y - radius)) + N.z*(P.z - pos->z)) / temp; break;
          case 3: temp = (N.x*(P.x - pos->x) + N.y*(P.y - (pos->y + radius)) + N.z*(P.z - pos->z)) / temp; break;
          case 4: temp = (N.x*(P.x - (pos->x - radius)) + N.y*(P.y - pos->y) + N.z*(P.z - pos->z)) / temp; break;
          case 5: temp = (N.x*(P.x - pos->x) + N.y*(P.y - pos->y) + N.z*(P.z - (pos->z - radius))) / temp; break;
        }
      }

      if (count == 0) {
        collisionTime = temp;
        hits |= 0x01 << side;
        count++;
      } else if (count > 0 && temp < collisionTime) {
        collisionTime = temp;
        hits = 0x01 << side;
      } else if (count > 0 && temp == collisionTime) {
        hits |= 0x01 << side;
      } else {
        // It's greater
      }
    } else {
      // Line and plane do not intersect
    }

    side++;
    sides >>= 1;
  }

  // The future particle position
  nextPos.x = pos->x + vel->x*collisionTime;
  nextPos.y = pos->y + vel->y*collisionTime;
  nextPos.z = pos->z + vel->z*collisionTime;

  msg->setPosition(nextPos);

  double boundaryCollisionTime = collisionTime;
  boundaryCollisionTime += lastCollisionTime;

  // Compute the next velocity vector
  nextVel.x = (hits & 0x02 || hits & 0x10) ? -vel->x : vel->x;
  nextVel.y = (hits & 0x04 || hits & 0x08) ? -vel->y : vel->y;
  nextVel.z = (hits & 0x01 || hits & 0x20) ? -vel->z : vel->z;

  msg->setVelocity(nextVel);
  msg->setHits(hits);

  return boundaryCollisionTime;

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
double SphereMobility::brownianMotion(BrownianMotionMessage *msg, Sphere *p) {

  double nextDeltaTime = NO_TIME;
  double BMStdDev;
  double dt;

  point3_t *pos = NULL;
  point3_t nextPos;

  vector3_t vel;

  dt = msg->getManager()->getDeltaTime();

  if (dt > 0) {

    nextDeltaTime = simTime().dbl() + dt;

    pos = p->getPosition();
    BMStdDev = p->getBrownianMotionStdDev();

    nextPos.x = normal(pos->x, BMStdDev);
    nextPos.y = normal(pos->y, BMStdDev);
    nextPos.z = normal(pos->z, BMStdDev);

    vel.x = (nextPos.x - pos->x)/dt;
    vel.y = (nextPos.y - pos->y)/dt;
    vel.z = (nextPos.z - pos->z)/dt;

    msg->setVelocity(vel);

    msg->setBrownianMotionTime(nextDeltaTime);

  }

  return nextDeltaTime;
}

/*
 *
 */
void SphereMobility::handleBoundaryCollision(CollisionMessage *msg, Sphere *s) {

  unsigned int boundariesMode = s->getBoundariesMode();

  switch (boundariesMode) {
    case BM_EXPIRE:
      s->expire();
      break;
    case BM_PERIODIC:
      handlePeriodicBoundary(msg, s);
      break;
    default:
    case BM_ELASTIC:
      handleWallCollision(msg, s);
      break;
  }

}

/*
 * Updates the particle position, velocity and the last collision time values
 * and does the same for the partner sphere.
 *
 * @param {CollisionMessage *} msg
 * @param {Sphere *} s
 */
void SphereMobility::handleCollision(CollisionMessage *msg, Sphere *s) {

  unsigned int imageIdx = 0;

  double collisionTime;
  double m1, m2, tmp, itmp;
  double v1n, v1e1, v1e2, v2n;

  point3_t c1, c2;
  vector3_t v1, n, e1, e2;

  Manager *manager = s->getManager();

  Sphere *partner = (Sphere *)msg->getPartner();
  collisionTime = msg->getCollisionTime();

  point3_t *pos = s->getPosition();
  vector3_t *vel = s->getVelocity();

  point3_t *ppos = partner->getPosition();
  vector3_t *pvel = partner->getVelocity();

  vector3_t *spaceSize = manager->getSpaceSize();

  double lastCollisionTime = s->getLastCollisionTime();
  double partnerLastCollisionTime = partner->getLastCollisionTime();

  // Find the center position of the spheres
  c1.x = pos->x + vel->x*(collisionTime - lastCollisionTime);
  c1.y = pos->y + vel->y*(collisionTime - lastCollisionTime);
  c1.z = pos->z + vel->z*(collisionTime - lastCollisionTime);

  c2.x = ppos->x + pvel->x*(collisionTime - partnerLastCollisionTime);
  c2.y = ppos->y + pvel->y*(collisionTime - partnerLastCollisionTime);
  c2.z = ppos->z + pvel->z*(collisionTime - partnerLastCollisionTime);

  // Periodic boundary conditions
  imageIdx = partner->getImageIdx();

  // Image transformation
  if (imageIdx > 0) {

    point3_t c2p = c2;

    c2p.x += (imageIdx & 0x20) ? -spaceSize->x : (imageIdx & 0x10 ? spaceSize->x : 0);
    c2p.y += (imageIdx & 0x08) ? -spaceSize->y : (imageIdx & 0x04 ? spaceSize->y : 0);
    c2p.z += (imageIdx & 0x02) ? -spaceSize->z : (imageIdx & 0x01 ? spaceSize->z : 0);

    // Find the normal vector of the plane of collision
    n.x = c2p.x - c1.x;
    n.y = c2p.y - c1.y;
    n.z = c2p.z - c1.z;
  } else {

    // Find the normal vector of the plane of collision
    n.x = c2.x - c1.x;
    n.y = c2.y - c1.y;
    n.z = c2.z - c1.z;
  }

  m1 = s->getMass();
  m2 = partner->getMass();

  // Change frame of reference of the system to one of the spheres
  v1.x = vel->x - pvel->x;
  v1.y = vel->y - pvel->y;
  v1.z = vel->z - pvel->z;

  itmp = 1/sqrt(n.x*n.x + n.y*n.y + n.z*n.z);

  n.x *= itmp;
  n.y *= itmp;
  n.z *= itmp;

  // Find e1 as the perpendicular vector to both n and v, and then e2 as the 
  // one perpendicular to n and e1
  e1.x = n.y*v1.z - n.z*v1.y;
  e1.y = n.z*v1.x - n.x*v1.z;
  e1.z = n.x*v1.y - n.y*v1.x;

  // Find the velocity vectors in the new basis and if ...
  v1n  = v1.x*n.x  + v1.y*n.y  + v1.z*n.z;

  if (e1.x == 0.0 && e1.y == 0.0 && e1.z == 0.0) {
    // n and v are parallel, we can solve directly
    tmp = (m1 - m2)*v1n/(m1 + m2);
    v2n = 2*m1*v1n/(m1 + m2);
    v1n = tmp;

    // Revert the frame of reference, the velocity vectors and set the new
    // velocity
    vel->x = v1n*n.x + pvel->x;
    vel->y = v1n*n.y + pvel->y;
    vel->z = v1n*n.z + pvel->z;

    pvel->x += v2n*n.x;
    pvel->y += v2n*n.y;
    pvel->z += v2n*n.z;

  } else {

    // Normalize the vector found, e1
    itmp = 1/sqrt(e1.x*e1.x + e1.y*e1.y + e1.z*e1.z);

    if (itmp > 0) {
      e1.x *= itmp;
      e1.y *= itmp;
      e1.z *= itmp;
    }
  
    e2.x = e1.y*n.z - e1.z*n.y;
    e2.y = e1.z*n.x - e1.x*n.z;
    e2.z = e1.x*n.y - e1.y*n.x;

    // Find the rest of the components
    v1e1 = v1.x*e1.x + v1.y*e1.y + v1.z*e1.z;
    v1e2 = v1.x*e2.x + v1.y*e2.y + v1.z*e2.z;

    v2n  = 0.0;

    // Find the new velocity in the normal component (remember that v2n 
    // initially is 0.0)
    tmp = (m1 - m2)*v1n/(m1 + m2);
    v2n = 2*m1*v1n/(m1 + m2);
    v1n = tmp;

    // Revert the frame of reference, the velocity vectors and set the new 
    // velocity
    vel->x = (v1n*n.x + v1e1*e1.x + v1e2*e2.x + pvel->x);
    vel->y = (v1n*n.y + v1e1*e1.y + v1e2*e2.y + pvel->y);
    vel->z = (v1n*n.z + v1e1*e1.z + v1e2*e2.z + pvel->z);

    pvel->x += v2n*n.x; 
    pvel->y += v2n*n.y;
    pvel->z += v2n*n.z;
  }

  // Update the particles position
  s->setPosition(c1);
  partner->setPosition(c2);

  // Update the last collision times
  s->setLastCollisionTime(collisionTime);
  partner->setLastCollisionTime(collisionTime);

  // Statistics
  s->logCollisionTime(collisionTime);
  partner->logCollisionTime(collisionTime);

  manager->registerCollision();
}

/*
 * Updates the position, velocity and the last collision time when the particle
 * collides with a wall.
 *
 * @param {CollisionMessage *} msg
 */
void SphereMobility::handleWallCollision(CollisionMessage *msg, Sphere *s) {

  point3_t pos = msg->getPosition();
  s->setPosition(pos);

  vector3_t vel = msg->getVelocity();
  s->setVelocity(vel);

  double collisionTime = msg->getCollisionTime();
  s->setLastCollisionTime(collisionTime);

  // Statistics
  msg->getManager()->registerWallCollision();
}

/*
 *
 */
void SphereMobility::handlePeriodicBoundary(CollisionMessage *msg, Sphere *s) {

  // TODO There is a problem with typedef uint8_t Hits in Collision.msg
  // TODO trying with int...
  int hits = msg->getHits();

  // Velocity goes unchanged

  // Update last collision time
  double collisionTime = msg->getCollisionTime();
  s->setLastCollisionTime(collisionTime);

  // Move particle position to the other side and update index3_t
  Manager *manager = s->getManager();

  point3_t pos = msg->getPosition();
  index3_t *idx = s->getSpaceCellIdx();
  index3_t nextIdx;

  vector3_t *spaceSize = manager->getSpaceSize();

  double maxSpaceSize = manager->getMaxSpaceSize();
  double sideLength = (maxSpaceSize / (1 << idx->depth));

  nextIdx.i = idx->i;
  nextIdx.j = idx->j;
  nextIdx.k = idx->k;
  nextIdx.depth = idx->depth;

  if (hits & 0x01) { pos.z -= spaceSize->z; } // side 0
  if (hits & 0x02) { pos.x -= spaceSize->x; } // side 1
  if (hits & 0x04) { pos.y += spaceSize->y; } // side 2
  if (hits & 0x08) { pos.y -= spaceSize->y; } // side 3
  if (hits & 0x10) { pos.x += spaceSize->x; } // side 4
  if (hits & 0x20) { pos.z += spaceSize->z; } // side 5

  nextIdx.i = floor(pos.x / sideLength);
  nextIdx.j = floor(pos.y / sideLength);
  nextIdx.k = floor(pos.z / sideLength);

  s->setPosition(pos);
  manager->transferParticle(s, idx, &nextIdx);

  idx->i = nextIdx.i;
  idx->j = nextIdx.j;
  idx->k = nextIdx.k;
  idx->depth = nextIdx.depth;
}

/*
 *
 */
void SphereMobility::handleBrownianMotion(BrownianMotionMessage *msg, Sphere *s) {

  Manager *manager = msg->getManager();

  double dt = manager->getDeltaTime();

  point3_t *pos = s->getPosition();
  vector3_t *vel = s->getVelocity();

  pos->x += vel->x*dt;
  pos->y += vel->y*dt;
  pos->z += vel->z*dt;

  // New velocity for the next brownian motion step
  s->setVelocity(msg->getVelocity());

  s->setLastCollisionTime(msg->getBrownianMotionTime());

}

/*
 * Solves the Sphere-Sphere collision problem and returns the collision time
 * 
 * @param {Particle *} pa
 * @param {Particle *} pb
 * @return {double} the collision time
 */
double SphereMobility::solveCollision(Sphere *sa, Sphere *sb) {
// Distance between centers A and B when t = tc (time of collision):
  //                 ______________________________________________
  //           \    / ( Ax + Avx*(tc-ta) - (Bx + Bvx*(tc-tb) )² +
  //  D(A, B) = \  /  ( Ay + Avy*(tc-ta) - (By + Bvy*(tc-tb) )² +   = Ra + Rb
  //             \/   ( Az + Avz*(tc-ta) - (Bz + Bvz*(tc-tb) )²
  //
  // ta: time when the previous collision took place for particle A
  // tb: same for particle B

  // (dxi + dvx*tc)² + (dyi + dvy*tc)² + (dyi + dvy*tc)²= (A.r + B.r)²
  double ta = sa->getLastCollisionTime();
  double tb = sb->getLastCollisionTime();

  point3_t *posa = sa->getPosition();
  point3_t *posb = sb->getPosition();

  vector3_t *vela = sa->getVelocity();
  vector3_t *velb = sb->getVelocity();

  Manager *manager = sa->getManager();
  vector3_t *spaceSize = manager->getSpaceSize();

  double dxi, dyi, dzi;

  unsigned int imageIdx = sb->getImageIdx();

  // Image transformation
  if (imageIdx > 0) {

    point3_t temp = *posb;

    temp.x += (imageIdx & 0x20) ? -spaceSize->x : (imageIdx & 0x10 ? spaceSize->x : 0);
    temp.y += (imageIdx & 0x08) ? -spaceSize->y : (imageIdx & 0x04 ? spaceSize->y : 0);
    temp.z += (imageIdx & 0x02) ? -spaceSize->z : (imageIdx & 0x01 ? spaceSize->z : 0);
  
    dxi = posa->x - temp.x - vela->x*ta + velb->x*tb;
    dyi = posa->y - temp.y - vela->y*ta + velb->y*tb;
    dzi = posa->z - temp.z - vela->z*ta + velb->z*tb;
  } else {

    dxi = posa->x - posb->x - vela->x*ta + velb->x*tb;
    dyi = posa->y - posb->y - vela->y*ta + velb->y*tb;
    dzi = posa->z - posb->z - vela->z*ta + velb->z*tb;
  }

  double dvx = vela->x - velb->x;
  double dvy = vela->y - velb->y;
  double dvz = vela->z - velb->z;

  double radd = sa->getRadius() + sb->getRadius();

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