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
// the transfer time. Defined in src/base/Common.h file.
uint8_t Mobility::sides_[] = { 
//0, 1, 2, 3, 4, 5 <- dice side
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
double Mobility::nextTransferTime(TransferMessage *msg, Particle *p) {

  index3_t *idx = p->getSpaceCellIdx();
  index3_t nextIdx;

  point3_t *pos = p->getPosition();
  vector3_t *vel = p->getVelocity();

  point3_t cPos;

  point3_t P, Q, R;
  vector3_t V, W, N;

  Manager *manager = NULL;

  double currentTime = simTime().dbl();
  double transferTime = NO_TIME;

  uint8_t side, sides, count;
  uint8_t hits;

  double spaceCellSize = 0;
  double temp = 0;

  memset(&nextIdx, 0, sizeof(index3_t));

  manager = msg->getManager();

  sides = hits = count = 0;

  spaceCellSize = manager->getMaxSpaceSize();

  spaceCellSize /= (1 << idx->depth);

  double diff = currentTime - p->getLastCollisionTime();

  cPos.x = pos->x + vel->x*diff;
  cPos.y = pos->y + vel->y*diff;
  cPos.z = pos->z + vel->z*diff;

  // In a space cell we have 6 possible sides but since we know the direction 
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

    P.x = (idx->i + sides_[0*6 + side])*spaceCellSize;
    P.y = (idx->j + sides_[1*6 + side])*spaceCellSize;
    P.z = (idx->k + sides_[2*6 + side])*spaceCellSize;

    Q.x = (idx->i + sides_[3*6 + side])*spaceCellSize;
    Q.y = (idx->j + sides_[4*6 + side])*spaceCellSize;
    Q.z = (idx->k + sides_[5*6 + side])*spaceCellSize;

    R.x = (idx->i + sides_[6*6 + side])*spaceCellSize;
    R.y = (idx->j + sides_[7*6 + side])*spaceCellSize;
    R.z = (idx->k + sides_[8*6 + side])*spaceCellSize;

    V.x = Q.x - P.x; V.y = Q.y - P.y; V.z = Q.z - P.z;
    W.x = R.x - P.x; W.y = R.y - P.y; W.z = R.z - P.z;

    // Cross product to find the Normal vector
    N.x = V.y * W.z - V.z * W.y;
    N.y = V.z * W.x - V.x * W.z;
    N.z = V.x * W.y - V.y * W.x;

    // The plane equation is: N.x*(x-P.x)+N.y*(y-P.y)+N.z*(z-P.z) = 0
    // Replace by: x = xi + vx*t, y = yi + vy*t, z = zi + vz*t and 
    // solve for t.
    temp = (N.x*vel->x + N.y*vel->y + N.z*vel->z);

    if (temp != 0) {

      temp = (N.x*(P.x-cPos.x) + N.y*(P.y-cPos.y) + N.z*(P.z-cPos.z))/temp;

      // Solution found. "temp" is the amount of time the centroid of the
      // particle takes to go from its current position to the space cell 
      // side where it is bounded.
      if (temp == 0) {
        // The centroid of the particle is on a space cell side and
        // needs to update its space cell value and recalculate the
        // transfer time.
        transferTime = 0;
        hits |= 0x01 << side;
        break;
      } else if (count == 0 && 0 < temp) {
        transferTime = temp;
        hits |= 0x01 << side;
        count++;
      } else if (count > 0 && 0 < temp && temp < transferTime) {
        transferTime = temp;
        hits = 0x01 << side;
      } else if (count > 0 && 0 < temp && temp == transferTime) {
        hits |= 0x01 << side;
      } else {
        // We don't want it
      }
    }

    side++;
    sides >>= 1;
  }

  if (transferTime > 0) {

    transferTime += currentTime;

  } else if (transferTime == 0) {
    // If transferTime equals the simuation time sTime means that temp 
    // equals 0, thus the centroid of the particle belongs to the plane it 
    // is crossing. Set an event transfer at the same simulation time so it
    // will update the NextSpaceCell value and compute again the next 
    // transfer time.
    transferTime = currentTime;
  } else {
    // transfer time not found (NO_TIME)
  }

  nextIdx = *idx;

  if (hits & 0x01) nextIdx.k = idx->k+1;
  if (hits & 0x02) nextIdx.i = idx->i+1;
  if (hits & 0x04) nextIdx.j = idx->j-1;
  if (hits & 0x08) nextIdx.j = idx->j+1;
  if (hits & 0x10) nextIdx.i = idx->i-1;
  if (hits & 0x20) nextIdx.k = idx->k-1;

  msg->setNextSpaceCell(nextIdx);
  msg->setPrevSpaceCell(*idx);

  msg->setTransferTime(transferTime);

  return transferTime;

}

/*
 *
 */
void Mobility::handleTransfer(TransferMessage *msg, Particle *p) {

  Manager *manager = p->getManager();
  index3_t prev, next;

  prev = msg->getPrevSpaceCell();
  next = msg->getNextSpaceCell();

  manager->transferParticle(p, &prev, &next);

  p->setSpaceCellIdx(next);
}


/*
 * Set default values for the collision message.
 *
 * @param {CollisionMessage *} msg
 */
void Mobility::resetCollisionMessage(CollisionMessage *msg) {

  msg->setCollisionTime(NO_TIME);

  point3_t pos;
  vector3_t vel;

  memset(&pos, 0, sizeof(point3_t));
  memset(&vel, 0, sizeof(vector3_t));

  msg->setPosition(pos);
  msg->setVelocity(vel);

  msg->setPartner(NULL);
  msg->setPrevPartner(NULL);

  msg->setHits(0);
}
