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

#include "ParticleDistributionHelper.h"

using namespace std;
/*
 * Place each particle at a random position.
 *
 * @param {vector3_t} spaceSize: the simulation space (x, y and z)
 * @param {std::list<Particle *> *} particles
 */
void uniformDistribution(vector3_t spaceSize, std::list<Particle *> *particles) {

  bool overlap;

  point3_t pos;
  vector3_t range;

  std::list<Particle *>::iterator p, q;

  p = particles->begin();

  while (p != particles->end()) {

    overlap = false;
    
    range.x = spaceSize.x - 2*(*p)->getRadius();
    range.y = spaceSize.y - 2*(*p)->getRadius();
    range.z = spaceSize.z - 2*(*p)->getRadius();

    // We are using the Omnet internal random number generator dblrand()
    pos.x = (*p)->getRadius() + range.x*dblrand();
    pos.y = (*p)->getRadius() + range.y*dblrand();
    pos.z = (*p)->getRadius() + range.z*dblrand();

    for (q = particles->begin(); q != p; ++q) {

      if (checkOverlap(&pos, (*p)->getRadius(), (*q)->getPosition(), (*q)->getRadius())) {
        overlap = true;
        break;
      }

    }

    if ( ! overlap) {
      (*p)->setPosition(pos);
      ++p;
    }

  }

}

void uniformDistribution2(vector3_t spaceSize, std::list<Particle *> *particles) {

  bool overlap;

  int Nx, Ny, Nz;
  int i, j, k;
  int a, b, c;

  int count;

  double maxRadius;

  point3_t pos;
  vector3_t range;

  std::list<Particle *>::iterator p, q;
  std::vector<std::list<Particle*> > spaceCellLists;
  std::vector<int> spaceCells;
  std::vector<int>::iterator sci;

  count = 0;
  maxRadius = 0;

  for (p = particles->begin(); p != particles->end(); ++p) {
    if (maxRadius < (*p)->getRadius()) maxRadius = (*p)->getRadius();
  }

  Nx = ceil(spaceSize.x/2/maxRadius);
  Ny = ceil(spaceSize.y/2/maxRadius);
  Nz = ceil(spaceSize.z/2/maxRadius);

  spaceCellLists.resize(Nx*Ny*Nz);

  p = particles->begin();

  while (p != particles->end()) {

    overlap = false;

    range.x = spaceSize.x - 2*(*p)->getRadius();
    range.y = spaceSize.y - 2*(*p)->getRadius();
    range.z = spaceSize.z - 2*(*p)->getRadius();

    // We are using the Omnet internal random number generator dblrand()
    pos.x = (*p)->getRadius() + range.x*dblrand();
    pos.y = (*p)->getRadius() + range.y*dblrand();
    pos.z = (*p)->getRadius() + range.z*dblrand();

    i = floor(pos.x/2/maxRadius);
    j = floor(pos.y/2/maxRadius);
    k = floor(pos.z/2/maxRadius);

    for (a = -1; a <= 1; a++)
    for (b = -1; b <= 1; b++)
    for (c = -1; c <= 1; c++) {
      if (CELLBELONGSTOSIMSPACE(i+a, j+b, k+c, Nx, Ny, Nz)) {
        spaceCells.push_back((i+a)*Ny*Nz + (j+b)*Nz + k+c);
      }
    }
    
    for (sci = spaceCells.begin(); sci != spaceCells.end(); ++sci) {
      for (q = spaceCellLists.at(*sci).begin(); q != spaceCellLists.at(*sci).end(); ++q) {
        if (checkOverlap(&pos, (*p)->getRadius(), (*q)->getPosition(), (*q)->getRadius())) {
          overlap = true;
          break;
        }
      }

      if (overlap) break;
    }

    if (overlap) {
      continue;
    } else {
        spaceCellLists.at(i*Ny*Nz + j*Nz + k).push_back(*p);
      (*p)->setPosition(pos);
      ++p;
      count++;
      std::cout << "Particles placed: " << count;
      std::cout << ". Remaining: " << (particles->size()-count) << std::endl;
    }
  }

}

void uniformDistribution3(vector3_t spaceSize, std::list<Particle *> *particles) {

  unsigned int Nx, Ny, Nz;
  unsigned int i, j, k;
  unsigned int idx;
  int count;
  unsigned int available;
  int cell;
  double maxRadius;
  double d;

  point3_t pos;
  vector3_t offset;
  std::vector<int> cellIndexes;

  std::list<Particle *>::iterator p;

  count = 0;
  maxRadius = 0;
  d = 1.001;

  for (p = particles->begin(); p != particles->end(); ++p) {
    if (maxRadius < (*p)->getRadius()) maxRadius = (*p)->getRadius();
  }

  Nx = floor(spaceSize.x/(2*maxRadius*d));
  Ny = floor(spaceSize.y/(2*maxRadius*d));
  Nz = floor(spaceSize.z/(2*maxRadius*d));

  std::cout << "Space cells: " << Nx << " " << Ny << " " << Nz << std::endl;

  offset.x = spaceSize.x - Nx*(2*maxRadius*d);
  offset.y = spaceSize.y - Ny*(2*maxRadius*d);
  offset.z = spaceSize.z - Nz*(2*maxRadius*d);

  available = Nx*Ny*Nz;

  if (available < particles->size()) {
    std::cout << "Unable to place particles." << std::endl;
    exit(-1);
  }

  for (i = 0; i < available; i++) {
    cellIndexes.push_back(i);
  }

  p = particles->begin();

  while (p != particles->end()) {

    cell = intrand(cellIndexes.size());
    available--;

    idx = cellIndexes.at(cell);
    cellIndexes.erase(cellIndexes.begin() + cell);

    i = idx/Nz/Ny;
    j = (idx%(Nz*Ny))/Nz;
    k = (idx%(Nz*Ny))%Nz;

    pos.x = offset.x/2 + maxRadius*d*(1 + 2*i);
    pos.y = offset.y/2 + maxRadius*d*(1 + 2*j);
    pos.z = offset.z/2 + maxRadius*d*(1 + 2*k);

    std::cout << "Particle " << count << " position:";
    std::cout << " x=" << pos.x;
    std::cout << " y=" << pos.y;
    std::cout << " z=" << pos.z << std::endl;

    (*p)->setPosition(pos);
    ++p;

    count++;
  }

}

/*
 * Place each particle following a cube pattern.
 *
 * @param {vector3_t} spaceSize: the simulation space (x, y and z)
 * @param {std::list<Particle *> *} particles
 */
void cubeDistribution(vector3_t spaceSize, std::list<Particle *> *particles) {

  int m;
  int i, j, k;

  std::list<Particle *>::iterator p;

  m = (int)round(pow(particles->size(), 1/3.0));

  i = 0;
  j = 0;
  k = 0;

  point3_t pos;

  for (p = particles->begin(); p != particles->end(); ++p) {

    pos.x = (0.5 + i)*spaceSize.x/m;
    pos.y = (0.5 + i)*spaceSize.y/m;
    pos.z = (0.5 + i)*spaceSize.z/m;

    (*p)->setPosition(pos);

    k++;

    if (k%m == 0) {
      k = 0; j++;
      if (j%m == 0) {
        j = 0; i++;
        if (i%m == 0) {
          i = 0;
        }
      }
    }
  }

}

/* 
 * Place each particle following a sphere surface. The particles are uniformly 
 * distributed over a sphere surface. More info here:
 *
 * - http://mathworld.wolfram.com/SpherePointPicking.html
 *
 * @param {vector3_t} spaceSize: the simulation space (x, y and z)
 * @param {std::list<Particle *> *} particles
 * @param {point3_t} c: the center of the sphere
 * @param {double} radius: the sphere radius containing the center of the particles
 */
void sphereDistribution(vector3_t spaceSize, std::list<Particle *> *particles, point3_t c, double radius) {

  bool overlap;

  double u, v;
  double theta, phi;

  double minx, miny, minz;
  double maxRadius;

  point3_t pos;

  std::list<Particle *>::iterator p, q;

  u = 0;
  v = 0;

  maxRadius = 0;

  if (c.x == 0) c.x = spaceSize.x/2;
  if (c.y == 0) c.y = spaceSize.y/2;
  if (c.z == 0) c.z = spaceSize.z/2;

  if (radius == 0) {
    // Use the maximum radius possible given the space size
    minx = std::min(spaceSize.x - c.x, c.x - 0);
    miny = std::min(spaceSize.y - c.y, c.y - 0);
    minz = std::min(spaceSize.z - c.z, c.z - 0);

    radius = std::min(std::min(minx, miny), minz);

    // Find the max radius from the particles and subtract half of it (and
    // a little more).
    for (p = particles->begin(); p != particles->end(); ++p) {
      maxRadius = std::max((*p)->getRadius(), maxRadius);
    }

    radius -= maxRadius*(1.0 + 0.0001);

  } else if (radius < 0) {
    radius = -radius;
  }

  p = particles->begin();

  // WARNING: if the space size is not big enough this function may hang the
  // simulation.
  while (p != particles->end()) {

    overlap = false;

    u = 0;
    v = 0;

    // dblrand() returns values in [0,1) while u and v must be in (0, 1)
    while (u == 0) u = dblrand();
    while (v == 0) v = dblrand();

    theta = 2*M_PI*u;
    phi = acos(2*v - 1);

    pos.x = c.x + radius*sin(theta)*cos(phi);
    pos.y = c.y + radius*sin(theta)*sin(phi);
    pos.z = c.z + radius*cos(theta);

    for (q = particles->begin(); q != p; ++q) {

      if (checkOverlap(&pos, (*p)->getRadius(), (*q)->getPosition(), (*q)->getRadius())) {
        overlap = true;
        break;
      }

    }

    if ( ! overlap) {
      (*p)->setPosition(pos);
      ++p;
    }

  }

}

/*
 *
 */
void sphereEquallyDistributed(vector3_t spaceSize, std::list<Particle *> *particles, point3_t c, double radius) {

  int m, n;
  int Ncount;
  int M_theta;
  int M_phi;
  int N;

  double minx, miny, minz;
  double maxRadius;

  double a, d;
  double theta, phi;
  double d_theta;
  double d_phi;

  point3_t pos;

  std::list<Particle *>::iterator p, q;

  if (c.x == 0) c.x = spaceSize.x/2;
  if (c.y == 0) c.y = spaceSize.y/2;
  if (c.z == 0) c.z = spaceSize.z/2;

  if (radius == 0) {
    // Use the maximum radius possible given the space size.
    minx = std::min(spaceSize.x - c.x, c.x - 0);
    miny = std::min(spaceSize.y - c.y, c.y - 0);
    minz = std::min(spaceSize.z - c.z, c.z - 0);

    radius = std::min(std::min(minx, miny), minz);

    // Find the max radius from the particles and subtract half of it (and
    // a little more).
    for (p = particles->begin(); p != particles->end(); ++p) {
      maxRadius = std::max((*p)->getRadius(), maxRadius);
    }

    radius -= maxRadius*(1.0 + 0.0001);

  } else if (radius < 0) {
    radius = -radius;
  }

  N = particles->size();
  Ncount = 0;
  a = 4*M_PI*radius*radius/N;
  d = sqrt(a);
  M_theta = round(M_PI/d);
  d_theta = M_PI/M_theta;
  d_phi = a/d_theta;

  p = particles->begin();

  for (m = 0; m < M_theta; m++) {

    theta = M_PI*(m + 0.5)/M_theta;
    M_phi = round(2*M_PI*sin(theta)/d_phi);

    for (n = 0; n < M_phi; n++) {
      phi = 2*M_PI*n/M_phi;

      pos.x = radius*sin(theta)*cos(phi);
      pos.y = radius*sin(theta)*sin(phi);
      pos.z = radius*cos(theta);

      (*p)->setPosition(pos);
      ++p;

      Ncount++;
    }
  }

  // Remove the particles that could not made it to the surface-
  for (q = particles->end(); q != p; --q) {
    particles->erase(q);
  }

}

/*
 * Put all the particles one near the other at a certain point.
 *
 * @param {vector3_t} spaceSize: the simulation space (x, y and z)
 * @param {std::list<Particle *> *} particles
 * @param {point3_t} c: center
 */
void highDensityDistribution(vector3_t spaceSize, std::list<Particle *> *particles, point3_t c) {

  bool overlap;

  uint8_t overlapCount;
  double movingVariance;

  point3_t pos;

  std::list<Particle *>::iterator p, q;

  overlapCount = 100;

  p = particles->begin();

  while (p != particles->end()) {

    overlap = false;

    pos.x = 0;
    pos.y = 0;
    pos.z = 0;

    while (pos.x - (*p)->getRadius() <= 0 || pos.x + (*p)->getRadius() >= spaceSize.x ||
      pos.y - (*p)->getRadius() <= 0 || pos.y + (*p)->getRadius() >= spaceSize.y ||
      pos.z - (*p)->getRadius() <= 0 || pos.z + (*p)->getRadius() >= spaceSize.z) {
      
      pos.x = dblRandNormal(c.x, movingVariance);
      pos.y = dblRandNormal(c.y, movingVariance);
      pos.z = dblRandNormal(c.z, movingVariance);
    }

    for (q = particles->begin(); q != p; ++q) {

      if (checkOverlap(&pos, (*p)->getRadius(), (*q)->getPosition(), (*q)->getRadius())) {
        overlap = true;
        break;
      }

    }

    if ( ! overlap) {

      (*p)->setPosition(pos);
      ++p;

    } else {

      overlapCount--;
      if (overlapCount == 0) {
        overlapCount = 100;
        movingVariance++;
      }

    }

  }

}

void densepacked(vector3_t spaceSize, std::list<Particle *> *particles, point3_t c) {

  //int cx, cy;
  int plane, row;
  unsigned int count;
  //double srt;
  double r, R;
  
  point3_t pos;
  std::list<Particle *>::iterator p;

  //cx = 0; cy = 0;
  //srt = sqrt(3);
  count = 0;
  r = 0;

  for (p = particles->begin(); p != particles->end(); ++p) {
    if (r < (*p)->getRadius()) r = (*p)->getRadius();
  }

  r = r*1.01;
  // R = pow(particles->size()*1.098/0.7404804897, 1.0/3.0)*r; // N = 1000
  R = pow(particles->size()*1.065/0.7404804897, 1.0/3.0)*r; // N = 10000

  pos.x = r; pos.y = r; pos.z = r;

  p = particles->begin();

  plane = 0;
  row = 0;

  pos.z = r;

  while (pos.z + r < spaceSize.z) {

    if (plane % 2 == 0) { // A plane
      pos.y = r;
    } else {
      pos.y = r + sqrt(3.0)*r/3.0; // B plane
    }

    while (pos.y + r < spaceSize.y) {

      if (plane % 2 == 0) {
        pos.x = row % 2 == 0 ? 2*r : r;
      } else {
        pos.x = row % 2 == 0 ? r : 2*r;
      }

      while (pos.x + r < spaceSize.x) {

        if ((c.x-pos.x)*(c.x-pos.x) + (c.y-pos.y)*(c.y-pos.y) + 
          (c.z-pos.z)*(c.z-pos.z) < R*R) {

          if (count < particles->size()) {
            std::cout << "Particle " << count << " position: ";
            std::cout << " x=" << pos.x;
            std::cout << " y=" << pos.y;
            std::cout << " z=" << pos.z << std::endl;
            (*p)->setPosition(pos);
            ++p;
            count++;
          }

        }

        pos.x += 2*r;

      }

      row++;

      pos.y += sqrt(3.0)*r;

    }

    // Finished a plane, add z
    pos.z += 2.0*sqrt(6.0)*r/3.0;

    plane++;
    row = 0;

  }

}

/*
 * Detect whether two sphere particles are overlaping or not.
 *
 * @param {point3_t} ca: the center coordinates of the first particle
 * @param {double} ra: the radius of the first particle
 * @param {point3_t} cb: the center coordinates of the second particle
 * @param {double} rb: the radius of the second particle
 * @return {bool}
 */
bool checkOverlap(point3_t *ca, double ra, point3_t *cb, double rb) {

  double dx = (ca->x-cb->x);
  double dy = (ca->y-cb->y);
  double dz = (ca->z-cb->z);

  return  sqrt(dx*dx + dy*dy + dz*dz) < ra + rb ? true : false;

}

/*
 * Returns a random number following a Normal distribution N(mean, var)
 * 
 * @param {double} mean
 * @param {double} vari
 */
double dblRandNormal(double mean, double var) {

  double x1, x2, w, y1;
  // double y2;

  do {
    x1 = 2.0*dblrand() - 1.0;
    x2 = 2.0*dblrand() - 1.0;
    w = x1*x1 + x2*x2;
  } while (w >= 1.0);

  w = sqrt((-2.0*log(w))/w);
  y1 = x1*w;
  //y2 = x2*w;

  return mean + y1*sqrt(var);
}
