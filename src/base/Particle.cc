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

#include "Particle.h"

using namespace std;

/*
 * Constructor.
 */
Particle::Particle() {}

/*
 * Constructor overload.
 */
Particle::Particle(
	double x,
	double y,
	double vx,
	double vy,
	double m) {

	position.x = x;
	position.y = y;
	position.z = 0;

	velocity.x = vx;
	velocity.y = vy;
	velocity.z = 0;

	spaceCell = -1;

	spaceCellIdx.im = 0x00;
	spaceCellIdx.i = 0;
	spaceCellIdx.j = 0;
	spaceCellIdx.k = 0;

	prevSpaceCellIdx.im = 0x00;
	prevSpaceCellIdx.i = 0;
	prevSpaceCellIdx.j = 0;
	prevSpaceCellIdx.k = 0;

	mass = m;
	lastCollisionTime = 0;

	listRadius = 1;

}