//  Omnet++ project to simulate cell signaling communications 
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

cplusplus {{
	#include "Manager.h"
	#include "Particle.h"

	typedef Manager *ManagerPtr;
	typedef Particle *ParticlePtr;
}}

class noncobject ManagerPtr;
class noncobject ParticlePtr;

message CollisionMessage {

	double collisionTime;

	double x; 	// x coordinate of the centroid after a wall collision occurs
	double y; 	// y coordinate of the centroid after a wall collision occurs
	double z; 	// z coordinate of the centroid after a wall collision occurs

	double vx; 	// vx component of the velocity after a wall collision occurs
	double vy;	// vy component of the velocity after a wall collision occurs
	double vz;	// vz component of the velocity after a wall collision occurs

	ManagerPtr manager;

	ParticlePtr partner;
	ParticlePtr prevPartner;

}