cplusplus {{
	#include "../base/Manager.h"
	#include "../base/Particle.h"

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