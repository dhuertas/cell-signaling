cplusplus {{
	#include "../base/Manager.h"
	#include "../base/Particle.h"

	typedef Manager *ManagerPtr;
	typedef Particle *ParticlePtr;
}}

class noncobject ManagerPtr;
class noncobject ParticlePtr;

message MobilityMessage {

	int eventType;			// Not in use
	int prevSpaceCell;
	int nextSpaceCell;

	double eventTime;

	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;

	ManagerPtr manager;
	ParticlePtr partner;

}