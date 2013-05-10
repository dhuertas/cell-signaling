cplusplus {{
	#include "../base/Particle.h"
	typedef Particle *ParticlePtr;
}}

class noncobject ParticlePtr;

message MobilityMessage {

	int eventType;
	int prevSpaceCell;
	int nextSpaceCell;

	double eventTime;

	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;

	ParticlePtr partner;

}