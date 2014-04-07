cplusplus {{
	#include "Manager.h"

	typedef Manager *ManagerPtr;
}}

class noncobject ManagerPtr;

message BrownianMotionMessage {

	double brownianMotionTime;

	double x;
	double y;
	double z;

	double vx;
	double vy;
	double vz;

	ManagerPtr manager;
}