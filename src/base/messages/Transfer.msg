cplusplus {{
	#include "Manager.h"
	#include "Particle.h"

	typedef Manager *ManagerPtr;
}}

class noncobject ManagerPtr;

message TransferMessage {

	int prevSpaceCell;
	int nextSpaceCell;

	double transferTime;

	ManagerPtr manager;

}