cplusplus {{
	#include "Manager.h"

	typedef Manager *ManagerPtr;
}}

class noncobject ManagerPtr;

message TimeToLiveMessage {

	int particleId;
	double timeToLive;

	ManagerPtr manager;

}