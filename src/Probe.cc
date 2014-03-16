#include "Probe.h"

using namespace std;

Define_Module(Probe);

Probe::~Probe() {}

/*
 *
 */
void Probe::initialize(int stage) {
	
}

/*
 *
 */
int Probe::numInitStages() const {

	return 3;

}

/*
 *
 */
void Probe::handleMessage(cMessage *msg) {

}

/*
 *
 */
void Probe::finish() {

}