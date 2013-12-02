#include "MoleculeEmitter.h"

Define_Module(MoleculeEmitter);

MoleculeEmitter::MoleculeEmitter() : Emitter() {

}

MoleculeEmitter::~MoleculeEmitter() {

}

void MoleculeEmitter::initialize() {

	emissionStart = par("emissionStart");

	emissionDuration = par("emissionDuration");

	emissionRate = par("emissionRate");

}

void MoleculeEmitter::handleMessage(cMessage *msg) {

}

void MoleculeEmitter::finish() {

}