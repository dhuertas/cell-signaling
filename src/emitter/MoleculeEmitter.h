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

#ifndef MOLECULEEMITTER_H
#define MOLECULEEMITTER_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <coutvector.h>
#include <omnetpp.h>
#include <random.h>

#include <list>
#include <vector>

#include "Emitter.h"
#include "SimpleCell.h"
#include "Manager.h"

class MoleculeEmitter : public Emitter, public cSimpleModule {

	private:

		Manager *manager;

		SimpleCell *mobility;

		bool checkOverlap(point_t, double);

	protected:

		double emissionStart;

		double emissionDuration;

		double emissionRate;

		double emissionParticleRadius;

		double emissionParticleMass;

		double emissionTimeToLive;

		uint8_t emissionFunction;

		int emissionBoundariesMode;

		double emissionVelocity;

		double emissionListRadius;

		double emissionRefreshListRadius;

		double emissionDiffusion;

	public:

		MoleculeEmitter();
		~MoleculeEmitter();

		virtual void initialize(int);
		virtual int numInitStages() const;
		virtual void handleMessage(cMessage *);
		virtual void finish();

		Molecule * createMolecule();

		void setManager(std::string);
};

#endif
