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

#ifndef MOLECULE_H
#define MOLECULE_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <coutvector.h>
#include <omnetpp.h>
#include <string.h>

#include <iostream>

#include "base/Sphere.h"
#include "receiver/MoleculeReceiver.h"
#include "messages/TimeToLive_m.h"

class Molecule : public Sphere {

	private:

		// Amount of time that will pass before it expires
		double timeToLive;

		// Update the cOutVectors periodically
		double statsRefreshRate;

		// Self messages
		TimeToLiveMessage *timeToLiveMsg;

	protected:

		// Vectors to track particle position over time
		cOutVector *xPositionVector;
		cOutVector *yPositionVector;
		cOutVector *zPositionVector;

	public:

		~Molecule();

		void expire();

		void scheduleExpire(double);

		bool isSignaling(cMessage *msg);

		void handleSignaling(cMessage *msg);

		//
		// cSimpleModule inheritance
		//
		virtual void initialize(int);

		virtual int numInitStages(void) const;

		virtual void handleMessage(cMessage *);

		virtual void finish();

};

#endif
