//  This file is part of the Cell-Signaling project. Cell-Signaling is an
//  Omnet++ project to simulate cell signaling communications.
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

#ifndef MOLECULERECEIVER_H
#define MOLECULERECEIVER_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <coutvector.h>
#include <omnetpp.h>

#include "Receiver.h"
#include "SimpleCell.h"
#include "Manager.h"

// Forward declarations
class SimpleCell;

class MoleculeReceiver : public Receiver, public cSimpleModule {

	private:

		Manager *manager;

		SimpleCell *mobility;

		// Update the cOutVectors periodically.
		double statsRefreshRate;

	protected:
	
		unsigned long long received;

		cOutVector particlesReceivedVector; 

	public:

		MoleculeReceiver();

		~MoleculeReceiver();

		void setManager(std::string);

		void registerReception(int particleType);

		//
		// cSimpleModule inheritance
		//
		virtual void initialize(int);

		virtual int numInitStages() const;

		virtual void handleMessage(cMessage *);

		virtual void finish();
};

#endif
