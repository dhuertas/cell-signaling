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

#ifndef PROBE_H
#define PROBE_H

#include <cmessage.h>
#include <csimplemodule.h>
#include <cqueue.h>
#include <omnetpp.h>

#include <cmath>
#include <vector>
#include <list>

#include "Manager.h"
#include "Common.h"

class Probe : public cSimpleModule {

	private:

		Manager *manager_;

		// Probe name
		// Also used to name the cOutVector
		std::string name_;

		int type_;

		point3_t position_;

		double radius_;

		// Update the cOutVectors periodically.
		double statsRefreshRate_;

	protected:

		cOutVector moleculeDensityVector;

	public:

		~Probe();

		//
		// cSimpleModule inheritance
		//
		virtual void initialize(int stage);

		virtual int numInitStages() const;

		virtual void handleMessage(cMessage *);

		virtual void finish();

		double getMoleculeDensity(void);

		//
		// Gets and sets
		//
		std::string getName(void) { return name_; };

		point3_t getPosition(void) { return position_; };

		double getRadius(void) { return radius_; };

		int getType(void) { return type_; };

		void setName(std::string n) { name_ = n; };

		void setPosition(point3_t p) { position_ = p; };

		void setRadius(double r) { radius_ = r; };

		void setType(int t) { type_ = t; };

		void setManager(std::string);

};

#endif
