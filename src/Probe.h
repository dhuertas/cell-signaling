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
#include "base/Defines.h"

class Probe : public cSimpleModule {

	private:

		Manager *manager;

		// Probe name
		// Also used to name the cOutVector
		std::string name;

		int type;

		point_t position;

		double radius;

		// Update the cOutVectors periodically.
		double statsRefreshRate;

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
		std::string getName(void) { return name; };

		point_t getPosition(void) { return position; };

		double getRadius(void) { return radius; };

		int getType(void) { return type; };

		void setName(std::string n) { name = n; };

		void setPosition(point_t p) { position = p; };

		void setRadius(double r) { radius = r; };

		void setType(int t) { type = t; };

		void setManager(std::string);

};

#endif
