#ifndef CELL_H
#define CELL_H

#include <string.h>
#include <iostream>

#include "Manager.h"
#include "base/Sphere.h"

class Cell : public Sphere {

	private:

		int identifier;

		double emitEvery;   // Dumb var to test simulation
		int emitCount;

		Manager *manager;

	protected:

	public:

		~Cell();

		int getIdentifier(void) { return identifier; };
		void setIdentifer(int id) { identifier = id; };

		// cSimpleModule inheritance
		virtual void initialize(int);
		virtual int numInitStages() const;
		virtual void handleMessage(cMessage *);
		virtual void finish();

};

#endif
