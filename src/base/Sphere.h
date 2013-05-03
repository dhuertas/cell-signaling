#ifndef SPHERE_H
#define SPHERE_H

#include <csimplemodule.h>
#include <cmessage.h>
#include <cqueue.h>
#include <omnetpp.h>

#include "Circle.h"

class Sphere : public Circle, public cSimpleModule {

	private:

	protected:

	public:

		Sphere() : Circle() {};
		Sphere(double, double, double, double , double, double, double, double);

		// tk Environment related methods
		void tkEnvDrawShape(void);
		void tkEnvUpdatePosition(void);
		void tkEnvUpdatePosition(double);
};

#endif
