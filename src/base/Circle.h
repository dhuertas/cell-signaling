#ifndef CIRCLE_H
#define CIRCLE_H

#include "Particle.h"

class Circle : public Particle {

	private:

		double radius;

	protected:

	public:

		Circle() : Particle () {};
		Circle(double, double, double, double, double, double);
		double getRadius(void) { return radius; };
		void setRadius(double r) { radius = r; };

};

#endif
