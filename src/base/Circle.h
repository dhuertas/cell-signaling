#ifndef CIRCLE_H
#define CIRCLE_H

#include "Particle.h"

class Circle : public Particle {

	private:

	protected:

        double radius;

	public:

		Circle() : Particle () {};

		Circle(double, double, double, double, double, double);

		inline double getRadius(void) { return radius; };

		void setRadius(double r) { radius = r; };

};

#endif
