#include "Sphere.h"

using namespace std;

Sphere::Sphere(
	double x, 
	double y, 
	double z, 
	double vx, 
	double vy, 
	double vz, 
	double radius, 
	double mass) : Circle(x, y, vx, vy, radius, mass) {

	position.z = z;
	velocity.z = vz;

}

void Sphere::tkEnvDrawShape() {

	std::stringstream buffer;

	// We will use the shape drawing tool to draw a circle around the particle
	// center instead of using
	buffer << 2*getRadius();
	getDisplayString().setTagArg("b", 0, buffer.str().c_str());
	getDisplayString().setTagArg("b", 1, buffer.str().c_str());

	getDisplayString().setTagArg("b", 2, "oval");
	getDisplayString().setTagArg("b", 3, "white");
	getDisplayString().setTagArg("b", 4, "black");
	getDisplayString().setTagArg("b", 5, 1);
}

// Update the module position in the tk environment. This method is used when
// the particle is initialized.
void Sphere::tkEnvUpdatePosition() {

	std::stringstream buffer;

	// Set position string for tkenv
	buffer << getY();

	getDisplayString().setTagArg("p", 0, buffer.str().c_str());
	buffer.str(std::string()); // clear buffer

	buffer << getX();

	getDisplayString().setTagArg("p", 1, buffer.str().c_str());
	buffer.str(std::string());

}

// Update the module position in the tk environment. This method is used 
// during the simulation.
void Sphere::tkEnvUpdatePosition(double t) {

	std::stringstream buffer;

	// Set position string for tkenv
	buffer << getY() + getVy()*(t-getLastCollisionTime());
	EV << "y: " << buffer.str().c_str() << "\n";

	getDisplayString().setTagArg("p", 0, buffer.str().c_str());
	buffer.str(std::string()); // clear buffer

	buffer << getX() + getVx()*(t-getLastCollisionTime());
	EV << "x: " << buffer.str().c_str() << "\n";
	getDisplayString().setTagArg("p", 1, buffer.str().c_str());
	buffer.str(std::string());

}
