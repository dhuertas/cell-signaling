#ifndef PARTICLE_H
#define PARTICLE_H

#include "Defines.h"
#include <list>

class Particle {

	private:

	protected:

		point_t position;
		vect_t velocity;

		// Cell List attributes
		int spaceCell;

		// Verlet List attributes
		std::list<Particle*> neighbourParticles;
		double cutOffRadius;
		double listRadius; // listRadius > cutOffRadius

		double mass;
		double lastCollisionTime;

	public:

		Particle();
		Particle(double, double, double, double, double);

		// Verlet methods
		void createVerletList(std::list<Particle*> *);
		void updateVerletList(std::list<Particle*> *);
		void addParticleToVerletList(Particle *);

		// Gets and sets
		double getX(void) { return position.x; };
		double getY(void) { return position.y; };
		double getZ(void) { return position.z; };
		point_t getPosition(void) { return position; };
		double getVx(void) { return velocity.x; };
		double getVy(void) { return velocity.y; };
		double getVz(void) { return velocity.z; };
		vect_t getVelocity(void) { return velocity; };
		double getMass(void) { return mass; };
		double getSpaceCell(void) { return spaceCell; }
		virtual double getRadius(void) = 0;
		double getLastCollisionTime(void) { return lastCollisionTime; };

		void setX(double x) { position.x = x; };
		void setY(double y) { position.y = y; };
		void setZ(double z) { position.z = z; };
		void setPosition(point_t p) { position = p; };
		void setVx(double vx) { velocity.x = vx; };
		void setVy(double vy) { velocity.y = vy; };
		void setVz(double vz) { velocity.z = vz; };
		void setVelocity(vect_t v) { velocity = v; };
		void setMass(double m) { mass = m; };
		void setSpaceCell(int c) { spaceCell = c; };
		virtual void setRadius(double) = 0;
		void setLastCollisionTime(double tc) { lastCollisionTime = tc; };

		// tk Environment methods
		virtual void tkEnvDrawShape(void) = 0;
		virtual void tkEnvUpdatePosition(void) = 0;
		virtual void tkEnvUpdatePosition(double) = 0;
		virtual void firstEventTime() = 0;
		virtual void nextEventTime() = 0;
		virtual double scheduledCollisionTime() = 0;

};

#endif
