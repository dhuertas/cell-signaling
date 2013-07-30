#ifndef PARTICLE_H
#define PARTICLE_H

#include "Defines.h"
#include <list>

class Particle {

	private:

	protected:

		point_t position;
		vect_t velocity;

		double mass;
		double lastCollisionTime;

		int identifier;

		int mode;

// Cell List attributes
		int spaceCell;

// Near Neighbor List attributes
		std::list<Particle*> neighborParticles;
		double listRadius;

	public:

		Particle();
		Particle(double, double, double, double, double);

		// Gets and sets
		int getIdentifier(void) { return identifier; };
		double getX(void) { return position.x; };
		double getY(void) { return position.y; };
		double getZ(void) { return position.z; };
		point_t getPosition(void) { return position; };
		double getVx(void) { return velocity.x; };
		double getVy(void) { return velocity.y; };
		double getVz(void) { return velocity.z; };
		vect_t getVelocity(void) { return velocity; };
		double getMass(void) { return mass; };
		int getMode(void) { return mode; };
		double getSpaceCell(void) { return spaceCell; };
		double getListRadius(void) { return listRadius; };
		virtual double getRadius(void) = 0;
		double getLastCollisionTime(void) { return lastCollisionTime; };
		std::list<Particle*> getNeighborParticles(void) { return neighborParticles; };
		
		void setIdentifier(int id) { identifier = id; };
		void setX(double x) { position.x = x; };
		void setY(double y) { position.y = y; };
		void setZ(double z) { position.z = z; };
		void setPosition(point_t p) { position = p; };
		void setVx(double vx) { velocity.x = vx; };
		void setVy(double vy) { velocity.y = vy; };
		void setVz(double vz) { velocity.z = vz; };
		void setVelocity(vect_t v) { velocity = v; };
		void setMass(double m) { mass = m; };
		void setMode(int m) { mode = m; };
		void setSpaceCell(int c) { spaceCell = c; };
		void setListRadius(double r) { listRadius = r; };
		virtual void setRadius(double) = 0;
		void setLastCollisionTime(double tc) { lastCollisionTime = tc; };

		// tk Environment methods
		virtual void tkEnvDrawShape(void) = 0;
		virtual void tkEnvUpdatePosition(void) = 0;
		virtual void tkEnvUpdatePosition(double) = 0;
		virtual void firstEventTime() = 0;
		virtual void nextEventTime() = 0;
		virtual double scheduledCollisionTime() = 0;

		virtual void createNearNeighborList() = 0;
		virtual void updateNearNeighborList() = 0;
};

#endif
