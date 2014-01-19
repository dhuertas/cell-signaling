#ifndef PARTICLE_H
#define PARTICLE_H

#include "Defines.h"
#include <vector>

class Particle {

	private:

	protected:

		point_t position;

		vect_t velocity;

		double mass;

		double lastCollisionTime;

		int particleId;

		int particleType;

		// Cell List attributes
		int spaceCell;	// TODO remove this and start using spaceCellIdx instead.

		index_t spaceCellIdx;

		int prevSpaceCell;

		index_t prevSpaceCellIdx;

		// Near Neighbor List of particles
		std::vector<Particle*> neighborParticles;

		// Near Neighbor List radius
		double listRadius;

		// Boundaries mode. Whether the particle must bounce, expire, etc.
		int boundariesMode;

	public:

		Particle();

		Particle(double, double, double, double, double);

		// Gets and sets
		inline int getParticleId(void) { return particleId; };

		inline int getParticleType(void) { return particleType; };

		inline double getX(void) { return position.x; };

		inline double getY(void) { return position.y; };

		inline double getZ(void) { return position.z; };

		inline point_t getPosition(void) { return position; };

		inline double getVx(void) { return velocity.x; };

		inline double getVy(void) { return velocity.y; };

		inline double getVz(void) { return velocity.z; };

		inline vect_t getVelocity(void) { return velocity; };

		inline double getMass(void) { return mass; };

		inline int getSpaceCell(void) { return spaceCell; }; // TODO remove this. Only use spaceCellIdx

		inline index_t getSpaceCellIdx(void) { return spaceCellIdx; };

		inline int getPrevSpaceCell(void) { return prevSpaceCell; }; // TODO remove this. Only use prevSpaceCellIdx

		inline index_t getPrevSpaceCellIdx(void) {return prevSpaceCellIdx; };

		inline double getListRadius(void) { return listRadius; };

		virtual double getRadius(void) = 0;

		inline double getLastCollisionTime(void) { return lastCollisionTime; };

		inline std::vector<Particle*> getNeighborParticles(void) { return neighborParticles; };

		inline int getBoundariesMode(void) const { return boundariesMode; };

		void setParticleId(int id) { particleId = id; };

		void setParticleType(int t) { particleType = t; };

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

		void setSpaceCellIdx(index_t sci) { spaceCellIdx = sci; };

		void setPrevSpaceCell(int c) { prevSpaceCell = c; };

		void setListRadius(double r) { listRadius = r; };

		virtual void setRadius(double) = 0;

		void setLastCollisionTime(double tc) { lastCollisionTime = tc; };

		void setBoundariesMode(int bm) { boundariesMode = bm; };

		//
		// tk Environment methods
		//
		virtual void tkEnvDrawShape(void) = 0;

		virtual void tkEnvUpdatePosition(void) = 0;

		virtual void tkEnvUpdatePosition(double) = 0;

		//
		// Event related methods
		//
		virtual void initMobilityMessages() = 0;

		virtual void initializeMobility() = 0;

		virtual void finishMobility() = 0;

		virtual void createNearNeighborList() = 0;

		virtual void updateNearNeighborList() = 0;

		virtual void expire() = 0;
};

#endif
